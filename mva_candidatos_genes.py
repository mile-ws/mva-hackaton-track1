#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Candidatos en genes asociados a MVA (BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5, BUB1, BUB3).

Qué hace:
  1. Pregunta a Ensembl (base pública, hg38) dónde están los genes y sus exones.
     No se envía ningún dato del paciente: solo nombres de genes.
  2. Consulta tu VCF solo en esas regiones (usa el índice .tbi, es rápido).
  3. Se queda con variantes que caen en la parte codificante de un exón o a ±8 pb de un
     exón (zona de splicing) y que pasan un control de calidad básico. Las que caen en
     las UTR (extremos del gen que no se traducen a proteína) se cuentan aparte.
  4. Busca genes con 2 o más variantes candidatas (posible compuesto heterocigota)
     y, si GATK dejó información de fase (PGT/PID), dice si están en la misma copia
     del cromosoma (cis) o en copias distintas (trans). Una pareja en trans es la
     que sugiere herencia recesiva compuesta.

Qué NO hace (limitaciones importantes):
  - No anota el efecto en la proteína (stop, missense, etc.). Solo marca indels que
    "podrían" desplazar el marco de lectura. Para anotar de verdad hace falta VEP.
  - No conoce la frecuencia poblacional (gnomAD): una variante común puede aparecer
    como candidata. Hay que revisar las pocas que queden.
  - Solo mira los genes de la lista GENES (los asociados a MVA). Alrededor de la mitad
    de los pacientes con MVA no tiene un gen identificado, así que el gen causante
    podría no estar en la lista. Si ninguno sirve, hay que abrir la búsqueda.

Salidas (carpeta Escritorio/mva_hackathon_salida):
  - resumen_genes_COMPARTIBLE.txt : solo conteos por gen. Es lo único que podés pasar por chat.
  - LOCAL_NO_COMPARTIR_variantes_genes.xlsx : variantes del niño. Queda SOLO en tu compu.
    No lo subas a GitHub ni lo pegues en chats, y borralo al cerrar el hackatón.

Uso en Spyder: abrilo y apretá F5.
Requisitos: pip install openpyxl   (requests ya viene con huggingface_hub)
"""

import json
from itertools import combinations
from pathlib import Path

import pandas as pd
import pysam
import requests
from huggingface_hub import HfApi, hf_hub_download

# --------------------------------------------------------------------------- #
# Configuración
# --------------------------------------------------------------------------- #
REPO_ID = "SageBio/mva-hackathon-2026-data"
REPO_TYPE = "dataset"
_VCF_ENCONTRADO = None


def nombre_vcf() -> str:
    """Busca el VCF en el dataset. Su nombre no se escribe en el código porque lleva
    identificadores de muestra del laboratorio, y el código se publica."""
    global _VCF_ENCONTRADO
    if _VCF_ENCONTRADO is None:
        vcfs = [f for f in HfApi().list_repo_files(REPO_ID, repo_type=REPO_TYPE)
                if f.endswith(".vcf.gz")]
        if len(vcfs) != 1:
            raise RuntimeError(f"Esperaba un solo VCF en el dataset y encontré {len(vcfs)}. "
                               "Revisá el login y el acceso aprobado.")
        _VCF_ENCONTRADO = vcfs[0]
    return _VCF_ENCONTRADO


def descargar_vcf(con_indice: bool = True) -> str:
    """Descarga (o toma de la caché) el VCF y, si se pide, su índice .tbi."""
    ruta = hf_hub_download(REPO_ID, nombre_vcf(), repo_type=REPO_TYPE)
    if con_indice:
        hf_hub_download(REPO_ID, nombre_vcf() + ".tbi", repo_type=REPO_TYPE)  # queda junto al VCF
    return ruta

# nombre a mostrar -> posibles nombres en Ensembl (CENATAC antes se llamaba CCDC84)
GENES = {
    "BUB1B": ["BUB1B"],                 # MVA1
    "CEP57": ["CEP57"],                 # MVA2
    "TRIP13": ["TRIP13"],               # MVA3
    "CENATAC": ["CENATAC", "CCDC84"],   # MVA4
    "SLF2": ["SLF2", "FAM178A"],        # MVA5
    "SMC5": ["SMC5"],                   # MVA6
    "BUB1": ["BUB1"],                   # asociado a MVA, sin subtipo
    "BUB3": ["BUB3"],                   # asociado a MVA, sin subtipo
}

MARGEN_SPLICE = 8      # pares de bases alrededor de cada exón que se consideran
MIN_QUAL = 30          # calidad mínima de la variante
MIN_DP = 8             # profundidad mínima de lectura en la muestra
MIN_GQ = 20            # calidad mínima del genotipo
SOLO_TRANSCRITO_CANONICO = True  # False = unir exones de todos los transcritos

ENSEMBL = "https://rest.ensembl.org"
CARPETA_SALIDA = Path.home() / "Desktop" / "mva_hackathon_salida"


# --------------------------------------------------------------------------- #
# Ensembl: coordenadas de genes y exones (información pública, no del paciente)
# --------------------------------------------------------------------------- #
def buscar_gen_ensembl(nombres: list[str]) -> dict:
    ultimo_error = None
    for nombre in nombres:
        r = requests.get(
            f"{ENSEMBL}/lookup/symbol/homo_sapiens/{nombre}",
            params={"expand": 1},
            headers={"Content-Type": "application/json"},
            timeout=60,
        )
        if r.status_code == 200:
            datos = r.json()
            assembly = datos.get("assembly_name")
            if assembly and assembly != "GRCh38":
                print(f"  AVISO: Ensembl devolvió {assembly}, se esperaba GRCh38.")
            return datos
        ultimo_error = f"{nombre}: HTTP {r.status_code}"
    raise RuntimeError(f"No pude obtener el gen en Ensembl ({ultimo_error}).")


def cargar_genes(carpeta: Path) -> dict:
    """Usa un caché local; solo consulta Ensembl por los genes que todavía no están."""
    cache = carpeta / "genes_ensembl_cache.json"
    genes = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    faltan = [n for n in GENES if n not in genes]
    for nombre in faltan:
        print(f"Buscando {nombre} en Ensembl...")
        genes[nombre] = buscar_gen_ensembl(GENES[nombre])
    if faltan:
        cache.write_text(json.dumps(genes), encoding="utf-8")
    return {n: genes[n] for n in GENES}


def exones_del_gen(info_gen: dict, solo_canonico: bool = SOLO_TRANSCRITO_CANONICO) -> list[dict]:
    transcritos = info_gen.get("Transcript", [])
    if solo_canonico:
        canonicos = [t for t in transcritos if t.get("is_canonical")]
        if canonicos:
            transcritos = canonicos[:1]
    vistos: dict[tuple, dict] = {}
    for t in transcritos:
        exones = sorted(t.get("Exon", []), key=lambda e: e["start"])
        if t.get("strand", info_gen.get("strand", 1)) == -1:
            exones = exones[::-1]  # numerar en el sentido de lectura del gen
        for i, ex in enumerate(exones, start=1):
            clave = (ex["start"], ex["end"])
            if clave not in vistos:
                vistos[clave] = {
                    "start": ex["start"],
                    "end": ex["end"],
                    "exon_n": i,
                    "transcrito": t.get("id"),
                }
    return sorted(vistos.values(), key=lambda e: e["start"])


# --------------------------------------------------------------------------- #
# Lectura del VCF
# --------------------------------------------------------------------------- #
def nombre_contig(vcf: pysam.VariantFile, cromosoma: str) -> str:
    contigs = set(vcf.header.contigs)
    if cromosoma in contigs:
        return cromosoma
    if "chr" + cromosoma in contigs:
        return "chr" + cromosoma
    raise KeyError(f"El cromosoma {cromosoma} no está en el VCF.")


def _get(muestra, campo):
    try:
        return muestra.get(campo)
    except (KeyError, ValueError):
        return None


def distancia_a_exones(ini: int, fin: int, exones: list[dict]) -> tuple[int | None, int | None]:
    """Devuelve (distancia al exón más cercano, número de exón). 0 = dentro del exón."""
    mejor_d, mejor_n = None, None
    for ex in exones:
        if fin >= ex["start"] and ini <= ex["end"]:
            d = 0
        elif fin < ex["start"]:
            d = ex["start"] - fin
        else:
            d = ini - ex["end"]
        if mejor_d is None or d < mejor_d:
            mejor_d, mejor_n = d, ex["exon_n"]
    return mejor_d, mejor_n


def cds_del_gen(info_gen: dict) -> tuple[int, int] | None:
    """Inicio y fin (genómicos) de la parte que codifica proteína en el transcrito canónico."""
    transcritos = info_gen.get("Transcript", [])
    canonicos = [t for t in transcritos if t.get("is_canonical")] or transcritos[:1]
    if not canonicos:
        return None
    tr = canonicos[0].get("Translation")
    if not tr or "start" not in tr or "end" not in tr:
        return None
    return min(tr["start"], tr["end"]), max(tr["start"], tr["end"])


ZONAS_CANDIDATAS = ("codificante", "borde de exón (splicing)", "exónica (sin CDS)")


def zona_variante(ini: int, fin: int, exones: list[dict], cds: tuple[int, int] | None,
                  dist: int | None) -> str:
    """Clasifica dónde cae la variante respecto del transcrito canónico."""
    if not exones or dist is None:
        return "sin exones"
    tx_ini = min(e["start"] for e in exones)
    tx_fin = max(e["end"] for e in exones)
    if dist == 0:
        if cds is None:
            return "exónica (sin CDS)"
        return "codificante" if (fin >= cds[0] and ini <= cds[1]) else "UTR"
    if fin < tx_ini or ini > tx_fin:
        return "fuera del transcrito"
    if dist <= MARGEN_SPLICE:
        return "borde de exón (splicing)"
    return "intrónica"


def leer_variantes_gen(vcf: pysam.VariantFile, nombre_gen: str, info_gen: dict,
                       exones: list[dict], cds: tuple[int, int] | None = None) -> list[dict]:
    muestra_id = list(vcf.header.samples)[0]
    contig = nombre_contig(vcf, str(info_gen["seq_region_name"]))
    filas = []
    for rec in vcf.fetch(contig, info_gen["start"] - 1, info_gen["end"]):
        s = rec.samples[muestra_id]
        gt = _get(s, "GT")
        if not gt:
            continue
        alelos_alt = [a for a in gt if a not in (None, 0)]
        if not alelos_alt:
            continue  # llamada de referencia
        alts = [rec.alts[a - 1] for a in alelos_alt]
        ref, alt = rec.ref, alts[0]
        dlen = len(alt) - len(ref)
        if alt.startswith("<") or alt == "*":
            tipo = "otro"
        elif len(ref) == 1 and len(alt) == 1:
            tipo = "SNV"
        elif dlen == 0:
            tipo = "MNV"
        elif dlen > 0:
            tipo = "inserción"
        else:
            tipo = "deleción"

        ini, fin = rec.start + 1, rec.stop
        dist, exon_n = distancia_a_exones(ini, fin, exones)
        en_exon = dist == 0
        cerca_exon = dist is not None and dist <= MARGEN_SPLICE
        zona = zona_variante(ini, fin, exones, cds, dist)
        indel_no_mult3 = tipo in ("inserción", "deleción") and abs(dlen) % 3 != 0

        filtros = list(rec.filter.keys())
        filtro_ok = (not filtros) or filtros == ["PASS"]
        dp, gq = _get(s, "DP"), _get(s, "GQ")
        pasa_qc = bool(
            filtro_ok
            and (rec.qual is None or rec.qual >= MIN_QUAL)
            and (dp is None or dp >= MIN_DP)
            and (gq is None or gq >= MIN_GQ)
        )
        sep = "|" if s.phased else "/"
        gt_txt = sep.join("." if a is None else str(a) for a in gt)
        es_het = len(set(gt)) > 1  # alelos distintos (0/1, 1/2); 1/1 y hemicigota = False

        filas.append(
            {
                "gen": nombre_gen,
                "CHROM": rec.chrom,
                "POS": rec.pos,
                "REF": ref,
                "ALT": ",".join(rec.alts),
                "tipo": tipo,
                "GT": gt_txt,
                "heterocigota": bool(es_het),
                "QUAL": rec.qual,
                "FILTER": ",".join(filtros) or ".",
                "DP": dp,
                "GQ": gq,
                "AD": ",".join(str(x) for x in _get(s, "AD") or []),
                "en_dbSNP": "DB" in rec.info,
                "PGT": _get(s, "PGT"),
                "PID": _get(s, "PID"),
                "dist_a_exon": dist,
                "exon_n": exon_n if cerca_exon else None,
                "exonica": en_exon,
                "exon_o_borde": cerca_exon,
                "zona": zona,
                "candidata": pasa_qc and zona in ZONAS_CANDIDATAS,
                "indel_no_mult3": indel_no_mult3,
                "posible_frameshift": indel_no_mult3 and zona in ("codificante", "exónica (sin CDS)"),
                "pasa_QC": pasa_qc,
            }
        )
    return filas


# --------------------------------------------------------------------------- #
# Parejas (compuesto heterocigota) y fase
# --------------------------------------------------------------------------- #
def relacion_fase(a: dict, b: dict) -> str:
    """Usa PGT/PID de GATK: mismo PID y PGT distinto = copias distintas (trans)."""
    if a["PID"] and b["PID"] and a["PID"] == b["PID"] and a["PGT"] and b["PGT"]:
        return "cis (misma copia)" if a["PGT"] == b["PGT"] else "trans (copias distintas)"
    return "sin información de fase"


def parejas_candidatas(filas_gen: list[dict]) -> list[dict]:
    cands = [f for f in filas_gen if f["candidata"] and f["heterocigota"]]
    parejas = []
    for a, b in combinations(cands, 2):
        parejas.append(
            {
                "gen": a["gen"],
                "variante_1": f'{a["CHROM"]}:{a["POS"]} {a["REF"]}>{a["ALT"]}',
                "variante_2": f'{b["CHROM"]}:{b["POS"]} {b["REF"]}>{b["ALT"]}',
                "distancia_pb": abs(b["POS"] - a["POS"]),
                "fase": relacion_fase(a, b),
                "alguna_posible_frameshift": bool(a["posible_frameshift"] or b["posible_frameshift"]),
                "alguna_en_dbSNP": bool(a["en_dbSNP"] and b["en_dbSNP"]),
            }
        )
    return parejas


def resumen_gen(nombre: str, info_gen: dict, filas: list[dict], parejas: list[dict]) -> dict:
    cand = [f for f in filas if f["candidata"]]
    largo_kb = (info_gen["end"] - info_gen["start"] + 1) / 1000
    return {
        "gen": nombre,
        "region_hg38": f'{info_gen["seq_region_name"]}:{info_gen["start"]}-{info_gen["end"]}',
        "largo_gen_kb": round(largo_kb, 1),
        "variantes_en_el_gen": len(filas),
        "variantes_por_kb": round(len(filas) / largo_kb, 2) if largo_kb else None,
        "pasan_QC": sum(f["pasa_QC"] for f in filas),
        "en_UTR_excluidas": sum(f["pasa_QC"] and f["zona"] == "UTR" for f in filas),
        "candidatas_codif_o_splice": len(cand),
        "heterocigotas_candidatas": sum(f["heterocigota"] for f in cand),
        "homocigotas_candidatas": sum(not f["heterocigota"] for f in cand),
        "posibles_frameshift": sum(f["posible_frameshift"] for f in cand),
        "parejas_posibles": len(parejas),
        "parejas_trans_faseadas": sum(p["fase"].startswith("trans") for p in parejas),
    }


# --------------------------------------------------------------------------- #
# Salidas
# --------------------------------------------------------------------------- #
def guardar_excel(ruta: Path, hojas: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(ruta, engine="openpyxl") as writer:
        for nombre, tabla in hojas.items():
            tabla.to_excel(writer, sheet_name=nombre, index=False)
            ws = writer.sheets[nombre]
            for col in ws.columns:
                ancho = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(ancho + 2, 60)


def texto_compartible(resumen: pd.DataFrame) -> str:
    lineas = [
        "Resumen por gen (solo conteos, sin coordenadas de variantes del paciente)",
        "Filtros: PASS/sin filtro, QUAL>=%d, DP>=%d, GQ>=%d; candidatas = parte codificante"
        " o borde de exón ±%d pb (las UTR se cuentan aparte)"
        % (MIN_QUAL, MIN_DP, MIN_GQ, MARGEN_SPLICE),
        "",
    ]
    cols = [c for c in resumen.columns if c != "region_hg38"]
    lineas.append(resumen[cols].to_string(index=False))
    return "\n".join(lineas)


def main() -> None:
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)

    ruta_vcf = descargar_vcf()
    vcf = pysam.VariantFile(ruta_vcf)

    genes = cargar_genes(CARPETA_SALIDA)

    todas, todas_parejas, resumenes = [], [], []
    for nombre, info_gen in genes.items():
        exones = exones_del_gen(info_gen)
        cds = cds_del_gen(info_gen)
        if cds is None:
            print(f"  AVISO: {nombre} sin parte codificante en Ensembl; se usan exones completos.")
        filas = leer_variantes_gen(vcf, nombre, info_gen, exones, cds)
        parejas = parejas_candidatas(filas)
        todas += filas
        todas_parejas += parejas
        resumenes.append(resumen_gen(nombre, info_gen, filas, parejas))
        print(f"{nombre}: {len(filas)} variantes en el gen, {len(parejas)} parejas posibles")

    resumen = pd.DataFrame(resumenes)
    variantes = pd.DataFrame(todas)
    parejas = pd.DataFrame(todas_parejas)

    (CARPETA_SALIDA / "resumen_genes_COMPARTIBLE.txt").write_text(
        texto_compartible(resumen), encoding="utf-8"
    )
    guardar_excel(
        CARPETA_SALIDA / "LOCAL_NO_COMPARTIR_variantes_genes.xlsx",
        {
            "resumen": resumen,
            "parejas": parejas if len(parejas) else pd.DataFrame({"info": ["sin parejas"]}),
            "variantes": variantes if len(variantes) else pd.DataFrame({"info": ["sin variantes"]}),
        },
    )
    print("\nListo. Archivos en:", CARPETA_SALIDA)
    print("Compartible por chat: resumen_genes_COMPARTIBLE.txt")
    print("NO compartir: LOCAL_NO_COMPARTIR_variantes_genes.xlsx")


if __name__ == "__main__":
    main()
