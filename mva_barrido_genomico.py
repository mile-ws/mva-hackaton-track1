#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Barrido genómico con modelo recesivo (corre en tu computadora).

Para qué sirve:
  Revisar TODOS los genes que codifican proteína, no solo los 8 de MVA, y quedarse con
  los que tienen el patrón de una enfermedad recesiva:
    - 2 o más variantes heterocigotas candidatas en el mismo gen (posible compuesto), o
    - 1 variante homocigota candidata.
  "Candidata" = pasa el control de calidad y cae en la parte codificante de un exón o en
  el borde de un exón (±8 pb, zona de splicing). Mismos criterios que mva_candidatos_genes.py.
  La lista de genes que sale sirve para restringir el ranking por síntomas de Cata
  (--only-genes).

Qué hace:
  1. Descarga la anotación pública de genes de GENCODE (GRCh38, ~50 MB, una sola vez).
  2. Para cada gen codificante toma el transcrito canónico (exones y parte codificante).
  3. Recorre el VCF completo y ubica cada variante respecto de esos exones.
  4. Agrupa por gen y arma categorías y listas de genes.
  5. CONTROL POSITIVO: BUB1B tiene que aparecer como "2+ heterocigotas".

Limitaciones (importante para el informe):
  - No hay frecuencia poblacional: va a haber MILES de genes con 2 variantes comunes e
    inofensivas. El barrido NO encuentra el gen: achica la lista para cruzarla con los
    síntomas, y después se revisan en gnomAD los pocos que queden arriba.
  - "Fuera de dbSNP" (campo DB del VCF) es una pista débil de que una variante es rara,
    pero muchas variantes patogénicas SÍ están en dbSNP. Por eso hay lista estricta y amplia.
  - No distingue sinónimas de missense (para eso hace falta la secuencia de referencia).
  - No considera el cromosoma X en varones (hemicigotas) ni variantes estructurales.

Salidas (carpeta Escritorio/mva_hackathon_salida):
  - resumen_barrido_COMPARTIBLE.txt   : conteos y genes, sin coordenadas.
  - genes_barrido_estricto.txt        : un gen por línea (para --only-genes).
  - genes_barrido_amplio.txt          : un gen por línea (para --only-genes).
  - LOCAL_NO_COMPARTIR_barrido.xlsx   : variantes del paciente. Queda SOLO en tu compu.

Uso en Spyder: misma carpeta que mva_candidatos_genes.py, y F5.
Tarda unos minutos: la primera vez se suma la descarga de GENCODE (~50 MB) y leer
la anotación lleva un par de minutos. La consola muestra el avance.
"""

import gzip
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
import requests
from huggingface_hub import hf_hub_download

try:
    _CARPETA_SCRIPT = Path(__file__).resolve().parent
except NameError:
    _CARPETA_SCRIPT = Path.cwd()
sys.path.insert(0, str(_CARPETA_SCRIPT))
import mva_candidatos_genes as base  # noqa: E402

# --------------------------------------------------------------------------- #
# Configuración
# --------------------------------------------------------------------------- #
GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
               "gencode.v50.basic.annotation.gtf.gz")
CARPETA_SALIDA = base.CARPETA_SALIDA
ARCHIVO_GTF = "gencode_basic_GRCh38.gtf.gz"
GEN_CONTROL = "BUB1B"
TOP_EN_RESUMEN = 60  # cuántos genes de cada categoría se listan en el resumen

CATEGORIAS = [  # en orden de prioridad
    "2+ heterocigotas fuera de dbSNP",
    "homocigota fuera de dbSNP",
    "2+ heterocigotas (alguna en dbSNP)",
    "homocigota en dbSNP",
]
ESTRICTAS = CATEGORIAS[:2]


# --------------------------------------------------------------------------- #
# Anotación de genes (GENCODE)
# --------------------------------------------------------------------------- #
def normalizar_crom(c) -> str:
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "m") else c


_ATRIB = re.compile(r'(\w+) "([^"]*)"')


def descargar_gtf(carpeta: Path) -> Path:
    destino = carpeta / ARCHIVO_GTF
    if destino.exists() and destino.stat().st_size > 1_000_000:
        print("Anotación GENCODE ya descargada:", destino)
        return destino
    print("Descargando anotación de genes GENCODE (pública, ~50 MB)...")
    tmp = destino.with_suffix(".parcial")
    with requests.get(GENCODE_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for trozo in r.iter_content(chunk_size=1 << 20):
                f.write(trozo)
    tmp.rename(destino)
    print("Descargada:", destino)
    return destino


def leer_gtf(ruta: Path) -> dict:
    """Devuelve {clave_gen: {nombre, crom, strand, exones, cds}} con el transcrito canónico."""
    transcritos: dict[str, dict] = {}
    abrir = gzip.open if str(ruta).endswith(".gz") else open
    with abrir(ruta, "rt", encoding="utf-8") as f:
        for linea in f:
            if linea.startswith("#"):
                continue
            campos = linea.rstrip("\n").split("\t")
            if len(campos) < 9 or campos[2] not in ("exon", "CDS"):
                continue
            atr = campos[8]
            if 'gene_type "protein_coding"' not in atr:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', atr)
            if not tid:
                continue
            t = transcritos.get(tid.group(1))
            if t is None:
                d = dict(_ATRIB.findall(atr))
                t = transcritos[tid.group(1)] = {
                    "tid": tid.group(1),
                    "gene_id": d.get("gene_id", ""),
                    "nombre": d.get("gene_name", d.get("gene_id", "")),
                    "crom": normalizar_crom(campos[0]),
                    "strand": -1 if campos[6] == "-" else 1,
                    "canonico": 'tag "Ensembl_canonical"' in atr,
                    "exones": [],
                    "cds": None,
                }
            ini, fin = int(campos[3]), int(campos[4])
            if 'tag "Ensembl_canonical"' in atr:
                t["canonico"] = True
            if campos[2] == "exon":
                t["exones"].append((ini, fin))
            else:
                c = t["cds"]
                t["cds"] = (ini, fin) if c is None else (min(c[0], ini), max(c[1], fin))

    por_gen: dict[str, list] = defaultdict(list)
    for t in transcritos.values():
        if t["exones"]:
            por_gen[t["gene_id"]].append(t)

    genes = {}
    for gid, ts in por_gen.items():
        canon = [t for t in ts if t["canonico"]]
        elegido = canon[0] if canon else max(
            ts, key=lambda t: (t["cds"][1] - t["cds"][0]) if t["cds"] else 0)
        ex = sorted(set(elegido["exones"]))
        orden = ex if elegido["strand"] == 1 else ex[::-1]
        numero = {e: i for i, e in enumerate(orden, start=1)}
        nombre = elegido["nombre"]
        clave = nombre if not gid.endswith("_PAR_Y") else f"{nombre}_PAR_Y"
        if clave in genes:  # mismo nombre en dos genes distintos: no pisar el anterior
            clave = f"{nombre}_{gid}"
        genes[clave] = {
            "nombre": nombre,
            "crom": elegido["crom"],
            "strand": elegido["strand"],
            "transcrito": elegido["tid"],
            "exones": [{"start": s, "end": e, "exon_n": numero[(s, e)], "transcrito": elegido["tid"]}
                       for s, e in ex],
            "cds": elegido["cds"],
        }
    return genes


class IndiceExones:
    """Índice por cromosoma de exones ±margen, para ubicar variantes rápido."""

    def __init__(self, genes: dict, margen: int):
        por_crom = defaultdict(list)
        for clave, g in genes.items():
            for e in g["exones"]:
                por_crom[g["crom"]].append((e["start"] - margen, e["end"] + margen, clave))
        self.datos = {}
        for c, ivs in por_crom.items():
            ivs.sort()
            pmax, m = [], -1
            for _, fin, _ in ivs:
                m = max(m, fin)
                pmax.append(m)
            self.datos[c] = ([i[0] for i in ivs], [i[1] for i in ivs], [i[2] for i in ivs], pmax)

    def genes_en(self, crom: str, ini: int, fin: int) -> set:
        d = self.datos.get(crom)
        if d is None:
            return set()
        starts, ends, claves, pmax = d
        i = bisect_right(starts, fin) - 1
        res = set()
        while i >= 0 and pmax[i] >= ini:
            if ends[i] >= ini:
                res.add(claves[i])
            i -= 1
        return res


# --------------------------------------------------------------------------- #
# Recorrido del VCF
# --------------------------------------------------------------------------- #
def datos_variante(rec, muestra) -> dict | None:
    gt = base._get(muestra, "GT")
    if not gt or not rec.alts:
        return None
    alelos = [a for a in gt if a not in (None, 0)]
    if not alelos:
        return None
    ref, alt = rec.ref, rec.alts[alelos[0] - 1]
    dlen = len(alt) - len(ref)
    if alt.startswith("<") or alt == "*":
        tipo = "otro"
    elif len(ref) == 1 and len(alt) == 1:
        tipo = "SNV"
    elif dlen == 0:
        tipo = "MNV"
    else:
        tipo = "inserción" if dlen > 0 else "deleción"
    filtros = list(rec.filter.keys())
    dp, gq = base._get(muestra, "DP"), base._get(muestra, "GQ")
    pasa = bool(((not filtros) or filtros == ["PASS"])
                and (rec.qual is None or rec.qual >= base.MIN_QUAL)
                and (dp is None or dp >= base.MIN_DP)
                and (gq is None or gq >= base.MIN_GQ))
    hom = len(set(gt)) == 1 and gt[0] not in (None, 0)
    return {
        "CHROM": rec.chrom, "POS": rec.pos, "REF": ref, "ALT": ",".join(rec.alts),
        "tipo": tipo,
        "GT": ("|" if muestra.phased else "/").join("." if a is None else str(a) for a in gt),
        "heterocigota": len(set(gt)) > 1, "homocigota": hom,
        "QUAL": rec.qual, "DP": dp, "GQ": gq,
        "AD": ",".join(str(x) for x in base._get(muestra, "AD") or []),
        "en_dbSNP": "DB" in rec.info,
        "indel_no_mult3": tipo in ("inserción", "deleción") and abs(dlen) % 3 != 0,
        "pasa_QC": pasa,
    }


def barrer(vcf: pysam.VariantFile, genes: dict, indice: IndiceExones, cada: int = 500_000) -> list[dict]:
    muestra_id = list(vcf.header.samples)[0]
    candidatas, n = [], 0
    for rec in vcf:
        n += 1
        if cada and n % cada == 0:
            print(f"  {n:,} variantes revisadas, {len(candidatas):,} candidatas")
        ini, fin = rec.start + 1, rec.stop
        claves = indice.genes_en(normalizar_crom(rec.chrom), ini, fin)
        if not claves:
            continue  # intergénica o intrónica profunda: se descarta rápido
        d = datos_variante(rec, rec.samples[muestra_id])
        if d is None or not d["pasa_QC"]:
            continue
        for clave in claves:
            g = genes[clave]
            dist, exon_n = base.distancia_a_exones(ini, fin, g["exones"])
            zona = base.zona_variante(ini, fin, g["exones"], g["cds"], dist)
            if zona not in base.ZONAS_CANDIDATAS:
                continue
            candidatas.append({
                "gen": clave, **d, "zona": zona, "exon_n": exon_n,
                "posible_frameshift": d["indel_no_mult3"] and zona in ("codificante", "exónica (sin CDS)"),
            })
    print(f"Listo: {n:,} variantes revisadas, {len(candidatas):,} candidatas en genes codificantes.")
    return candidatas


def categorizar(candidatas: list[dict]) -> pd.DataFrame:
    por_gen = defaultdict(list)
    for c in candidatas:
        por_gen[c["gen"]].append(c)
    filas = []
    for g, cs in por_gen.items():
        het = [c for c in cs if c["heterocigota"]]
        hom = [c for c in cs if c["homocigota"]]
        het_nodb = [c for c in het if not c["en_dbSNP"]]
        hom_nodb = [c for c in hom if not c["en_dbSNP"]]
        if len(het_nodb) >= 2:
            cat = CATEGORIAS[0]
        elif hom_nodb:
            cat = CATEGORIAS[1]
        elif len(het) >= 2:
            cat = CATEGORIAS[2]
        elif hom:
            cat = CATEGORIAS[3]
        else:
            cat = "1 heterocigota (no encaja en recesivo)"
        filas.append({
            "gen": g, "categoria": cat,
            "heterocigotas": len(het), "het_fuera_dbSNP": len(het_nodb),
            "homocigotas": len(hom), "hom_fuera_dbSNP": len(hom_nodb),
            "posibles_frameshift": sum(c["posible_frameshift"] for c in cs),
            "en_splicing": sum(c["zona"] == "borde de exón (splicing)" for c in cs),
        })
    df = pd.DataFrame(filas)
    if len(df):
        orden = {c: i for i, c in enumerate(CATEGORIAS)}
        df["_p"] = df["categoria"].map(orden).fillna(len(CATEGORIAS))
        df = (df.sort_values(["_p", "posibles_frameshift", "het_fuera_dbSNP", "gen"],
                             ascending=[True, False, False, True])
              .drop(columns="_p").reset_index(drop=True))
    return df


def texto_compartible(tabla: pd.DataFrame, n_genes_anotados: int) -> str:
    cuenta = tabla["categoria"].value_counts() if len(tabla) else {}
    lineas = [
        "Barrido genómico, modelo recesivo. Solo genes y conteos, sin coordenadas del paciente.",
        f"Genes codificantes anotados (GENCODE, transcrito canónico): {n_genes_anotados:,}",
        f"Candidata = QC (QUAL>={base.MIN_QUAL}, DP>={base.MIN_DP}, GQ>={base.MIN_GQ}, PASS) y "
        f"parte codificante o borde de exón ±{base.MARGEN_SPLICE} pb.",
        "",
    ]
    for c in CATEGORIAS:
        lineas.append(f"  {c}: {cuenta.get(c, 0):,} genes")
    estricto = sum(cuenta.get(c, 0) for c in ESTRICTAS)
    amplio = sum(cuenta.get(c, 0) for c in CATEGORIAS)
    lineas += ["", f"Lista estricta: {estricto:,} genes | lista amplia: {amplio:,} genes", ""]

    control = tabla[tabla["gen"] == GEN_CONTROL] if len(tabla) else tabla
    if len(control):
        f = control.iloc[0]
        lineas.append(f"CONTROL POSITIVO {GEN_CONTROL}: aparece como '{f['categoria']}' "
                      f"(het={f['heterocigotas']}, het fuera de dbSNP={f['het_fuera_dbSNP']}). "
                      f"En lista estricta: {'sí' if f['categoria'] in ESTRICTAS else 'no'}.")
    else:
        lineas.append(f"CONTROL POSITIVO {GEN_CONTROL}: NO APARECE (revisar antes de sacar conclusiones)")

    for c in ESTRICTAS:
        sub = tabla[tabla["categoria"] == c].head(TOP_EN_RESUMEN)
        if len(sub):
            lineas += ["", f"--- {c} (primeros {len(sub)}) ---", sub.to_string(index=False)]
    return "\n".join(lineas)


def main() -> None:
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)
    ruta_vcf = base.descargar_vcf()

    print("Leyendo anotación de genes (tarda un par de minutos)...")
    genes = leer_gtf(descargar_gtf(CARPETA_SALIDA))
    print(f"  {len(genes):,} genes codificantes.")
    indice = IndiceExones(genes, base.MARGEN_SPLICE)

    print("Recorriendo el VCF completo (esto es lo largo)...")
    candidatas = barrer(pysam.VariantFile(ruta_vcf), genes, indice)
    tabla = categorizar(candidatas)

    texto = texto_compartible(tabla, len(genes))
    (CARPETA_SALIDA / "resumen_barrido_COMPARTIBLE.txt").write_text(texto, encoding="utf-8")
    for nombre, cats in (("estricto", ESTRICTAS), ("amplio", CATEGORIAS)):
        lista = tabla[tabla["categoria"].isin(cats)]["gen"] if len(tabla) else []
        (CARPETA_SALIDA / f"genes_barrido_{nombre}.txt").write_text(
            "\n".join(lista) + "\n", encoding="utf-8")
    base.guardar_excel(
        CARPETA_SALIDA / "LOCAL_NO_COMPARTIR_barrido.xlsx",
        {"genes": tabla if len(tabla) else pd.DataFrame({"info": ["sin candidatas"]}),
         "candidatas": pd.DataFrame(candidatas) if candidatas else pd.DataFrame({"info": ["sin candidatas"]})},
    )
    print("\n" + "\n".join(texto.split("\n")[:12]))
    print("\nArchivos en:", CARPETA_SALIDA)


if __name__ == "__main__":
    main()
