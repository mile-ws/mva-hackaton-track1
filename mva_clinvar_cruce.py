#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cruce del VCF del paciente con ClinVar en TODO el genoma (corre en tu computadora).

Para qué sirve:
  Chequear que no se nos esté escapando una variante ya conocida como patogénica
  fuera de los 8 genes de MVA que revisamos.

Qué hace:
  1. Descarga el archivo público de ClinVar para GRCh38 (unos 200 MB, una sola vez).
  2. Se queda con las variantes que ClinVar clasifica como patogénicas, probablemente
     patogénicas o "conflictivas" con algún voto patogénico.
  3. Recorre el VCF del paciente y busca coincidencias exactas (cromosoma, posición,
     letra de referencia y letra alternativa). No sale ningún dato del paciente:
     todo se compara localmente.
  4. Agrupa por gen y clasifica:
       - homocigota patogénica           -> candidata recesiva directa
       - 2 o más heterocigotas           -> posible compuesto heterocigota
       - 1 heterocigota                  -> busca un "segundo golpe" en ese gen
                                            (otra variante codificante o de splicing)
  5. CONTROL POSITIVO: tiene que encontrar la variante stop de BUB1B (patogénica en
     ClinVar) y, como segundo golpe, la missense de BUB1B. Si no, algo anda mal.

Limitaciones:
  - Solo encuentra variantes que alguien ya cargó en ClinVar. Una variante nueva no aparece
    como primer golpe (sí puede aparecer como "segundo golpe").
  - Los indels pueden estar escritos distinto en los dos archivos y no coincidir.
  - Para el segundo golpe no hay frecuencia poblacional: se marca si está en dbSNP
    y hay que revisar en gnomAD las pocas que queden.

Salidas (carpeta Escritorio/mva_hackathon_salida):
  - resumen_clinvar_COMPARTIBLE.txt : solo nombres de genes y conteos.
  - LOCAL_NO_COMPARTIR_clinvar.xlsx : variantes del paciente. Queda SOLO en tu compu.

Uso en Spyder: poné este archivo en la MISMA carpeta que mva_candidatos_genes.py y F5.
Tarda unos minutos (recorre millones de variantes). La consola muestra el avance.
"""

import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
import requests
from huggingface_hub import hf_hub_download

try:
    _CARPETA_SCRIPT = Path(__file__).resolve().parent
except NameError:  # por si se corre pegado en la consola
    _CARPETA_SCRIPT = Path.cwd()
sys.path.insert(0, str(_CARPETA_SCRIPT))
import mva_candidatos_genes as base  # noqa: E402  (reutiliza Ensembl, exones y filtros)

# --------------------------------------------------------------------------- #
# Configuración
# --------------------------------------------------------------------------- #
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CARPETA_SALIDA = base.CARPETA_SALIDA
INCLUIR_CONFLICTIVAS = True   # conflictivas con al menos un voto patogénico
UMBRAL_COMUN = 0.01           # frecuencia (ExAC/1000G/ESP de ClinVar) para marcar "común"
BUSCAR_SEGUNDO_GOLPE = True
MAX_GENES_SEGUNDO_GOLPE = 150  # tope de consultas a Ensembl
GEN_CONTROL = "BUB1B"

PRIORIDAD = {
    "homocigota": 0,
    "2+ heterocigotas (posible compuesto)": 1,
    "1 heterocigota + segundo golpe": 2,
    "1 heterocigota, sin segundo golpe": 3,
    "1 heterocigota (segundo golpe no revisado)": 4,
}


# --------------------------------------------------------------------------- #
# Utilidades de ClinVar
# --------------------------------------------------------------------------- #
def normalizar_crom(c) -> str:
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "m") else c


def _texto(valor) -> str:
    """pysam devuelve tuplas en los campos Number=. (separa por comas); las unimos."""
    if valor is None:
        return ""
    if isinstance(valor, (tuple, list)):
        return ",".join(str(v) for v in valor)
    return str(valor)


def clasificar_clnsig(clnsig, clnsigconf=None) -> str | None:
    tokens = [t.strip().strip("_").lower() for t in re.split(r"[|/,;]", _texto(clnsig)) if t.strip()]
    if any(t.startswith("conflicting") for t in tokens):
        conf = _texto(clnsigconf).lower()
        return "conflictiva" if (not conf or "pathogenic" in conf) else None
    if "pathogenic" in tokens:
        return "patogénica"
    if "likely_pathogenic" in tokens:
        return "probablemente patogénica"
    return None


def estrellas(revstat) -> int:
    t = _texto(revstat).lower()
    if "practice_guideline" in t:
        return 4
    if "reviewed_by_expert_panel" in t:
        return 3
    if "multiple_submitters" in t and "no_conflicts" in t:
        return 2
    if "criteria_provided" in t:
        return 1
    return 0


def genes_de(geneinfo) -> list[str]:
    return [g.split(":")[0] for g in _texto(geneinfo).split("|") if g and g.split(":")[0]]


def consecuencias_de(mc) -> str:
    return ", ".join(sorted({p.split("|")[-1] for p in _texto(mc).split(",") if "|" in p}))


def descargar_clinvar(carpeta: Path) -> Path:
    destino = carpeta / "clinvar_GRCh38.vcf.gz"
    if destino.exists() and destino.stat().st_size > 1_000_000:
        print("ClinVar ya descargado:", destino)
        return destino
    print("Descargando ClinVar (público, ~200 MB)...")
    tmp = destino.with_suffix(".parcial")
    with requests.get(CLINVAR_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        bajado, ultimo = 0, 0
        with open(tmp, "wb") as f:
            for trozo in r.iter_content(chunk_size=1 << 20):
                f.write(trozo)
                bajado += len(trozo)
                if bajado - ultimo >= 20 << 20:
                    print(f"  {bajado / 1e6:.0f} MB")
                    ultimo = bajado
    tmp.rename(destino)
    print("ClinVar descargado:", destino)
    return destino


def _info(rec, campo):
    """Lee un campo INFO sin romperse si no está definido en el encabezado."""
    try:
        return rec.info.get(campo)
    except (KeyError, ValueError):
        return None


def cargar_clinvar(ruta: Path) -> dict:
    """(crom, pos, ref, alt) -> datos de ClinVar, solo para las clases que nos interesan.

    Si aparece el aviso 'Could not retrieve index file', se puede ignorar: acá el
    archivo se lee entero y no hace falta índice.
    """
    tabla = {}
    for rec in pysam.VariantFile(str(ruta)):
        if not rec.alts:
            continue
        clase = clasificar_clnsig(_info(rec, "CLNSIG"), _info(rec, "CLNSIGCONF"))
        if clase is None or (clase == "conflictiva" and not INCLUIR_CONFLICTIVAS):
            continue
        afs = [_info(rec, k) for k in ("AF_EXAC", "AF_TGP", "AF_ESP")]
        afs = [float(a) for a in afs if a is not None]
        datos = {
            "clinvar_id": rec.id,
            "clase_clinvar": clase,
            "clnsig": _texto(_info(rec, "CLNSIG")),
            "estrellas": estrellas(_info(rec, "CLNREVSTAT")),
            "genes": genes_de(_info(rec, "GENEINFO")),
            "consecuencia_clinvar": consecuencias_de(_info(rec, "MC")),
            "enfermedad": _texto(_info(rec, "CLNDN")).replace("_", " "),
            "af_max_clinvar": max(afs) if afs else None,
        }
        crom = normalizar_crom(rec.chrom)
        for alt in rec.alts:
            tabla[(crom, rec.pos, rec.ref.upper(), alt.upper())] = datos
    return tabla


# --------------------------------------------------------------------------- #
# Cruce con el paciente
# --------------------------------------------------------------------------- #
def pasa_qc(rec, muestra) -> bool:
    filtros = list(rec.filter.keys())
    filtro_ok = (not filtros) or filtros == ["PASS"]
    dp, gq = base._get(muestra, "DP"), base._get(muestra, "GQ")
    return bool(
        filtro_ok
        and (rec.qual is None or rec.qual >= base.MIN_QUAL)
        and (dp is None or dp >= base.MIN_DP)
        and (gq is None or gq >= base.MIN_GQ)
    )


def cruzar(vcf: pysam.VariantFile, clinvar: dict, cada: int = 500_000) -> list[dict]:
    posiciones = {(c, p) for (c, p, _, _) in clinvar}
    muestra_id = list(vcf.header.samples)[0]
    hits, n = [], 0
    for rec in vcf:
        n += 1
        if cada and n % cada == 0:
            print(f"  {n:,} variantes revisadas, {len(hits)} coincidencias")
        crom = normalizar_crom(rec.chrom)
        if (crom, rec.pos) not in posiciones:
            continue
        s = rec.samples[muestra_id]
        gt = base._get(s, "GT")
        if not gt or not rec.alts:
            continue
        idx_alt = sorted({a for a in gt if a not in (None, 0)})
        for a in idx_alt:
            alt = rec.alts[a - 1]
            datos = clinvar.get((crom, rec.pos, rec.ref.upper(), alt.upper()))
            if datos is None:
                continue
            copias = sum(1 for x in gt if x == a)
            hits.append(
                {
                    **{k: v for k, v in datos.items() if k != "genes"},
                    "genes": ",".join(datos["genes"]),
                    "CHROM": rec.chrom,
                    "POS": rec.pos,
                    "REF": rec.ref,
                    "ALT": alt,
                    "GT": ("|" if s.phased else "/").join("." if x is None else str(x) for x in gt),
                    "cigosidad": "homocigota" if copias >= 2 and len(set(gt)) == 1 else "heterocigota",
                    "QUAL": rec.qual,
                    "DP": base._get(s, "DP"),
                    "GQ": base._get(s, "GQ"),
                    "AD": ",".join(str(x) for x in base._get(s, "AD") or []),
                    "pasa_QC": pasa_qc(rec, s),
                    "comun": (datos["af_max_clinvar"] or 0) > UMBRAL_COMUN,
                }
            )
    print(f"Listo: {n:,} variantes del paciente revisadas, {len(hits)} coincidencias con ClinVar.")
    return hits


def agrupar_por_gen(hits: list[dict]) -> dict:
    genes = defaultdict(list)
    for h in hits:
        if not h["pasa_QC"]:
            continue
        for g in (h["genes"].split(",") if h["genes"] else ["(sin gen)"]):
            genes[g].append(h)
    return dict(genes)


def cargar_genes_lista(nombres: list[str], carpeta: Path) -> tuple[dict, dict]:
    """Comparte el caché con mva_candidatos_genes.py; devuelve (genes, errores)."""
    cache = carpeta / "genes_ensembl_cache.json"
    datos = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    errores, nuevos = {}, False
    for n in nombres:
        if n in datos:
            continue
        try:
            print(f"Buscando {n} en Ensembl...")
            datos[n] = base.buscar_gen_ensembl([n])
            nuevos = True
            time.sleep(0.1)
        except Exception as e:  # noqa: BLE001
            errores[n] = str(e)[:120]
    if nuevos:
        cache.write_text(json.dumps(datos), encoding="utf-8")
    return {n: datos[n] for n in nombres if n in datos}, errores


def segundo_golpe(vcf_indexado, gen: str, info_gen: dict, hit: dict) -> list[dict]:
    exones = base.exones_del_gen(info_gen)
    filas = base.leer_variantes_gen(vcf_indexado, gen, info_gen, exones, base.cds_del_gen(info_gen))
    return [
        {**f, "acompaña_a": f'{hit["CHROM"]}:{hit["POS"]}'}
        for f in filas
        if f["candidata"] and f["heterocigota"] and f["POS"] != hit["POS"]
    ]


def resumir(genes: dict, vcf_indexado, carpeta: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    filas_resumen, filas_segundas = [], []
    unicos = [g for g, hs in genes.items()
              if g != "(sin gen)" and len(hs) == 1 and hs[0]["cigosidad"] == "heterocigota"]
    a_revisar = unicos[:MAX_GENES_SEGUNDO_GOLPE] if BUSCAR_SEGUNDO_GOLPE else []
    info_genes, errores = cargar_genes_lista(a_revisar, carpeta) if a_revisar else ({}, {})

    for g, hs in genes.items():
        hom = [h for h in hs if h["cigosidad"] == "homocigota"]
        het = [h for h in hs if h["cigosidad"] == "heterocigota"]
        segundas = None
        if hom:
            cat = "homocigota"
        elif len(het) >= 2:
            cat = "2+ heterocigotas (posible compuesto)"
        elif g in info_genes:
            segundas = segundo_golpe(vcf_indexado, g, info_genes[g], het[0])
            filas_segundas += segundas
            cat = "1 heterocigota + segundo golpe" if segundas else "1 heterocigota, sin segundo golpe"
        else:
            cat = "1 heterocigota (segundo golpe no revisado)"
        filas_resumen.append(
            {
                "gen": g,
                "categoria": cat,
                "patogenicas": sum(h["clase_clinvar"] == "patogénica" for h in hs),
                "prob_patogenicas": sum(h["clase_clinvar"] == "probablemente patogénica" for h in hs),
                "conflictivas": sum(h["clase_clinvar"] == "conflictiva" for h in hs),
                "homocigotas": len(hom),
                "heterocigotas": len(het),
                "max_estrellas": max(h["estrellas"] for h in hs),
                "alguna_comun": any(h["comun"] for h in hs),
                "segundas_candidatas": None if segundas is None else len(segundas),
                "segundas_fuera_de_dbSNP": None if segundas is None else sum(not s["en_dbSNP"] for s in segundas),
                "nota": f"Ensembl: {errores[g]}" if g in errores else "",
            }
        )
    resumen = pd.DataFrame(filas_resumen)
    if len(resumen):
        for col in ("segundas_candidatas", "segundas_fuera_de_dbSNP"):
            resumen[col] = resumen[col].astype("Int64")  # enteros, vacío si no se revisó
        resumen["_p"] = resumen["categoria"].map(PRIORIDAD).fillna(9)
        resumen = (resumen.sort_values(["_p", "alguna_comun", "max_estrellas"], ascending=[True, True, False])
                   .drop(columns="_p").reset_index(drop=True))
    return resumen, pd.DataFrame(filas_segundas)


def texto_compartible(resumen: pd.DataFrame, n_hits: int, n_qc: int) -> str:
    lineas = [
        "Cruce con ClinVar (todo el genoma). Solo genes y conteos, sin coordenadas del paciente.",
        f"Clases incluidas: patogénica, probablemente patogénica"
        f"{', conflictiva con voto patogénico' if INCLUIR_CONFLICTIVAS else ''}.",
        f"Coincidencias con ClinVar: {n_hits} (pasan QC: {n_qc}). 'Común' = frecuencia > {UMBRAL_COMUN:.0%}.",
        "",
    ]
    if len(resumen):
        cat = resumen["categoria"].value_counts()
        for c in PRIORIDAD:
            if c in cat:
                lineas.append(f"  {c}: {cat[c]} genes")
        lineas.append("")
        control = resumen[resumen["gen"] == GEN_CONTROL]
        if len(control):
            fila = control.iloc[0]
            lineas.append(f"CONTROL POSITIVO {GEN_CONTROL}: OK, categoría '{fila['categoria']}', "
                          f"segundas candidatas = {fila['segundas_candidatas']}")
        else:
            lineas.append(f"CONTROL POSITIVO {GEN_CONTROL}: NO APARECE (revisar antes de sacar conclusiones)")
        lineas.append("")
        lineas.append(resumen.to_string(index=False))
    else:
        lineas.append("Sin coincidencias que pasen QC.")
    return "\n".join(lineas)


def main() -> None:
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)
    ruta_vcf = base.descargar_vcf()

    ruta_clinvar = descargar_clinvar(CARPETA_SALIDA)
    print("Leyendo ClinVar...")
    clinvar = cargar_clinvar(ruta_clinvar)
    print(f"  {len(clinvar):,} variantes de ClinVar en las clases elegidas.")

    print("Cruzando con el VCF del paciente (tarda unos minutos)...")
    hits = cruzar(pysam.VariantFile(ruta_vcf), clinvar)
    genes = agrupar_por_gen(hits)
    resumen, segundas = resumir(genes, pysam.VariantFile(ruta_vcf), CARPETA_SALIDA)

    n_qc = sum(h["pasa_QC"] for h in hits)
    texto = texto_compartible(resumen, len(hits), n_qc)
    (CARPETA_SALIDA / "resumen_clinvar_COMPARTIBLE.txt").write_text(texto, encoding="utf-8")
    base.guardar_excel(
        CARPETA_SALIDA / "LOCAL_NO_COMPARTIR_clinvar.xlsx",
        {
            "genes": resumen if len(resumen) else pd.DataFrame({"info": ["sin coincidencias"]}),
            "coincidencias": pd.DataFrame(hits) if hits else pd.DataFrame({"info": ["sin coincidencias"]}),
            "segundo_golpe": segundas if len(segundas) else pd.DataFrame({"info": ["sin segundas"]}),
        },
    )
    print("\n" + texto.split("\n\n")[0])
    print("\nArchivos en:", CARPETA_SALIDA)
    print("Compartible: resumen_clinvar_COMPARTIBLE.txt | NO compartir: LOCAL_NO_COMPARTIR_clinvar.xlsx")


if __name__ == "__main__":
    main()
