#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv
from pathlib import Path

import pandas as pd

try:
    import mva_candidatos_genes as base
    CARPETA_SALIDA = base.CARPETA_SALIDA
except ImportError:  # si se corre suelto
    CARPETA_SALIDA = Path.home() / "Desktop" / "mva_hackathon_salida"

# --------------------------------------------------------------------------- #
# Completar acá
# --------------------------------------------------------------------------- #
EQUIPO = "nombre-del-equipo"          # el MISMO en todos los envíos
MODELO = "panel-mva-clinvar-barrido"  # descripción corta del método

FILAS = [
    {
        "gen": "BUB1B",
        "epcr": 0.9,
        "finding_type": "primary",
        "notes": "Stop gained P/LP en ClinVar (2 estrellas) + missense ultra rara; "
                 "compuesto heterocigota, fase no confirmada",
        "fuente": "candidatos",
    },
    # Ejemplo de alternativa (descomentar y completar si corresponde):
    # {"gen": "PEX5", "epcr": 0.05, "finding_type": "primary",
    #  "notes": "Alternativa por fenotipo; variantes a confirmar en gnomAD", "fuente": "barrido"},
]

PROBAND_ID = "PROBAND01"
COLUMNAS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
            "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
ARCHIVOS = {
    "candidatos": ("LOCAL_NO_COMPARTIR_variantes_genes.xlsx", "variantes"),
    "barrido": ("LOCAL_NO_COMPARTIR_barrido.xlsx", "candidatas"),
}


# --------------------------------------------------------------------------- #
# Armado
# --------------------------------------------------------------------------- #
def con_chr(c) -> str:
    c = str(c)
    return c if c.lower().startswith("chr") else "chr" + c


def alt_real(alt: str, gt: str) -> str:
    """Si el sitio es multialélico, devuelve el alelo que realmente tiene el paciente."""
    alts = str(alt).split(",")
    if len(alts) == 1:
        return alts[0]
    idx = sorted({int(a) for a in str(gt).replace("|", "/").split("/") if a.isdigit() and a != "0"})
    if len(idx) != 1:
        raise ValueError(f"Sitio multialélico con dos alelos alternativos ({gt}); "
                         "elegí la variante a mano con 'posiciones'.")
    return alts[idx[0] - 1]


def candidatas_gen(tabla: pd.DataFrame, gen: str, posiciones=None) -> pd.DataFrame:
    sub = tabla[tabla["gen"] == gen]
    if "candidata" in sub.columns:
        sub = sub[sub["candidata"] == True]  # noqa: E712  (viene de Excel)
    if "pasa_QC" in sub.columns:
        sub = sub[sub["pasa_QC"] == True]  # noqa: E712
    if posiciones:
        sub = sub[sub["POS"].isin(posiciones)]
    return sub.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"]).sort_values("POS")


def fila_para(entrada: dict, tabla: pd.DataFrame) -> dict:
    gen = entrada["gen"]
    sub = candidatas_gen(tabla, gen, entrada.get("posiciones"))
    het = sub[sub["heterocigota"] == True]  # noqa: E712
    hom = sub[sub["heterocigota"] == False]  # noqa: E712

    if len(het) == 2 and len(hom) == 0:
        v1, v2 = het.iloc[0], het.iloc[1]
    elif len(hom) == 1 and len(het) == 0:
        v1, v2 = hom.iloc[0], None
    else:
        raise ValueError(
            f"{gen}: hay {len(het)} heterocigotas y {len(hom)} homocigotas candidatas. "
            "Una fila necesita exactamente 2 heterocigotas (pareja) o 1 homocigota. "
            "Usá 'posiciones' para elegir cuáles."
        )
    fila = {
        "proband_id": PROBAND_ID,
        "chrom_1": con_chr(v1["CHROM"]), "pos_1": int(v1["POS"]),
        "ref_1": str(v1["REF"]).upper(), "alt_1": alt_real(v1["ALT"], v1["GT"]).upper(),
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "epcr": entrada["epcr"],
        "finding_type": entrada.get("finding_type", "primary"),
        "notes": entrada.get("notes", ""),
    }
    if v2 is not None:
        fila.update({
            "chrom_2": con_chr(v2["CHROM"]), "pos_2": int(v2["POS"]),
            "ref_2": str(v2["REF"]).upper(), "alt_2": alt_real(v2["ALT"], v2["GT"]).upper(),
        })
    return fila


def armar(filas_cfg: list[dict], carpeta: Path) -> list[dict]:
    tablas = {}
    filas = []
    for entrada in filas_cfg:
        fuente = entrada.get("fuente", "candidatos")
        if fuente not in tablas:
            archivo, hoja = ARCHIVOS[fuente]
            ruta = carpeta / archivo
            if not ruta.exists():
                raise FileNotFoundError(f"No encuentro {ruta}. Corré antes el script que lo genera.")
            tablas[fuente] = pd.read_excel(ruta, sheet_name=hoja)
        filas.append(fila_para(entrada, tablas[fuente]))
    return sorted(filas, key=lambda f: -float(f["epcr"]))


# --------------------------------------------------------------------------- #
# Validación (también sirve para revisar un CSV hecho a mano)
# --------------------------------------------------------------------------- #
def validar(filas: list[dict]) -> list[str]:
    errores = []
    if not filas:
        errores.append("El CSV no tiene filas.")
    if len(filas) > 10:
        errores.append(f"Tiene {len(filas)} filas; el máximo es 10.")
    vistas = set()
    for i, f in enumerate(filas, start=1):
        pre = f"Fila {i}:"
        if f.get("proband_id") != PROBAND_ID:
            errores.append(f"{pre} proband_id debe ser {PROBAND_ID}.")
        for k in ("chrom_1", "pos_1", "ref_1", "alt_1"):
            if f.get(k) in (None, ""):
                errores.append(f"{pre} falta {k}.")
        for n in ("1", "2"):
            c = f.get(f"chrom_{n}")
            if c and not str(c).startswith("chr"):
                errores.append(f"{pre} chrom_{n} debe llevar prefijo 'chr'.")
            for k in (f"ref_{n}", f"alt_{n}"):
                v = str(f.get(k) or "")
                if v and (not set(v) <= set("ACGTN") or "," in v):
                    errores.append(f"{pre} {k}='{v}' no es un alelo simple (solo A, C, G, T, N).")
        segunda = [f.get(k) for k in ("chrom_2", "pos_2", "ref_2", "alt_2")]
        if any(x not in (None, "") for x in segunda) and not all(x not in (None, "") for x in segunda):
            errores.append(f"{pre} la segunda variante está incompleta.")
        try:
            e = float(f.get("epcr"))
            if not 0 < e <= 1:
                errores.append(f"{pre} epcr={e} fuera de (0, 1].")
        except (TypeError, ValueError):
            errores.append(f"{pre} epcr no es un número.")
        if f.get("finding_type") not in ("primary", "secondary"):
            errores.append(f"{pre} finding_type debe ser 'primary' o 'secondary'.")
        clave = (f.get("chrom_1"), f.get("pos_1"), f.get("chrom_2"), f.get("pos_2"))
        if clave in vistas:
            errores.append(f"{pre} está repetida.")
        vistas.add(clave)
    epcrs = [float(f["epcr"]) for f in filas if str(f.get("epcr", "")).replace(".", "", 1).isdigit()]
    if epcrs != sorted(epcrs, reverse=True):
        errores.append("Las filas no están ordenadas por epcr de mayor a menor.")
    return errores


def escribir(filas: list[dict], ruta: Path) -> None:
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS)
        w.writeheader()
        for fila in filas:
            w.writerow({k: fila.get(k, "") for k in COLUMNAS})


def main() -> None:
    if EQUIPO == "nombre-del-equipo":
        print("AVISO: completá EQUIPO con el nombre real del equipo antes del envío final.")
    filas = armar(FILAS, CARPETA_SALIDA)
    errores = validar(filas)
    if errores:
        print("NO se generó el CSV. Errores:")
        for e in errores:
            print("  -", e)
        return
    ruta = CARPETA_SALIDA / f"{EQUIPO}_{MODELO}.csv"
    escribir(filas, ruta)
    print(f"CSV válido guardado en: {ruta}\n")
    print(pd.DataFrame(filas)[COLUMNAS].to_string(index=False))
    print("\nRevisá que las variantes coincidan con las que miraste en gnomAD antes de subirlo.")


if __name__ == "__main__":
    main()
