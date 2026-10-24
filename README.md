# MVA Hackathon 2026 — Track 1 (Variant Prediction)

**Team:** EMBS Equipo C1 (Milena W, Catalina W)

Pipeline used to propose the causal variant(s) for the proband in the
*Rare Disease, Real Kid: MVA Hackathon 2026* (Sage Bionetworks, MVA Society, Hugging Face).

> **No patient data is stored in this repository.** All scripts read the gated dataset
> from Hugging Face at runtime and write their outputs to local folders that are excluded
> by `.gitignore`. See *Data policy* below.

## Result in one sentence

A compound heterozygous pair in **BUB1B** (MVA1): a stop-gained variant classified as
pathogenic/likely pathogenic in ClinVar, plus an ultra-rare missense variant. Phase could not
be confirmed (single sample, no parents). Full reasoning in the methods report.

## Approach in one paragraph

The VCF (single sample, GRCh38, GATK HaplotypeCaller, unannotated) is filtered for call
quality and restricted to coding sequence or exon boundaries (±8 bp, splice region),
excluding UTRs. Under a recessive model (two heterozygous variants in the same gene, or one
homozygous variant) we (1) screened the 8 genes associated with MVA, (2) intersected the whole
VCF with ClinVar locally, (3) ran a genome-wide recessive sweep over 20,107 protein-coding
genes, and (4) ranked the resulting genes by phenotype (HPO) overlap, with a sensitivity
analysis. Every gene ranked at or above BUB1B was then reviewed at variant level in gnomAD v4
and ClinVar.

## Scripts (run in this order)

| Step | Script | What it does | Local outputs |
| --- | --- | --- | --- |
| 0 | `mva_exploracion.py` | Lists dataset files and reads the VCF header (genome build, INFO/FORMAT fields) | header summary (`.xlsx`) |
| 1 | `mva_candidatos_genes.py` | Gene panel (8 MVA genes: BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5, BUB1, BUB3). Ensembl GRCh38 coordinates, canonical transcript, CDS vs UTR, QC, pair detection with GATK phasing (PGT/PID) | per-gene counts, candidate variants |
| 2 | `mva_clinvar_cruce.py` | Downloads the ClinVar GRCh38 VCF and intersects it locally. Groups pathogenic/likely pathogenic/conflicting hits by gene and searches for a second candidate allele in single-hit genes. Includes a positive control | per-gene categories |
| 3 | `mva_barrido_genomico.py` | Genome-wide recessive sweep using GENCODE v50 (basic, Ensembl canonical transcripts) | strict and broad gene lists |
| 4 | `mva_cruce_hpo_genes.py` | Phenotype ranking: HPO terms × `phenotype_to_genes.txt`, weighted by term specificity, restricted to the genes from step 3 (`--only-genes`), with tie-aware ranks and sensitivity analysis (`--sensitivity-gene`) | gene ranking, sensitivity table |
| 5 | `mva_armar_csv.py` | Builds and validates the submission CSV directly from the local outputs (no manual copying of coordinates) | submission CSV |

### Main parameters

- QC: `FILTER` = PASS (or empty), `QUAL` ≥ 30, `DP` ≥ 8, `GQ` ≥ 20.
- Candidate region: coding sequence of the canonical transcript, or within 8 bp of an exon boundary. UTRs excluded.
- ClinVar classes kept: Pathogenic, Likely pathogenic, Conflicting (with at least one pathogenic vote).
- "Common" flag in ClinVar hits: max(AF_EXAC, AF_TGP, AF_ESP) > 1 %.
- Sweep, strict list: ≥2 heterozygous or ≥1 homozygous candidates absent from dbSNP (VCF `DB` flag). Broad list: any.
- Phenotype ranking, primary configuration: IDF weighting, no per-term cap, strict list.
- Manual review: variants above 5 % frequency in any gnomAD population considered benign (ACMG BA1).

## Reproducing

Requirements: Python 3.12+, an approved request for the gated dataset, and a Hugging Face read token.

```bash
python3 -m pip install -r requirements.txt
hf auth login                      # or: huggingface-cli login
python3 mva_exploracion.py
python3 mva_candidatos_genes.py
python3 mva_clinvar_cruce.py
python3 mva_barrido_genomico.py
python3 mva_cruce_hpo_genes.py --hpo-ids <local HPO file> --annotations phenotype_to_genes.txt \
    --only-genes genes_barrido_estricto.txt      # see the section below
python3 mva_armar_csv.py           # edit EQUIPO and FILAS at the top of the script first
```

Steps 0–3 and 5 write to `~/Desktop/mva_hackathon_salida/` (change `CARPETA_SALIDA` in
`mva_candidatos_genes.py`). Files whose name starts with `LOCAL_NO_COMPARTIR_` contain
patient-level variants and must never be shared. Files ending in `_COMPARTIBLE` contain only
gene names and counts. Steps 1–3 must be run from the same folder, because steps 2 and 3
import functions from `mva_candidatos_genes.py`.

### Input for the phenotype ranking (`mva_cruce_hpo_genes.py`)

The script needs a local file with the patient's HPO IDs. It is **not included** in this
repository because it describes the patient; it can be rebuilt from the phenotype table
provided by the organisers. Accepted formats:

- a CSV/TXT with a column named `HPO ID` (other columns are ignored), or
- a plain list with one HPO ID per line, or
- a spreadsheet (`.xlsx`) with a sheet `Síntomas HPO` (or the first sheet) and a column `HPO ID`.

Only the HPO IDs are read; any other column (clinical notes, gene lists) is ignored and should
not be shared. IDs must have the format `HP:` followed by 7 digits (e.g. `HP:0000001`).

The public file `phenotype_to_genes.txt` must be downloaded from HPO (see *External resources*).
Outputs are written to `salida_hpo/`, which is excluded by `.gitignore`.

To reproduce the sensitivity table of the report (rank of BUB1B under 3 weightings × 3 caps),
using the gene lists produced in step 3:

```bash
python3 mva_cruce_hpo_genes.py --hpo-ids <local HPO file> --annotations phenotype_to_genes.txt \
    --only-genes genes_barrido_estricto.txt --sensitivity-gene BUB1B

python3 mva_cruce_hpo_genes.py --hpo-ids <local HPO file> --annotations phenotype_to_genes.txt \
    --only-genes genes_barrido_amplio.txt --sensitivity-gene BUB1B
```

Run `python3 mva_cruce_hpo_genes.py --help` for all options.

## External resources

| Resource | Version | Accessed |
| --- | --- | --- |
| Ensembl REST API (gene/transcript coordinates) | GRCh38 | [DATE] |
| ClinVar `vcf_GRCh38/clinvar.vcf.gz` | release [CLINVAR RELEASE DATE] | [DATE] |
| GENCODE `gencode.v50.basic.annotation.gtf.gz` | v50 | [DATE] |
| HPO `phenotype_to_genes.txt` | — | 2026-10-20 |
| gnomAD browser (manual annotation) | v4 | 2026-10-19 to 2026-10-21 |
| OMIM / Orphanet (gene–disease relationships) | online | [DATE] |

## Limitations

- Single sample, no parents: the phase of the two candidate variants (cis vs trans) cannot be confirmed.
- The compound heterozygous hypothesis was inferred from the challenge's evaluation code, not from the clinical data.
- No population frequency is available locally; dbSNP membership is only a weak proxy for rarity, and is unreliable for indels.
- Synonymous and missense SNVs are not distinguished by the scripts (manual annotation does it for the shortlist).
- One transcript per gene (canonical); structural variants, deep intronic variants and X-chromosome hemizygosity are not modelled automatically.
- Phenotype ranking is consistent with, but not independent from, the gene-panel hypothesis (the phenotype is MVA-like).

## Data policy

- The dataset is gated and may not be redistributed. Nothing derived from it at genotype
  or phenotype level is committed here (see `.gitignore`).
- All patient data and derived files are deleted within 30 days of the hackathon close,
  and each team member confirms the deletion to the organisers.

## Use of AI assistants

Anthropic Claude via claude.ai — Pro plan (Mile) and [PLAN] plan (Cata), Consumer Terms,
model-training setting [ON/OFF]. Used for code, drafting, editing and translating the report
and this README, public literature, and interpreting public annotations (gnomAD, ClinVar) of a
small number of named candidate variants. All numbers were verified against our own outputs.
No dataset files, variant tables, HPO terms or clinical descriptions were shared.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in
partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation,
and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We
are deeply grateful to the child and their family who generously contributed their data
and their story to advance research into this rare disease. We acknowledge their trust in
making this Hackathon possible.

## License

Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as required for hackathon submissions.
