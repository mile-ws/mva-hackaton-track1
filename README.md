MVA Hackathon 2026 — Track 1 (Variant Prediction)

Team: EMBS Equipo C1 (Milena W, Catalina W)

Pipeline used to propose the causal variant(s) for the proband in the Rare Disease, Real Kid: MVA Hackathon 2026 (Sage Bionetworks, MVA Society, Hugging Face).

No patient data is stored in this repository. All scripts read the gated dataset from Hugging Face at runtime and write their outputs to a local folder that is excluded by .gitignore. See Data policy below.

Approach in one paragraph

The VCF (single sample, GRCh38, GATK HaplotypeCaller, unannotated) is filtered for call quality and restricted to coding sequence or exon boundaries (±8 bp, splice region), excluding UTRs. Under a recessive model (biallelic: two heterozygous variants in the same gene, or one homozygous variant) we (1) screened the 8 genes associated with MVA, (2) cross-referenced the whole VCF with ClinVar locally, (3) ran a genome-wide recessive sweep over 20,107 protein-coding genes, and (4) ranked the resulting genes by phenotype (HPO) overlap. Candidates were then annotated manually with gnomAD v4, ClinVar, CADD, REVEL and SpliceAI.

Scripts (run in this order)
Step	Script	What it does	Local outputs
0	mva_exploracion.py	Lists dataset files, reads VCF header (genome build, INFO/FORMAT fields)	header summary (.xlsx)
1	mva_candidatos_genes.py	Gene panel (8 MVA genes: BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5, BUB1, BUB3). Ensembl GRCh38 coordinates, canonical transcript, CDS vs UTR, QC, pair detection with GATK phasing (PGT/PID)	per-gene counts, candidate variants
2	mva_clinvar_cruce.py	Downloads ClinVar GRCh38 VCF and intersects locally. Groups P/LP/conflicting hits by gene and searches for a second hit in single-hit genes. Includes a positive control	per-gene categories
3	mva_barrido_genomico.py	Genome-wide recessive sweep using GENCODE v50 (basic, Ensembl canonical transcripts)	strict and broad gene lists
4	mva_cruce_hpo_genes.py	Phenotype-driven ranking: HPO terms × phenotype_to_genes.txt, IDF-weighted, restricted to genes from step 3 (--only-genes), with sensitivity analysis	gene ranking
5	mva_armar_csv.py	Builds and validates the submission CSV directly from the local outputs (no manual copying of coordinates)	submission CSV
Main parameters
QC: FILTER = PASS (or empty), QUAL ≥ 30, DP ≥ 8, GQ ≥ 20.
Candidate region: coding sequence of the canonical transcript, or within 8 bp of an exon boundary.
ClinVar classes kept: Pathogenic, Likely pathogenic, Conflicting (with at least one pathogenic vote).
"Common" flag in ClinVar hits: max(AF_EXAC, AF_TGP, AF_ESP) > 1 %.
Strict sweep list: ≥2 heterozygous or ≥1 homozygous candidates not in dbSNP (VCF DB flag); broad list: any.
Reproducing
bash
python3 -m pip install -r requirements.txt
hf auth login            # a read token; the dataset is gated and needs approved access
python3 mva_candidatos_genes.py
python3 mva_clinvar_cruce.py
python3 mva_barrido_genomico.py
python3 mva_cruce_hpo_genes.py --only-genes <strict list> ...   # see script help
python3 mva_armar_csv.py

Outputs are written to ~/Desktop/mva_hackathon_salida/ (change CARPETA_SALIDA in mva_candidatos_genes.py). Files whose name starts with LOCAL_NO_COMPARTIR_ contain patient-level variants and must never be shared. Files ending in _COMPARTIBLE contain only gene names and counts.

External resources (record the download date for each run)
Resource	Version	Downloaded on
Ensembl REST (gene/transcript coordinates)	GRCh38	YYYY-MM-DD
ClinVar vcf_GRCh38/clinvar.vcf.gz	weekly release	YYYY-MM-DD
GENCODE gencode.v50.basic.annotation.gtf.gz	v50	YYYY-MM-DD
HPO phenotype_to_genes.txt	—	YYYY-MM-DD
gnomAD browser (manual annotation)	v4	YYYY-MM-DD
Limitations
Single sample, no parents: the phase of the two candidate variants (cis vs trans) cannot be confirmed.
No population frequency is available locally; dbSNP membership is only a weak proxy for rarity.
Synonymous and missense SNVs are not distinguished by the scripts (manual annotation does it for the shortlist).
Structural variants, deep intronic variants and X-chromosome hemizygosity are not modelled.
Phenotype ranking is consistent with, but not independent from, the gene-panel hypothesis (the phenotype is MVA-like).
Data policy
The dataset is gated and may not be redistributed. Nothing derived from it at genotype scale is committed here (see .gitignore).
All patient data and derived files are deleted within 30 days of the hackathon close, and deletion is confirmed to the organisers.
Use of AI assistants

Anthropic Claude (claude.ai, Pro plan / Free plan for Cata, Consumer Terms, model-training setting: OFF). Used for writing code, searching public literature and interpreting public annotations (gnomAD, ClinVar) of a small number of named candidate variants. No dataset files, variant tables, or clinical descriptions were shared.

Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

License

Released under CC BY 4.0, as required for hackathon submissions.
