<div align="center">

<img src="./images/graphabstract.png" width="520" alt="xMetaVar graphical abstract">

<!-- TODO[image]: replace ./images/graphabstract.png with the final Fig. 1 overview from the revised manuscript -->

# xMetaVar

### Scalable harmonization and interpretation of multi-layer microbial genomic variation across metagenomic cohorts

[![Backend Image](https://img.shields.io/badge/backend%20image-ghcr.io%2Fldearlistm%2Fxmetavar%3A1.0.0-2496ED?logo=docker)](https://github.com/ldearlistm/xMetaVar/pkgs/container/xmetavar)
[![Frontend Image](https://img.shields.io/badge/frontend%20image-ghcr.io%2Fldearlistm%2Fxmetavar--frontend%3A1.0-2496ED?logo=docker)](https://github.com/ldearlistm/xMetaVar/pkgs/container/xmetavar-frontend)
[![Web Server](https://img.shields.io/badge/web%20server-biosino.org%2FiMAC-success)](https://www.biosino.org/iMAC/xmetavar)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Database](https://img.shields.io/badge/reference%20database-Figshare-orange)](https://doi.org/10.6084/m9.figshare.30846347)

</div>

---

xMetaVar is a **metagenome-specific, reference-aware framework** for profiling and harmonizing single-nucleotide variants (SNVs), predefined SNP genotypes, short insertions/deletions (InDels), structural-variation-related signals (dSV/vSV), inversions, and gene-level copy-number variations (CNVs) into **standardized, locus-traceable cohort matrices**, while preserving variant-specific representations.

xMetaVar adopts a **two-stage local–web design** that separates computationally intensive read-level variant calling from reusable, interactive matrix-level interpretation:

- **Local workflow (this repository):** a containerized [Snakemake](https://snakemake.github.io/) pipeline that turns raw or host-depleted metagenomic reads into harmonized sample-by-feature matrices.
- **Web interpretation layer:** an interactive platform for cohort-level analysis and locus-level genomic inspection — available as a public web server or deployable as a self-hosted frontend.

> **Raw sequencing reads never leave your machine.** Only the compact standardized matrices, feature annotations and sample metadata are uploaded to the web layer, so read-level calling is performed once while downstream analyses can be repeated freely.

---

## Table of Contents

- [Features](#features)
- [Which way should I use xMetaVar?](#which-way-should-i-use-xmetavar)
- [Part I — Local variant-calling workflow](#part-i--local-variant-calling-workflow)
  - [1. Installation](#1-installation)
  - [2. Prepare input data](#2-prepare-input-data)
  - [3. Run the workflow](#3-run-the-workflow)
  - [4. Output files](#4-output-files)
  - [5. Usage tips & troubleshooting](#5-usage-tips--troubleshooting)
- [Part II — Web-based interpretation](#part-ii--web-based-interpretation)
- [Part III — Self-host the web frontend](#part-iii--self-host-the-web-frontend)
- [Reference database](#reference-database)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

- **Multi-layer variant profiling in one run** — integrates five complementary callers behind a unified reference framework:

  | Variant layer | Caller | Representation |
  | :-- | :-- | :-- |
  | SNVs (de novo) | MIDAS v3 | Minor-allele-frequency matrix |
  | Predefined SNP genotypes | GT-Pro | Reference/alternative allele counts |
  | Short InDels | QuickVariants (BWA-MEM) | Binary presence/absence matrix |
  | dSV / vSV signals | SGVFinder2 (ICRA) | Deletion states / variable-segment scores |
  | Inversions | PhaseFinder | Orientation states |
  | Gene-level CNVs | MIDAS gene module | Normalized gene copy number |

- **Harmonized, type-aware outputs** — native matrices preserve module-specific evidence; standardized matrices use variant-class-specific encoding for integrated cohort analysis.
- **Locus traceability** — every feature keeps its variant class, species, reference coordinate/interval, and gene/product annotation, so a cohort-level hit can always be traced back to its genomic context.
- **Reliability-aware** — simulation benchmarking defines variant-class-specific reliability boundaries by sequencing depth and community complexity.
- **Modular execution** — run the full pipeline or any combination of variant layers; regenerate deliverables from existing results without re-calling.
- **Reproducible by default** — one Docker image, fixed Conda environments, retained logs and QC summaries.

<!-- TODO[figure]: insert a compact variant-layer / architecture schematic here (Fig. 1B–D from the manuscript). Keep images under ./images/ and reference them with relative paths. -->

---

## Which way should I use xMetaVar?

| You want to… | Use this | Setup needed |
| :-- | :-- | :-- |
| Try analysis & visualization with **no installation** | [Public web server](https://www.biosino.org/iMAC/xmetavar) (use bundled demo data) | None |
| Call variants from **your own FASTQ reads** on a server/HPC | [Part I: local Docker workflow](#part-i--local-variant-calling-workflow) | Backend image + database |
| Explore the matrices you produced locally | Upload them to the [public web server](https://www.biosino.org/iMAC/xmetavar) | None |
| Host your **own private instance** of the web interface | [Part III: self-hosted frontend](#part-iii--self-host-the-web-frontend) | Frontend image + `public/` + `database/` |

---

# Part I — Local variant-calling workflow

## 1. Installation

### 1.1 Hardware requirements

The workflow runs on Linux (or any host with a Docker engine). Recommended resources:

| Resource | Minimum | Recommended |
| :-- | :-- | :-- |
| CPU cores | 8 | 16–32 (embarrassingly parallel across samples/rules) |
| RAM | 32 GB | 64 GB+ for large cohorts / many parallel species |
| Disk | ~60 GB | Input FASTQ + database (~7 GB unpacked) + intermediate/results; reserve headroom |

> Benchmark reference: with 16 cores, wall-clock time was ~3.0 h for 2 samples and ~12.5 h for 8 samples; peak memory ranged 22–52 GB. Adding cores shortens wall time (4 samples: 13.9 h @ 8 cores → 2.95 h @ 32 cores).

### 1.2 Install Docker

Install the Docker Engine (and the Compose plugin if you plan to self-host the frontend). No Conda, Snakemake or bioinformatics tool needs to be installed on the host — everything is inside the image.

### 1.3 Pull the backend image

```bash
docker pull ghcr.io/ldearlistm/xmetavar:1.0.0
```

The image bundles Snakemake v8.25.3 and all module Conda environments. Verify it runs:

```bash
docker run --rm ghcr.io/ldearlistm/xmetavar:1.0.0 --help
```

### 1.4 Obtain the reference database

See [Reference database](#reference-database) below. Unpack it so that a single `database/` directory contains the `kneaddata/`, `midasv3/`, `gt-pro/`, `QuickVariant/`, `SGVFinder2/`, `PhaseFinder/` subfolders and `43_species_features.tsv`.

### 1.5 Quick test with bundled example data

A minimal paired-end example (sample sheets, `config.yaml` and demo FASTQs) is provided under [`test/`](./test):

```bash
# test/example_pe_samples.tsv, test/example_se_samples.tsv, test/config.yaml, test/rawdata/
```

<!-- TODO: add a verified copy-paste command that runs the test/ example end-to-end once paths are finalized (dry-run first with -n). -->

---

## 2. Prepare input data

### 2.1 Project layout

Organize a working directory on the host:

```text
my_project/
├── rawdata/          # input .fastq / .fastq.gz (PE or SE), gzip or uncompressed
├── database/         # unpacked xMetaVar reference database (see 1.4)
├── results/          # output directory (created automatically)
├── samples.tsv       # sample sheet (see 2.2)
└── config.yaml       # pipeline configuration (see 2.3)
```

### 2.2 The sample sheet (`samples.tsv`)

A **tab-separated** file. Paths inside the sheet use **in-container** paths (`/pipeline/rawdata/...`), because the sheet is read inside the container.

**Paired-end (PE)** — columns `sample`, `r1`, `r2`:

```tsv
sample	r1	r2
CCMD19168690ST-11-0	/pipeline/rawdata/CCMD19168690ST-11-0_1.fastq.gz	/pipeline/rawdata/CCMD19168690ST-11-0_2.fastq.gz
CCMD46727384ST-11-0	/pipeline/rawdata/CCMD46727384ST-11-0_1.fastq.gz	/pipeline/rawdata/CCMD46727384ST-11-0_2.fastq.gz
```

**Single-end (SE)** — columns `sample`, `se`:

```tsv
sample	se
SampleA	/pipeline/rawdata/SampleA.fastq.gz
```

> One row per sample. Match `sequencing_type` in `config.yaml` to the sheet layout (`PE` or `SE`).

### 2.3 The configuration file (`config.yaml`)

A ready-to-use template is at [`test/config.yaml`](./test/config.yaml). Key fields:

| Field | Meaning |
| :-- | :-- |
| `raw_data_dir` | In-container input dir (`/pipeline/rawdata/`) |
| `results_dir` | In-container output dir (`/pipeline/results`) |
| `kneaddata_db` / `contaminant_db_prefix` | Host-depletion reference (e.g. `hg_38`); host build GRCh38 or T2T-CHM13 |
| `trimmomatic_path` | Trimmomatic adapters/tools path |
| `midasv3_db` | MIDAS v3 local database |
| `GT_Pro_db` / `GT_dict_path` | GT-Pro catalog and SNP dictionary |
| `QuickVariant_db` | Merged reference FASTA for BWA/QuickVariants |
| `SGVFinder2_db` | SGVFinder2 reference |
| `PhaseFinder_db` | PhaseFinder invertible-region reference |
| `annotation_file` | `43_species_features.tsv` feature annotation table |
| `samples_tsv` | In-container path to the sample sheet |
| `sequencing_type` | `"PE"` or `"SE"` |
| `skip_qc` | `true` to bypass KneadData when reads are already host-depleted |

Read preprocessing defaults (KneadData v0.12.0): Trimmomatic sliding-window `SLIDINGWINDOW:4:20`, `LEADING:3`, `TRAILING:3`, discard reads < 50 bp; host reads removed by Bowtie 2 alignment; tandem-repeat filtering bypassed. Set `skip_qc: true` only when supplying pre-cleaned, host-depleted reads.

---

## 3. Run the workflow

### 3.1 Mount mapping

Host directories are mounted to fixed in-container paths:

| Host path | Container path | Mode | Purpose |
| :-- | :-- | :-- | :-- |
| `./rawdata` | `/pipeline/rawdata` | `ro` | Input reads |
| `./results` | `/pipeline/results` | `rw` | All outputs |
| `./database` | `/pipeline/database` | `rw` | Reference database |
| `./config.yaml` | `/pipeline/config.yaml` | `ro` | Configuration |
| `./samples.tsv` | `/pipeline/samples.tsv` | `ro` | Sample sheet |

### 3.2 Run everything (`all`)

```bash
docker run --rm \
  -v "$(pwd)/rawdata:/pipeline/rawdata:ro" \
  -v "$(pwd)/results:/pipeline/results:rw" \
  -v "$(pwd)/database:/pipeline/database:rw" \
  -v "$(pwd)/config.yaml:/pipeline/config.yaml:ro" \
  -v "$(pwd)/samples.tsv:/pipeline/samples.tsv:ro" \
  -e XMETAVAR_CORES=16 \
  -e XMETAVAR_CHOWN_TO="$(id -u):$(id -g)" \
  ghcr.io/ldearlistm/xmetavar:1.0.0 all
```

`all` runs every variant module and then prepares the full deliverable set.

### 3.3 Modular target keys

Replace `all` with one or more target keys to run only selected layers:

| Target key | Produces |
| :-- | :-- |
| `snp_gtpro` | GT-Pro predefined SNP table |
| `snp_midas` | MIDAS SNV outputs + done-flag |
| `indel` | QuickVariants InDel annotation matrix |
| `sv_sgvfinder` | SGVFinder2 dSV **and** vSV annotations |
| `sv_inversion` | PhaseFinder inversion annotation |
| `sv_midas` | MIDAS gene CNV output + done-flag |
| `deliverables` | Re-build web-ready files, summaries and HTML report **from existing results** (allows missing variant types) |
| `all` | Full workflow + full deliverables (default when no key is given) |

Example — run only InDel + SV layers on 8 cores:

```bash
docker run --rm [mounts...] -e XMETAVAR_CORES=8 \
  ghcr.io/ldearlistm/xmetavar:1.0.0 indel sv_sgvfinder sv_inversion sv_midas
```

> `all` cannot be combined with specific keys. Keys may be combined freely with each other.

### 3.4 Environment variables (the wrapper API)

The entrypoint is a wrapper around Snakemake. Control it through environment variables rather than Snakemake flags:

| Variable | Default | Meaning |
| :-- | :-- | :-- |
| `XMETAVAR_CORES` | `16` | Snakemake `--cores` (**do not pass `--cores/-c/--jobs/-j` yourself**) |
| `XMETAVAR_CONFIGFILE` | `/pipeline/config.yaml` | Config path |
| `XMETAVAR_RESULTS_DIR` | `/pipeline/results` | Results path |
| `XMETAVAR_CONDA_PREFIX` | `/pipeline/.snakemake/conda` | Conda env cache |
| `XMETAVAR_FIX_PERMISSIONS` | `1` | `chmod/chown` results on exit |
| `XMETAVAR_RESULT_CHMOD` | `a+rwX` | Permission mode applied on exit |
| `XMETAVAR_CHOWN_TO` | *(unset)* | Optional `uid:gid` owner for results |

Pass raw Snakemake arguments after `--`, e.g.:

```bash
docker run --rm [mounts...] -e XMETAVAR_CORES=16 \
  ghcr.io/ldearlistm/xmetavar:1.0.0 all -- --printshellcmds --rerun-incomplete --keep-going
```

> The wrapper rejects `--cores/-c/--jobs/-j`, `--configfile`, `--use-conda` and `--conda-prefix` because it manages them for you — use the environment variables instead.

### 3.5 HPC / Slurm

Wrap the `docker run` command in a batch script (template: [`test/example_pipeline.job`](./test/example_pipeline.job)). A common pattern is to (1) run the chosen variant modules, then (2) run `deliverables` to assemble web-ready outputs:

```bash
# step 1: variant modules
docker run --rm [mounts...] -e XMETAVAR_CORES=16 -e XMETAVAR_CHOWN_TO="$(id -u):$(id -g)" \
  ghcr.io/ldearlistm/xmetavar:1.0.0 snp_midas indel sv_sgvfinder sv_inversion sv_midas \
  -- --printshellcmds --rerun-incomplete --keep-going

# step 2: assemble deliverables from the results above
docker run --rm [same mounts...] -e XMETAVAR_CORES=16 \
  ghcr.io/ldearlistm/xmetavar:1.0.0 deliverables -- --rerun-incomplete --keep-going
```

---

## 4. Output files

Results are written under `results/02-variant-calling/` (plus `results/logs/`). Key harmonized outputs:

| Layer | File (relative to `results/`) |
| :-- | :-- |
| SNP – GT-Pro | `02-variant-calling/SNP/GT-Pro/across-sample/snp.tsv` |
| SNP – MIDAS | `02-variant-calling/SNP/MIDAS/across_sample/snps/snps_summary.tsv` (+ `logs/MIDASv3/snp_done.txt`) |
| InDel | `02-variant-calling/INDEL/across-sample/indel_anno.tsv` |
| dSV / vSV | `02-variant-calling/SV/dSV-vSV/across-sample/dsgv_anno.tsv`, `vsgv_anno.tsv` |
| Inversion | `02-variant-calling/SV/Inversion/across-sample/inversion_anno.tsv` |
| CNV | `02-variant-calling/SV/CNV/across_sample/merge.genes_copynum.tsv` (+ `logs/MIDASv3/cnv_done.txt`) |

Each layer ships both a **native/near-native matrix** (module-specific values) and a **standardized matrix** (variant-type-specific binary or scaled encoding), together with an output manifest and feature annotation linking every feature to its module, variant class, species, reference locus and gene/product.

<!-- TODO: document the exact deliverables/ folder layout and the web-upload file set once finalized (which matrices + annotation + metadata the web server expects). -->

To visualize these matrices (landscape summaries, group comparisons, UpSet cross-layer overlap, phenotype association, Boruta/SHAP prioritization, JBrowse locus inspection), upload them to the web layer — see Part II.

---

## 5. Usage tips & troubleshooting

- **Always dry-run first** for a new config: append `-- -n` to print the planned jobs without executing.
- **Resume after interruption:** pass `-- --rerun-incomplete --keep-going`; completed jobs are cached and skipped.
- **Already host-depleted reads?** Set `skip_qc: true` to skip KneadData and save time.
- **PE vs SE mismatch:** the sample-sheet columns and `sequencing_type` must agree; a mismatch fails at input validation, not mid-run.
- **Permission-denied on `results/`:** set `-e XMETAVAR_CHOWN_TO="$(id -u):$(id -g)"` so outputs are owned by the host user.
- **Depth matters by layer:** SNVs/inversions stay robust at low depth; predefined SNPs/InDels need moderate depth; dSV/vSV require higher coverage. Treat vSV as an exploratory, coverage-sensitive signal.
- **"Missing" vs "absent":** a non-evaluable feature (species absent/too shallow) is tracked separately from a true, evaluable non-event — do not equate the two when filtering.
- **Cores have no effect?** You likely passed `--cores` directly; the wrapper blocks it. Use `XMETAVAR_CORES`.

<!-- TODO: expand with real recurring issues from your users as they come up. -->

---

# Part II — Web-based interpretation

The fastest way to explore results is the **public, login-free web server**:

:point_right: **https://www.biosino.org/iMAC/xmetavar**

It offers two modules:

1. **Cohort-level interpretation** — variant landscape, sample/group comparisons, variant-type selection, prevalence filtering, cross-layer UpSet overlap, phenotype association (microSLAM GLMM), random-forest + Boruta prioritization and SHAP interpretation.
2. **Locus-level genomic inspection** — query features by id/species/gene/region and inspect them on a JBrowse 2 genome view with GFF gene annotations.

### Step-by-step website usage

For the full click-by-click guide (upload format, module walkthrough, demo dataset, visualization screenshots), see the dedicated tutorial:

:point_right: **https://www.biosino.org/iMAC/xmetavar/tutorial**

> You upload only the standardized matrices + feature annotations + a sample-metadata table — never raw reads. A lightweight demo dataset is built in so you can try every module without running the local pipeline.

---

# Part III — Self-host the web frontend

The web interface is distributed as a separate, lightweight image (it does **not** perform read-level variant calling). Two deployment modes are supported.

### Prerequisites

- Pull the frontend image:
  ```bash
  docker pull ghcr.io/ldearlistm/xmetavar-frontend:1.0
  ```
- Prepare two data folders on the host:
  - `public/` — static web assets, JBrowse reference data and the bundled example set (tracked in this repository).
  - `database/` — the reference/annotation data needed for locus display (download from [Figshare](https://doi.org/10.6084/m9.figshare.30846347); see [Reference database](#reference-database)).
- An empty `shared_data/` is created at runtime for user sessions/results.

### Option A — Personal / local use (single container)

No `cleaner` sidecar is needed for a single local user:

```bash
docker run -d \
  --name xmetavar-frontend \
  --restart always \
  -p 3000:3000 \
  -e NODE_ENV=production \
  -e PORT=3000 \
  -v "$(pwd)/public:/app/public:ro" \
  -v "$(pwd)/database:/app/database:ro" \
  -v "$(pwd)/shared_data:/app/shared_data:rw" \
  ghcr.io/ldearlistm/xmetavar-frontend:1.0
```

Open <http://localhost:3000>.

### Option B — Multi-user server (frontend + periodic cleaner)

For shared deployments, add the `cleaner` sidecar that periodically removes stale files from `shared_data` (default: delete files/empty dirs older than 3 days every 24 h). Use Docker Compose:

```yaml
services:
  frontend:
    image: ghcr.io/ldearlistm/xmetavar-frontend:1.0
    container_name: xmetavar-frontend
    restart: always
    ports:
      - "3000:3000"
    environment:
      NODE_ENV: production
      HOSTNAME: 0.0.0.0
      PORT: 3000
      MAX_CONCURRENT_ANALYSES: "6"
    volumes:
      - ./frontend/shared_data:/app/shared_data
      - ./frontend/public:/app/public
      # TODO[deploy]: confirm the exact in-container database mount path for the frontend
      - ./database:/app/database:ro
    cpus: 48
    mem_limit: 44g
    command: ["node", "server.js"]
    logging:
      driver: json-file
      options: { max-size: "50m", max-file: "5" }

  cleaner:
    # Build locally from the cleaner Dockerfile (tiny Alpine image, not published to the registry)
    image: xmetavar-cleaner:v0
    container_name: xmetavar-cleaner
    restart: always
    volumes:
      - ./frontend/shared_data:/data_to_clean:rw
    logging:
      driver: json-file
      options: { max-size: "20m", max-file: "3" }
```

```bash
docker compose up -d
```

> The `cleaner` image is intentionally tiny (Alpine + a shell loop) and is only relevant for shared servers; build it locally — it does not need to be published. Tune `cpus`/`mem_limit`/`MAX_CONCURRENT_ANALYSES` to your host.

<!-- TODO[deploy]: (1) confirm exact frontend in-container mount paths for public/database/shared_data against server.js; (2) add the cleaner Dockerfile/build command; (3) state which database subset the frontend actually needs (full vs annotation-only). -->

---

## Reference database

The default framework covers **43 prevalent human-gut bacterial species** (selected from curatedMetagenomicData at mean relative abundance > 0.5% and prevalence > 50%), with representative genomes and gene annotations from BV-BRC cross-checked against NCBI RefSeq. Module-specific indices (Bowtie 2, GT-Pro, MIDAS, BWA, SGVFinder2, PhaseFinder) are pre-built on this common framework.

- **Download (pre-compiled):** [Figshare — 10.6084/m9.figshare.30846347](https://doi.org/10.6084/m9.figshare.30846347)
- **Archive:** `xMetaVar_database.tar.xz`
- **Unpack:**
  ```bash
  tar -xvf xMetaVar_database.tar.xz
  ```
- After unpacking, point every `*_db` / `*_path` field in `config.yaml` at the matching subfolder under `/pipeline/database/...` (in-container paths).

### Custom reference panel

xMetaVar also supports a custom panel: provide representative genomes + GFF annotations and build module-specific indices following each tool's official procedure, then update `config.yaml`. Catalog-based modules (GT-Pro, PhaseFinder) keep their predefined locus definitions and act as complementary layers.

<!-- TODO[database]: verify the Figshare archive's top-level layout and checksums (md5/sha256), and record the exact unpacked size here. -->

---

## Citation

If you use xMetaVar in your research, please cite:

> Cao Y., Li J., Chen W., et al. *xMetaVar enables scalable harmonization and interpretation of multi-layer microbial genomic variation across metagenomic cohorts.* **[Journal], [Year].** DOI: [pending]

<!-- TODO[citation]: fill in journal / volume / pages / DOI upon publication. -->

Key methods build on: Snakemake; MIDAS; GT-Pro; QuickVariants; SGVFinder2/ICRA; PhaseFinder; KneadData/Trimmomatic/Bowtie 2; microSLAM; JBrowse 2. Full reference list is provided in the manuscript.

---

## License

Released under the [MIT License](./LICENSE) — Copyright (c) 2026 NajiaoLab.

---

## Contact

- **Bug reports & feature requests:** [GitHub Issues](https://github.com/NajiaoLab/xMetaVar/issues)
- **Web server:** <https://www.biosino.org/iMAC/xmetavar>
- **Correspondence:** najiao@fudan.edu.cn

<!-- TODO[contact]: add a contributors/acknowledgements block and any funding numbers if desired. -->
