<p align="center">
  <img src="https://i.imgur.com/lMrgrJG.png" alt="Treeforge" width="1000">
</p>

**Treeforge** is a pipeline for building phylogenetic trees from genomic data. This tool automates the phylogenomic workflow developed by Yang Y. and S.A. Smith. Treeforge combines several bioinformatics tools to do sequence alignment, clustering, tree building, and refinement. The whole process is automated and runs in steps, so, all you need to enter is a directory of .fasta files containing your genomic or transcriptomic sequences.
To learn more about the original workflow behind this pipeline, see: https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/
<p align="right">EC Bretz</p>

[TreeForge Metrics Reference](https://github.com/ericbretz/TreeForge/wiki/TreeForge-Metrics-Reference)

> [!CAUTION]
> $\Huge\textcolor[RGB]{248, 82, 73}{\textsf{TreeForge is currently under development}}$

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Dependencies</h2>

- **MAFFT** - Multiple sequence alignment
- **MCL** - Markov Clustering for sequence clustering
- **BLAST+** - Sequence similarity search (`blastn`, `blastx`, `makeblastdb`)
- **DIAMOND** - Fast protein alignment for BUSCO extraction (`diamond`)
- **PRANK** - Phylogeny-aware alignment
- **IQ-TREE2** - Maximum likelihood tree inference (`iqtree2`)
- **ASTRAL** - Coalescent-based species tree estimation
- **PHYX** - Phylogenetic tools (`pxclsq`)
- **MMseqs2** - Fast sequence clustering for HCluster mode (`mmseqs`) *(default HCluster tool)*
- **vsearch** - Alternative clustering tool for HCluster mode (`vsearch`) *(optional)*
- **Python** (>=3.7) with:
  - numpy
  - pysam
  - biopython
  - ete3
  - pyyaml

#### Ubuntu/Debian:
```bash
sudo apt-get install mafft mcl ncbi-blast+ prank iqtree phyx vsearch
```

#### DIAMOND:
```bash
# Download from https://github.com/bbuchfink/diamond/releases
# or via conda:
conda install -c bioconda diamond
```

#### MMseqs2:
```bash
# Download from https://github.com/soedinglab/MMseqs2/releases
# or via conda:
conda install -c conda-forge -c bioconda mmseqs2
```

#### Python dependencies:
```bash
pip install -r requirements.txt
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Installation</h2>

1. Clone the repository:
```bash
git clone https://github.com/ericbretz/TreeForge.git
cd TreeForge
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Usage</h2>

### Basic Usage

```bash
python treeforge.py -d /path/to/fasta/files
```

### Command Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| **Basic Parameters** | | | |
| `--input-dir` | `-d` | Directory containing FASTA files | Current directory |
| `--output-dir` | `-o` | Output directory | Same as input |
| `--iter` | `-i` | Number of MAFFT/Tree iterations | 5 |
| `--threads` | `-t` | Number of threads for parallel processing | 2 |
| `--clutter` | `-c` | Remove intermediate files after completion | False |
| `--subprocess-logs` | `-sl` | Save subprocess stdout/stderr to log files | False |
| **BLAST Stage** | | | |
| `--blast-evalue` | `-be` | BLAST E-value threshold | 10.0 |
| `--blast-max-targets` | `-bm` | BLAST max target sequences | 1000 |
| **MCL Stage** | | | |
| `--mcl-hit-frac-cutoff` | `-mf` | Hit fraction cutoff for MCL clustering | 0.3 |
| `--mcl-minimum-taxa` | `-mt` | Minimum taxa for MCL clustering | 4 |
| `--mcl-inflation` | `-mi` | MCL inflation parameter | 1.4 |
| `--mcl-perfect-identity` | `-mp` | Perfect identity threshold | 100.0 |
| `--mcl-coverage-threshold` | `-mc` | Coverage threshold for identical sequences | 0.9 |
| `--mcl-min-seq-length` | `-ml` | Minimum sequence length | 300 |
| **MAFFT Stage** | | | |
| `--mafft-maxiter` | `-mm` | MAFFT max iterations | 1000 |
| `--mafft-pxclsq-threshold` | `-mx` | pxclsq probability threshold (MAFFT) | 0.1 |
| `--mafft-thread-divisor` | `-md` | Thread division factor | 4 |
| **Tree Stage** | | | |
| `--tree-start-from-prev` | `-tg` | Use previous iteration's trees as IQ-TREE starting trees | False |
| `--tree-relative-cutoff` | `-tr` | Relative cutoff for trimming tips | 0.2 |
| `--tree-absolute-cutoff` | `-ta` | Absolute cutoff for trimming tips | 0.3 |
| `--tree-branch-cutoff` | `-tb` | Branch cutoff for cutting branches | 0.02 |
| `--tree-mask-paralogs` | `-tm` | Mask paraphyletic tips (y/n) | n |
| `--tree-outlier-ratio` | `-to` | Ratio threshold for outlier detection | 20.0 |
| `--tree-max-trim-iterations` | `-ti` | Maximum trimming iterations | 10 |
| `--tree-min-subtree-taxa` | `-ts` | Minimum taxa for valid subtrees | 4 |
| `--tree-min-leaves` | `-tl` | Minimum leaves for valid tree | 4 |
| **Prune Stage** | | | |
| `--prune-orthocutoff` | `-po` | Ortholog minimum taxa cutoff | 4 |
| `--prune-relative-cutoff` | `-pr` | Relative tip cutoff for pruning | 0.2 |
| `--prune-absolute-cutoff` | `-pa` | Absolute tip cutoff for pruning | 0.3 |
| `--prune-outlier-ratio` | `-por` | Outlier ratio for pruning | 20.0 |
| `--prune-max-trim-iterations` | `-pmt` | Max trim iterations for pruning | 10 |
| `--prune-min-tree-leaves` | `-pml` | Min tree leaves for pruning | 3 |
| **PRANK Stage** | | | |
| `--prank-seqtype` | `-ps` | Sequence type for PRANK (dna/aa) | aa |
| `--prank-pxclsq-threshold` | `-pp` | pxclsq probability threshold (PRANK) | 0.3 |
| `--prank-bootstrap` | `-pb` | IQ-TREE bootstrap replicates | 1000 |
| **Super Stage** | | | |
| `--super-bootstrap` | `-sb` | Supermatrix bootstrap replicates | 1000 |
| `--bes-support` | `-bs` | Molecular distance support threshold | 0.0 |
| `--super-matrix` | `-sm` | Output supermatrix tree | False |
| **HCluster Stage** | | | |
| `--hcluster-enabled` | `-hc` | Run Hierarchical Clustering | False |
| `--hcluster-id` | `-hi` | Clustering identity threshold | 0.9 |
| `--hcluster-iddef` | `-hid` | Clustering identity definition | 2 |
| `--hcluster-tool` | `-ht` | Clustering tool (`vsearch` or `mmseqs2`) | mmseqs2 |
| `--hcluster-tree` | `-hg` | Path to a custom guide tree (skips guide tree generation) | - |
| `--hcluster-use-busco` | `-hub` | Use BUSCO-filtered files for clustering | False |
| **BUSCO Stage** *(used with HCluster)* | | | |
| `--busco-evalue` | `-bce` | DIAMOND E-value threshold for BUSCO extraction | 1e-5 |
| `--busco-max-targets` | `-bct` | DIAMOND max target sequences for BUSCO | 1 |
| `--busco-coverage-threshold` | `-bcc` | Coverage threshold for BUSCO hits | 0.5 |
| **Configuration** | | | |
| `--config` | | Path to configuration file | - |
| `--config-create` | | Create default configuration file | - |
| `--config-save` | | Save current arguments to config.yaml | - |
| **Utility** | | | |
| `--version` | `-v` | Print version | - |

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Examples</h2>

#### Basic run with default parameters:
```bash
python treeforge.py -d /path/to/sequences
```

#### Specify output directory and supermatrix:
```bash
python treeforge.py -d /path/to/sequences -o /path/to/output --super-matrix
```

#### Typical parameters for DNA sequences:
```bash
python treeforge.py -d /path/to/sequences -i 5 -t 8 -mf 0.4 -mt 8 -po 9
```

#### Clean Clutter (remove intermediate files):
```bash
python treeforge.py -d /path/to/sequences -c
```

#### Using configuration files:
```bash
# Create a default configuration file
python treeforge.py --config-create

# Load configuration file
python treeforge.py --config config.yaml

# Save current arguments to config
python treeforge.py -d /path/to/sequences -t 8 --config-save
```

#### Run with Hierarchical Clustering (default mmseqs2 tool):
```bash
python treeforge.py -d /path/to/sequences --hcluster-enabled -hi 0.9
```

#### Run HCluster with vsearch and BUSCO filtering:
```bash
python treeforge.py -d /path/to/sequences --hcluster-enabled -ht vsearch -hub
```

#### Run HCluster with a custom guide tree (skip guide tree generation):
```bash
python treeforge.py -d /path/to/sequences --hcluster-enabled -hg /path/to/guide.tre
```

#### Save subprocess logs for debugging:
```bash
python treeforge.py -d /path/to/sequences -sl
```

#### Use previous iteration trees as IQ-TREE starting trees:
```bash
python treeforge.py -d /path/to/sequences -tg
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Overview</h2>

#### Standard Pipeline

1. **BLAST** - All-by-all sequence similarity search
2. **MCL** - Markov Clustering to identify orthologous groups
3. **Iterative MAFFT/Tree**:
   - **MAFFT** - Multiple sequence alignment
   - **Tree** - Phylogenetic tree construction, tip trimming, masking, and refinement
4. **Prune** - Filter for 1-to-1 orthologs and prune paralogs
5. **PRANK** - Phylogeny-aware alignment
6. **ASTRAL** - Coalescent-based species tree estimation
7. **Supermatrix (optional)** - If `--super-matrix` is set, a supermatrix tree is produced

#### HCluster Pipeline (`--hcluster-enabled`)

1. **BUSCO Extraction** - Extract conserved eukaryotic genes using DIAMOND against a built-in BUSCO database
2. **BLAST** - All-by-all BLAST on the BUSCO-filtered sequences
3. **GeneCluster** - Group sequences by gene identity from BLAST results
4. **Guide Tree** - MAFFT alignment + ASTRAL to produce a species guide tree
5. **HCluster** - Hierarchical clustering guided by the species tree using `mmseqs2` or `vsearch`

A custom guide tree can be provided with `--hcluster-tree` to skip guide tree generation entirely.

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Output</h2>

TreeForge creates a `TreeForge/` directory in your output directory (or input directory by default) with the following structure:

```
TreeForge/
├── 01_input/                    # Input processing
│   ├── fai/                     # FASTA index files
│   └── renamed_fasta/           # Renamed FASTA files
│
├── 02_analysis/                 # Core analysis steps
│   ├── blast/                   # BLAST results
│   │   ├── concatenated.fasta
│   │   └── raw.blast
│   ├── mcl/                     # MCL clustering results
│   ├── iterations/              # Iterative refinement (iter_0, iter_1, etc.)
│   │   ├── mafft/               # MAFFT alignments
│   │   └── tree/                # Tree processing
│   │       └── trimmed/         # Trimmed trees
│   ├── prune/                   # Pruning results
│   │   └── orthologs/           # 1-to-1 ortholog trees
│   ├── prank/                   # PRANK alignments
│   └── super/                   # Supermatrix analysis
│       └── concat.tre
│
├── 03_results/                  # Final results
│   ├── species_trees/           # Species trees (standard pipeline)
│   │   ├── SpeciesTree.coalescent.tre
│   │   ├── SpeciesTree.molecular.tre
│   │   └── SuperMatrix.tre      # (if --super-matrix is set)
│   ├── gene_trees/              # Individual gene trees
│   │   └── individual/          # Per-gene tree files
│   └── hcluster/                # HCluster results (if --hcluster-enabled)
│       ├── blast/               # HCluster BLAST results
│       ├── busco/               # BUSCO extraction results
│       ├── genecluster/         # Gene cluster FASTA files
│       └── clusters/            # Final clustering output
│
└── 04_metrics/                  # Metrics and logs
    ├── run_metrics.json         # Full pipeline metrics
    ├── pipeline_summary.csv     # Per-stage summary
    ├── iteration_flow.csv       # Per-iteration statistics
    ├── hcluster_metrics.json    # HCluster metrics (if enabled)
    └── logs/                    # Log files
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Output Files</h2>

### **Results** (in `03_results/`)
- **`SpeciesTree.coalescent.tre`** - The final coalescent-based phylogenetic tree
- **`SpeciesTree.molecular.tre`** - The molecular distance-based phylogenetic tree  
- **`gene_trees/`** - Directory containing individual gene trees
- **`SuperMatrix.tre`** - The supermatrix tree (if `--super-matrix` is used)

### **Metrics and Analysis** (in `04_metrics/`)
- **`run_metrics.json`** - Full pipeline metrics in JSON format
- **`pipeline_summary.csv`** - Summary metrics for each pipeline step
- **`iteration_flow.csv`** - Per-iteration statistics
- **`mafft_summary.csv`** - MAFFT alignment metrics
- **`tree_processing.csv`** - Tree trimming and refinement metrics
- **`prune_metrics.csv`** - Pruning stage metrics
- **`final_alignment.csv`** - PRANK alignment metrics
- **`hcluster_metrics.json`** - HCluster run metrics (if `--hcluster-enabled` is used)
- **`hcluster_metrics.csv`** - HCluster summary CSV (if `--hcluster-enabled` is used)
- **`logs/`** - Detailed log files for debugging (including subprocess logs if `-sl` is used)

### **Additional Outputs**
- **`hcluster/`** - Hierarchical clustering results (if `--hcluster-enabled` is used)
- **`supermatrix/`** - Supermatrix alignment and model files (if `--super-matrix` is used)
- **`orthologs/*.tre`** - 1-to-1 ortholog trees (in `02_analysis/prune/`)
