<p align="center">
  <img src="https://i.imgur.com/lMrgrJG.png" alt="Treeforge" width="1000">
</p>

**Treeforge** is a pipeline for building phylogenetic trees from genomic data. This tool automates the phylogenomic workflow developed by Yang Y. and S.A. Smith. Treeforge combines several bioinformatics tools to do sequence alignment, clustering, tree building, and refinement. The whole process is automated and runs in steps, so, all you need to enter is a directory of .fasta files containing your genomic or transcriptomic sequences.
To learn more about the original workflow behind this pipeline, see: https://bitbucket.org/yanglab/phylogenomic_dataset_construction/
<p align="right">EC Bretz</p>

> [!CAUTION]
> $\Huge\textcolor[RGB]{248, 82, 73}{\textsf{TreeForge is currently under development}}$

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Dependencies</h2>

- **MAFFT** - Multiple sequence alignment
- **MCL** - Markov Clustering for sequence clustering
- **BLAST+** - Sequence similarity search
- **PRANK** - Phylogeny-aware alignment
- **IQ-TREE2** - Maximum likelihood tree inference
- **ASTRAL** - Coalescent-based species tree estimation
- **PHYX** - Phylogenetic tools (pxcat, pxclsq)

#### Ubuntu/Debian:
```bash
sudo apt-get install mafft mcl ncbi-blast+ prank iqtree phyx
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
| `--dir` | `-d` | Directory containing FASTA files | Current directory |
| `--iter` | `-i` | Number of MAFFT/Tree iterations | 3 |
| `--threads` | `-t` | Number of threads for parallel processing | 2 |
| `--clutter` | `-c` | Remove intermediate files after completion | False |
| **BLAST Stage** | | | |
| `--blast-evalue` | `-be` | BLAST E-value threshold | 10.0 |
| `--blast-max-targets` | `-bm` | BLAST max target sequences | 1000 |
| **MCL Stage** | | | |
| `--mcl-hit-frac-cutoff` | `-hf` | Hit fraction cutoff for MCL clustering | 0.3 |
| `--mcl-minimum-taxa` | `-mt` | Minimum taxa for MCL clustering | 10 |
| `--mcl-inflation` | `-mi` | MCL inflation parameter | 1.4 |
| `--mcl-perfect-identity` | `-mp` | Perfect identity threshold | 100.0 |
| `--mcl-coverage-threshold` | `-mc` | Coverage threshold for identical sequences | 0.9 |
| `--mcl-min-seq-length` | `-ml` | Minimum sequence length | 300 |
| **MAFFT Stage** | | | |
| `--mafft-maxiter` | `-mm` | MAFFT max iterations | 1000 |
| `--mafft-pxclsq-threshold` | `-mx` | pxclsq probability threshold (MAFFT) | 0.1 |
| `--mafft-thread-divisor` | `-md` | Thread division factor | 4 |
| **Tree Stage** | | | |
| `--tree-relative-cutoff` | `-tr` | Relative cutoff for trimming tips | 0.2 |
| `--tree-absolute-cutoff` | `-ta` | Absolute cutoff for trimming tips | 0.3 |
| `--tree-branch-cutoff` | `-tb` | Branch cutoff for cutting branches | 0.02 |
| `--tree-mask-paralogs` | `-tm` | Mask paraphyletic tips (y/n) | n |
| `--tree-outlier-ratio` | `-to` | Ratio threshold for outlier detection | 20.0 |
| `--tree-max-trim-iterations` | `-ti` | Maximum trimming iterations | 10 |
| `--tree-min-subtree-taxa` | `-ts` | Minimum taxa for valid subtrees | 4 |
| `--tree-min-leaves` | `-tl` | Minimum leaves for valid tree | 4 |
| **Prune Stage** | | | |
| `--prune-orthocutoff` | `-po` | Ortholog minimum taxa cutoff | 20 |
| `--prune-relative-cutoff` | `-pr` | Relative tip cutoff for pruning | 0.2 |
| `--prune-absolute-cutoff` | `-pa` | Absolute tip cutoff for pruning | 0.3 |
| **PRANK Stage** | | | |
| `--prank-seqtype` | `-ps` | Sequence type for PRANK (dna/aa) | aa |
| `--prank-pxclsq-threshold` | `-pp` | pxclsq probability threshold (PRANK) | 0.3 |
| `--prank-bootstrap` | `-pb` | IQ-TREE bootstrap replicates | 1000 |
| **Utility** | | | |
| `--version` | `-v` | Print version | - |

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Examples</h2>

#### Basic run with default parameters:
```bash
python treeforge.py -d /path/to/sequences
```

#### Recommended parameters to use (custom values):
```bash
python treeforge.py -d /path/to/sequences -i 5 -t 8 -hf 0.4 -mt 15 -po 25 -ps dna
```

#### Clean Clutter (remove intermediate files):
```bash
python treeforge.py -d /path/to/sequences -c
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Overview</h2>

1. **BLAST** - All-by-all sequence similarity search
2. **MCL** - Markov Clustering to identify orthologous groups
3. **Iterative MAFFT/Tree** (configurable iterations):
   - **MAFFT** - Multiple sequence alignment
   - **Tree** - Phylogenetic tree construction and refinement
4. **Prune** - Filter for 1-to-1 orthologs
5. **PRANK** - Phylogeny-aware alignment
6. **ASTRAL** - Coalescent-based species tree estimation

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Output</h2>

TreeForge creates a `TreeForge/` directory in your input directory with the following structure:

```
TreeForge/
├── blast/          # BLAST results
├── fai/            # Stores fai files
├── mcl/            # MCL clustering results
├── iter_0/         # First iteration results
│   ├── mafft/      # MAFFT alignments
│   └── tree/       # Tree files
├── iter_1/         # Second iteration results
│   ├── mafft/
│   └── tree/
├── prune/          # Pruned orthologs
├── prank/          # PRANK alignments
├── super/          # Supertree results
├── gene_trees/     # Gene tree files
├── logs/           # Log files
├── summary.csv     # Metrics summary
└── FinalTree.tree  # Final phylogenetic tree
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Output Files</h2>

- **`FinalTree.tre`**  - The final phylogenetic tree in Newick format
- **`cluster_x.tree`** - Collection of gene trees in the Newick format
- **`summary.csv`**  - Summary metrics for each pipeline step

<p align="center">
  <img src="https://i.imgur.com/v4OAUEE.png" alt="Treeforge Workflow" width="1000">
</p>
