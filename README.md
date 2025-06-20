<p align="center">
  <img src="https://i.imgur.com/m3zp2j1.png" alt="Treeforge" width="1000">
</p>

**Treeforge** is a pipeline for building phylogenetic trees from genomic data. This tool automates the phylogenomic workflow developed by Yang Y. and S.A. Smith. Treeforge combines several bioinformatics tools to do sequence alignment, clustering, tree building, and refinement. The whole process is automated and runs in steps, so, all you need to enter is a directory of .fasta files containing your genomic or transcriptomic sequences.
To learn more about the original workflow behind this pipeline, see: https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/
<p align="right">EC Bretz</p>

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Dependencies</h2>

- **MAFFT** - Multiple sequence alignment
- **MCL** - Markov Clustering for sequence clustering
- **BLAST+** - Sequence similarity search
- **PRANK** - Phylogeny-aware alignment
- **IQ-TREE2** - Maximum likelihood tree inference
- **ASTRAL** - Coalescent-based species tree estimation
- **PHYX** - Phylogenetic tools (pxcat, pxclsq)

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Installation of Dependencies</h2>

#### Ubuntu/Debian:
```bash
sudo apt-get install mafft mcl ncbi-blast+ prank iqtree phyx
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Installation</h2>

1. Clone the repository:
```bash
git clone https://github.com/yourusername/TreeForge.git
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
| `--dir` | `-d` | Directory containing FASTA files | Current directory |
| `--iter` | `-i` | Number of MAFFT/Tree iterations | 3 |
| `--threads` | `-t` | Number of threads for parallel processing | 2 |
| `--hit-frac-cutoff` | `-f` | Hit fraction cutoff for MCL clustering | 0.3 |
| `--minimum-taxa` | `-m` | Minimum taxa for MCL clustering | 10 |
| `--orthocutoff` | `-o` | Ortholog minimum taxa cutoff | 20 |
| `--seqtype` | `-s` | Sequence type (`dna` or `aa`) | `aa` |
| `--clutter` | `-c` | Remove intermediate files after completion | False |

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Examples</h2>

#### Basic run with default parameters:
```bash
python treeforge.py -d /path/to/sequences
```

#### Custom parameters for DNA sequences:
```bash
python treeforge.py -d /path/to/sequences -i 5 -t 8 -f 0.4 -m 15 -o 25 -s aa
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
├── base/           # Base files and indices
├── blast/          # BLAST results
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
├── logs/           # Log files
├── treeforge.csv   # Metrics summary
└── FinalTree.tree  # Final phylogenetic tree
```

<h2><img src="https://i.imgur.com/kEuy7Sd.png" width="20" align="top">&ensp;Output Files</h2>

- **`FinalTree.tree`** - The final phylogenetic tree in Newick format
- **`treeforge.csv`** - Summary metrics for each pipeline step
