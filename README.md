# TreeForge

A comprehensive phylogenetic tree construction pipeline for comparative genomics.

## Overview

TreeForge is an automated pipeline for constructing phylogenetic trees from genomic or transcriptomic data. It integrates multiple bioinformatics tools (BLAST, MCL, MAFFT, IQ-TREE, PRANK, ASTRAL) into a streamlined workflow.

## Features

- **Automated ortholog detection** using BLAST and MCL clustering
- **Multiple sequence alignment** with MAFFT
- **Gene tree construction** with IQ-TREE
- **Tree refinement** through iterative pruning and realignment
- **Species tree inference** using both coalescent (ASTRAL) and concatenation methods
- **Hierarchical clustering** support for large datasets
- **BUSCO-based gene extraction** for high-quality phylogenomics
- **Configurable parameters** via command-line or YAML files

## Installation

### Dependencies

TreeForge requires Python 3.7+ and the following Python packages:

```bash
pip install -r requirements.txt
```

### External Tools

Install the following bioinformatics tools and ensure they're in your PATH:

- **BLAST+** (blastn, makeblastdb)
- **MCL** (mcl)
- **MAFFT** (mafft)
- **IQ-TREE** (iqtree2)
- **PRANK** (prank)
- **ASTRAL** (astral)
- **phyx tools** (pxclsq) - for sequence cleaning
- **vsearch** or **mmseqs2** - for hierarchical clustering

## Quick Start

### Basic Usage

```bash
python treeforge.py -d /path/to/fasta_files -i 3 -t 8
```

### Common Options

```bash
# Run with custom parameters
python treeforge.py \\
  -d input_sequences/ \\
  -o output_results/ \\
  -i 5 \\
  -t 16 \\
  -mt 10 \\
  -po 20

# Use a configuration file
python treeforge.py --config my_config.yaml

# Create a configuration template
python treeforge.py --config-create my_config.yaml
```

## Configuration

TreeForge supports configuration via:
1. Command-line arguments
2. YAML configuration files
3. Combination of both (command-line overrides config file)

### Creating a Config File

```bash
python treeforge.py --config-create my_config.yaml
```

Edit the generated file to customize parameters:

```yaml
blast_evalue: 10.0
mcl_minimum_taxa: 10
mafft_maxiter: 1000
prune_orthocutoff: 20
# ... more parameters
```

## Pipeline Stages

1. **BLAST**: All-vs-all sequence comparison
2. **MCL**: Clustering of homologous sequences
3. **MAFFT**: Multiple sequence alignment (iterative)
4. **Tree Building**: Gene tree construction with IQ-TREE
5. **Pruning**: Ortholog identification and paralog removal
6. **PRANK**: Codon-aware alignment refinement
7. **ASTRAL**: Species tree inference from gene trees
8. **Supermatrix**: Concatenation-based species tree

## Testing

### Run Unit Tests

```bash
# Test formatters
python tests/unit/test_formatters.py

# Test constants
python tests/unit/test_constants.py
```

### Run Integration Tests

```bash
python tests/integration/test_stages.py
```

### Generate Test Data

```bash
python tests/generate_test_data.py
```

### Validation Tests

```bash
python tests/test_refactoring_equivalence.py
```

## Project Structure

```
TreeForge/
├── core/
│   ├── config/          # Configuration dataclasses
│   ├── stages/          # Pipeline stage implementations
│   ├── treeutils/       # Tree manipulation utilities
│   └── utils/           # Shared utilities
├── tests/               # Test suite
│   ├── unit/           # Unit tests
│   ├── integration/    # Integration tests
│   └── generate_test_data.py
├── treeforge.py        # Main entry point
└── requirements.txt    # Python dependencies
```

## Advanced Features

### Hierarchical Clustering

For large datasets with many taxa:

```bash
python treeforge.py -d input/ -hc --hcluster-tool vsearch
```

### Output Supermatrix

Generate concatenated supermatrix alignment:

```bash
python treeforge.py -d input/ --super-matrix
```

## Output Files

TreeForge creates the following directory structure:

```
output_dir/TreeForge/
├── 01_input/          # Renamed FASTA files
├── 02_analysis/       # Intermediate analysis files
│   ├── blast/
│   ├── mcl/
│   ├── iterations/
│   ├── prune/
│   └── prank/
├── 03_results/        # Final results
│   ├── species_trees/
│   │   ├── SpeciesTree.coalescent.tre
│   │   ├── SpeciesTree.molecular.tre
│   │   └── SuperMatrix.tre
│   └── gene_trees/
└── 04_metrics/        # Statistics and logs
    ├── summary.csv
    └── logs/
```

## Contributing

TreeForge has been refactored for improved:
- **Code organization**: Base classes and configuration management
- **Efficiency**: Optimized path operations and data structures
- **Maintainability**: Centralized constants and utilities
- **Testing**: Comprehensive test suite

## Citation

If you use TreeForge in your research, please cite:
[Citation information to be added]

## License

[License information to be added]

## Contact

For questions or issues, please open an issue on GitHub or contact [contact info].

## Acknowledgments

TreeForge integrates and builds upon many excellent bioinformatics tools.
We thank the developers of BLAST, MCL, MAFFT, IQ-TREE, PRANK, and ASTRAL.
