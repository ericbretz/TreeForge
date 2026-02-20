from typing import Dict, Set, List
from pathlib import Path

FASTA_EXTENSIONS    = ['.fasta', '.fa', '.fas', '.fna']

CENTROIDS_SUFFIX    = '.centroids.fasta'
CLUSTERED_SUFFIX    = '.clustered.fasta'
UNCLUSTERED_SUFFIX  = '.unclustered.fasta'

CLUSTER_ID_FORMAT   = "{:06d}"

ITERATIVE_STEPS     = {'mafft', 'tree'}

STEP_METRICS_CONFIG: Dict[str, Dict[str, Set]] = {
    'blast' : {'exclude_keys': {'files'}},
    'mcl'   : {'exclude_keys': set()},
    'prune' : {'nested_key'  : 'prune'},
    'prank' : {'exclude_keys': set()},
    'astral': {'exclude_keys': set()}
}

FUNCTION_PARAMS_MAP = {
    'blast': [
        ('blast_evalue',              'E-value threshold'),
        ('blast_max_targets',         'Max target sequences'),
        ('threads',                   'Threads'),
    ],
    'mcl': [
        ('hit_frac_cutoff',           'Hit fraction cutoff'),
        ('minimum_taxa',              'Minimum taxa'),
        ('mcl_inflation',             'MCL inflation'),
        ('perfect_identity',          'Perfect identity threshold'),
        ('coverage_threshold',        'Coverage threshold'),
        ('min_seq_length',            'Min sequence length'),
        ('threads',                   'Threads'),
    ],
    'mafft': [
        ('mafft_maxiter',             'MAFFT max iterations'),
        ('pxclsq_threshold',          'pxclsq threshold'),
        ('thread_divisor',            'Thread divisor'),
        ('threads',                   'Threads'),
        ('current_iter',              'Current iteration'),
    ],
    'tree': [
        ('relative_cutoff',           'Relative tip cutoff'),
        ('absolute_cutoff',           'Absolute tip cutoff'),
        ('branch_cutoff',             'Branch cutoff'),
        ('mask_paralogs',             'Mask paralogs'),
        ('outlier_ratio',             'Outlier ratio'),
        ('max_trim_iterations',       'Max trim iterations'),
        ('min_subtree_taxa',          'Min subtree taxa'),
        ('min_tree_leaves',           'Min tree leaves'),
        ('current_iter',              'Current iteration'),
    ],
    'prune': [
        ('orthocutoff',               'Ortholog minimum taxa cutoff'),
        ('prune_relative_cutoff',     'Relative tip cutoff'),
        ('prune_absolute_cutoff',     'Absolute tip cutoff'),
        ('prune_outlier_ratio',       'Outlier ratio'),
        ('prune_max_trim_iterations', 'Max trim iterations'),
        ('prune_min_tree_leaves',     'Min tree leaves'),
    ],
    'prank': [
        ('seqtype',                   'Sequence type'),
        ('prank_pxclsq_threshold',    'pxclsq threshold'),
        ('bootstrap_replicates',      'Bootstrap replicates'),
        ('threads',                   'Threads'),
    ],
    'astral': [
        ('super_bootstrap',           'Supermatrix bootstrap'),
        ('bes_support',               'BES support'),
    ],
    'hcluster': [
        ('hcluster_id',               'Identity threshold'),
        ('hcluster_iddef',            'Identity definition'),
        ('hcluster_tool',             'Clustering tool (vsearch or mmseqs2)'),
        ('threads',                   'Threads'),
    ],
    'busco': [
        ('busco_evalue',              'E-value threshold'),
        ('busco_max_targets',         'Max target sequences'),
        ('busco_coverage_threshold',  'Coverage threshold'),
    ],
    'busco_extract': [
        ('busco_evalue',              'BUSCO E-value'),
        ('busco_max_targets',         'BUSCO max targets'),
        ('busco_coverage_threshold',  'BUSCO coverage threshold'),
        ('threads',                   'Threads'),
    ],
    'hcluster_guide_tree': [
        ('hcluster_id',               'HCluster identity'),
        ('hcluster_iddef',            'HCluster identity definition'),
        ('hcluster_tool',             'HCluster tool'),
        ('threads',                   'Threads'),
    ],
    'genecluster_for_hcluster': [
        ('minimum_taxa',              'Minimum taxa'),
        ('threads',                   'Threads'),
    ],
}

DEFAULT_CONFIG_FILENAME = 'config.yaml'

def get_fasta_files(directory: Path) -> List[Path]:
    files = []
    for ext in FASTA_EXTENSIONS:
        files.extend(directory.glob(f'*{ext}'))
    return sorted(files)
