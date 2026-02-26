from dataclasses import dataclass
from pathlib import Path


@dataclass
class BlastConfig:
    evalue              : float
    max_targets         : int


@dataclass
class MCLConfig:
    hit_frac_cutoff     : float
    minimum_taxa        : int
    inflation           : float
    perfect_identity    : float
    coverage_threshold  : float
    min_seq_length      : int


@dataclass
class MafftConfig:
    maxiter             : int
    pxclsq_threshold    : float
    thread_divisor      : int


@dataclass
class TreeConfig:
    relative_cutoff     : float
    absolute_cutoff     : float
    branch_cutoff       : float
    mask_paralogs       : str
    outlier_ratio       : float
    max_trim_iterations : int
    min_subtree_taxa    : int
    min_leaves          : int
    start_from_prev     : bool = False


@dataclass
class PruneConfig:
    orthocutoff         : int
    relative_cutoff     : float
    absolute_cutoff     : float
    outlier_ratio       : float
    max_trim_iterations : int
    min_tree_leaves     : int


@dataclass
class PrankConfig:
    seqtype             : str
    pxclsq_threshold    : float
    bootstrap           : int


@dataclass
class SuperConfig: 
      bootstrap         : int
      bes_support       : float
      output_matrix     : bool = False


@dataclass
class HClusterConfig: 
      enabled           : bool
      id                : float
      iddef             : int
      tool              : str
      tree              : str = ''
      use_busco         : bool = False

@dataclass
class BuscoConfig       : 
      evalue            : float
      max_targets       : int
      coverage_threshold: float


@dataclass
class PipelineConfig:
    iterations     : int
    threads        : int
    clutter        : bool
    seqtype        : str
    subprocess_logs: bool

    dir_base       : Path
    output_dir     : Path
    dir_treeforge  : Path
    files_fasta    : list

    blast          : BlastConfig
    mcl            : MCLConfig
    mafft          : MafftConfig
    tree           : TreeConfig
    prune          : PruneConfig
    prank          : PrankConfig
    super          : SuperConfig
    hcluster       : HClusterConfig
    busco          : BuscoConfig
    
    log: int
    hc : str
    bc : str
    
    hcluster_tree_override: str = ''
    
    current_iter: int = 0
    
    @classmethod
    def from_args(cls, args):
        blast_config = BlastConfig(
            evalue              = args.blast_evalue,
            max_targets         = args.blast_max_targets
        )
        
        mcl_config = MCLConfig(
            hit_frac_cutoff     = args.mcl_hit_frac_cutoff,
            minimum_taxa        = args.mcl_minimum_taxa,
            inflation           = args.mcl_inflation,
            perfect_identity    = args.mcl_perfect_identity,
            coverage_threshold  = args.mcl_coverage_threshold,
            min_seq_length      = args.mcl_min_seq_length
        )
        
        mafft_config = MafftConfig(
            maxiter             = args.mafft_maxiter,
            pxclsq_threshold    = args.mafft_pxclsq_threshold,
            thread_divisor      = args.mafft_thread_divisor
        )
        
        tree_config = TreeConfig(
            relative_cutoff     = args.tree_relative_cutoff,
            absolute_cutoff     = args.tree_absolute_cutoff,
            branch_cutoff       = args.tree_branch_cutoff,
            mask_paralogs       = args.tree_mask_paralogs,
            outlier_ratio       = args.tree_outlier_ratio,
            max_trim_iterations = args.tree_max_trim_iterations,
            min_subtree_taxa    = args.tree_min_subtree_taxa,
            min_leaves          = args.tree_min_leaves,
            start_from_prev     = getattr(args, 'tree_start_from_prev', False)
        )
        
        prune_config = PruneConfig(
            orthocutoff         = args.prune_orthocutoff,
            relative_cutoff     = args.prune_relative_cutoff,
            absolute_cutoff     = args.prune_absolute_cutoff,
            outlier_ratio       = args.prune_outlier_ratio,
            max_trim_iterations = args.prune_max_trim_iterations,
            min_tree_leaves     = args.prune_min_tree_leaves
        )
        
        prank_config = PrankConfig(
            seqtype             = args.prank_seqtype,
            pxclsq_threshold    = args.prank_pxclsq_threshold,
            bootstrap           = args.prank_bootstrap
        )
        
        super_config = SuperConfig(
            bootstrap           = args.super_bootstrap,
            bes_support         = args.bes_support,
            output_matrix       = args.super_matrix
        )
        
        hcluster_config = HClusterConfig(
            enabled             = args.hcluster_enabled,
            id                  = args.hcluster_id,
            iddef               = args.hcluster_iddef,
            tool                = args.hcluster_tool,
            tree                = getattr(args, 'hcluster_tree', ''),
            use_busco           = getattr(args, 'hcluster_use_busco', False)
        )
        
        busco_config = BuscoConfig(
            evalue              = args.busco_evalue,
            max_targets         = args.busco_max_targets,
            coverage_threshold  = args.busco_coverage_threshold
        )
        
        from core.utils.constants import FASTA_EXTENSIONS
        
        dir_base      = Path(args.input_dir)
        output_dir    = Path(args.output_dir) if args.output_dir else dir_base
        dir_treeforge = output_dir / 'TreeForge'
        files_fasta   = [f for f in Path(dir_base).iterdir() if f.suffix in FASTA_EXTENSIONS]
        
        return cls(
            iterations             = args.iter,
            threads                = args.threads,
            clutter                = args.clutter,
            seqtype                = args.prank_seqtype,
            subprocess_logs        = getattr(args, 'subprocess_logs', False),
            dir_base               = dir_base,
            output_dir             = output_dir,
            dir_treeforge          = dir_treeforge,
            files_fasta            = files_fasta,
            blast                  = blast_config,
            mcl                    = mcl_config,
            mafft                  = mafft_config,
            tree                   = tree_config,
            prune                  = prune_config,
            prank                  = prank_config,
            super                  = super_config,
            hcluster               = hcluster_config,
            busco                  = busco_config,
            log                    = args.log,
            hc                     = args.highlight_color,
            bc                     = args.background_color,
            hcluster_tree_override = getattr(args, 'hcluster_tree', '')
        )
