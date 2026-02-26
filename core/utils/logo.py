import os
import numpy as np

COLORS = {
    'red'   : ('\033[101m', '\033[91m'),
    'green' : ('\033[102m', '\033[92m'),
    'yellow': ('\033[103m', '\033[93m'),
    'blue'  : ('\033[104m', '\033[94m'),
    'purple': ('\033[105m', '\033[95m'),
}

def get_random_colors():
    colors = list(COLORS.keys())
    if hasattr(get_random_colors, '_previous_color'):
        colors.remove(get_random_colors._previous_color)
    if len(colors) < 3:
        colors = list(COLORS.keys())
    selected = np.random.choice(colors, 3, replace=False)
    get_random_colors._previous_color = selected[-1]
    return [COLORS[color] for color in selected]

def draw_box(lines, hcolor, log_func=None, nocolor=False):
    reset = '' if nocolor else '\033[0m'
    width  = 78
    top    = f'{hcolor}╭{"─" * width}╮{reset}'
    bottom = f'{hcolor}╰{"─" * width}╯{reset}'
    print(top)
    if log_func:
        log_func(f"+{'-' * width}+")
    for line in lines:
        if len(line) > width - 1:
            line_to_print = '...' + line[-((width - 1) - 3):]
        else:
            line_to_print = line
        output = f'{hcolor}│ {reset}{line_to_print:<{width-1}}{hcolor}│{reset}'
        print(output)
        if log_func:
            import re
            clean_text = re.sub(r'\033\[[0-9;]+m', '', output)
            log_func(clean_text)
    print(bottom)
    if log_func:
        log_func(f"+{'-' * width}+")

def print_logo(version, log_func=None, nocolor=False):
    version = f'v{version}'
    reset = '' if nocolor else '\033[0m'
    logo = [
        '              ╭──────────────────╮   ',
        '      ╭───────┴──────╮           │   ',
        '  ╭───┴───╮       ╭──┴──╮     ╭──┴──╮',
        '╭─┴─╮   ╭─┴─╮   ╭─┴─╮   │   ╭─┴─╮   │',
        'T   r   e   e   -   F   o   r   g   e',
        version.rjust(68)
    ]
    
    for line in logo[:-1]:
        centered = line.center(80)
        print(centered)
        if log_func:
            log_func(centered)
    
    outer, inner, title = (('', ''), ('', ''), ('', '')) if nocolor else get_random_colors()
    title_line = f'{outer[0]}   {reset}{inner[0]}  {reset}{title[0]}{logo[-1].center(70)}{reset}{inner[0]}  {reset}{outer[0]}   {reset}'
    print(title_line)
    if log_func:
        import re
        clean_text = re.sub(r'\033\[[0-9;]+m', '', title_line)
        log_func(clean_text)
    return title[1], title[0]

def print_help(hcolor, defaults, log_func=None, nocolor=False):
    help_lines = [
        'BASIC:',
        f'--input-dir                 -d    DIR     Directory of FASTA files   (./)',
        f'--iter                      -i    INT     Number of iterations       ({defaults["iter"]})',
        f'--threads                   -t    INT     Number of threads          ({defaults["threads"]})',
        f'--clutter                   -c    BOOL    Remove Intermediate Files  ({defaults["clutter"]})',
        f'--output-dir                -o    DIR     Output directory           (./)',
        f'--subprocess-logs           -sl   BOOL    Save subprocess logs       ({defaults["subprocess_logs"]})',
        f'--nocolor                         BOOL    Disable color in terminal',

        '',
        'BLAST:',
        f'--blast-evalue              -be   FLOAT   BLAST E-value threshold    ({defaults["blast_evalue"]})',
        f'--blast-max-targets         -bm   INT     BLAST max target seqs      ({defaults["blast_max_targets"]})',
        '',
        'MCL:',
        f'--mcl-hit-frac-cutoff       -mf   FLOAT   Hit fraction cutoff        ({defaults["mcl_hit_frac_cutoff"]})',
        f'--mcl-minimum-taxa          -mt   INT     Min taxa for MCL           ({defaults["mcl_minimum_taxa"]})',
        f'--mcl-inflation             -mi   FLOAT   MCL inflation param        ({defaults["mcl_inflation"]})',
        f'--mcl-perfect-identity      -mp   FLOAT   Perfect identity           ({defaults["mcl_perfect_identity"]})',
        f'--mcl-coverage-threshold    -mc   FLOAT   Coverage threshold         ({defaults["mcl_coverage_threshold"]})',
        f'--mcl-min-seq-length        -ml   INT     Min seq length             ({defaults["mcl_min_seq_length"]})',
        '',
        'MAFFT:',
        f'--mafft-maxiter             -mm   INT     MAFFT max iterations       ({defaults["mafft_maxiter"]})',
        f'--mafft-pxclsq-threshold    -mx   FLOAT   pxclsq prob threshold      ({defaults["mafft_pxclsq_threshold"]})',
        f'--mafft-thread-divisor      -md   INT     Thread division factor     ({defaults["mafft_thread_divisor"]})',
        '',
        'TREE:',
        f'--tree-relative-cutoff      -tr   FLOAT   Rel cutoff for trim tips   ({defaults["tree_relative_cutoff"]})',
        f'--tree-absolute-cutoff      -ta   FLOAT   Abs cutoff for trim tips   ({defaults["tree_absolute_cutoff"]})',
        f'--tree-branch-cutoff        -tb   FLOAT   Branch cutoff              ({defaults["tree_branch_cutoff"]})',
        f'--tree-start-from-prev      -tg   BOOL    Use prev trees as start    ({defaults.get("tree_start_from_prev", False)})',
        f'--tree-mask-paralogs        -tm   STR     Mask paraphyletic tips     ({defaults["tree_mask_paralogs"]})',
        f'--tree-outlier-ratio        -to   FLOAT   Outlier ratio threshold    ({defaults["tree_outlier_ratio"]})',
        f'--tree-max-trim-iterations  -ti   INT     Max trim iterations        ({defaults["tree_max_trim_iterations"]})',
        f'--tree-min-subtree-taxa     -ts   INT     Min taxa for subtrees      ({defaults["tree_min_subtree_taxa"]})',
        f'--tree-min-leaves           -tl   INT     Min leaves for tree        ({defaults["tree_min_leaves"]})',
        '',
        'PRUNE:',
        f'--prune-orthocutoff         -po   INT     Ortho min taxa cutoff      ({defaults["prune_orthocutoff"]})',
        f'--prune-relative-cutoff     -pr   FLOAT   Rel tip cutoff for prune   ({defaults["prune_relative_cutoff"]})',
        f'--prune-absolute-cutoff     -pa   FLOAT   Abs tip cutoff for prune   ({defaults["prune_absolute_cutoff"]})',
        f'--prune-outlier-ratio       -por  FLOAT   Outlier ratio for prune    ({defaults["prune_outlier_ratio"]})',
        f'--prune-max-trim-iter       -pmt  INT     Max trim iter for prune    ({defaults["prune_max_trim_iterations"]})',
        f'--prune-min-tree-leaves     -pml  INT     Min tree leaves for prune  ({defaults["prune_min_tree_leaves"]})',
        '',
        'PRANK:',
        # f'--prank-seqtype             -ps   STR     Sequence type              ({defaults["prank_seqtype"]})',
        f'--prank-pxclsq-threshold    -pp   FLOAT   pxclsq prob threshold      ({defaults["prank_pxclsq_threshold"]})',
        f'--prank-bootstrap           -pb   INT     IQ-TREE bootstraps         ({defaults["prank_bootstrap"]})',
        '',
        'SUPER:',
        f'--super-matrix              -sm   BOOL    Output supermatrix         ({defaults["super_matrix"]})',
        f'--super-bootstrap           -sb   INT     Supermatrix bootstraps     ({defaults["super_bootstrap"]})',
        f'--bes-support               -bs   FLOAT   BES support                ({defaults["bes_support"]})',
        '',
        'HCLUSTER:',
        f'--hcluster-enabled          -hc   BOOL    Hierarchical Clustering    ({defaults["hcluster_enabled"]})',
        f'--hcluster-id               -hi   FLOAT   Vsearch id                 ({defaults["hcluster_id"]})',
        f'--hcluster-iddef            -hid  INT     Vsearch iddef              ({defaults["hcluster_iddef"]})',
        f'--hcluster-tree             -hg   PATH    Custom guide tree path     ({defaults["hcluster_tree"]})',
        # f'--hcluster-tool             -ht   STR     Clustering tool            ({defaults["hcluster_tool"]})',
        f'--hcluster-tool             -ht   STR     vsearch or mmseqs2         ()',
        f'--hcluster-use-busco        -hub  BOOL    Use BUSCO-filtered files   ({defaults["hcluster_use_busco"]})',
        '',
        'BUSCO:',
        f'--busco-evalue              -bce  FLOAT   BUSCO BLAST E-value        ({defaults["busco_evalue"]})',
        f'--busco-max-targets         -bct  INT     BUSCO BLAST max targets    ({defaults["busco_max_targets"]})',
        f'--busco-coverage-threshold  -bcc  FLOAT   BUSCO coverage threshold   ({defaults["busco_coverage_threshold"]})',
        '',
        'CONFIG:',
        '--config                          PATH    Path to config file',
        '--config-create [NAME]            BOOL    Create config template',
        '--config-save [NAME]              BOOL    Save args to config',
    ]
    draw_box(help_lines, hcolor, log_func, nocolor=nocolor)

def print_args(args, hcolor, passed_args, log_func=None, nocolor=False):
    """Print command line args in formatted box."""
    dir_path = str(os.getcwd())[-30:] if not args.input_dir else ('..' + args.input_dir[-38:] if len(args.input_dir) > 40 else args.input_dir)
    
    arg_mappings = {
        'BASIC:': {
            'input_dir'                 : [f'Directory:', dir_path],
            'iter'                      : [f'Iterations:', args.iter],
            'threads'                   : [f'Threads:', args.threads],
            'clutter'                   : [f'Clutter:', args.clutter],
            'output_dir'                : [f'Output directory:', args.output_dir if args.output_dir else 'Same as input'],
            'subprocess_logs'           : [f'Subprocess logs:', args.subprocess_logs],
        },
        'BLAST:': {
            'blast_evalue'              : [f'E-value:', args.blast_evalue],
            'blast_max_targets'         : [f'Max targets:', args.blast_max_targets],
        },
        'MCL:': {
            'mcl_hit_frac_cutoff'       : [f'Hit fraction cutoff:', args.mcl_hit_frac_cutoff],
            'mcl_minimum_taxa'          : [f'Minimum taxa:', args.mcl_minimum_taxa],
            'mcl_inflation'             : [f'Inflation:', args.mcl_inflation],
            'mcl_perfect_identity'      : [f'Perfect identity:', args.mcl_perfect_identity],
            'mcl_coverage_threshold'    : [f'Coverage threshold:', args.mcl_coverage_threshold],
            'mcl_min_seq_length'        : [f'Min sequence length:', args.mcl_min_seq_length],
        },
        'MAFFT:': {
            'mafft_maxiter'             : [f'Max iterations:', args.mafft_maxiter],
            'mafft_pxclsq_threshold'    : [f'pxclsq threshold:', args.mafft_pxclsq_threshold],
            'mafft_thread_divisor'      : [f'Thread divisor:', args.mafft_thread_divisor],
        },
        'TREE:': {
            'tree_relative_cutoff'      : [f'Relative cutoff:', args.tree_relative_cutoff],
            'tree_absolute_cutoff'      : [f'Absolute cutoff:', args.tree_absolute_cutoff],
            'tree_branch_cutoff'        : [f'Branch cutoff:', args.tree_branch_cutoff],
            'tree_start_from_prev'      : [f'Start from previous trees:', args.tree_start_from_prev],
            'tree_mask_paralogs'        : [f'Mask paralogs:', 'true' if args.tree_mask_paralogs == 'y' else 'false'],
            'tree_outlier_ratio'        : [f'Outlier ratio:', args.tree_outlier_ratio],
            'tree_max_trim_iterations'  : [f'Max trim iterations:', args.tree_max_trim_iterations],
            'tree_min_subtree_taxa'     : [f'Min subtree taxa:', args.tree_min_subtree_taxa],
            'tree_min_leaves'           : [f'Min leaves:', args.tree_min_leaves],
        },
        'PRUNE:': {
            'prune_orthocutoff'         : [f'Orthocutoff:', args.prune_orthocutoff],
            'prune_relative_cutoff'     : [f'Relative cutoff:', args.prune_relative_cutoff],
            'prune_absolute_cutoff'     : [f'Absolute cutoff:', args.prune_absolute_cutoff],
            'prune_outlier_ratio'       : [f'Outlier ratio:', args.prune_outlier_ratio],
            'prune_max_trim_iterations' : [f'Max trim iterations:', args.prune_max_trim_iterations],
            'prune_min_tree_leaves'     : [f'Min tree leaves:', args.prune_min_tree_leaves],
        },
        'PRANK:': {
            # 'prank_seqtype'             : [f'Sequence type:', args.prank_seqtype],
            'prank_pxclsq_threshold'    : [f'pxclsq threshold:', args.prank_pxclsq_threshold],
            'prank_bootstrap'           : [f'Bootstrap replicates:', args.prank_bootstrap],
        },
        'SUPER:': {
            'super_bootstrap'           : [f'Bootstrap replicates:', args.super_bootstrap],
            'super_matrix'              : [f'Output supermatrix:', args.super_matrix],
        },
        'HCLUSTER:': {
            'hcluster_enabled'          : [f'Hierarchical Clustering:', args.hcluster_enabled],
            'hcluster_id'               : [f'Vsearch id:', args.hcluster_id],
            'hcluster_iddef'            : [f'Vsearch iddef:', args.hcluster_iddef],
            'hcluster_tree'             : [f'Custom guide tree:', args.hcluster_tree],
            'hcluster_tool'             : [f'Clustering tool:', args.hcluster_tool],
        },
        'BUSCO:': {
            'busco_evalue'              : [f'BUSCO E-value:', args.busco_evalue],
            'busco_max_targets'         : [f'BUSCO max targets:', args.busco_max_targets],
            'busco_coverage_threshold'  : [f'BUSCO coverage threshold:', args.busco_coverage_threshold],
        },
    }
    
    args_list = []
    for stage, arg_dict in arg_mappings.items():
        stage_args = [f'  {v[0]:<33} {v[1]}' for k, v in arg_dict.items() if k in passed_args]
        
        if stage_args:
            args_list.append(stage)
            args_list.extend(stage_args)
    
    draw_box(args_list, hcolor, log_func, nocolor=nocolor)

def print_hcluster_warning(hcolor, log_func=None, nocolor=False):
    warning          = 'WARNING: HIERARCHICAL CLUSTERING MODE'
    visible_width    = 77
    warning_centered = warning.center(visible_width)
    reset = '' if nocolor else '\033[0m'
    red_bg = '' if nocolor else '\033[41m'
    warning_with_bg  = f'{red_bg}{warning_centered}{reset}'
    red              = '' if nocolor else '\033[31m'
    width            = 78
    
    lines = [
        f'{red}╭{"─" * width}╮{reset}',
        f'{red}│ {reset}{warning_with_bg}{red}│{reset}',
        f'{red}│ {reset}{"Hierarchical clustering is still under construction and untested.":<77}{red}│{reset}',
        f'{red}│ {reset}{"Please use this feature with caution and verify results carefully.":<77}{red}│{reset}',
        f'{red}╰{"─" * width}╯{reset}'
    ]
    
    for line in lines:
        print(line)
        if log_func:
            import re
            clean_text = re.sub(r'\033\[[0-9;]+m', '', line)
            log_func(clean_text)

if __name__ == "__main__":
    from core.utils.config import Config
    hcolor, bcolor = print_logo("11110.1.0")
    defaults = Config().get_defaults_dict()
    print_help(hcolor, defaults)