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
    """Get three random colors, avoiding repetition."""
    colors = list(COLORS.keys())
    if hasattr(get_random_colors, '_previous_color'):
        colors.remove(get_random_colors._previous_color)
    if len(colors) < 3:
        colors = list(COLORS.keys())
    selected = np.random.choice(colors, 3, replace=False)
    get_random_colors._previous_color = selected[-1]
    return [COLORS[color] for color in selected]

def draw_box(lines, hcolor):
    """Draw box around lines with highlight color."""
    width  = 78
    top    = f'{hcolor}╭{"─" * width}╮\033[0m'
    bottom = f'{hcolor}╰{"─" * width}╯\033[0m'
    print(top)
    for line in lines:
        print(f'{hcolor}│ \033[0m{line:<{width-1}}{hcolor}│\033[0m')
    print(bottom)

def print_logo(version):
    """Print TreeForge logo with random colors."""
    version = f'v{version}'
    logo = [
        '              ╭──────────────────╮   ',
        '      ╭───────┴──────╮           │   ',
        '  ╭───┴───╮       ╭──┴──╮     ╭──┴──╮',
        '╭─┴─╮   ╭─┴─╮   ╭─┴─╮   │   ╭─┴─╮   │',
        'T   r   e   e   -   F   o   r   g   e',
        version.rjust(68)
    ]
    
    for line in logo[:-1]:
        print(line.center(80))
    
    outer, inner, title = get_random_colors()
    title_line = f'{outer[0]}   \033[0m{inner[0]}  \033[0m{title[0]}{logo[-1].center(70)}\033[0m{inner[0]}  \033[0m{outer[0]}   \033[0m'
    print(title_line)
    return title[1], title[0]

def print_help(hcolor, defaults):
    """Print help info in formatted box, using defaults from config."""
    help_lines = [
        'BASIC:',
        f'--input-dir                 -d    DIR     Directory of FASTA files   (./)',
        f'--iter                      -i    INT     Number of iterations       ({defaults["iter"]})',
        f'--threads                   -t    INT     Number of threads          ({defaults["threads"]})',
        f'--clutter                   -c    BOOL    Remove Intermediate Files  ({defaults["clutter"]})',
        f'--output-dir                -o    DIR     Output directory           ({defaults["output_dir"] or "./"})',
        '',
        'BLAST:',
        f'--blast-evalue              -be   FLOAT   BLAST E-value threshold    ({defaults["blast_evalue"]})',
        f'--blast-max-targets         -bm   INT     BLAST max target seqs      ({defaults["blast_max_targets"]})',
        '',
        'MCL:',
        f'--mcl-hit-frac-cutoff       -hf   FLOAT   Hit fraction cutoff        ({defaults["mcl_hit_frac_cutoff"]})',
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
        # f'--prank-seqtype/-ps            STR     PRANK seq type            (dna/aa, {defaults["prank_seqtype"]})',
        f'--prank-pxclsq-threshold    -pp   FLOAT   pxclsq prob threshold      ({defaults["prank_pxclsq_threshold"]})',
        f'--prank-bootstrap           -pb   INT     IQ-TREE bootstraps         ({defaults["prank_bootstrap"]})',
        '',
        'SUPER:',
        f'--super-bootstrap           -sb   INT     Supermatrix bootstraps     ({defaults["super_bootstrap"]})',
        '',
        'CONFIG:',
        '--config                          PATH    Path to config file',
        '--config-create [NAME]            BOOL    Create config template',
        '--config-save [NAME]              BOOL    Save  args to config',
        '',
        'OUTPUT:',
        f'--output-super-matrix       -om   BOOL    Output supermatrix         ({defaults["output_super_matrix"]})',
    ]
    draw_box(help_lines, hcolor)

def print_args(args, hcolor, passed_args):
    """Print command line args in formatted box."""
    skip_map = {
        'b': 'BLAST', 'm': 'MCL', 'a': 'MAFFT', 't': 'Tree',
        'p': 'Prune', 'r': 'Prank', 's': 'Astral'
    }
    
    dir_path = str(os.getcwd())[-30:] if not args.input_dir else ('..' + args.input_dir[-38:] if len(args.input_dir) > 40 else args.input_dir)
    
    arg_mappings = {
        'BASIC:': {
            'input_dir'                 : [f'Directory:', dir_path],
            'iter'                      : [f'Iterations:', args.iter],
            'threads'                   : [f'Threads:', args.threads],
            'clutter'                   : [f'Clutter:', args.clutter],
            'output_dir'                : [f'Output directory:', args.output_dir if args.output_dir else 'Same as input'],
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
            'tree_mask_paralogs'        : [f'Mask paralogs:', args.tree_mask_paralogs],
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
            'prank_seqtype'             : [f'Sequence type:', args.prank_seqtype],
            'prank_pxclsq_threshold'    : [f'pxclsq threshold:', args.prank_pxclsq_threshold],
            'prank_bootstrap'           : [f'Bootstrap replicates:', args.prank_bootstrap],
        },
        'SUPER:': {
            'super_bootstrap'           : [f'Bootstrap replicates:', args.super_bootstrap],
        },
    }
    
    args_list = []
    for stage, arg_dict in arg_mappings.items():
        stage_args = [f'  {v[0]:<33} {v[1]}' for k, v in arg_dict.items() if k in passed_args]
        
        if stage_args:
            args_list.append(stage)
            args_list.extend(stage_args)

    if args.SKIP != '-':
        skipped_stages = [skip_map[x] for x in args.SKIP]
        if len(skipped_stages) > 4:
            mid_point  = (len(skipped_stages) + 1)              // 2
            first_row  = ', '.join(skipped_stages[:mid_point])
            second_row = ', '.join(skipped_stages[mid_point:])
            args_list.append('SKIPPED STAGES:')
            args_list.append(f'                                 {first_row}')
            args_list.append(f'                                 {second_row}')
        else:
            skipped = ', '.join(skipped_stages)
            args_list.append('SKIPPED STAGES:')
            args_list.append(f'                                 {skipped}')
    
    draw_box(args_list, hcolor)

if __name__ == "__main__":
    from core.utils.config import Config
    hcolor, bcolor = print_logo("11110.1.0")
    defaults = Config().get_defaults_dict()
    print_help(hcolor, defaults)