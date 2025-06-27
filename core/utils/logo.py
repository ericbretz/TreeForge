import os
import numpy as np

COLORS = {
    'red': ('\033[101m', '\033[91m'),
    'green': ('\033[102m', '\033[92m'),
    'yellow': ('\033[103m', '\033[93m'),
    'blue': ('\033[104m', '\033[94m'),
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
    width = 78
    top = f'{hcolor}╭{"─" * width}╮\033[0m'
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

def print_help(hcolor):
    """Print help info in formatted box."""
    help_lines = [
        'BASIC:',
        '--dir/-d                       DIR     Directory of FASTA files    (./)',
        '--iter/-i                      INT     Number of iterations        (3)',
        '--threads/-t                   INT     Number of threads           (2)',
        '--clutter/-c                   BOOL    Remove Intermediate Files   (False)',
        '',
        'BLAST:',
        '--blast-evalue/-be             FLOAT   BLAST E-value threshold     (10.0)',
        '--blast-max-targets/-bm        INT     BLAST max target seqs       (1000)',
        '',
        'MCL:',
        '--mcl-hit-frac-cutoff/-hf      FLOAT   Hit fraction cutoff         (0.3)',
        '--mcl-minimum-taxa/-mt         INT     Min taxa for MCL            (10)',
        '--mcl-inflation/-mi            FLOAT   MCL inflation param         (1.4)',
        '--mcl-perfect-identity/-mp     FLOAT   Perfect identity            (100.0)',
        '--mcl-coverage-threshold/-mc   FLOAT   Coverage threshold          (0.9)',
        '--mcl-min-seq-length/-ml       INT     Min seq length              (300)',
        '',
        'MAFFT:',
        '--mafft-maxiter/-mm            INT     MAFFT max iterations        (1000)',
        '--mafft-pxclsq-threshold/-mx   FLOAT   pxclsq prob threshold       (0.1)',
        '--mafft-thread-divisor/-md     INT     Thread division factor      (4)',
        '',
        'TREE:',
        '--tree-relative-cutoff/-tr     FLOAT   Rel cutoff for trim tips    (0.2)',
        '--tree-absolute-cutoff/-ta     FLOAT   Abs cutoff for trim tips    (0.3)',
        '--tree-branch-cutoff/-tb       FLOAT   Branch cutoff               (0.02)',
        '--tree-mask-paralogs/-tm       STR     Mask paraphyletic tips      (y/n, n)',
        '--tree-outlier-ratio/-to       FLOAT   Outlier ratio threshold     (20.0)',
        '--tree-max-trim-iterations/-ti INT     Max trim iterations         (10)',
        '--tree-min-subtree-taxa/-ts    INT     Min taxa for subtrees       (4)',
        '--tree-min-leaves/-tl          INT     Min leaves for tree         (4)',
        '',
        'PRUNE:',
        '--prune-orthocutoff/-po        INT     Ortho min taxa cutoff       (20)',
        '--prune-relative-cutoff/-pr    FLOAT   Rel tip cutoff for prune    (0.2)',
        '--prune-absolute-cutoff/-pa    FLOAT   Abs tip cutoff for prune    (0.3)',
        '',
        'PRANK:',
        '--prank-seqtype/-ps            STR     PRANK seq type            (dna/aa, aa)',
        '--prank-pxclsq-threshold/-pp   FLOAT   pxclsq prob threshold       (0.3)',
        '--prank-bootstrap/-pb          INT     IQ-TREE bootstraps          (1000)',
        # '',
        # 'SUPER STAGE:',
        # '--super-bootstrap/-sb    INT         Supermatrix bootstraps (1000)',
    ]
    draw_box(help_lines, hcolor)

def print_args(args, hcolor, passed_args):
    """Print command line args in formatted box."""
    skip_map = {
        'b': 'BLAST', 'm': 'MCL', 'a': 'MAFFT', 't': 'Tree',
        'p': 'Prune', 'r': 'Prank', 's': 'Astral'
    }
    
    dir_path = str(os.getcwd())[-30:] if not args.dir else ('..' + args.dir[-38:] if len(args.dir) > 40 else args.dir)
    
    # Define argument mappings for each stage
    arg_mappings = {
        'BASIC:': {
            'dir': [f'Directory:', dir_path],
            'iter': [f'Iterations:', args.iter],
            'threads': [f'Threads:', args.threads],
            'clutter': [f'Clutter:', args.clutter],
        },
        'BLAST:': {
            'blast_evalue': [f'E-value:', args.blast_evalue],
            'blast_max_targets': [f'Max targets:', args.blast_max_targets],
        },
        'MCL:': {
            'mcl_hit_frac_cutoff': [f'Hit fraction cutoff:', args.mcl_hit_frac_cutoff],
            'mcl_minimum_taxa': [f'Minimum taxa:', args.mcl_minimum_taxa],
            'mcl_inflation': [f'Inflation:', args.mcl_inflation],
            'mcl_perfect_identity': [f'Perfect identity:', args.mcl_perfect_identity],
            'mcl_coverage_threshold': [f'Coverage threshold:', args.mcl_coverage_threshold],
            'mcl_min_seq_length': [f'Min sequence length:', args.mcl_min_seq_length],
        },
        'MAFFT:': {
            'mafft_maxiter': [f'Max iterations:', args.mafft_maxiter],
            'mafft_pxclsq_threshold': [f'pxclsq threshold:', args.mafft_pxclsq_threshold],
            'mafft_thread_divisor': [f'Thread divisor:', args.mafft_thread_divisor],
        },
        'TREE:': {
            'tree_relative_cutoff': [f'Relative cutoff:', args.tree_relative_cutoff],
            'tree_absolute_cutoff': [f'Absolute cutoff:', args.tree_absolute_cutoff],
            'tree_branch_cutoff': [f'Branch cutoff:', args.tree_branch_cutoff],
            'tree_mask_paralogs': [f'Mask paralogs:', args.tree_mask_paralogs],
            'tree_outlier_ratio': [f'Outlier ratio:', args.tree_outlier_ratio],
            'tree_max_trim_iterations': [f'Max trim iterations:', args.tree_max_trim_iterations],
            'tree_min_subtree_taxa': [f'Min subtree taxa:', args.tree_min_subtree_taxa],
            'tree_min_leaves': [f'Min leaves:', args.tree_min_leaves],
        },
        'PRUNE:': {
            'prune_orthocutoff': [f'Orthocutoff:', args.prune_orthocutoff],
            'prune_relative_cutoff': [f'Relative cutoff:', args.prune_relative_cutoff],
            'prune_absolute_cutoff': [f'Absolute cutoff:', args.prune_absolute_cutoff],
        },
        'PRANK:': {
            'prank_seqtype': [f'Sequence type:', args.prank_seqtype],
            'prank_pxclsq_threshold': [f'pxclsq threshold:', args.prank_pxclsq_threshold],
            'prank_bootstrap': [f'Bootstrap replicates:', args.prank_bootstrap],
        },
    }
    
    # Build args_list based on passed_args
    args_list = []
    for stage, arg_dict in arg_mappings.items():
        stage_args = [f'  {v[0]:<33} {v[1]}' for k, v in arg_dict.items() if k in passed_args]
        
        if stage_args:  # Only add stage if it has arguments to display
            # args_list.append('')
            args_list.append(stage)
            args_list.extend(stage_args)

    if args.SKIP != '-':
        skipped_stages = [skip_map[x] for x in args.SKIP]
        if len(skipped_stages) > 4:
            mid_point = (len(skipped_stages) + 1) // 2
            first_row = ', '.join(skipped_stages[:mid_point])
            second_row = ', '.join(skipped_stages[mid_point:])
            # args_list.append('')
            args_list.append('SKIPPED STAGES:')
            args_list.append(f'                                 {first_row}')
            args_list.append(f'                                 {second_row}')
        else:
            skipped = ', '.join(skipped_stages)
            # args_list.append('')
            args_list.append('SKIPPED STAGES:')
            args_list.append(f'                                 {skipped}')
    
    draw_box(args_list, hcolor)

if __name__ == "__main__":
    hcolor, bcolor = print_logo("11110.1.0")
    print_help(hcolor)