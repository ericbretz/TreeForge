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
        '--dir/-d                 DIR         Directory of FASTA files',
        '--iter/-i                INT         Number of iterations',
        '--threads/-t             INT         Number of threads',
        '--hit-frac-cutoff/-f     FLOAT       Hit fraction cutoff',
        '--minimum-taxa/-m        INT         Minimum taxa for MCL clustering',
        '--orthocutoff/-o         FLOAT       Ortholog Minimum Taxa Cutoff',
        '--seqtype/-s             STR         Sequence type',
        '--clutter/-c             STR         Remove Intermediate Files'
    ]
    draw_box(help_lines, hcolor)

def print_args(args, hcolor):
    """Print command line args in formatted box."""
    skip_map = {
        'b': 'BLAST', 'm': 'MCL', 'a': 'MAFFT', 't': 'Tree',
        'p': 'Prune', 'r': 'Prank', 's': 'Astral'
    }
    
    dir_path = str(os.getcwd())[-30:] if not args.dir else ('..' + args.dir[-38:] if len(args.dir) > 40 else args.dir)
    
    arg_lines = [
        f'Directory:                       {dir_path}',
        f'Iterations:                      {args.iter}',
        f'Threads:                         {args.threads}',
        f'Hit fraction cutoff:             {args.hit_frac_cutoff}',
        f'Minimum taxa:                    {args.minimum_taxa}',
        f'Ortholog Minimum Taxa Cutoff:    {args.orthocutoff}',
        f'Sequence type:                   {args.seqtype}'
    ]
    
    if args.SKIP != '-':
        skipped_stages = [skip_map[x] for x in args.SKIP]
        if len(skipped_stages) > 4:
            mid_point = (len(skipped_stages) + 1) // 2
            first_row = ', '.join(skipped_stages[:mid_point])
            second_row = ', '.join(skipped_stages[mid_point:])
            arg_lines.append(f'Skipped:                         {first_row}')
            arg_lines.append(f'                                 {second_row}')
        else:
            skipped = ', '.join(skipped_stages)
            arg_lines.append(f'Skipped:                         {skipped}')
    
    draw_box(arg_lines, hcolor)

if __name__ == "__main__":
    hcolor, bcolor = print_logo("11110.1.0")
    print_help(hcolor)