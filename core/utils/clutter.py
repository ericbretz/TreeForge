import os
from pathlib import Path

"""If the --clutter flag is passed, this is going to try and delete files on the fly. So that a massive tree's intermediate files don't take up too much space.
Convoluted? Yes. Works? Sometimes."""

def cleanup_files(master_dict, current_iter, stage, clutter, printout, hit_frac_cutoff=None, mcl_inflation=None, cleanup_type='intermediate'):
    if not clutter:
        return
    cleanup_rules = {
        'intermediate': {
            'blast': {
                'dir'     : master_dict['blast']['dir'],
                'patterns': ['concatenated.fasta.nhr', 'concatenated.fasta.nin', 'concatenated.fasta.nsq',
                             'concatenated.fasta.ndb', 'concatenated.fasta.nos', 'concatenated.fasta.ntf',
                             'concatenated.fasta.nto', 'concatenated.fasta.not', 'concatenated.fasta.nog']
            },
            'mcl': {
                'dir'     : master_dict['mcl']['dir'],
                'patterns': [
                    f'raw.hit-frac{hit_frac_cutoff}.minusLogEvalue' if hit_frac_cutoff is not None else 'raw.hit-frac.minusLogEvalue',
                    f'raw.hit-frac{hit_frac_cutoff}_I{mcl_inflation}_e5' if hit_frac_cutoff is not None and mcl_inflation is not None else 'raw.hit-frac_I_e5',
                    'raw.ident'],
                'check_iteration': True
            },
            'mafft': {
                'patterns': ['*.aln', '*.cln.tt']
            },
            'tree': {
                'patterns': ['*.subtree']
            },
            'prank': {
                'dir'     : master_dict['prank']['dir'],
                'patterns': ['*.fa', '*.fas', '*.iqtree', '*.log', '*.ckp.gz', '*.splits.nex', '*.bionj', '*.mldist']
            }
        },
        'preserved': {
            'mcl': {
                'dir'     : master_dict['blast']['dir'],
                'patterns': ['concatenated.fasta', 'raw.blast']
            },
            'prune': {
                'patterns': ['*.cln']
            },
            'prank': {
                'patterns': ['*.subtree']
            },
            'astral': {
                'dir'     : master_dict['prank']['dir'],
                'patterns': ['*.treefile', '*-cln']
            }
        }
    }
    rules = cleanup_rules.get(cleanup_type, {})
    if stage not in rules:
        return
    rule     = rules[stage]
    dir_path = get_cleanup_directory(master_dict, current_iter, stage, rule)
    if not dir_path or not dir_path.exists():
        return
    if (cleanup_type == 'intermediate' and stage == 'mcl' and rule.get('check_iteration', False)):
        if not should_cleanup_mcl_files(master_dict, current_iter):
            return
    if not has_files_to_cleanup(dir_path, rule['patterns']):
        return
    remove_files_by_patterns(dir_path, rule['patterns'], cleanup_type, printout)

def get_cleanup_directory(master_dict, current_iter, stage, rule):
    if stage in ['mafft', 'tree']:
        if f'iter_{current_iter}' not in master_dict:
            return None
        if stage == 'mafft':
            return master_dict[f'iter_{current_iter}']['mafft']['dir']
        else:
            return master_dict[f'iter_{current_iter}']['tree']['trimmed']
    elif stage == 'prune':
        if f'iter_{current_iter}' not in master_dict:
            return None
        return master_dict[f'iter_{current_iter}']['mafft']['dir']
    elif stage == 'prank':
        if f'iter_{current_iter}' not in master_dict:
            return None
        return master_dict[f'iter_{current_iter}']['tree']['trimmed']
    else:
        return rule.get('dir')

def should_cleanup_mcl_files(master_dict, current_iter):
    if f'iter_{current_iter}' not in master_dict:
        return False
    iter_dir = master_dict[f'iter_{current_iter}']['mafft']['dir']
    if not iter_dir.exists():
        return False
    fa_files = list(iter_dir.glob('*.fa'))
    return len(fa_files) > 0

def has_files_to_cleanup(dir_path, patterns):
    for pattern in patterns:
        if list(dir_path.glob(pattern)):
            return True
    return False

def remove_files_by_patterns(dir_path, patterns, cleanup_type, printout):
    # file_type = "intermediate" if cleanup_type == "intermediate" else "preserved"
    for pattern in patterns:
        for file_path in dir_path.glob(pattern):
            if file_path.exists():
                try:
                    file_path.unlink()
                    # printout('metric', f'Removed {file_type} file: {file_path.name}')
                except OSError as e:
                    pass
                    # printout('error', f'Failed to remove {file_path.name}: {e}')

def remove_empty_dirs(root_dir, protected_dirs, printout):
    if protected_dirs is None:
        protected_dirs = set()
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        dirpath = Path(dirpath)
        if dirpath.name in protected_dirs:
            continue
        if not any(dirpath.iterdir()):
            try:
                dirpath.rmdir()
                # if printout:
                    # printout('metric', f'Removed empty directory: {dirpath}')
            except OSError as e:
                # if printout:
                    pass
                    # printout('error', f'Failed to remove directory {dirpath}: {e}')

def final_cleanup(root_dir, logs_dir, printout, protected_files, protected_dirs, protected_exts):
    if protected_files is None:
        protected_files = {'SpeciesTree.tre', 'FinalTree.tre', 'summary.csv', 'SuperMatrix.tre'}
    if protected_dirs is None:
        protected_dirs  = {'logs', 'gene_trees'}
    if protected_exts is None:
        protected_exts  = {'.csv', '.yaml'}
    if not root_dir.exists():
        return
    if printout:
        # printout('metric', 'Performing targeted final cleanup')
        pass
    for item in root_dir.rglob('*'):
        if item.is_file():
            if item.name in protected_files:
                continue
            if any(prot in item.parts for prot in protected_dirs):
                if item.suffix == '.log':
                    try:
                        item.relative_to(logs_dir)
                        continue
                    except ValueError:
                        pass
                else:
                    continue
            if item.suffix in protected_exts:
                continue
            try:
                item.unlink()
                if printout:
                    # printout('metric', f'Removed file: {item}')
                    pass
            except OSError as e:
                if printout:
                    # printout('error', f'Failed to remove {item}: {e}')
                    pass
    remove_empty_dirs(root_dir, protected_dirs=protected_dirs, printout=printout)
    for protected_dir in protected_dirs:
        protected_path = root_dir / protected_dir
        if not protected_path.exists():
            protected_path.mkdir(parents=True, exist_ok=True)
            if printout:
                # printout('metric', f'Created protected directory: {protected_dir}')
                pass
    if printout:
        # printout('metric', 'Final cleanup complete')
        pass