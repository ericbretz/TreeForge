import os
import pysam
import shutil
from core.treeutils.trimtips    import TrimTips
from core.treeutils.masktips    import MaskTips
from core.treeutils.cutbranches import CutBranches
from core.treeutils.utils       import get_clusterID, get_front_labels
from core.treeutils.newick      import parse
from core.utils.printout        import PrintOut

class Tree:
    def __init__(self,
                 dir_base,
                 dir_treeforge,
                 dir_tree,
                 dir_mafft,
                 dir_trimmed,
                 dir_prune,
                 iter_total,
                 iter_current,
                 minimum_taxa,
                 concatenated_fasta,
                 threads,
                 log,
                 hc,
                 bc,
                 relative_cutoff,
                 absolute_cutoff,
                 branch_cutoff,
                 mask_paralogs,
                 outlier_ratio,
                 max_trim_iterations,
                 min_subtree_taxa,
                 min_tree_leaves,
                 clutter=False):
        
        self.dir_base            = dir_base
        self.dir_treeforge       = dir_treeforge
        self.dir_tree            = dir_tree
        self.dir_mafft           = dir_mafft
        self.dir_trimmed         = dir_trimmed
        self.dir_prune           = dir_prune
        self.iter_total          = iter_total
        self.iter_current        = iter_current
        self.minimum_taxa        = minimum_taxa
        self.concatenated_fasta  = concatenated_fasta
        self.threads             = threads
        self.log                 = log
        self.hc                  = hc
        self.bc                  = bc

        self.tree_ending         = '.treefile'          # Tree ending
        self.relative_cutoff     = relative_cutoff      # Relative cutoff
        self.absolute_cutoff     = absolute_cutoff      # Absolute cutoff
        self.branch_cutoff       = branch_cutoff        # Branch cutoff
        self.para                = mask_paralogs        # Run Paralogs
        self.contig_dct          = {}                   # Contig dictionary
        self.outlier_ratio       = outlier_ratio        # Ratio threshold for outlier detection
        self.max_trim_iterations = max_trim_iterations  # Maximum trimming iterations
        self.min_subtree_taxa    = min_subtree_taxa     # Minimum taxa for valid subtrees
        self.min_tree_leaves     = min_tree_leaves      # Minimum leaves for valid tree

        self.clutter             = clutter              # Clutter flag
        self.printClass          = PrintOut(log, hc, bc)
        self.printout            = self.printClass.printout
        self.return_dict         = {}
        
        self.modified_files      = set()
        self.final_files         = []
        
    def run(self):
        self.trim_tips()
        self.mask_tips()
        self.cut_branches()
        self.write_tree(clutter=self.clutter)
        return self.return_dict

    def trim_tips(self):
        """Remove outlier sequences from phylogenetic trees by trimming tips that fall outside our relative and absolute cutoff thresholds."""
        self.printout('metric', 'Trimming tips')
        trim_tips = TrimTips(self.dir_tree,
                             self.dir_mafft,
                             self.tree_ending,
                             self.relative_cutoff,
                             self.absolute_cutoff,
                             self.contig_dct,
                             self.hc,
                             self.outlier_ratio,
                             self.max_trim_iterations,
                             self.min_tree_leaves)
        result = trim_tips.run()
        self.printout('metric', result)
        self.return_dict['trim_tips'] = result
        
        if 'file_details' in result:
            for file_detail in result['file_details']:
                if file_detail.get('status') == 'processed' and file_detail.get('modified', False):
                    cluster_id = get_clusterID(file_detail['filename'])
                    self.modified_files.add(cluster_id)

    def mask_tips(self):
        """Mask tips that are not in the tree."""
        self.printout('metric', 'Masking tips')
        mask_tips = MaskTips(self.dir_tree,
                             self.dir_mafft,
                             self.para,
                             self.tree_ending,
                             self.hc,
                             self.min_tree_leaves)
        result = mask_tips.run()
        self.printout('metric', result)
        self.return_dict['mask_tips'] = result
        
        if 'file_details' in result:
            for file_detail in result['file_details']:
                if file_detail.get('status') == 'processed' and file_detail.get('modified', False):
                    cluster_id = get_clusterID(file_detail['filename'])
                    self.modified_files.add(cluster_id)

    def cut_branches(self):
        """Remove branches that are too short or too long."""
        self.printout('metric', 'Cutting branches')
        cut_branches = CutBranches(self.dir_tree,
                                   self.dir_trimmed,
                                   '.mm',
                                   self.minimum_taxa,
                                   self.branch_cutoff,
                                   self.hc,
                                   self.min_subtree_taxa)
        result = cut_branches.run()
        self.printout('metric', result)
        self.return_dict['cut_branches'] = result
        
        if 'file_details' in result:
            for file_detail in result['file_details']:
                if file_detail.get('status') == 'cut':
                    cluster_id = get_clusterID(file_detail['filename'])
                    self.modified_files.add(cluster_id)

    def write_tree(self, clutter=None):
        """Write trees to file."""
        if clutter is None:
            clutter = self.clutter
        self.printout('metric', 'Writing trees')
        current_iter_dir = self.dir_trimmed
        
        if self.iter_current == self.iter_total - 1 or len(self.modified_files) == 0:
            self.printout('info', f'No files modified in iteration {int(self.iter_current) + 1}')
            self.move_to_prune_stage(clutter=clutter)
            return
        else:
            next_iter = self.iter_current + 1
            out_dir   = os.path.join(self.dir_treeforge, f'iter_{next_iter}', 'mafft')

        os.makedirs(out_dir, exist_ok=True)
        
        seq_dict = {}
        seq_dict_alt = {}
        with pysam.FastxFile(self.concatenated_fasta) as f:
            for entry in f:
                seq_dict[entry.name] = entry.sequence
                seq_dict_alt[entry.name.replace('@', '_')] = entry.sequence

        written_files   = []
        subtree_files   = [f for f in os.listdir(current_iter_dir) if f.endswith('.subtree')]
        total_sequences = 0
        files_written   = 0
        files_skipped   = 0
        
        for filename in subtree_files:
            cluster_id = get_clusterID(filename)
            
            if cluster_id in self.modified_files:
                with open(os.path.join(current_iter_dir, filename), "r") as infile:
                    intree = parse(infile.readline())
                outname           = os.path.join(out_dir, f"{cluster_id}.fa")
                labels            = get_front_labels(intree)
                sequences_written = 0
                with open(outname, "w") as outfile:
                    for label in labels:
                        if label in seq_dict:
                            outfile.write(f">{label}\n{seq_dict[label]}\n")
                            sequences_written += 1
                        elif label in seq_dict_alt:
                            outfile.write(f">{label}\n{seq_dict_alt[label]}\n")
                            sequences_written += 1
                if sequences_written > 1:
                    written_files.append(outname)
                    total_sequences += sequences_written
                    files_written   += 1
            else:
                self.final_files.append(os.path.join(current_iter_dir, filename))
                files_skipped += 1
        
        if self.final_files:
            self.move_to_prune_stage(clutter=clutter)
        
        metric = {
            'written_files'       : written_files,
            'total_sequences'     : total_sequences,
            'files_written'       : files_written,
            'files_skipped'       : files_skipped,
            'modified_files_count': len(self.modified_files),
            'files_moved_to_prune': len(self.final_files)
        }
        
        self.printout('metric', metric)
        self.return_dict['write_tree'] = metric

    def move_to_prune_stage(self, clutter=None):
        """Move all remaining files to the prune stage and remove iteration directories if clutter is set."""
        if clutter is None:
            clutter = self.clutter
        current_iter_dir = self.dir_trimmed
        subtree_files    = [f for f in os.listdir(current_iter_dir) if f.endswith('.subtree')]
        
        self.final_files = []
        for filename in subtree_files:
            cluster_id = get_clusterID(filename)
            if cluster_id not in self.modified_files:
                self.final_files.append(os.path.join(current_iter_dir, filename))
        
        os.makedirs(self.dir_prune, exist_ok=True)
        for file_path in self.final_files:
            filename  = os.path.basename(file_path)
            dest_path = os.path.join(self.dir_prune, filename)
            shutil.copy2(file_path, dest_path)
        
        self.printout('metric', f'Moved {len(self.final_files)} unmodified files to prune stage')
        if 'write_tree' not in self.return_dict:
            self.return_dict['write_tree'] = {}
        self.return_dict['write_tree'].update({
            'files_moved_to_prune': len(self.final_files),
            'final_files': self.final_files
        })

        if clutter:
            try:
                if not any(os.scandir(current_iter_dir)):
                    os.rmdir(current_iter_dir)
                    self.printout('metric', f'Removed empty directory: {current_iter_dir}')
            except Exception:
                pass
            iteration_dir = os.path.dirname(current_iter_dir)
            try:
                if all(
                    not any(os.scandir(os.path.join(iteration_dir, subdir)))
                    for subdir in os.listdir(iteration_dir)
                    if os.path.isdir(os.path.join(iteration_dir, subdir))
                ):
                    os.rmdir(iteration_dir)
                    self.printout('metric', f'Removed empty iteration directory: {iteration_dir}')
            except Exception:
                pass