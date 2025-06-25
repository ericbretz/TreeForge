import os
import pysam
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
                 bc):
        
        self.dir_base           = dir_base
        self.dir_treeforge      = dir_treeforge
        self.dir_tree           = dir_tree
        self.dir_mafft          = dir_mafft
        self.dir_trimmed        = dir_trimmed
        self.dir_prune          = dir_prune
        self.iter_total         = iter_total
        self.iter_current       = iter_current
        self.minimum_taxa       = minimum_taxa
        self.concatenated_fasta = concatenated_fasta
        self.threads            = threads
        self.log                = log
        self.hc                 = hc
        self.bc                 = bc

        self.tree_ending     = '.treefile'   # Tree ending
        self.relative_cutoff = 0.2           # Relative cutoff
        self.absolute_cutoff = 0.3           # Absolute cutoff
        self.branch_cutoff   = 0.02          # Branch cutoff
        self.para            = 'n'           # Run Paralogs
        self.contig_dct      = {}            # Contig dictionary

        self.printClass         = PrintOut(log, hc, bc)
        self.printout           = self.printClass.printout
        self.return_dict        = {}
        
    def run(self):
        self.trim_tips()
        self.mask_tips()
        self.cut_branches()
        self.write_tree()
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
                             self.hc)
        result = trim_tips.run()
        self.printout('metric', result)
        self.return_dict['trim_tips'] = result

    def mask_tips(self):
        """Mask tips that are not in the tree."""
        self.printout('metric', 'Masking tips')
        mask_tips = MaskTips(self.dir_tree,
                             self.dir_mafft,
                             self.para,
                             self.tree_ending,
                             self.hc)
        result = mask_tips.run()
        self.printout('metric', result)
        self.return_dict['mask_tips'] = result

    def cut_branches(self):
        """Remove branches that are too short or too long."""
        self.printout('metric', 'Cutting branches')
        cut_branches = CutBranches(self.dir_tree,
                                   self.dir_trimmed,
                                   '.mm',
                                   self.minimum_taxa,
                                   self.branch_cutoff,
                                   self.hc)
        result = cut_branches.run()
        self.printout('metric', result)
        self.return_dict['cut_branches'] = result

    def write_tree(self):
        """Write trees to file."""
        self.printout('metric', 'Writing trees')
        current_iter_dir = self.dir_trimmed        
        if self.iter_current == self.iter_total - 1:
            return
        else:
            next_iter = self.iter_current + 1
            out_dir = os.path.join(self.dir_treeforge, f'iter_{next_iter}', 'mafft')

        os.makedirs(out_dir, exist_ok=True)
        
        seq_dict = {}
        seq_dict_alt = {}
        with pysam.FastxFile(self.concatenated_fasta) as f:
            for entry in f:
                seq_dict[entry.name] = entry.sequence
                seq_dict_alt[entry.name.replace('@', '_')] = entry.sequence

        written_files = []
        subtree_files = [f for f in os.listdir(current_iter_dir) if f.endswith('.subtree')]
        total_sequences = 0
        files_written = 0
        
        for filename in subtree_files:
            with open(os.path.join(current_iter_dir, filename), "r") as infile:
                intree = parse(infile.readline())
            cluster_id = get_clusterID(filename)
            outname = os.path.join(out_dir, f"{cluster_id}.fa")
            labels = get_front_labels(intree)
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
                files_written += 1
            metric = {'written_files': written_files, 'total_sequences': total_sequences, 'files_written': files_written}
            self.printout('metric', metric)
            self.return_dict['write_tree'] = metric