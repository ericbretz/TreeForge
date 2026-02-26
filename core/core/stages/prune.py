import os
import shutil
from ete3 import Tree
from pathlib                    import Path
from core.stages.base_stage     import BaseStage
from core.treeutils.trimtips    import TrimTips
from core.treeutils.newick      import parse, tostring
from core.treeutils.utils       import get_front_labels, get_name, get_back_labels, remove_kink

class Prune(BaseStage):
    def __init__(self,
                 dir_base,
                 dir_prune,
                 dir_trimmed,
                 dir_ortho,
                 orthocutoff,
                 threads,
                 log,
                 hc,
                 bc,
                 prune_relative_cutoff,
                 prune_absolute_cutoff,
                 prune_outlier_ratio=20.0,
                 prune_max_trim_iterations=10,
                 prune_min_tree_leaves=3,
                 shared_printClass=None):
        
        super().__init__(log, hc, bc, threads, shared_printClass=shared_printClass)
        self.dir_base            = Path(dir_base)
        self.dir_prune           = Path(dir_prune)
        self.dir_trimmed         = Path(dir_trimmed)
        self.dir_ortho           = Path(dir_ortho)
        self.orthocutoff         = orthocutoff

        self.subtree_files       = [self.dir_prune / f for f in os.listdir(self.dir_prune) if f.endswith('.subtree')]
        self.relative_tip_cutoff = prune_relative_cutoff
        self.absolute_tip_cutoff = prune_absolute_cutoff
        self.minimum_taxa        = orthocutoff
        self.outlier_ratio       = prune_outlier_ratio
        self.max_trim_iterations = prune_max_trim_iterations
        self.min_tree_leaves     = prune_min_tree_leaves
        self.outpath             = self.dir_ortho
        self.ortho1to1_files     = []
        self.ortho_mi_files      = []
        self.outpath.mkdir(parents=True, exist_ok=True)

    def run(self):
        self.printout('metric', 'Finding 1-to-1 orthologs and pruning paralogs')
        
        metrics = {
            'total_trees'    : len(self.subtree_files),
            'ortho1to1_count': 0,
            'ortho_mi_count' : 0,
            'pruned_trees'   : 0,
            'skipped_trees'  : 0,
            'file_details'   : []
        }
        
        for st in self.subtree_files:
            file_metrics = {
                'filename'      : st.name,
                'stem'          : st.stem,
                'status'        : '',
                'ortho1to1'     : False,
                'ortho_mi_count': 0,
                'error'         : None
            }
            
            try:
                with open(st, 'r') as infile:
                    intree = parse(infile.readline())
                
                curroot = intree
                pp_trees = []
                
                if self.get_front_score(curroot) >= self.minimum_taxa:
                    self.ortho1to1_files.append(st)
                    metrics['ortho1to1_count'] += 1
                    file_metrics['status'] = '1to1_ortholog'
                    file_metrics['ortho1to1'] = True
                else:
                    # self.printout('metric', f'Pruning paralogs in {st.name}')
                    pruned_trees = self.prune_paralogs(curroot, st.stem)
                    
                    if pruned_trees:
                        metrics['pruned_trees']       += 1
                        file_metrics['ortho_mi_count'] = len(pruned_trees)
                        metrics['ortho_mi_count']     += len(pruned_trees)
                        file_metrics['status']         = 'pruned'
                        self.ortho_mi_files.extend(pruned_trees)
                    else:
                        metrics['skipped_trees'] += 1
                        file_metrics['status']    = 'skipped'
                        file_metrics['error']     = 'No valid orthologs after pruning'
                        # self.printout('metric', f'Skipping {st.name}: No valid orthologs found after pruning.')
                        # self.printout('error', f'Minimum taxa required: {self.minimum_taxa}')
                        # self.printout('error', f'Taxa Present: {len(curroot.leaves())}')
                
                metrics['file_details'].append(file_metrics)
                
            except Exception as e:
                metrics['skipped_trees'] += 1
                file_metrics['status']    = 'error'
                file_metrics['error']     = str(e)
                metrics['file_details'].append(file_metrics)
                self.printout('error', f'Error processing {st.name}: {str(e)}')

        self.write_ortho1to1_files()
        self.write_ortho_mi_files()
        
        self.return_dict['prune'] = {
            'ortho1to1'      : self.ortho1to1_files,
            'ortho1to1_dir'  : self.outpath,
            'ortho1to1_count': len(self.ortho1to1_files),
            'ortho_mi_files' : self.ortho_mi_files,
            'ortho_mi_count' : len(self.ortho_mi_files),
            'metrics'        : metrics
        }
        return self.return_dict
    
    def write_ortho1to1_files(self):
        """Write 1-to-1 ortholog files that don't need pruning."""
        self.printout('metric', 'Writing 1-to-1 ortholog trees')
        for st in self.ortho1to1_files:
            outfile = self.outpath / (Path(st).stem + '_1to1ortho.tre')
            shutil.copy(st, outfile)

    def write_ortho_mi_files(self):
        """Write multiple ortholog files from pruned trees."""
        self.printout('metric', 'Writing multiple ortholog trees')
        for file_path in self.ortho_mi_files:
            name = file_path if len(file_path) < 35 else '...' + file_path[-35:]
            self.printout('metric', f'Written: {name}')

    def get_front_score(self, node):
        """Number of unique taxa in front."""
        front_labels = get_front_labels(node)
        num_labels   = len(front_labels)
        num_taxa     = len(set([get_name(i) for i in front_labels]))
        if num_taxa == num_labels:
            return num_taxa
        return -1

    def get_back_score(self, node, root):
        """Number of unique taxa in back."""
        back_labels = get_back_labels(node, root)
        num_labels  = len(back_labels)
        num_taxa    = len(set([get_name(i) for i in back_labels]))
        if num_taxa == num_labels:
            return num_taxa
        return -1

    def prune_paralogs(self, curroot, cluster_id):
        """Split paralogous trees into ortholog groups."""
        pp_trees = []
        going    = True
        
        while going:
            highest      = 0
            highest_node = None
            score_hashes = {}
            
            for node in curroot.iternodes():
                front_score        = self.get_front_score(node)
                back_score         = self.get_back_score(node, curroot)
                score_hashes[node] = (front_score, back_score)
                if front_score > highest or back_score > highest:
                    highest_node = node
                    highest      = max(front_score, back_score)
            
            if highest >= self.minimum_taxa:
                curroot, done = self.prune_node(score_hashes[highest_node], highest_node, curroot, pp_trees)
                if done or len(curroot.leaves()) < self.minimum_taxa:
                    going = False
                    break
            else:
                going = False
                break
        
        if len(pp_trees) > 0:
            return self.process_pruned_trees(pp_trees, cluster_id)
        
        return []

    def prune_node(self, score_tuple, node, root, pp_trees):
        """Prune a node either front or back based on scores."""
        if score_tuple[0] > score_tuple[1]:
            self.printout('metric', f'Pruning front at node with {score_tuple[0]} taxa')
            pp_trees.append(node)
            par = node.prune()
            if par is not None and len(root.leaves()) >= 3:
                par, root = remove_kink(par, root)
            return root, node == root
        else:
            self.printout('metric', f'Pruning back at node with {score_tuple[1]} taxa')
            if node != root:
                par = node.parent
                par.remove_child(node)
                if par.parent is not None:
                    par, root = remove_kink(par, root)
            node.prune()
            pp_trees.append(root)
            if len(node.leaves()) >= 3:
                node, newroot = remove_kink(node, node)
            else:
                newroot = node
            return newroot, False

    def process_pruned_trees(self, pp_trees, cluster_id):
        """Process and write the pruned trees as individual ortholog files."""
        written_files = []
        count         = 1
        
        for tree in pp_trees:
            if tree.nchildren == 2:
                _, tree = remove_kink(tree, tree)
    
            trimmed_tree = self.trim_tree(tree)
            
            if trimmed_tree is not None and len(trimmed_tree.leaves()) >= self.minimum_taxa:
                outfile = self.outpath / f"{cluster_id}_MIortho{count}.tre"
                with open(outfile, "w") as outfile_handle:
                    outfile_handle.write(tostring(trimmed_tree) + ";\n")
                written_files.append(outfile)
                count += 1
        
        return written_files

    def trim_tree(self, tree):
        try:
            trimmer = TrimTips(
                dir_tree            = ".",
                dir_mafft           = ".",
                tree_ending         = ".tre",
                relative_cutoff     = self.relative_tip_cutoff,
                absolute_cutoff     = self.absolute_tip_cutoff,
                contig_dct          = {},
                hcolor              = self.hc,
                outlier_ratio       = self.outlier_ratio,
                max_trim_iterations = self.max_trim_iterations,
                min_tree_leaves     = self.min_tree_leaves
            )
            
            ete_tree = Tree(tostring(tree) + ";")
            
            trimmed_ete_tree = trimmer.trim(ete_tree, self.relative_tip_cutoff, self.absolute_tip_cutoff)
            
            if trimmed_ete_tree is not None:
                trimmed_str = trimmed_ete_tree.write(format=1)
                return parse(trimmed_str)
            
        except Exception as e:
            self.printout('error', f'Error trimming tree: {str(e)}')
        
        return tree