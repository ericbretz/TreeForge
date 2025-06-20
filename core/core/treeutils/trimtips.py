import os
from pathlib import Path
from core.treeutils.newick import parse, tostring
from itertools import combinations
from core.treeutils.utils import remove_kink
from core.treeutils.phylo import Node
from ete3 import Tree

class TrimTips:
    def __init__(self,
                 dir_tree,          # Directory of the tree files
                 dir_mafft,         # Directory of the mafft files
                 tree_ending,       # Ending of the tree files
                 relative_cutoff,   # Relative cutoff for trimming tips
                 absolute_cutoff,   # Absolute cutoff for trimming tips
                 contig_dct,        # Dictionary of contigs
                 hcolor,            # Color for logging
                 ):
        
        self.dir_tree        = Path(dir_tree)
        self.dir_mafft       = Path(dir_mafft)
        self.tree_ending     = tree_ending
        self.treefile_files  = list(Path(os.path.join(self.dir_mafft, f)) for f in os.listdir(self.dir_mafft) if f.endswith(self.tree_ending))
        self.relative_cutoff = relative_cutoff
        self.absolute_cutoff = absolute_cutoff
        self.contig_dct      = contig_dct
        self.hcolor          = hcolor

    def check_contrast_outlier(self, node0, node1, above0, above1, relative_cutoff):
        if node0.is_leaf() and above0 > relative_cutoff:
            if above1 == 0.0 or above0/above1 > 20:  # Increased ratio threshold
                return node0
        if node1.is_leaf() and above1 > relative_cutoff:
            if above0 == 0.0 or above1/above0 > 20:  # Increased ratio threshold
                return node1
        return None

    def remove_a_tip(self, root, tip_node):
        if not tip_node.is_leaf():
            return root

        parent = tip_node.parent
        if parent is None:
            return root
        parent.remove_child(tip_node)

        if len(parent.children) == 1:
            grandparent = parent.parent
            if grandparent is not None:
                child = parent.children[0]
                child.length += parent.length
                grandparent.remove_child(parent)
                grandparent.add_child(child)

        if len(parent.children) == 0:
            if parent.parent is not None:
                parent.parent.remove_child(parent)

        remaining_leaves = len(root.get_leaves())
        if remaining_leaves > 3:
            return root
        else:
            return None

    def get_node_length(self, node):
        return node.get_feature('len') if node.get_feature('len') is not None else 0

    def set_node_length(self, node, length):
        node.add_feature('len', length)

    def handle_single_child(self, node):
        if len(node.children) != 1:
            return False
            
        parent = node.parent
        if parent is None:
            return False
            
        child = node.children[0]
        child.length += node.length
        parent.remove_child(node)
        parent.add_child(child)
        return True

    def trim(self, tree, relative_cutoff, absolute_cutoff):
        if tree is None:
            return None
        
        if len(tree.children) == 2:
            root_node = Node(tree)
            node_node = Node(tree)
            node_node, root_node = remove_kink(node_node, root_node)
            tree = Tree(root_node._ete_node.write(format=1))

        going = True
        iteration = 0
        max_iterations = 10  # Limit the number of iterations
        while going and tree is not None and len(tree.get_leaves()) > 3 and iteration < max_iterations:
            iteration += 1
            going = False
            
            for node in tree.traverse():
                if len(node.children) == 0:
                    self.set_node_length(node, node.length)
                    if node.length > absolute_cutoff:
                        tree = self.remove_a_tip(tree, node)
                        if tree is None:
                            return None
                        going = True
                        break
                elif len(node.children) == 1:
                    if self.handle_single_child(node):
                        going = True
                        break
                elif len(node.children) == 2:
                    child0, child1 = node.children[0], node.children[1]
                    above0 = self.get_node_length(child0)
                    above1 = self.get_node_length(child1)
                    self.set_node_length(node, ((above0 + above1) / 2.0) + node.length)
                    outlier = self.check_contrast_outlier(child0, child1, above0, above1, relative_cutoff)
                    if outlier is not None:
                        tree = self.remove_a_tip(tree, outlier)
                        if tree is None:
                            return None
                        going = True
                        break
                else:
                    children = node.children
                    total_len = sum(self.get_node_length(child) for child in children)
                    self.set_node_length(node, total_len / len(children))
                    
                    for child1, child2 in combinations(children, 2):
                        above1 = self.get_node_length(child1)
                        above2 = self.get_node_length(child2)
                        outlier = self.check_contrast_outlier(child1, child2, above1, above2, relative_cutoff)
                        if outlier is not None:
                            tree = self.remove_a_tip(tree, outlier)
                            if tree is None:
                                return None
                            going = True
                            break
                    if going:
                        break

        if tree is None:
            return None

        return tree

    def run(self):
        file_count = len(self.treefile_files)

        # Initialize metrics dictionary
        metrics = {
            'trimming_thresholds': {
                'relative_cutoff': self.relative_cutoff,
                'absolute_cutoff': self.absolute_cutoff
            },
            'total_trees': file_count,
            'processed_count': 0,
            'error_count': 0,
            'skipped_count': 0,
            'file_details': []
        }

        for tree_file in self.treefile_files:
            file_metrics = {
                'filename': tree_file.name,
                'stem': tree_file.stem,
                'status': '',
                'error': None,
                'leaves_before': 0,
                'leaves_after': 0
            }
            
            try:
                with open(tree_file, 'r') as infile:
                    tree_str = infile.readline().strip()
                    if not tree_str:
                        metrics['skipped_count'] += 1
                        file_metrics['error'] = 'empty tree string'
                        file_metrics['status'] = 'skipped'
                        metrics['file_details'].append(file_metrics)
                        continue
                    
                    try:
                        intree = parse(tree_str)
                        file_metrics['leaves_before'] = len(intree.get_leaves())
                    except Exception as e:
                        metrics['error_count'] += 1
                        file_metrics['error'] = f'parse error: {str(e)}'
                        file_metrics['status'] = 'error'
                        metrics['file_details'].append(file_metrics)
                        continue

                outtree = self.trim(intree, float(self.relative_cutoff), float(self.absolute_cutoff))
                if outtree is None:
                    metrics['skipped_count'] += 1
                    file_metrics['error'] = 'trimming resulted in null tree'
                    file_metrics['status'] = 'skipped'
                    metrics['file_details'].append(file_metrics)
                    continue

                file_metrics['leaves_after'] = len(outtree.get_leaves())
                if len(outtree.get_leaves()) < 4:
                    metrics['skipped_count'] += 1
                    file_metrics['error'] = f'insufficient leaves after trimming ({len(outtree.get_leaves())} < 4)'
                    file_metrics['status'] = 'skipped'
                    metrics['file_details'].append(file_metrics)
                    continue

                tree_name = tree_file.with_suffix('.tt')
                with open(tree_name, 'w') as outfile:
                    outfile.write(tostring(outtree) + ';\n')
                metrics['processed_count'] += 1
                file_metrics['status'] = 'processed'
                metrics['file_details'].append(file_metrics)

            except Exception as e:
                metrics['error_count'] += 1
                file_metrics['error'] = f'unexpected error: {str(e)}'
                file_metrics['status'] = 'error'
                metrics['file_details'].append(file_metrics)
                continue
        return metrics