import os
from ete3 import Tree
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
from core.treeutils.newick import parse, tostring
from core.treeutils.phylo  import Node
from core.treeutils.seq    import read_fasta_file
from core.treeutils.utils  import get_name, remove_kink

@dataclass
class TreeData:
    tree: Node
    cluster_id: str
    char_dict: Dict[str, int]

class MaskTips:
    def __init__(self,
                 dir_tree,
                 dir_cln,
                 para,
                 tree_ending,
                 hcolor
                 ):
        
        self.dir_tree    = Path(dir_tree)
        self.dir_cln     = Path(dir_cln)
        self.cln_files   = list(Path(os.path.join(self.dir_cln, f)) for f in os.listdir(self.dir_cln) if f.endswith('.cln'))
        self.para        = para
        self.tree_ending = tree_ending
        self.tree_files  = list(Path(os.path.join(self.dir_cln, f)) for f in os.listdir(self.dir_cln) if f.endswith(self.tree_ending))
        self.hcolor      = hcolor
        self.return_dict = {'tree': {}}
        
        self.metrics = {
            'trees_processed': 0,
            'trees_skipped': 0,
            'monophyletic_tips_masked': 0,
            'paraphyletic_tips_masked': 0,
            'clusters_processed': 0,
            'clusters_skipped': 0,
            'file_read_errors': 0,
            'tree_parse_errors': 0,
            'invalid_tree_errors': 0,
            'cln_read_errors': 0,
            'duplicate_cluster_ids': 0,
            'missing_cln_files': 0
        }
        
        cln_ids = {self.get_cluster_id(os.path.basename(f)) for f in self.cln_files}
        tree_ids = {self.get_cluster_id(os.path.basename(f)) for f in self.tree_files}
        matching_ids = cln_ids.intersection(tree_ids)
        
    def run(self):
        # Initialize comprehensive metrics dictionary
        metrics = {
            'masking_settings': {
                'mask_paraphyletic': self.para.lower() == 'y',
                'tree_ending': self.tree_ending
            },
            'file_counts': {
                'total_cln_files': len(self.cln_files),
                'total_tree_files': len(self.tree_files)
            },
            'processing_results': {
                'trees_processed': 0,
                'trees_skipped': 0,
                'clusters_processed': 0,
                'clusters_skipped': 0
            },
            'masking_results': {
                'monophyletic_tips_masked': 0,
                'paraphyletic_tips_masked': 0
            },
            'errors': {
                'file_read_errors': 0,
                'tree_parse_errors': 0,
                'invalid_tree_errors': 0,
                'cln_read_errors': 0,
                'duplicate_cluster_ids': 0,
                'missing_cln_files': 0
            },
            'file_details': []
        }
        
        mask_para = self.para.lower() == 'y'    #Mask paraphyletic tips?
        file_match = {}

        for cf in self.cln_files:
            cluster_id = self.get_cluster_id(os.path.basename(cf))
            if cluster_id in file_match:
                metrics['errors']['duplicate_cluster_ids'] += 1
                continue
            file_match[cluster_id] = cf
        
        for tf in self.tree_files:
            file_metrics = {
                'filename': tf.name,
                'cluster_id': self.get_cluster_id(os.path.basename(tf)),
                'status': '',
                'error': None,
                'leaves_before': 0,
                'leaves_after': 0,
                'monophyletic_masked': 0,
                'paraphyletic_masked': 0
            }
            
            cluster_id = self.get_cluster_id(os.path.basename(tf))
            if cluster_id not in file_match:
                metrics['processing_results']['clusters_skipped'] += 1
                metrics['errors']['missing_cln_files'] += 1
                file_metrics['error'] = 'no matching CLN file found'
                file_metrics['status'] = 'skipped'
                metrics['file_details'].append(file_metrics)
                continue

            metrics['processing_results']['clusters_processed'] += 1
            cln_path = Path(file_match[cluster_id])
            tree_data = self.process_tree_file(tf, cln_path, mask_para)
            if tree_data is None:
                metrics['processing_results']['trees_skipped'] += 1
                file_metrics['error'] = 'failed to process tree file'
                file_metrics['status'] = 'skipped'
                metrics['file_details'].append(file_metrics)
                continue
            
            file_metrics['leaves_before'] = len(tree_data.tree.leaves())
            metrics['processing_results']['trees_processed'] += 1
            curroot = self.mask_monophyletic_tips(tree_data.tree, tree_data.char_dict)
            if curroot is None:
                metrics['processing_results']['trees_skipped'] += 1
                file_metrics['error'] = 'monophyletic masking failed'
                file_metrics['status'] = 'skipped'
                metrics['file_details'].append(file_metrics)
                continue

            if mask_para:
                curroot = self.mask_paraphyletic_tips(curroot, tree_data.char_dict)
                if curroot is None:
                    metrics['processing_results']['trees_skipped'] += 1
                    file_metrics['error'] = 'paraphyletic masking failed'
                    file_metrics['status'] = 'skipped'
                    metrics['file_details'].append(file_metrics)
                    continue
                
            file_metrics['leaves_after'] = len(curroot.leaves())
            file_metrics['monophyletic_masked'] = self.metrics['monophyletic_tips_masked']
            file_metrics['paraphyletic_masked'] = self.metrics['paraphyletic_tips_masked']
            file_metrics['status'] = 'processed'
            metrics['file_details'].append(file_metrics)
            
            self.write_tree(curroot, os.path.basename(tf))

        # Update final metrics from the instance metrics
        metrics['masking_results']['monophyletic_tips_masked'] = self.metrics['monophyletic_tips_masked']
        metrics['masking_results']['paraphyletic_tips_masked'] = self.metrics['paraphyletic_tips_masked']
        metrics['errors']['file_read_errors'] = self.metrics['file_read_errors']
        metrics['errors']['tree_parse_errors'] = self.metrics['tree_parse_errors']
        metrics['errors']['invalid_tree_errors'] = self.metrics['invalid_tree_errors']
        metrics['errors']['cln_read_errors'] = self.metrics['cln_read_errors']
        metrics['errors']['duplicate_cluster_ids'] = self.metrics['duplicate_cluster_ids']
        metrics['errors']['missing_cln_files'] = self.metrics['missing_cln_files']

        return metrics

    def get_cluster_id(self, filename):
        return os.path.basename(filename).split('.')[0]
    
    def process_tree_file(self, tree_path, cln_path, mask_para):
        try:
            with open(tree_path, "r") as infile:
                tree_str = infile.readline().strip()
                try:
                    parsed_tree = parse(tree_str)
                    tree = Node(parsed_tree)
                except Exception as e:
                    self.metrics['tree_parse_errors'] += 1
                    return None
        except Exception as e:
            self.metrics['file_read_errors'] += 1
            return None
            
        if not self.is_valid_tree(tree):
            self.metrics['invalid_tree_errors'] += 1
            return None
            
        cluster_id = self.get_cluster_id(os.path.basename(tree_path))
        char_dict = {}
        
        try:
            for seq in read_fasta_file(cln_path):
                seq.seq = ''.join(c for c in seq.seq if c not in ['-', 'X', 'x', '?', '*'])
                char_dict[seq.name] = len(seq.seq)
        except Exception as e:
            self.metrics['cln_read_errors'] += 1
            return None
            
        return TreeData(tree, cluster_id, char_dict)
    
    def is_valid_tree(self, node: Node):
        if node is None:
            return False
        try:
            leaves = node.leaves()
            num_leaves = len(leaves)
            return num_leaves >= 4
        except Exception as e:
            return False
    
    def mask_monophyletic_tips(self, curroot, char_dict):
        if not self.is_valid_tree(curroot):
            return None
            
        while True:
            modified = False
            for node in curroot.iternodes():
                if not node.istip:
                    continue
                    
                for sister in node.get_sisters():
                    if sister.istip and get_name(node.label) == get_name(sister.label):
                        pruned_node = self.prune_node(node, sister, char_dict)
                        if pruned_node is None:
                            return None
                            
                        if self.is_valid_tree(curroot):
                            if (pruned_node == curroot and pruned_node.nchildren == 2) or \
                            (pruned_node != curroot and pruned_node.nchildren == 1):
                                pruned_node, curroot = remove_kink(pruned_node, curroot)
                                if curroot is None:
                                    return None
                        self.metrics['monophyletic_tips_masked'] += 1
                        modified = True
                        break
                        
                if modified:
                    break
                    
            if not modified:
                break
                
        return curroot

    def mask_paraphyletic_tips(self, curroot: Node, char_dict: Dict[str, int]) -> Optional[Node]:
        if not self.is_valid_tree(curroot):
            return None
            
        while True:
            modified = False
            for node in curroot.iternodes():
                if not node.istip:
                    continue
                    
                parent = node.parent
                if node == curroot or parent == curroot:
                    continue
                    
                for para in parent.get_sisters():
                    if para.istip and get_name(node.label) == get_name(para.label):
                        pruned_node = self.prune_node(node, para, char_dict)
                        if pruned_node is None:
                            return None
                            
                        if self.is_valid_tree(curroot):
                            if (pruned_node == curroot and pruned_node.nchildren == 2) or \
                            (pruned_node != curroot and pruned_node.nchildren == 1):
                                pruned_node, curroot = remove_kink(pruned_node, curroot)
                                if curroot is None:
                                    return None
                        self.metrics['paraphyletic_tips_masked'] += 1
                        modified = True
                        break
                        
                if modified:
                    break
                    
            if not modified:
                break
                
        return curroot

    def prune_node(self, node, sister, char_dict):
        if char_dict[node.label] > char_dict[sister.label]:
            return sister.prune()
        return node.prune()
    
    def write_tree(self, curroot, filename):
        outpath = Path(os.path.join(self.dir_tree, filename))
        outfile = outpath.with_suffix('').with_suffix('.mm')
        with open(outfile, "w") as outfile:
            outfile.write(tostring(curroot))