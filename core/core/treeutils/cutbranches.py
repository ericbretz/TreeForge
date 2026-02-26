import os
from pathlib import Path
from core.treeutils.newick  import parse, tostring
from core.treeutils.utils 	import get_front_names, remove_kink

class CutBranches:
	def __init__(self,
				dir_tree,
				dir_trimmed,
				tree_ending,
				min_taxa,
				cutoff,
				hcolor,
				min_subtree_taxa):

		self.dir_tree         = Path(dir_tree)
		self.dir_trimmed      = Path(dir_trimmed)
		self.tree_ending      = tree_ending
		self.tree_files       = list(self.dir_tree / f for f in os.listdir(self.dir_tree) if f.endswith(self.tree_ending))
		self.min_taxa         = min_taxa
		self.cutoff           = cutoff
		self.hcolor           = hcolor
		self.min_subtree_taxa = min_subtree_taxa
		self.error_count      = 0
		self.error_dict       = {'insufficient_taxa': 0, 'no_branches_cut': 0}
		self.total_counts     = {
			'total_taxa'    : 0,
			'total_tips'    : 0,
			'total_subtrees': 0,
			'large_subtrees': 0,
			'small_subtrees': 0
		}
		self.processed_count = 0

	def run(self):
		"""
		Cut branches that have too many tips.
		"""
		self.dir_trimmed.mkdir(parents=True, exist_ok=True)
        
		self.processed_count = 0
		self.error_dict      = {
			'insufficient_taxa': 0,
			'no_valid_subtrees': 0,
			'no_branches_cut'  : 0
		}
		
		metrics = {
			'cutting_threshold': self.cutoff,
			'total_trees'      : len(self.tree_files),
			'processed_count'  : 0,
			'total_counts'     : {
				'total_taxa'    : 0,
				'total_tips'    : 0,
				'total_subtrees': 0,
				'large_subtrees': 0,
				'small_subtrees': 0
			},
			'errors': {
				'insufficient_taxa': 0,
				'no_valid_subtrees': 0,
				'no_branches_cut'  : 0
			},
			'file_details': []
		}
        
		for tree_file in self.tree_files: 
			file_metrics = {
				'filename'      : tree_file.name,
				'stem'          : tree_file.stem,
				'taxa'          : 0,
				'tips'          : 0,
				'total_subtrees': 0,
				'large_subtrees': 0,
				'small_subtrees': 0,
				'subtree_sizes' : [],
				'status'        : '',
				'error'         : None
			}
			
			with open(tree_file, 'r') as infile:
				intree        = parse(infile.readline())
				raw_tree_size = len(get_front_names(intree))
				num_taxa      = self.count_taxa(intree)
				metrics['total_counts']['total_taxa'] += num_taxa
				metrics['total_counts']['total_tips'] += raw_tree_size
				
				file_metrics['taxa'] = num_taxa
				file_metrics['tips'] = raw_tree_size

				if num_taxa < self.min_taxa:
					metrics['errors']['insufficient_taxa'] += 1
					file_metrics['error']  = f'insufficient taxa ({num_taxa} < {self.min_taxa})'
					file_metrics['status'] = 'skipped'
					metrics['file_details'].append(file_metrics)
					continue

				subtrees = self.cut_long_internal_branches(intree, float(self.cutoff))
				file_metrics['total_subtrees'] = len(subtrees)
            
				if len(subtrees) == 0:
					metrics['errors']['no_valid_subtrees'] += 1
					file_metrics['error']  = 'no valid subtrees'
					file_metrics['status'] = 'skipped'
					metrics['file_details'].append(file_metrics)
					continue

				if len(subtrees) == 1 and self.count_taxa(subtrees[0]) == num_taxa:
					output_path = self.dir_trimmed / f"{tree_file.stem}.subtree"
					with open(tree_file, "r") as infile:
						with open(output_path, "w") as outfile:
							outfile.write(infile.read())
					self.processed_count += 1
					metrics['processed_count'] += 1
					file_metrics['status'] = 'unchanged'
					metrics['file_details'].append(file_metrics)
					continue

				count = 0
				valid_subtrees = 0
				for i, subtree in enumerate(subtrees, 1):
					subtree_taxa = self.count_taxa(subtree)
					if subtree_taxa >= self.min_subtree_taxa:
						valid_subtrees += 1
						if subtree.nchildren == 2:
							_, subtree = remove_kink(subtree, subtree)
						count += 1
						# output_path = self.dir_trimmed / f"{tree_file.stem}.subtree"  # whoops..
						output_path = self.dir_trimmed / f"{tree_file.stem}_{i}.subtree"
						with open(output_path, "w") as outfile:
							outfile.write(tostring(subtree)+";\n")
						file_metrics['subtree_sizes'].append(len(subtree.leaves()))
                            
				if count > 0:
					self.processed_count  			                += 1
					metrics['processed_count']                      += 1
					metrics['total_counts']['large_subtrees']       += count
					metrics['total_counts']['small_subtrees']       += len(subtrees) - count
					file_metrics['large_subtrees']                   = count
					file_metrics['small_subtrees']                   = len(subtrees) - count
					file_metrics['status']                           = 'cut'
				else:
					metrics['total_counts']['small_subtrees'] += len(subtrees)
					output_path = self.dir_trimmed / f"{tree_file.stem}.subtree"
					with open(tree_file, "r") as infile:
						with open(output_path, "w") as outfile:
							outfile.write(infile.read())
					self.processed_count 			  += 1
					metrics['processed_count']  	  += 1
					file_metrics['small_subtrees']     = len(subtrees)
					file_metrics['status']             = 'copied_original'
				
				metrics['file_details'].append(file_metrics)

		metrics['total_counts']['total_subtrees'] = metrics['total_counts']['large_subtrees'] + metrics['total_counts']['small_subtrees']
		
		return metrics

	def count_taxa(self, node):
		return len(set(get_front_names(node)))
	
	def get_clean_filename(self, filename):
		return filename.stem.replace('cluster', '')
	
	def cut_long_internal_branches(self, curroot, cutoff):
		going = True
		subtrees = []
		while going:
			going = False
			for node in curroot.iternodes():
				if node.istip or node == curroot: continue
				if node.nchildren == 1:
					node,curroot = remove_kink(node, curroot)
					going = True
					break
				child0,child1 = node.children[0],node.children[1]
				if node.length > cutoff:
					if not child0.istip and not child1.istip and child0.length+child1.length>cutoff:
						if self.count_taxa(child0) >= self.min_subtree_taxa:
							subtrees.append(child0)
						if self.count_taxa(child1) >= self.min_subtree_taxa:
							subtrees.append(child1)						
					else: subtrees.append(node)
					node = node.prune()
					if len(curroot.leaves()) > 2:
						node,curroot = remove_kink(node,curroot)
						going = True
					break
		if self.count_taxa(curroot) >= self.min_subtree_taxa:
			subtrees.append(curroot)
		return subtrees