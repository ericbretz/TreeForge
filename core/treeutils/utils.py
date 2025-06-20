from __future__ import annotations
import os
from pathlib import Path
from typing import List, Set, Dict, Optional, Iterator
from dataclasses import dataclass
from functools import lru_cache
from Bio.SeqIO import parse as seq_parse
from core.treeutils.newick import parse
from core.treeutils.phylo import reroot, Node

@dataclass
class TreeData:
	tree: Node
	cluster_id: str
	seq_dict: Dict[str, str]

@lru_cache(maxsize=128)
def get_name(label: str) -> str:
	return label.split("@")[0]
	
def get_clusterID(filename: str) -> str:
	return Path(filename).stem

def get_front_labels(node: Node) -> List[str]:
	return [i.label for i in node.leaves()]

@lru_cache(maxsize=128)
def get_back_labels(node: Node, root: Node) -> Set[str]:
	all_labels = set(get_front_labels(root))
	front_labels = set(get_front_labels(node))
	return all_labels - front_labels
	
def get_front_names(node: Node) -> List[str]:
	return [get_name(i) for i in get_front_labels(node)]

@lru_cache(maxsize=128)
def get_back_names(node: Node, root: Node) -> List[str]:
	return [get_name(i) for i in get_back_labels(node,root)]


def remove_kink(node, curroot):
	if node == curroot and len(curroot.children) == 2:
		if curroot.children[0].is_leaf():
			curroot = reroot(curroot, curroot.children[1])
		else: 
			curroot = reroot(curroot, curroot.children[0])
	
	if not node.children:
		return None, curroot
		
	# Can't remove kink if node has no parent
	if node.parent is None:
		return None, curroot
		
	length = node.dist + (node.children[0]).dist
	par = node.parent
	kink = node
	node = node.children[0]
	par.remove_child(kink)
	par.add_child(node)
	node.dist = length
	return None, curroot

def pass_boot_filter(node: Node, min_ave_boot: float) -> bool:
	total = 0.0
	count = 0
	for i in node.iternodes():
		if not i.istip and i.parent != None:
			total += float(i.label)
			count += 1
	if count == 0:
		return True
	boot_average = total / float(count)
	print(boot_average)
	return boot_average >= float(min_ave_boot)

def get_ortho_from_rooted_inclade(inclade: Node) -> List[Node]:
	assert inclade.nchildren == 2, "input clade not properly rooted"
	orthologs = []
	clades = [inclade]
	while True:
		newclades = []
		for clade in clades:
			num_taxa = len(set(get_front_names(clade)))
			num_tips = len((get_front_labels(clade)))
			if num_taxa == num_tips:
				orthologs.append(clade)
			else:
				for node in clade.iternodes(order=0):
					if node.istip: continue
					child0,child1 = node.children[0],node.children[1]
					name_set0 = set(get_front_names(child0))
					name_set1 = set(get_front_names(child1))
					if len(name_set0.intersection(name_set1)) > 0:
						if node == clade:
							newclades += [child0,child1]
						elif len(name_set0) > len(name_set1):
							node.remove_child(child1)
							child1.prune()
							node,clade = remove_kink(node,clade)
							newclades += [clade,child1]
						else:
							node.remove_child(child0)
							child0.prune()
							node,clade = remove_kink(node,clade)
							newclades += [clade,child0]
						break
		if newclades == []: break
		clades = newclades
	return orthologs

def extract_rooted_ingroup_clades(root: Node, ingroups: List[str], outgroups: List[str], min_ingroup_taxa: int) -> List[Node]:
	inclades = []
	ingroups_set = set(ingroups)
	outgroups_set = set(outgroups)
	
	while True:
		max_score, direction, max_node = 0, "", None
		
		node_data = {}
		for node in root.iternodes():
			front_names = set(get_front_names(node))
			back_names = set(get_back_names(node, root))
			node_data[node] = (front_names, back_names)
		
		for node in root.iternodes():
			front_names, back_names = node_data[node]
			
			if any(name in outgroups_set for name in front_names):
				front = -1
			else:
				front = sum(1 for name in front_names if name in ingroups_set)
			
			if any(name in outgroups_set for name in back_names):
				back = -1
			else:
				back = sum(1 for name in back_names if name in ingroups_set)
			
			if front > max_score:
				max_score, direction, max_node = front, "front", node
			if back > max_score:
				max_score, direction, max_node = back, "back", node
		
		if max_score >= min_ingroup_taxa:
			if direction == "front":
				inclades.append(max_node)
				kink = max_node.prune()
				if len(root.leaves()) > 3:
					newnode, root = remove_kink(kink, root)
				else: 
					break
			elif direction == "back":
				par = max_node.parent
				par.remove_child(max_node)
				max_node.prune()
				inclades.append(reroot(root, par))
				if len(max_node.leaves()) > 3:
					max_node, root = remove_kink(max_node, max_node)
				else: 
					break
		else: 
			break
	return inclades

def write_tree(all_fasta: str, trim_dir: str, iter_dir: str, tree_file_ending: str) -> None:
	fasta = all_fasta
	treDIR = Path(trim_dir)
	outIncr = int(iter_dir.split('_')[-1]) + 1
	outDIR = Path(iter_dir.split('_')[:-1]) / str(outIncr)
	tree_file_ending = '.subtree'

	print("Reading fasta file",fasta)
	seqDICT = {}
	for s in seq_parse(fasta, 'fasta'):
		seqDICT[s.id] = str(s.seq)
	print("Writing fasta files")
	filecount = 0
	for tree_file in treDIR.glob(f'*{tree_file_ending}'):
		print(tree_file.name)
		filecount += 1
		with open(tree_file, "r") as infile:
			intree = parse(infile.readline())
		clusterID = get_clusterID(tree_file.name)
		if clusterID.endswith("rr"):
			outname = outDIR / f"{clusterID}_rr.fa"
		else: outname = outDIR / f"{clusterID}rr.fa"
		with open(outname,"w") as outfile:
			for label in get_front_labels(intree):
				outfile.write(f">{label}\n{seqDICT[label]}\n")
	assert filecount > 0, f"No file ends with {tree_file_ending} found in {treDIR}"