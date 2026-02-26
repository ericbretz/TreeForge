import csv
from typing import List, Dict, Set, Optional, Tuple, Any
from pathlib import Path

class Node:
    def __init__(self):
        self.label      : str              = ""
        self.length     : float            = 0.0
        self.time_length: float            = 0.0
        self.parent     : Optional['Node'] = None
        self.children   : List['Node']     = []
        self.data       : Dict[str, Any]   = {}
        self.istip      : bool             = False
        self.height     : int              = 0

    def add_child(self, child: 'Node') -> None:
        assert child not in self.children
        self.children.append(child)
        child.parent = self

    def remove_child(self, child: 'Node') -> None:
        assert child in self.children
        self.children.remove(child)
        child.parent = None

    def leaves(self, v: Optional[List['Node']] = None) -> List['Node']:
        if v is None:
            v = []
        if not self.children:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self) -> List['Node']:
        return [n for n in self.iternodes() if n.istip]

    def lvsnms(self) -> List[str]:
        return [n.label for n in self.iternodes() if n.istip]

    def iternodes(self, order: str = "preorder") -> 'Node':
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            yield from child.iternodes(order)
        if order.lower() == "postorder":
            yield self

    def prune(self) -> Optional['Node']:
        p = self.parent
        if p is not None:
            p.remove_child(self)
        return p

    def get_newick_repr(self, showbl: bool = True) -> str:
        ret = ""
        for i, child in enumerate(self.children):
            if i == 0:
                ret += "("
            ret += child.get_newick_repr(showbl)
            ret += ")" if i == len(self.children) - 1 else ","
        if self.label:
            ret += self.label
        if showbl:
            ret += f":{self.length}"
        return ret

    def set_height(self) -> None:
        if not self.children:
            self.height = 0
        else:
            tnode = self
            h = 0
            while tnode.children:
                if len(tnode.children) > 1 and tnode.children[1].length < tnode.children[0].length:
                    tnode = tnode.children[1]
                else:
                    tnode = tnode.children[0]
                h += tnode.length
            self.height = h

    def get_tip(self, label: str) -> Optional['Node']:
        for i in self.leaves():
            if label == i.label:
                return i
        return None

    def get_newick_otherlen(self, data: str) -> str:
        ret = ""
        for i, child in enumerate(self.children):
            if i == 0:
                ret += "("
            ret += child.get_newick_otherlen(data)
            ret += ")" if i == len(self.children) - 1 else ","
        if self.label:
            ret += self.label
        if data in self.data:
            ret += f":{self.data[data]}"
        return ret

class Quartet:
    def __init__(self, lfs: List[Set[str]], rts: List[Set[str]], length: float):
        self.left  : Set[str]       = set().union(*lfs)
        self.right : Set[str]       = set().union(*rts)
        self.lefts : List[Set[str]] = lfs
        self.rights: List[Set[str]] = rts
        self.length: float          = length

    def conflict(self, inbp: 'Quartet') -> bool:
        if (self.right & inbp.right and self.right & inbp.left and
            self.left & inbp.right and self.left & inbp.left):
            return True
        if (self.left & inbp.left and self.left & inbp.right and
            self.right & inbp.left and self.right & inbp.right):
            return True
        return False

    def same(self, inbp: 'Quartet') -> bool:
        if len(inbp.right) != len(self.right) and len(inbp.right) != len(self.left):
            return False
        return ((inbp.right == self.right and inbp.left == self.left) or
                (inbp.right == self.left and inbp.left == self.right))

    def match(self, inq: 'Quartet') -> bool:
        return self.same(inq)

def get_quartet(nd: Node, rt: Node) -> Optional[Quartet]:
    if not nd.children or nd == rt:
        return None
    rights = [set(child.lvsnms()) for child in nd.children]
    right  = set().union(*rights)
    p      = nd.parent
    lefts  = []
    if p:
        for sib in p.children:
            if set(sib.lvsnms()) & right:
                continue
            lefts.append(set(sib.lvsnms()))
        if p != rt:
            out = set(rt.lvsnms())
            lefts.append(out - set(p.lvsnms()))
    return Quartet(lefts, rights, nd.length)

def get_quartets(rt: Node) -> List[Quartet]:
    return [q for i in rt.iternodes() for q in [get_quartet(i, rt)] if q is not None]

def build(instr: str) -> Node:
    root         = None
    index        = 0
    n            = len(instr)
    nextchar     = instr[index]
    begining     = True
    keepgoing    = True
    current_node = None
    name         = ""
    branch       = ""
    while keepgoing:
        if nextchar == "(" and begining:
            root = Node()
            current_node = root
            begining = False
        elif nextchar == "(" and not begining:
            if current_node is None:
                raise ValueError("Unexpected '(' - no parent node available")
            newnode = Node()
            current_node.add_child(newnode)
            current_node = newnode
        elif nextchar == ',':
            if current_node is None:
                raise ValueError("Unexpected ',' - no parent node available")
            current_node = current_node.parent
        elif nextchar == ")":
            if current_node is None:
                raise ValueError("Unexpected ')' - no parent node available")
            current_node = current_node.parent
            index += 1
            if index < n:
                nextchar = instr[index]
            else:
                break
            while nextchar not in ',):;[' and not nextchar.isspace():
                name += nextchar
                index += 1
                if index >= n:
                    break
                nextchar = instr[index]
            if current_node is not None:
                current_node.label = name
            index -= 1
        elif nextchar == ';':
            break
        elif nextchar == ':':
            index += 1
            if index < n:
                nextchar = instr[index]
            else:
                break
            while nextchar not in ',):;[' and not nextchar.isspace():
                branch += nextchar
                index += 1
                if index >= n:
                    break
                nextchar = instr[index]
            try:
                if current_node is not None:
                    current_node.length = float(branch)
            except Exception:
                if current_node is not None:
                    current_node.length = 0.0
            index -= 1
        elif nextchar == ' ':
            index += 1
            if index < n:
                nextchar = instr[index]
            else:
                break
        else:
            if current_node is None:
                raise ValueError("Invalid tree structure: current_node is None when encountering leaf node")
            newnode = Node()
            current_node.add_child(newnode)
            current_node       = newnode
            current_node.istip = True
            while nextchar not in ',):;[' and not nextchar.isspace():
                name += nextchar
                index += 1
                if index >= n:
                    break
                nextchar = instr[index]
            current_node.label  = name
            index              -= 1
        if index < n - 1:
            index    += 1
            nextchar  = instr[index]
        name = ""
        branch = ""
    return root

def mean(array: List[float]) -> float:
    return sum(array) / len(array) if array else 0.0

def median(array: List[float]) -> float:
    n = len(array)
    if n == 0:
        return 0.0
    array = sorted(array)
    mid = n // 2
    return (array[mid] + array[mid - 1]) / 2.0 if n % 2 == 0 else array[mid]

def minv(array: List[float]) -> float:
    return min(array) if array else 0.0

def maxv(array: List[float]) -> float:
    return max(array) if array else 0.0

def ssd(array: List[float]) -> float:
    c = mean(array)
    return sum((x - c) ** 2 for x in array)

def stddev(array: List[float], ddof: float = 0.0) -> float:
    n = len(array)
    if n <= 1:
        return 0.0
    ss   = ssd(array)
    pvar = ss / (n - ddof)
    return pvar ** 0.5

def ci(array: List[float], zval: float) -> Tuple[float, float]:
    n = len(array)
    if n == 0:
        return 0.0, 0.0
    std      = stddev(array)
    sqrt_pop = n ** 0.5
    ci_l     = mean(array) - zval * (std / sqrt_pop)
    ci_h     = mean(array) + zval * (std / sqrt_pop)
    return ci_l, ci_h

def prepare_sp_tree(sp_tree: Node) -> Tuple[Dict[Quartet, List[float]], Dict[Quartet, List[Node]]]:
    sp_tree_quartets: Dict[Quartet, List[float]] = {}
    q_to_n: Dict[Quartet, List[Node]] = {}
    for i in sp_tree.iternodes():
        i.data['q'] = get_quartet(i, sp_tree)
        if i.istip:
            i.data['qln']      = []
            i.data['qlntrees'] = set()
        elif i.data['q'] is not None:
            sp_tree_quartets[i.data['q']] = []
            q_to_n[i.data['q']]           = []
    return sp_tree_quartets, q_to_n

def process_gene_trees(g_tree: List[str], sp_quartet_tree: Dict[Quartet, List[float]], sp_tree: Node, supval: float) -> Tuple[Node, Dict[Quartet, List[float]]]:
    processed_trees  = 0
    matched_quartets = 0
    
    for count, gene in enumerate(g_tree):
        try:
            tree = build(gene.strip())
            processed_trees += 1
            
            for i in sp_quartet_tree:
                for j in tree.iternodes():
                    gqu = get_quartet(j, tree)
                    if gqu and i.match(gqu):
                        sp_quartet_tree[i].append(j.length)
                        matched_quartets += 1
                        check             = {list(k)[0] for k in gqu.lefts + gqu.rights if len(k) == 1}
                        keep              = {list(k)[0] for k in i.lefts + i.rights if len(k) == 1 and list(k)[0] in check}
                        for k in keep:
                            ln = tree.get_tip(k).length
                            st = sp_tree.get_tip(k)
                            if tree not in st.data["qlntrees"]:
                                st.data["qlntrees"].add(tree)
                                st.data["qln"].append(ln)
        except Exception as e:
            continue
    
    return sp_tree, sp_quartet_tree

def summarizer(sp_tree: Node, sp_quartet_tree: Dict[Quartet, List[float]], output_dir: Optional[str]) -> None:
    output_dir     = Path(output_dir)
    
    bes_trees_dir = output_dir / "BES_Trees"
    bes_trees_dir.mkdir(exist_ok=True)
    
    molecular_path = output_dir / "SpeciesTree.molecular.tre"
    verbose_path   = output_dir / "SpeciesTree.verbose.csv"
    v_out          = open(verbose_path, "w", newline='')
    outf           = open(molecular_path, "w")
    v_writer       = csv.writer(v_out) if v_out else None
    
    for i in sp_tree.iternodes(): 
        if i == sp_tree:
            continue
        holder = i.data["qln"] if not i.children else sp_quartet_tree[i.data["q"]]
        if v_writer:
            temp = list(map(str, holder))
            v_writer.writerow([(i.label if not i.children else i.get_newick_repr(False))] + temp)
        if not holder:
            mean_ = median_ = min_ = max_ = CIL = CIH = 0.0
            i.data["concord"] = 0
        else:
            mean_             = mean(holder)
            median_           = median(holder)
            min_              = minv(holder)
            max_              = maxv(holder)
            i.data["concord"] = len(holder)
            if len(holder) > 1:
                CIL, CIH = ci(holder, 1.96)
            else:
                CIL = CIH = 0.0
        i.data["mean"]   = mean_
        i.data["median"] = median_
        i.data["min"]    = min_
        i.data["max"]    = max_
        i.data["cih"]    = CIH
        i.data["cil"]    = CIL
    
    stat_names = ["mean", "median", "min", "max", "cil", "cih", "concord"]
    for stat in stat_names:
        tree_string = sp_tree.get_newick_otherlen(stat) + ";"
        stat_file = bes_trees_dir / f"SpeciesTree.{stat}.tre"
        with open(stat_file, "w") as f:
            f.write(tree_string + "\n")
    
    molecular_tree = sp_tree.get_newick_otherlen("mean") + ";"
    if outf:
        outf.write(molecular_tree + "\n")
    
    if v_out:
        v_out.close()
    if outf:
        outf.close()

def validate_tree_structure(node: Node, visited: Optional[Set[int]] = None, depth: int = 0, max_depth: int = 10000) -> Tuple[bool, str, int]:
    if visited is None:
        visited = set()
    
    node_id = id(node)
    if node_id in visited:
        return False, f"Circular reference detected at depth {depth}", depth
    
    if depth > max_depth:
        return False, f"Tree depth exceeds maximum ({max_depth})", depth
    
    visited.add(node_id)
    
    max_child_depth = depth
    for child in node.children:
        is_valid, error_msg, child_depth = validate_tree_structure(child, visited.copy(), depth + 1, max_depth)
        if not is_valid:
            return False, error_msg, child_depth
        max_child_depth = max(max_child_depth, child_depth)
    
    return True, "", max_child_depth

def get_tree_stats(node: Node) -> Dict[str, Any]:
    """Get statistics about a tree structure."""
    try:
        leaves_count   = len(list(node.leaves()))
        internal_nodes = sum(1 for n in node.iternodes() if not n.istip)
        total_nodes    = sum(1 for _ in node.iternodes())
        max_depth      = 0
        stack          = [(node, 0)]
        while stack:
            current, depth = stack.pop()
            max_depth = max(max_depth, depth)
            for child in current.children:
                stack.append((child, depth + 1))
        
        return {
            'leaves'        : leaves_count,
            'internal_nodes': internal_nodes,
            'total_nodes'   : total_nodes,
            'max_depth'     : max_depth
        }
    except Exception as e:
        return {'error': str(e)}

class BES:
    def __init__(self, printout=None):
        self.printout = printout

    def run(self, species_tree_file, gene_trees_file, output_prefix):
        try:
            with open(species_tree_file, 'r') as f:
                species_tree_str = f.read().strip()
            
            with open(gene_trees_file, 'r') as f:
                gene_trees = [line.strip() for line in f if line.strip()]
            
            sp_tree = build(species_tree_str)
            
            is_valid, error_msg, depth = validate_tree_structure(sp_tree)
            if not is_valid:
                if self.printout:
                    self.printout('error', f'Invalid species tree structure: {error_msg}')
                return False
            
            sp_quartet_tree, q_to_n = prepare_sp_tree(sp_tree)
            sp_tree, sp_quartet_tree = self._process_gene_trees(gene_trees, sp_quartet_tree, sp_tree)
            summarizer(sp_tree, sp_quartet_tree, output_prefix)
            
        except RecursionError:
            if self.printout:
                self.printout('error', 'Recursion depth exceeded in BES')
                self.printout('error', 'Tree may have circular references or excessive depth')
            return False
        except Exception as e:
            if self.printout:
                self.printout('error', f'BES pipeline error: {e}')
            else:
                msg = f'BES pipeline error: {e}'
                print(msg if len(msg) <= 80 else '...' + msg[-77:])
            return False
        return True
    
    def _process_gene_trees(self, g_tree: List[str], 
                           sp_quartet_tree: Dict[Quartet, List[float]], 
                           sp_tree: Node) -> Tuple[Node, Dict[Quartet, List[float]]]:
        """Process gene trees and match quartets."""
        processed_trees = 0
        
        for count, gene in enumerate(g_tree):
            try:
                tree = build(gene.strip())
                processed_trees += 1
                
                for i in sp_quartet_tree:
                    for j in tree.iternodes():
                        gqu = get_quartet(j, tree)
                        if gqu and i.match(gqu):
                            sp_quartet_tree[i].append(j.length)
                            check = {list(k)[0] for k in gqu.lefts + gqu.rights if len(k) == 1}
                            keep  = {list(k)[0] for k in i.lefts + i.rights if len(k) == 1 and list(k)[0] in check}
                            for k in keep:
                                ln = tree.get_tip(k).length
                                st = sp_tree.get_tip(k)
                                if tree not in st.data["qlntrees"]:
                                    st.data["qlntrees"].add(tree)
                                    st.data["qln"].append(ln)
            except RecursionError:
                if self.printout:
                    self.printout('error', f'Recursion error at gene tree {count}')
                raise
            except Exception:
                continue
        
        return sp_tree, sp_quartet_tree