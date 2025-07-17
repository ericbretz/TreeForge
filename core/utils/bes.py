import sys
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
        if self.conflict(inq):
            return False
        lmatchedl, lmatchedr = set(), set()
        for i in self.lefts:
            tl = 0
            for j, l in enumerate(inq.lefts):
                if i & l:
                    tl += 1
                    lmatchedl.add(j)
            for j, r in enumerate(inq.rights):
                if i & r:
                    tl += 1
                    lmatchedr.add(j)
            if tl > 1:
                return False
        samed = True
        lmatched = lmatchedl
        if lmatchedl and lmatchedr:
            return False
        if lmatchedr:
            samed = False
            lmatched = lmatchedr
        if len(lmatched) < 2:
            return False
        rmatchedl, rmatchedr = set(), set()
        for i in self.rights:
            tl = 0
            for j, l in enumerate(inq.lefts):
                if i & l:
                    tl += 1
                    rmatchedl.add(j)
            for j, r in enumerate(inq.rights):
                if i & r:
                    tl += 1
                    rmatchedr.add(j)
            if tl > 1:
                return False
        if rmatchedl and rmatchedr:
            return False
        rmatched = rmatchedr
        if samed and rmatchedl:
            return False
        if not samed and rmatchedr:
            return False
        if not samed:
            rmatched = rmatchedl
        if len(rmatched) < 2:
            return False
        return True

def get_quartet(nd: Node, rt: Node) -> Optional[Quartet]:
    if not nd.children or nd == rt:
        return None
    rights = [set(child.lvsnms()) for child in nd.children]
    right = set().union(*rights)
    p = nd.parent
    lefts = []
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
    return [q for i in rt.iternodes() if (q := get_quartet(i, rt)) is not None]

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
                raise ValueError("Invalid tree structure: current_node is None when encountering '('")
            newnode = Node()
            current_node.add_child(newnode)
            current_node = newnode
        elif nextchar == ',':
            if current_node is None:
                raise ValueError("Invalid tree structure: current_node is None when encountering ','")
            current_node = current_node.parent
        elif nextchar == ")":
            if current_node is None:
                raise ValueError("Invalid tree structure: current_node is None when encountering ')'")
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
            current_node = newnode
            current_node.istip = True
            while nextchar not in ',):;[' and not nextchar.isspace():
                name += nextchar
                index += 1
                if index >= n:
                    break
                nextchar = instr[index]
            current_node.label = name
            index -= 1
        if index < n - 1:
            index += 1
            nextchar = instr[index]
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
            i.data['qln'] = []
            i.data['qlntrees'] = set()
        elif i.data['q'] is not None:
            sp_tree_quartets[i.data['q']] = []
            q_to_n[i.data['q']] = []
    return sp_tree_quartets, q_to_n

def process_gene_trees(g_tree: List[str], sp_quartet_tree: Dict[Quartet, List[float]], sp_tree: Node, supval: float) -> Tuple[Node, Dict[Quartet, List[float]]]:
    for count, gene in enumerate(g_tree):
        tree = build(gene.strip())
        for i in sp_quartet_tree:
            for j in tree.iternodes():
                gqu = get_quartet(j, tree)
                if gqu and i.match(gqu):
                    if not j.label:
                        j.label = 0.0
                    if supval <= float(j.label):
                        sp_quartet_tree[i].append(j.length)
                        check = {list(k)[0] for k in gqu.lefts + gqu.rights if len(k) == 1}
                        keep = {list(k)[0] for k in i.lefts + i.rights if len(k) == 1 and list(k)[0] in check}
                        for k in keep:
                            ln = tree.get_tip(k).length
                            st = sp_tree.get_tip(k)
                            if tree not in st.data["qlntrees"]:
                                st.data["qlntrees"].add(tree)
                                st.data["qln"].append(ln)
    return sp_tree, sp_quartet_tree

def summarizer(sp_tree: Node, sp_quartet_tree: Dict[Quartet, List[float]], output_dir: Optional[str]) -> None:
    output_dir     = Path(output_dir)
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
    array = ["mean", "median", "min", "max", "cil", "cih", "concord"]
    for stat in array:
        outstr = sp_tree.get_newick_otherlen(stat) + ";"
        if outf:
            outf.write(outstr + "\n")
    if v_out:
        v_out.close()
    if outf:
        outf.close()

class BES:
    def __init__(self, printout=None):
        self.printout = printout

    def run(self, species_tree_file, gene_trees_file, support, output_prefix):
        try:
            with open(species_tree_file, 'r') as f:
                species_tree_str = f.read().strip()
            
            with open(gene_trees_file, 'r') as f:
                gene_trees = [line.strip() for line in f if line.strip()]
            
            sp_tree                  = build(species_tree_str)
            sp_quartet_tree, q_to_n  = prepare_sp_tree(sp_tree)
            sp_tree, sp_quartet_tree = process_gene_trees(gene_trees, sp_quartet_tree, sp_tree, support)
            summarizer(sp_tree, sp_quartet_tree, output_prefix)
        except Exception as e:
            if self.printout:
                self.printout('error', f'BES pipeline error: {e}')
            else:
                print(f'BES pipeline error: {e}')
            return False
        return True

# if __name__ == "__main__":
#     bes = BES()
#     species_tree_str = '/home/eric/02_Tests/joeTest_4/TreeForge/SpeciesTree.tre'
#     gene_trees       = '/home/eric/02_Tests/joeTest_4/TreeForge/super/concat.tre'
#     support          = 0.0
#     output_prefix    = '/home/eric/02_Tests/besTest/output3'
#     bes.run(species_tree_str, gene_trees, support, output_prefix)