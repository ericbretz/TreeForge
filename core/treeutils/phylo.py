
from functools import lru_cache
from ete3      import TreeNode

PREORDER     = 0
POSTORDER    = 1
BRANCHLENGTH = 0
INTERNODES   = 1

class Node:
    def __init__(self, ete_node=None):
        self._ete_node       = ete_node if ete_node is not None else TreeNode()
        self.data            = {}
        self.isroot          = False
        self.istip           = self._ete_node.is_leaf()
        self._label          = self._ete_node.name
        self._length         = self._ete_node.dist
        self._parent         = None
        self._children       = []
        self.nchildren       = 0
        self.excluded_dists  = []
        self._leaf_distances = None
        self._path_to_root   = None
        self._size           = None
        self._features       = {}
        
        if not self.istip:
            self._children = [Node(child) for child in self._ete_node.children]
            self.nchildren = len(self._children)
            for child in self._children:
                child._parent = self

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label         = value
        self._ete_node.name = value

    @property
    def name(self):
        return self._label

    @name.setter
    def name(self, value):
        self._label         = value
        self._ete_node.name = value

    @property
    def dist(self):
        return self._length

    @dist.setter
    def dist(self, value):
        self._length        = value
        self._ete_node.dist = value

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        self._length        = value
        self._ete_node.dist = value

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value
        if value:
            self._ete_node.up = value._ete_node
        else:
            self._ete_node.up = None

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, value):
        self._children          = value
        self._ete_node.children = [child._ete_node for child in value]
        self.nchildren          = len(value)

    def is_leaf(self):
        return self._ete_node.is_leaf()

    def get_leaves(self):
        return [Node(leaf) for leaf in self._ete_node.get_leaves()]

    def traverse(self):
        """Return an iterator over all nodes in the tree."""
        for node in self._ete_node.traverse():
            yield Node(node)

    def add_feature(self, name, value):
        """Add a feature to the node."""
        self._features[name] = value
        self._ete_node.add_feature(name, value)

    def get_feature(self, name):
        """Get a feature from the node."""
        return self._features.get(name)

    def order_subtrees_by_size(self, n2s = None, recurse = False, reverse = False):
        if n2s is None:
            n2s = node2size(self)
        if not self.istip:
            v = [(n2s[c], c.label, c) for c in self._children]
            v.sort()
            if reverse:
                v.reverse()
            self._children = [x[-1] for x in v]
            if recurse:
                for c in self._children:
                    c.order_subtrees_by_size(n2s, recurse=True, reverse=reverse)

    def add_child(self, child):
        assert child not in self._children
        self._ete_node.add_child(child._ete_node)
        self._children.append(child)
        child._parent        = self
        self.nchildren       = len(self._children)
        self._leaf_distances = None
        self._path_to_root   = None
        self._size           = None

    def remove_child(self, child):
        assert child in self._children
        self._ete_node.remove_child(child._ete_node)
        self._children.remove(child)
        child._parent        = None
        self.nchildren       = len(self._children)
        self._leaf_distances = None
        self._path_to_root   = None
        self._size           = None
        
    def leaves(self):
        """Return a list of all leaf nodes in the tree."""
        if self.istip:
            return [self]
        leaves = []
        for child in self._children:
            leaves.extend(child.leaves())
        return leaves

    def iternodes(self, order = POSTORDER, v = None):
        if order == PREORDER:
            yield self
        for child in self._children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self
        
    def descendants(self, order = PREORDER, v = None):
        if v is None:
            v = []
        assert order in (PREORDER, POSTORDER)
        for child in self._children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child._children:
                child.descendants(order, v)
        return v
    
    def find_descendant(self, label):
        nodes = self._ete_node.search_nodes(name=label)
        return Node(nodes[0]) if nodes else None

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p

    def graft(self, node):
        parent = self.parent
        parent.remove_child(self)
        n = Node()
        n.add_child(self)
        n.add_child(node)
        parent.add_child(n)

    def leaf_distances(self, store=None, measure=BRANCHLENGTH):
        if self._leaf_distances is not None:
            return self._leaf_distances
            
        if store is None:
            store = {}
        leaf2len = {}
        if self._children:
            for child in self._children:
                if measure == BRANCHLENGTH:
                    assert child.length is not None
                    dist = child.length
                elif measure == INTERNODES:
                    dist = 1
                else:
                    raise "InvalidMeasure"
                child.leaf_distances(store, measure)
                if child.istip:
                    leaf2len[child.label] = dist
                else:
                    for k, v in store[child].items():
                        leaf2len[k] = v + dist
        else:
            leaf2len[self] = {self.label: 0}
        store[self] = leaf2len
        self._leaf_distances = store
        return store
        
    def rootpath(self):
        if self._path_to_root is not None:
            return self._path_to_root
            
        path = []
        n = self
        while 1:
            path.append(n)
            if n.parent:
                n = n.parent
            else:
                break
        self._path_to_root = path
        return path
            
    def subtree_mapping(self, labels, clean=False):
        d = {}
        oldtips = [ x for x in self.leaves() if x.label in labels ]
        for tip in oldtips:
            path = list(tip.rootpath())
            for node in path:
                if node not in d:
                    newnode          = Node()
                    newnode.istip    = node.istip
                    newnode.length   = node.length
                    newnode.label    = node.label
                    d[node]          = newnode
                    d[newnode]       = node
                else:
                    newnode = d[node]

                for child in node._children:
                    if child in d:
                        newchild = d[child]
                        if newchild not in newnode._children:
                            newnode.add_child(newchild)
        d["oldroot"] = self
        d["newroot"] = d[self]
        if clean:
            n = d["newroot"]
            while 1:
                if n.nchildren == 1:
                    oldnode = d[n]
                    del d[oldnode]
                    del d[n]
                    child           = n._children[0]
                    child._parent   = None
                    child.isroot    = True
                    d["newroot"]    = child
                    d["oldroot"]    = d[child]
                    n = child
                else:
                    break
                    
            for tip in oldtips:
                newnode = d[tip]
                while 1:
                    newnode = newnode.parent
                    oldnode = d[newnode]
                    if newnode.nchildren == 1:
                        child = newnode._children[0]
                        if newnode.length:
                            child.length += newnode.length
                        newnode.remove_child(child)
                        if newnode.parent:
                            parent = newnode.parent
                            parent.remove_child(newnode)
                            parent.add_child(child)
                        del d[oldnode]; del d[newnode]
                    if not newnode.parent:
                        break
        return d
        
    def get_sisters(self):
        if self.parent == None:
            return
        return [Node(sister) for sister in self._ete_node.get_sisters()]

    def write(self, format=1):
        """Write the tree to a string in Newick format."""
        return self._ete_node.write(format=format)

@lru_cache(maxsize=None)
def node2size(node, d=None):
    if d is None:
        d = {}
    size = int(node.istip)
    if not node.istip:
        for child in node._children:
            node2size(child, d)
    d[node] = size
    return d

def reroot(oldroot, newroot):
    oldroot.isroot = False
    newroot.isroot = True
    v = []
    n = newroot
    while 1:
        v.append(n)
        if not n.parent: break
        n = n.parent
    v.reverse()
    for i, cp in enumerate(v[:-1]):
        node = v[i+1]
        cp.remove_child(node)
        node.add_child(cp)
        cp.length = node.length
        cp.label  = node.label
    return newroot

def getMRCA(innames, tree):
    mrca = None
    if len(innames) == 1:
        return None
    else:
        outgroup = []
        for name in innames:
            for i in range(len(tree.leaves())):
                if tree.leaves()[i].label == name:
                    outgroup.append(tree.leaves()[i])
        cur2     = None
        tempmrca = None
        cur1     = outgroup.pop()
        while len(outgroup)>0:
            cur2     = outgroup.pop()
            tempmrca = getMRCATraverse(cur1,cur2)
            cur1     = tempmrca
        mrca = cur1
    return mrca

def getMRCATraverse(curn1, curn2):
    mrca   = None
    path1  = []
    parent = curn1
    path1.append(parent)
    while parent != None:
        path1.append(parent);
        if parent.parent != None:
            parent = parent.parent
        else:
            break
    parent = curn2
    x = True;
    while x == True:
        for i in range(len(path1)):
            if parent == path1[i]:
                mrca = parent
                x = False
                break
        parent = parent.parent
    return mrca
    
def getMRCATraverseFromPath(path1, curn2):
    mrca   = None
    parent = curn2
    x      = True;
    while x == True:
        for i in range(len(path1)):
            if parent == path1[i]:
                mrca = parent
                x    = False
                break
        parent = parent.parent
    return mrca