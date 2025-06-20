from ete3 import Tree
from core.treeutils.phylo import Node

def parse(input, ttable=None):
    if isinstance(input, str):
        return Node(Tree(input, format=1))
    else:
        return Node(Tree(input.readline(), format=1))

def tostring(tree):
    return tree._ete_node.write(format=1)
