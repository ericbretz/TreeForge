from functools import lru_cache
from ete3 import Tree
from core.treeutils.phylo import Node

@lru_cache(maxsize=128)
def _parse_cached(tree_string: str):
    return Node(Tree(tree_string, format=1))

def parse(input, ttable=None):
    if isinstance(input, str):
        return _parse_cached(input)
    else:
        tree_string = input.readline()
        return _parse_cached(tree_string)

def tostring(tree):
    return tree._ete_node.write(format=1)
