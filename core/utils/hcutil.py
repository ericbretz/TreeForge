import csv
import io
import re
import shutil
import subprocess
from collections     import OrderedDict, defaultdict
from pathlib         import Path
from tempfile        import mkstemp, mkdtemp
from typing          import Union, Optional, List, Tuple, Dict, Set, Any, Iterable
from Bio.Seq         import Seq as BioSeq
from Bio             import SeqIO

Node = Dict[str, Optional["Node"]]

RES_COLS_VSCH: List[str] = ["query", "target"]

KEY_COL_VSCH : str       = RES_COLS_VSCH[1]
VAL_COL_VSCH : str       = RES_COLS_VSCH[0]


def _truncate_for_terminal(message: str, width: int = 80) -> str:
    if len(message) <= width:
        return message
    return '...' + message[-(width - 3):]

def run(cmd: List[str], in_txt: Optional[str] = None, capture: bool = True, 
        cwd: Optional[str] = None, do_not_raise: bool = False, text: bool = True) -> subprocess.CompletedProcess:
    kwargs = {'input': in_txt, 'cwd': cwd, 'universal_newlines': text}
    if capture:
        kwargs['stdout'] = subprocess.PIPE
        kwargs['stderr'] = subprocess.PIPE
    out = subprocess.run(cmd, **kwargs)
    if not do_not_raise and out.returncode > 0:
        raise Exception(out.stderr)
    return out

def list_of_files_at_path(path: str, regexp: Optional[Union[str, bytes]] = None) -> Tuple[Optional[List[str]], Optional[str]]:
    path_obj = Path(path).resolve()
    
    if path_obj.is_file():
        return None, f'Item at path "{path}" is not a directory.'
    
    try:
        files = [str(path_obj / x) for x in path_obj.iterdir() if x.is_file()]
        files.sort()
    except FileNotFoundError:
        return None, f'Directory "{path}" does not exist.'
    
    if regexp is not None:
        regexp_new = f'.*{regexp}.*'
        pattern = re.compile(regexp_new)
        files = [f for f in files if pattern.search(Path(f).name)]
    
    return files, None

def combine_text_files(paths: Iterable[str], out_path: str) -> None:
    with open(out_path, 'w') as out_file:
        for p in paths:
            with open(p, 'r') as in_file:
                shutil.copyfileobj(in_file, out_file)

def read_fasta(f: Union[str, io.IOBase, io.StringIO]) -> List[Dict[str, Any]]:
    if isinstance(f, str):
        bio_records = list(SeqIO.parse(f, "fasta"))
    elif isinstance(f, (io.IOBase, io.StringIO)):
        bio_records = list(SeqIO.parse(f, "fasta"))
    else:
        raise ValueError("f must be a string (file path) or file-like object")

    return_object = []
    for bio_record in bio_records:
        seq        = str(bio_record.seq).upper()
        defn       = bio_record.id
        seq_record = {'definition': defn, 'seq': seq}
        return_object.append(seq_record)

    return return_object

def dict_to_fasta(d: Dict[str, Union[str, BioSeq]], max_line_len: Optional[int] = None) -> str:
    if not d:
        return ''

    def split_seq(seq: Union[str, BioSeq]) -> str:
        seq_str = str(seq)
        if max_line_len is None:
            return seq_str
        
        lines = []
        for i in range(0, len(seq_str), max_line_len):
            lines.append(seq_str[i:i + max_line_len])
        return '\n'.join(lines)

    fasta_lines = []
    for k, v in d.items():
        fasta_lines.append(f">{k}\n{split_seq(v)}")
    
    return '\n'.join(fasta_lines) + '\n'

def seq_records_to_dict(records: List[Dict[str, Any]], prepend_acc: bool = False, prepend_org: bool = False) -> Dict[str, str]:
    d = {}
    for rec in records:
        if isinstance(rec, dict):
            dfn = rec.get('definition', '')
            if prepend_org and rec.get('organism') is not None:
                dfn = rec['organism'].replace(' ', '_') + ' ' + dfn
            if prepend_acc and rec.get('accession_version') is not None:
                dfn = rec['accession_version'] + ' ' + dfn
            d[dfn] = rec.get('seq', '')
    return d

def write_fasta(data: Union[Dict[str, Union[str, BioSeq]], List[Dict[str, Any]]], f: Union[str, io.IOBase, io.StringIO], 
                max_line_len: Optional[int] = None, prepend_acc: bool = False, prepend_org: bool = False) -> None:
    if isinstance(f, (io.IOBase, io.StringIO)):
        if isinstance(data, (dict, OrderedDict)):
            text = dict_to_fasta(data, max_line_len)
            f.write(text)
        elif isinstance(data, (list, tuple)):
            text = seq_records_to_fasta(data, max_line_len, prepend_acc, prepend_org)
            f.write(text)
    else:
        if not isinstance(f, str):
            raise ValueError("f must be a string or file-like object")
        with open(f, 'w') as file_handle:
            if isinstance(data, (dict, OrderedDict)):
                text = dict_to_fasta(data, max_line_len)
                file_handle.write(text)
            elif isinstance(data, (list, tuple)):
                text = seq_records_to_fasta(data, max_line_len, prepend_acc, prepend_org)
                file_handle.write(text)

def seq_records_to_fasta(records: List[Dict[str, Any]], max_line_len: Optional[int] = None, 
                        prepend_acc: bool = False, prepend_org: bool = False) -> str:
    d = seq_records_to_dict(records, prepend_acc, prepend_org)
    return dict_to_fasta(d, max_line_len)

def parse_blast_results_file_raw(blast_results_file: str, col_names: List[str]) -> List[Dict[str, str]]:
    with Path(blast_results_file).open('r') as f:
        blast_results = f.read()
    return parse_blast_results_raw(blast_results, col_names)

def parse_blast_results_raw(blast_results: str, col_names: List[str]) -> List[Dict[str, str]]:
    rdr = csv.DictReader(blast_results.splitlines(), fieldnames=col_names, delimiter='\t')
    return list(rdr)

def noder(start: int = 0) -> Iterable[str]:
    """Generate unique node IDs starting from a given number."""
    i = start
    while True:
        yield f"N{i}"
        i += 1

def parse_newick(s: str) -> Node:
    if not s:
        raise ValueError("Empty Newick string")

    newick = s.strip()
    if not newick.endswith(";"):
        newick += ";"

    # Remove all whitespace
    newick = re.sub(r"\s+", "", newick)
    # Remove branch lengths
    newick = re.sub(r":[\d\.eE\-\+]+", "", newick)
    # Remove support values in square brackets
    newick = re.sub(r"\[[^\]]*\]", "", newick)
    # Remove numeric support values immediately following a closing parenthesis
    newick = re.sub(r"\)(\d+(\.\d+)?)", ")", newick)
    # Drop trailing semicolon
    newick = newick.rstrip(";")

    if not newick:
        raise ValueError("Newick string is empty after preprocessing")

    max_n_number = -1
    i = 0
    while i < len(newick):
        if newick[i] in "(),":
            i += 1
            continue
        label_start = i
        while i < len(newick) and newick[i] not in ",()":
            i += 1
        label = newick[label_start:i]
        if label.startswith('N') and label[1:].isdigit():
            max_n_number = max(max_n_number, int(label[1:]))
    
    n_id = noder(start=max_n_number + 1)
    idx = 0

    def parse_subtree() -> Tuple[str, Optional[Node]]:
        nonlocal idx
        if idx >= len(newick):
            raise ValueError("Unexpected end of Newick string while parsing subtree")

        c = newick[idx]
        if c == "(":
            idx += 1
            children: List[Tuple[str, Optional[Node]]] = []

            while True:
                name, subtree = parse_subtree()
                children.append((name, subtree))

                if idx >= len(newick):
                    raise ValueError("Unexpected end of Newick string inside internal node")

                if newick[idx] == ",":
                    idx += 1
                    continue
                elif newick[idx] == ")":
                    idx += 1
                    break
                else:
                    raise ValueError(f"Unexpected character '{newick[idx]}' in internal node")

            label_chars: List[str] = []
            while idx < len(newick) and newick[idx] not in ",()":
                label_chars.append(newick[idx])
                idx += 1
            label = "".join(label_chars)

            if label == "":
                label = next(n_id)

            children_dict: Node = {}
            for child_name, child_subtree in children:
                if child_name in children_dict:
                    raise ValueError(f"Duplicate node label '{child_name}' in Newick tree")
                children_dict[child_name] = child_subtree

            return label, children_dict

        label_chars: List[str] = []
        while idx < len(newick) and newick[idx] not in ",()":
            label_chars.append(newick[idx])
            idx += 1

        label = "".join(label_chars)
        if label == "":
            raise ValueError("Encountered empty leaf label in Newick string")

        return label, None

    root_name, root_children = parse_subtree()

    if idx != len(newick):
        raise ValueError(
            f"Unexpected trailing content in Newick string at position {idx}: {newick[idx:]}"
        )

    if root_children is None:
        return {root_name: None}

    return {root_name: root_children}

def parse_newick_file(path: str) -> Node:
    with open(path, "r", encoding="utf8") as f:
        newick: str = f.readline().strip()
        if not newick:
            raise ValueError("Empty Newick string in file")
        try:
            return parse_newick(newick)
        except Exception as e:
            raise ValueError(f"Failed to parse Newick string: {e}. String preview: {newick[:100]}...")

def parse_nodes_dict_to_steps(
    nodes    : Node,
    node_name: Optional[str] = None,
    steps    : Optional[Dict[str, List[str]]] = None,
) -> Tuple[Dict[str, List[str]], Set[str]]:
    if steps is None:
        steps = OrderedDict()
    node_set: Set[str] = set()
    for k, v in nodes.items():
        if v is not None:
            recurse: Set[str]
            steps, recurse = parse_nodes_dict_to_steps(v, k, steps)
            node_set.update(recurse)
        else:
            node_set.add(k)

    if node_name is not None:
        steps[node_name] = sorted(list(node_set))
        node_set = {node_name}

    return (steps, node_set)

def calculate_node_levels(
    nodes    : Node,
    node_name: Optional[str] = None,
    levels   : Optional[Dict[str, int]] = None,
    current_level: int = 0,
) -> Dict[str, int]:
    if levels is None:
        levels = {}
    
    for k, v in nodes.items():
        if v is not None:
            calculate_node_levels(v, k, levels, current_level)
            child_levels = [levels.get(child, 0) for child in v.keys()]
            levels[k] = max(child_levels) + 1 if child_levels else current_level
        else:
            levels[k] = 0
    
    if node_name is not None and node_name not in levels:
        child_levels = [levels.get(child, 0) for child in nodes.keys()]
        levels[node_name] = max(child_levels) + 1 if child_levels else current_level
    
    return levels

def parse_results_file(
    results_file: str,
    col_names   : Iterable[str],
    key_col_name: str,
    val_col_name: str,
) -> Dict[str, Set[str]]:
    parsed_raw: List[Dict[str, str]] = parse_blast_results_file_raw(
        results_file, col_names
    )
    parsed: Dict[str, Set[str]] = defaultdict(set)
    for rec in parsed_raw:
        parsed[rec[key_col_name]].add(rec[val_col_name])

    return parsed

def combine_results(mult_results: List[Dict[str, Set[str]]]) -> Dict[str, Set[str]]:
    combined: Dict[str, Set[str]] = defaultdict(set)
    for results in mult_results:
        for k, v in results.items():
            combined[k].update(v)

    tmp: Dict[str, Set[str]] = combined.copy()
    for k, vals in tmp.items():
        for v in list(vals):
            if v in tmp:
                del combined[v]
                combined[k].update(tmp[v])

    return combined

def make_cluster_sets(results: Dict[str, Set[str]]) -> List[Set[str]]:
    return [v | {k} for k, v in results.items()]

def run_cluster_fast(
    ident                   : float,
    in_file                 : str,
    out_file_centroids      : Optional[str] = None,
    out_file_user           : Optional[str] = None,
    user_fields             : Optional[Iterable[str]] = None,
    out_file_prefix_clusters: Optional[str] = None,
    clusterout_id           : bool = False,
    iddef                   : Optional[int] = None,
    strand                  : str = "both",
    threads                 : Optional[int] = None,
) -> None:
    cmd = ["vsearch", "--cluster_fast", in_file]

    cmd += ["--fasta_width", "0"]
    cmd += ["--id", str(ident)]
    cmd += ["--strand", strand]

    if threads is not None:
        cmd += ["--threads", str(threads)]

    if clusterout_id:
        cmd.append("--clusterout_id")

    if iddef is not None and 0 <= iddef <= 4:
        cmd += ["--iddef", str(iddef)]

    if out_file_prefix_clusters is not None:
        cmd += ["--clusters", out_file_prefix_clusters]

    if out_file_centroids is not None:
        cmd += ["--centroids", out_file_centroids]

    if out_file_user is not None and user_fields is not None:
        cmd += ["--userout", out_file_user]
        cmd += ["--userfields", "+".join(user_fields)]

    run(cmd=cmd, do_not_raise=False)

def run_sub_hc_vsearch(
    fastas          : List[str],
    out_dir         : str,
    ident           : float,
    iddef           : int,
    grp_name        : str,
    threads         : Optional[int] = None,
    precombined_file: Optional[str] = None,
) -> Dict[str, Set[str]]:
    results_id  : str = grp_name
    results_dir : str = out_dir
    results_fp  : str = str(Path(results_dir) / f"{results_id}.vsearch.tsv")
    centroids_fp: str = str(Path(results_dir) / f"{results_id}.centroids.fasta")

    if not Path(results_fp).exists():
        if precombined_file and Path(precombined_file).exists():
            input_file = precombined_file
            cleanup_temp = False
        else:
            _   : int
            tmpf: str
            try:
                _, tmpf = mkstemp(text=True)
            except TypeError:
                _, tmpf = mkstemp()
            combine_text_files(fastas, tmpf)
            input_file = tmpf
            cleanup_temp = True
        
        run_cluster_fast(
            iddef              = iddef,
            ident              = ident,
            in_file            = input_file,
            out_file_user      = results_fp,
            user_fields        = RES_COLS_VSCH,
            out_file_centroids = centroids_fp,
            clusterout_id      = False,
            threads            = threads,
        )
        
        if cleanup_temp:
            Path(input_file).unlink()

    results: Dict[str, Set[str]] = parse_results_file(
        results_fp,
        col_names    = RES_COLS_VSCH,
        key_col_name = KEY_COL_VSCH,
        val_col_name = VAL_COL_VSCH,
    )

    return results

def run_sub_hc_mmseqs2(
    fastas          : List[str],                # Input FASTA file paths
    out_dir         : str,                      # Output directory for results
    ident           : float,                    # Minimum sequence identity for clustering
    iddef           : int,                      # Identity definition mode
    grp_name        : str,                      # Group name for output files
    threads         : Optional[int] = None,     # Number of threads to use
    precombined_file: Optional[str] = None,     # Precombined FASTA file
) -> Dict[str, Set[str]]:
    results_id   : str = grp_name
    results_dir  : str = out_dir
    results_fp   : str = str(Path(results_dir) / f"{results_id}.mmseqs2.tsv")
    centroids_fp : str = str(Path(results_dir) / f"{results_id}.centroids.fasta")
    
    if not Path(results_fp).exists():
        if precombined_file and Path(precombined_file).exists():
            input_file = precombined_file
            cleanup_input = False
        else:
            try:
                _, tmpf = mkstemp(text=True, suffix='.fasta')
            except TypeError:
                _, tmpf = mkstemp(suffix='.fasta')
            combine_text_files(fastas, tmpf)
            input_file = tmpf
            cleanup_input = True
        
        tmp_dir = mkdtemp(prefix='mmseqs2_tmp_')
        
        try:
            cluster_prefix = str(Path(tmp_dir) / "clusterRes")
            
            cmd = [
                "mmseqs", "easy-cluster",
                input_file,
                cluster_prefix,
                tmp_dir,
                "--min-seq-id", str(ident),
                "-c", "0.8",
            ]
            
            if iddef <= 2:
                cmd.extend(["--cov-mode", "0"])
            else:
                cmd.extend(["--cov-mode", "1"])
            
            cmd.extend(["--cluster-mode", "0"])
            
            if threads is not None:
                cmd.extend(["--threads", str(threads)])
            
            run(cmd=cmd, do_not_raise=False)

            cluster_file = f"{cluster_prefix}_cluster.tsv"
            results: Dict[str, Set[str]] = defaultdict(set)
            
            with open(cluster_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        representative = parts[0]
                        member = parts[1]
                        if member != representative:
                            results[representative].add(member)
            
            rep_seq_file = f"{cluster_prefix}_rep_seq.fasta"
            if Path(rep_seq_file).exists():
                shutil.copy(rep_seq_file, centroids_fp)
            
            with open(results_fp, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                for rep_id, member_ids in results.items():
                    for member_id in member_ids:
                        writer.writerow([member_id, rep_id])
        
        finally:
            if Path(tmp_dir).exists():
                shutil.rmtree(tmp_dir)
            
            if cleanup_input and Path(input_file).exists():
                Path(input_file).unlink()
    
    results: Dict[str, Set[str]] = parse_results_file(
        results_fp,
        col_names    = RES_COLS_VSCH,
        key_col_name = KEY_COL_VSCH,
        val_col_name = VAL_COL_VSCH,
    )
    
    return results

def node_to_newick(node: Node) -> str:
    if not node:
        return ""
    
    if len(node) == 1:
        root_name, root_children = list(node.items())[0]
        if root_children is None:
            return f"{root_name};"
        else:
            children_strs = []
            for child_name, child_subtree in root_children.items():
                if child_subtree is None:
                    children_strs.append(child_name)
                else:
                    child_node = {child_name: child_subtree}
                    child_newick = node_to_newick(child_node).rstrip(';')
                    children_strs.append(child_newick)
            return f"({','.join(children_strs)}){root_name};"
    else:
        children_strs = []
        for child_name, child_subtree in node.items():
            if child_subtree is None:
                children_strs.append(child_name)
            else:
                child_node = {child_name: child_subtree}
                child_newick = node_to_newick(child_node).rstrip(';')
                children_strs.append(child_newick)
        return f"({','.join(children_strs)});"

def create_node_mapping_text(node: Node, node_name: Optional[str] = None, indent: int = 0) -> str:
    lines = []
    
    if not node:
        return ""
    
    def get_all_leaves(n: Node) -> List[str]:
        leaves = []
        for name, subtree in n.items():
            if subtree is None:
                leaves.append(name)
            else:
                leaves.extend(get_all_leaves(subtree))
        return leaves
    
    def traverse(n: Node, name: Optional[str], depth: int):
        if not n:
            return
        
        for node_label, children in n.items():
            if children is not None:
                leaves = get_all_leaves({node_label: children})
                indent_str = "  " * depth
                lines.append(f"{indent_str}Node: {node_label}")
                lines.append(f"{indent_str}  └─ {', '.join(sorted(leaves))}")
                lines.append("")
                traverse(children, node_label, depth + 1)
    
    traverse(node, node_name, indent)
    
    return "\n".join(lines)