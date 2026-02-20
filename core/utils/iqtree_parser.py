import re
from pathlib import Path
from typing import Dict, Any, Optional


def parse_iqtree_report(iqtree_file: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        'gene_id'                   : Path(iqtree_file).stem,
        'num_sequences'             : None,
        'num_sites'                 : None,
        'num_constant_sites'        : None,
        'pct_constant'              : None,
        'num_parsimony_informative' : None,
        'pct_parsimony_informative' : None,
        'num_distinct_patterns'     : None,
        'substitution_model'        : None,
        'log_likelihood'            : None,
        'AIC'                       : None,
        'AICc'                      : None,
        'BIC'                       : None,
        'gamma_alpha'               : None,
        'total_tree_length'         : None,
        'internal_branch_length'    : None,
        'pct_internal_length'       : None,
        'cpu_time_s'                : None,
        'wall_time_s'               : None,
    }

    try:
        with open(iqtree_file, 'r') as f:
            text = f.read()
    except (OSError, IOError):
        return result

    m = re.search(r'Input data:\s+(\d+)\s+sequences?\s+with\s+(\d+)\s+nucleotide sites?', text)
    if m:
        result['num_sequences'] = int(m.group(1))
        result['num_sites']     = int(m.group(2))

    m = re.search(r'Number of constant sites:\s+(\d+)\s+\(=\s*([\d.]+)%', text)
    if m:
        result['num_constant_sites'] = int(m.group(1))
        result['pct_constant']       = float(m.group(2))

    m = re.search(r'Number of parsimony informative sites:\s+(\d+)', text)
    if m:
        result['num_parsimony_informative'] = int(m.group(1))
        if result['num_sites']:
            result['pct_parsimony_informative'] = round(
                result['num_parsimony_informative'] / result['num_sites'] * 100, 4)

    m = re.search(r'Number of distinct site patterns:\s+(\d+)', text)
    if m:
        result['num_distinct_patterns'] = int(m.group(1))

    m = re.search(r'Model of substitution:\s+(\S+)', text)
    if m:
        result['substitution_model'] = m.group(1)

    m = re.search(r'Gamma shape alpha:\s+([\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['gamma_alpha'] = float(m.group(1))

    m = re.search(r'Log-likelihood of the tree:\s+([-\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['log_likelihood'] = float(m.group(1))

    m = re.search(r'Akaike information criterion \(AIC\) score:\s+([\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['AIC'] = float(m.group(1))

    m = re.search(r'Corrected Akaike information criterion \(AICc\) score:\s+([\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['AICc'] = float(m.group(1))

    m = re.search(r'Bayesian information criterion \(BIC\) score:\s+([\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['BIC'] = float(m.group(1))

    m = re.search(r'Total tree length \(sum of branch lengths\):\s+([\d.]+(?:e[+-]?\d+)?)', text)
    if m:
        result['total_tree_length'] = float(m.group(1))

    m = re.search(r'Sum of internal branch lengths:\s+([\d.]+(?:e[+-]?\d+)?)\s+\(([\d.]+)%', text)
    if m:
        result['internal_branch_length'] = float(m.group(1))
        result['pct_internal_length']    = float(m.group(2))

    m = re.search(r'Total CPU time used:\s+([\d.]+)\s+seconds', text)
    if m:
        result['cpu_time_s'] = float(m.group(1))

    m = re.search(r'Total wall-clock time used:\s+([\d.]+)\s+seconds', text)
    if m:
        result['wall_time_s'] = float(m.group(1))

    return result


def parse_iqtree_directory(directory: str, suffix: str = '.iqtree') -> list:
    dirpath = Path(directory)
    if not dirpath.is_dir():
        return []
    files = sorted(dirpath.glob(f'*{suffix}'))
    return [parse_iqtree_report(str(f)) for f in files]
