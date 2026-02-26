import csv
import re
from pathlib import Path
from typing import Dict, Any, List, Optional

_NEWICK_SUPPORT_RE = re.compile(r'\)([\d.]+(?:e[+-]?\d+)?)')


def _support_from_newick(topology: str) -> Optional[float]:
    m = _NEWICK_SUPPORT_RE.search(topology)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


def _detect_format(rows: List[List[str]]) -> str:
    for row in rows:
        if len(row) >= 2:
            try:
                float(row[1])
                return 'multi_column'
            except (ValueError, IndexError):
                pass
    return 'embedded_support'


def parse_verbose_csv(verbose_csv: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        'branches'            : [],
        'num_branches'        : 0,
        'mean_support'        : None,
        'median_support'      : None,
        'min_support'         : None,
        'max_support'         : None,
        'branches_gt_0_95'    : 0,
        'branches_lt_0_50'    : 0,
        'pct_branches_gt_0_95': None,
    }

    path = Path(verbose_csv)
    if not path.exists():
        return result

    try:
        with open(path, 'r', newline='') as f:
            all_rows = list(csv.reader(f))
    except (OSError, IOError):
        return result

    if not all_rows:
        return result

    fmt = _detect_format(all_rows)

    branch_means: List[float] = []
    branches: List[Dict[str, Any]] = []

    if fmt == 'multi_column':
        for row in all_rows:
            if len(row) < 2:
                continue
            clade = row[0].strip().strip('"')
            try:
                posteriors = [float(v) for v in row[1:] if v.strip()]
            except ValueError:
                continue
            if not posteriors:
                continue

            mean_pp   = sum(posteriors) / len(posteriors)
            sorted_pp = sorted(posteriors)
            n         = len(sorted_pp)
            median_pp = (sorted_pp[n // 2] if n % 2 == 1
                         else (sorted_pp[n // 2 - 1] + sorted_pp[n // 2]) / 2)

            branches.append({
                'branch_topology'       : clade,
                'num_gene_trees'        : n,
                'mean_local_posterior'  : round(mean_pp, 6),
                'median_local_posterior': round(median_pp, 6),
                'min_local_posterior'   : round(sorted_pp[0], 6),
                'max_local_posterior'   : round(sorted_pp[-1], 6),
            })
            branch_means.append(mean_pp)

    else:
        for row in all_rows:
            clade = row[0].strip().strip('"') if row else ''
            if not clade:
                continue
            support = _support_from_newick(clade)
            if support is None:
                continue

            branches.append({
                'branch_topology'       : clade,
                'num_gene_trees'        : '',
                'mean_local_posterior'  : support,
                'median_local_posterior': support,
                'min_local_posterior'   : support,
                'max_local_posterior'   : support,
            })
            branch_means.append(support)

    if not branch_means:
        return result

    n_branches    = len(branch_means)
    gt_095        = sum(1 for v in branch_means if v > 0.95)
    lt_050        = sum(1 for v in branch_means if v < 0.50)
    sorted_means  = sorted(branch_means)
    overall_median = (sorted_means[n_branches // 2] if n_branches % 2 == 1
                      else (sorted_means[n_branches // 2 - 1] + sorted_means[n_branches // 2]) / 2)

    result.update({
        'branches'            : branches,
        'num_branches'        : n_branches,
        'mean_support'        : round(sum(branch_means) / n_branches, 6),
        'median_support'      : round(overall_median, 6),
        'min_support'         : round(min(branch_means), 6),
        'max_support'         : round(max(branch_means), 6),
        'branches_gt_0_95'    : gt_095,
        'branches_lt_0_50'    : lt_050,
        'pct_branches_gt_0_95': round(gt_095 / n_branches * 100, 2),
    })
    return result
