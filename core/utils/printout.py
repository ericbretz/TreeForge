import os
import sys
from typing import Any, Dict, List, Union

class PrintOut:
    def __init__(self, level, hc, bc):
        self.level = level
        self.hc = hc  # Highlight color
        self.bc = bc  # Background color

        # Formatting styles
        self.styles = {
            'title': {
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'subtitle': {
                # 'color'   : '\033[30;47m',
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'metric': {
                'color'   : self.bc,
                'width'   : 80,
                'char'    : ' '
            },
            'info': {
                'color'   : '\033[43m',
                'width'   : 80,
                'char'    : ' '
            },
            'progress': {
                'color'   : '\033[47m',
                'width'   : 80,
                'char'    : ' '
            },
            'error': {
                'color'   : '\033[41m',
                'width'   : 80,
                'char'    : ' '
            },
            'debug': {
                'color'   : '\033[44m',
                'width'   : 80,
                'char'    : ' '
            },
            'success': {    
                'color'   : '\033[42m',
                'width'   : 80,
                'char'    : ' '
            },
        }
 
        self.key_translate = {
            'time'        : 'Total Time',
            'dir'         : 'Directory',
            # BLAST
            'contig_count': 'Contig Count',
            'raw_blast'   : 'Raw BLAST',
            'concatenated': 'Concatenated',
            # MCL
            'mcl_count'   : 'MCL Count',
            'mcl_out'     : 'MCL',
            'ident_out'   : 'Ident',
            'blast_mcl_out': 'BLAST to MCL',
            # MAFFT
            'aln_files'   : 'Alignment',
            'cln_files'   : 'Cleaned',
            'iqtree_files': 'IQ-TREE Files',
            # Tree
            'trim_tips'               : 'Trim Tips',
            'mask_tips'               : 'Mask Tips',
            'cut_branches'            : 'Cut Branches',
            'written_files'           : 'Written Files',
            'total_sequences'         : 'Total Sequences',
            'files_written'           : 'Files Written',
            'trees_processed'         : 'Trees Processed',
            'trees_skipped'           : 'Trees Skipped',
            'monophyletic_tips_masked': 'MonoTips Masked',
            'paraphyletic_tips_masked': 'ParaTips Masked',
            'clusters_processed'      : 'Clstr Processed',
            'clusters_skipped'        : 'Clstr Skipped',
            'file_read_errors'        : 'Read Errors',
            'tree_parse_errors'       : 'Parse Errors',
            'invalid_tree_errors'     : 'Invalid Tree',
            'cln_read_errors'         : 'CLN Read Errors',
            'duplicate_cluster_ids'   : 'Dup Clstr IDs',
            'missing_cln_files'       : 'Missing CLN',
            'cutting_threshold'       : 'Cut Threshold',
            'total_trees'             : 'Total Trees',
            'processed_count'         : 'Processed',
            'total_counts'            : 'Total',
            'large_subtrees'          : 'Large Subtrees',
            'small_subtrees'          : 'Small Subtrees',
            'subtree_sizes'           : 'Subtree Sizes',
            'insufficient_taxa'       : 'Low Taxa',
            'no_valid_subtrees'       : 'No Valid Sub',
            'no_branches_cut'         : 'No Branches Cut',
            'written_files'           : 'Written Files',
            'total_sequences'         : 'Total Sequences',
            'file_details'            : 'File Details',
            'total_tree_files'        : 'Total Files',
            'total_cln_files'         : 'Total CLN',
            'trees_processed'         : 'Trees Processed',
            'trees_skipped'           : 'Trees Skipped',
            'clusters_processed'      : 'Clstr Processed',
            'clusters_skipped'        : 'Clstr Skipped',
            'monophyletic_tips_masked': 'MonoTips Masked',
            'paraphyletic_tips_masked': 'ParaTips Masked',
            'mask_paraphyletic'       : 'Mask Para',
            'tree_ending'             : 'Tree Ending',
            'total_taxa'              : 'Total Taxa',
            'total_tips'              : 'Total Tips',
            'total_subtrees'          : 'Total Subtrees',
            'relative_cutoff'         : 'Relative Cutoff',
            'absolute_cutoff'         : 'Absolute Cutoff',
            'trimmed'                 : 'Trimmed Dir',
            # Prune
            'ortho1to1'               : 'Orthologs',
            'ortho1to1_dir'           : 'Ortho Directory',
            'ortho1to1_count'         : 'Ortho Count',
            # Prank
            'alignments'              : 'Alignments',
            'cleaned_alignments'      : 'Cleaned Aln',
            'trees'                   : 'Trees',
            # Super
            'super'                   : 'Super',
            'super_dir'               : 'Super Directory',
            'super_count'             : 'Super Count',
            'super_time'              : 'Super Time',
            # Astral
            'concat_tree'             : 'Concat Tree',
            'gene_trees'              : 'Gene Trees',
            'finaltree'               : 'Final Tree',
            'num_trees'               : '# Trees',
            # Contig
            'file'                    : 'File',
            'nucleotide_sequences'    : 'Nucleotide Seq',
            'protein_sequences'       : 'Protein Seq',
            'translated_sequences'    : 'Translated Seq',
            'files_processed'         : 'Files Processed',
            'name_mappings'           : 'Name Mappings',
        }

    def fmt_key(self, key: str) -> str:
        if key in self.key_translate:
            return self.key_translate[key]
        else:
            return key
    
    def fmt_lst_cnt(self, lst: List[Any]) -> str:
        return f"{len(lst)}"

    def fmt_dict(self, dct: Dict[Any, Any]) -> List[tuple]:
        dct_list = []
        for k, v in dct.items():
            if isinstance(v, list):
                dct_list.append((k, self.fmt_lst_cnt(v)))
            elif isinstance(v, dict):
                nested_items = self.fmt_dict(v)
                for nk, nv in nested_items:
                    dct_list.append((nk, nv))
            else:
                dct_list.append((k, str(v)))
        return dct_list
    
    def fmt_str(self, string: str, value: bool = False) -> str:
        if value is None:
            return string
        elif len(string) >= 20 and not value:
            return '...' + string[-17:]
        elif len(string) >= 58 and value:
            return '...' + string[-55:]
        else:
            return string

    def check_type(self, content: Any) -> None:
        if isinstance(content, dict):
            content = self.fmt_dict(content)
        elif isinstance(content, list):
            if content and isinstance(content[0], tuple):
                content = content
            else:
                content = self.fmt_lst_cnt(content)
        elif isinstance(content, str):
            content = content
        elif isinstance(content, int):
            content = str(content)
        else:
            raise ValueError(f"Invalid content type: {type(content)}")
        return content
    
    def check_style(self, style: str, content: Any) -> None:
        if style == 'title':
            self.p_title(content)
        elif style == 'subtitle':
            self.p_subtitle(content)
        elif style == 'metric':
            self.p_metric(content)
        elif style == 'info':
            self.p_info(content)
        elif style == 'error':
            self.p_error(content)
        elif style == 'debug':
            self.p_debug(content)
        elif style == 'success':
            self.p_success(content)
        elif style == 'progress':
            self.p_progress(content)
        elif style == 'final':
            self.p_final(content)
        else:
            raise ValueError(f"Invalid style: {style}")
        return content

    def p_title(self, title: str,) -> None:
        color    = self.styles['title']['color']
        width    = self.styles['title']['width']
        char     = self.styles['title']['char']
        string = f'{title:^{width}}'
        print(f"{color}{string}{'\033[0m'}")

    def p_subtitle(self, subtitle: str) -> None:
        color    = self.styles['subtitle']['color']
        width    = self.styles['subtitle']['width']
        char     = self.styles['subtitle']['char']
        string = f'{subtitle:^{width}}'
        print(f"{color}{string}{'\033[0m'}")

    def p_metric(self, metric: str) -> None:
        color    = self.styles['metric']['color']
        width    = self.styles['metric']['width']
        char     = self.styles['metric']['char']
        if isinstance(metric, list):
            for line in metric:
                if isinstance(line, tuple) and len(line) == 2:
                    key, value = line
                    metric_str = self.fmt_str(f'    {self.fmt_key(key)}', value=False)
                    metric_str = f'{metric_str:<20}│'
                    value_str  = self.fmt_str(str(value), value=True)
                    value_str  = f'{value_str:<{width - 21}}'
                    print(f"{metric_str}\033[0m {value_str}")
                else:
                    print(f"{color}{str(line):^{width}}\033[0m")
        elif isinstance(metric, dict):
            for key, value in metric.items():
                if key not in self.key_translate:
                    continue
                metric_str = self.fmt_str(f'    {self.fmt_key(key)}', value=False)
                metric_str = f'{metric_str:<20}│'
                value_str  = self.fmt_str(str(value), value=True)
                value_str  = f'{value_str:<{width - 21}}'
                print(f"{metric_str}\033[0m {value_str}")
        else:
            metric_str = self.fmt_str(self.fmt_key(str(metric)), value=None)
            print(f"\033[47m{metric_str:^{width}}\033[0m")

    def p_info(self, info: str) -> None:
        color    = self.styles['info']['color']
        width    = self.styles['info']['width']
        char     = self.styles['info']['char']
        title    = f'{"INFO":<20}'
        string = f'{info:<{width - 20}}'
        print(f"{color}{title}\033[0m {string}")

    def p_error(self, error: str) -> None:
        color    = self.styles['error']['color']
        width    = self.styles['error']['width']
        char     = self.styles['error']['char']
        title    = f'{"ERROR":<20}'
        string = f'{error:<{width - 20}}'
        print(f"{color}{title}\033[0m {string}")

    def p_debug(self, debug: str) -> None:
        color    = self.styles['debug']['color']
        width    = self.styles['debug']['width']
        char     = self.styles['debug']['char']
        title    = f'{"DEBUG":<20}'
        string = f'{debug:<{width - 20}}'
        print(f"{color}{title}\033[0m {string}")

    def p_success(self, success: str) -> None:
        color    = self.styles['success']['color']
        width    = self.styles['success']['width']
        char     = self.styles['success']['char']
        title    = f'{"SUCCESS":<20}'
        string = f'{success:<{width - 20}}'
        print(f"{color}{title}\033[0m {string}")

    def p_progress(self, progress: str) -> None:
        color = self.styles['progress']['color']
        width = self.styles['progress']['width']
        print(f"\033[47m{progress:^{width}}\033[0m", end='\r', flush=True)

    def p_final(self, files: dict) -> None:
        def draw_box(lines, hcolor):
            width = 78
            top = f'{hcolor}╭{"─" * width}╮\033[0m'
            bottom = f'{hcolor}╰{"─" * width}╯\033[0m'
            print(top)
            for line in lines:
                print(f'{hcolor}│ \033[0m{line:<{width-1}}{hcolor}│\033[0m')
            print(bottom)
        out_list = []
        for dir, file in files:
            if dir != 'Gene Tree Count':
                s = f'{dir:<18}│ {self.fmt_str(str(file), value=True)}'
                out_list.append(s)
            else:
                s = f'{dir:<18}│ {file}'
                out_list.append(s)
        draw_box(out_list, self.hc)

    def printout(self, style: str, content: Any) -> None:
        content = self.check_type(content)
        self.check_style(style, content)

if __name__ == "__main__":
    import time
    test_dict = {
        'a': 1,
        'b': 2,
        'ccccccccccccccccccccccccccccccccasddasasqweqweqasdasda': 3,
        'd': 44444444444444444444444444444444444444234242342344444444444444444444444,
        'e': 5,
        'f': {
            'g': 6,
            'h': 7,
            'i': 888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888,
            'j': [1, 2, 3, 4, 5],
        }
    }
    
    printout = PrintOut(level='info', hc='\033[94m', bc='\033[44m')
    
    printout.printout('title',      'Test Title')
    printout.printout('subtitle',   'Test Subtitle')
    printout.printout('metric',     'Test Metric Group')
    printout.printout('info',       'Test Info')
    printout.printout('error',      'Test Error')
    printout.printout('success',    'Test Success')
    
    printout.printout('title',      'Test Dictionary')
    printout.printout('subtitle',   'Test Dictionary')
    printout.printout('metric',     'Test Dictionary Group')
    printout.printout('metric',      test_dict)
    for i in range(20):
        printout.printout('progress', f"MAFFT progress: {i}/{20} files completed")
        time.sleep(0.2)
