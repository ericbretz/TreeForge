import re
from typing import Any, Dict, List, Optional

class PrintOut:
    def __init__(self, level, hc, bc, log_file: Optional[str] = None):
        self.level = level
        self.hc = hc or ''
        self.bc = bc or ''
        self.nocolor = False
        self.log_file = log_file

        if self.log_file:
            self.log_handle = open(self.log_file, 'w', encoding='utf-8', buffering=1)
        else:
            self.log_handle = None
        self._last_progress_msg = None
        self._progress_started = False

        self._init_styles()

    def _init_styles(self):
        _rev   = '\033[7m' if self.nocolor else None
        c      = _rev if self.nocolor else self.bc
        info_c = _rev if self.nocolor else '\033[43m'
        warn_c = _rev if self.nocolor else '\033[48;5;208m'
        prog_c = _rev if self.nocolor else '\033[47m'
        err_c  = _rev if self.nocolor else '\033[41m'
        dbg_c  = _rev if self.nocolor else '\033[44m'
        succ_c = _rev if self.nocolor else '\033[42m'

        self._reset = '\033[0m'
        self.styles = {
            'title':    {'color': c,   'width': 80, 'char': ' '},
            'subtitle': {'color': c,   'width': 80, 'char': ' '},
            'metric':   {'color': '\033[107m',   'width': 80, 'char': ' '},
            'info':     {'color': info_c, 'width': 80, 'char': ' '},
            'warning':  {'color': warn_c, 'width': 80, 'char': ' '},
            'progress': {'color': prog_c, 'width': 80, 'char': ' '},
            'error':    {'color': err_c,  'width': 80, 'char': ' '},
            'debug':    {'color': dbg_c,  'width': 80, 'char': ' '},
            'success':  {'color': succ_c, 'width': 80, 'char': ' '},
        }
 
        self.key_translate = {
            'time'                           : 'Total Time',
            'dir'                            : 'Directory',
            # BLAST
            'contig_count'                   : 'Contig Count',
            'raw_blast'                      : 'Raw BLAST',
            'concatenated'                   : 'Concatenated',
            # MCL
            'mcl_count'                      : 'MCL Count',
            'mcl_out'                        : 'MCL',
            'ident_out'                      : 'Ident',
            'blast_mcl_out'                  : 'BLAST to MCL',
            # MAFFT
            'aln_files'                      : 'Alignment',
            'cln_files'                      : 'Cleaned',
            'iqtree_files'                   : 'IQ-TREE Files',
            # Tree
            'trim_tips'                      : 'Trim Tips',
            'mask_tips'                      : 'Mask Tips',
            'cut_branches'                   : 'Cut Branches',
            'written_files'                  : 'Written Files',
            'total_sequences'                : 'Total Sequences',
            'files_written'                  : 'Files Written',
            'trees_processed'                : 'Trees Processed',
            'trees_skipped'                  : 'Trees Skipped',
            'monophyletic_tips_masked'       : 'MonoTips Masked',
            'paraphyletic_tips_masked'       : 'ParaTips Masked',
            'clusters_processed'             : 'Clstr Processed',
            'clusters_skipped'               : 'Clstr Skipped',
            'file_read_errors'               : 'Read Errors',
            'tree_parse_errors'              : 'Parse Errors',
            'invalid_tree_errors'            : 'Invalid Tree',
            'cln_read_errors'                : 'CLN Read Errors',
            'duplicate_cluster_ids'          : 'Dup Clstr IDs',
            'missing_cln_files'              : 'Missing CLN',
            'cutting_threshold'              : 'Cut Threshold',
            'total_trees'                    : 'Total Trees',
            'processed_count'                : 'Processed',
            'total_counts'                   : 'Total',
            'large_subtrees'                 : 'Large Subtrees',
            'small_subtrees'                 : 'Small Subtrees',
            'subtree_sizes'                  : 'Subtree Sizes',
            'insufficient_taxa'              : 'Low Taxa',
            'no_valid_subtrees'              : 'No Valid Sub',
            'no_branches_cut'                : 'No Branches Cut',
            'written_files'                  : 'Written Files',
            'total_sequences'                : 'Total Sequences',
            'file_details'                   : 'File Details',
            'total_tree_files'               : 'Total Files',
            'total_cln_files'                : 'Total CLN',
            'trees_processed'                : 'Trees Processed',
            'trees_skipped'                  : 'Trees Skipped',
            'clusters_processed'             : 'Clstr Processed',
            'clusters_skipped'               : 'Clstr Skipped',
            'monophyletic_tips_masked'       : 'MonoTips Masked',
            'paraphyletic_tips_masked'       : 'ParaTips Masked',
            'mask_paraphyletic'              : 'Mask Para',
            'tree_ending'                    : 'Tree Ending',
            'total_taxa'                     : 'Total Taxa',
            'total_tips'                     : 'Total Tips',
            'total_subtrees'                 : 'Total Subtrees',
            'relative_cutoff'                : 'Relative Cutoff',
            'absolute_cutoff'                : 'Absolute Cutoff',
            'trimmed'                        : 'Trimmed Dir',
            'modified_files_count'           : 'Modified Files',
            'files_moved_to_prune'           : 'Moved to Prune',
            'error_count'                    : 'Error Count',
            'skipped_count'                  : 'Skipped Count',
            # Prune
            'ortho1to1'                      : 'Orthologs',
            'ortho1to1_dir'                  : 'Ortho Directory',
            'ortho1to1_count'                : 'Ortho Count',
            'ortho_mi_count'                 : 'Ortho MI Count',
            'pruned_trees'                   : 'Pruned Trees',
            'skipped_trees'                  : 'Skipped Trees',
            # Prank
            'fasta_from_tree'                : 'FAs from Trees',
            'prank_trees_generated'          : 'Trees Generated',
            'alignments'                     : 'Alignments',
            'cleaned_alignments'             : 'Cleaned Aln',
            'trees'                          : 'Trees',
            # Super
            'super'                          : 'Super',
            'super_dir'                      : 'Super Directory',
            'super_count'                    : 'Super Count',
            'super_time'                     : 'Super Time',
            'concat_matrix'                  : 'Concat Matrix',
            'partition_file'                 : 'Partition File',
            'total_bp'                       : 'Total BP',
            # Astral
            'concat_tree'                    : 'Concat Tree',
            'gene_trees'                     : 'Gene Trees',
            'finaltree'                      : 'Final Tree',
            'num_trees'                      : '# Trees',
            # Contig
            'file'                           : 'File',
            'nucleotide_sequences'           : 'Nucleotide Seq',
            'protein_sequences'              : 'Protein Seq',
            'translated_sequences'           : 'Translated Seq',
            'files_processed'                : 'Files Processed',
            'name_mappings'                  : 'Name Mappings',
            # BUSCO
            'busco_count'                    : 'BUSCO Count',
            'busco_output_dir'               : 'BUSCO Output',
            'busco_files_processed'          : 'Files Processed',
            'species_name'                   : 'Species',
            'diamond_db_status'              : 'Database',
            'diamond_search_status'          : 'Search',
            'busco_files_for_pipeline'       : 'BUSCO Files',
            'busco_files_for_blast'          : 'BUSCO for BLAST',
            # Gene Clustering
            'gene_count'                     : 'Gene Count',
            'gene_files'                     : 'Gene Files',
            'gene_clusters_created'          : 'Gene Clusters',
            'gene_cluster_files_for_pipeline': 'GeneClstr Files',
            'blast_hits_parsed'              : 'BLAST Hits',
            'minimum_taxa'                   : 'Min Taxa',
            'unique_genes_found'             : 'Unique Genes',
            'genes_selected'                 : 'Genes Selected',
            'gene_fasta_files_written'       : 'Gene FASTA',
            'gene_trees_found'               : 'Gene Trees',
            # HCluster
            'clustering_groups'              : 'Cluster Groups',
            'cluster_count'                  : 'Cluster Count',
            'clustered_sequences'            : 'Cluster Seqs',
            'nodes_fasta_dir'                : 'Nodes Dir',
            'clustered_dir'                  : 'Cluster Dir',
            'unclustered_dir'                : 'Uncluster Dir',
            'group_name'                     : 'Group',
            'group_species'                  : 'Species',
            'files_matched'                  : 'Files Matched',
            'clusters_written'               : 'Clusters Written',
            'annotated_internal_nodes'       : 'Annotated Nodes',
            'guide_tree_file'                : 'Guide Tree',
            'gene_trees_count'               : 'Gene Trees',
            'hcluster_guide_tree_created'    : 'Guide Tree Path',
            'hcluster_groups_processing'     : 'Groups',
            'hcluster_levels'                : 'Levels',
            'hcluster_using_busco_files'     : 'Using BUSCO Files',
            'guide_tree_clades'              : 'Tree Clades',
            'guide_tree_depth'               : 'Tree Depth',
            
        }

    def _strip_ansi(self, text: str) -> str:
        """Remove ANSI color codes from text"""
        return re.sub(r'\033\[[0-9;]+m', '', text)
    
    def _flush_pending_progress(self) -> None:
        """Write the last progress message if there was one"""
        if self._last_progress_msg and self.log_handle:
            clean_text = self._strip_ansi(self._last_progress_msg)
            self.log_handle.write(f"PROGRESS END: {clean_text}\n")
            self.log_handle.flush()
            self._last_progress_msg = None
            self._progress_started = False
    
    def _write_to_log(self, text: str) -> None:
        """Write plain text (no ANSI codes) to log file"""
        if self.log_handle:
            self._flush_pending_progress()
            clean_text = self._strip_ansi(text)
            self.log_handle.write(clean_text + '\n')
            self.log_handle.flush()
    
    def close(self) -> None:
        """Close the log file if open"""
        if self.log_handle:
            self._flush_pending_progress()
            self.log_handle.close()
            self.log_handle = None
    
    def __del__(self):
        """Ensure log file is closed on cleanup"""
        self.close()

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
    
    def fmt_str(self, string: str, value: bool = False, max_len: int = None) -> str:
        if value is False and len(string) >= 20:
            return '...' + string[-17:]
        elif value is True and len(string) >= 56:
            return '...' + string[-53:]
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
        elif style == 'warning':
            self.p_warning(content)
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
        # print(f"{color}{string}{'\033[0m'}") #3.12
        print(str(color) + string + self._reset)
        self._write_to_log(string)

    def p_subtitle(self, subtitle: str) -> None:
        color    = self.styles['subtitle']['color']
        width    = self.styles['subtitle']['width']
        char     = self.styles['subtitle']['char']
        string = f'{subtitle:^{width}}'
        print(str(color) + string + self._reset)
        self._write_to_log(string)

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
                    output = f"{metric_str}{self._reset} {value_str}"
                    print(output)
                    self._write_to_log(f"{metric_str} {value_str}")
                else:
                    output = f"{color}{str(line):^{width}}{self._reset}"
                    print(output)
                    self._write_to_log(f"{str(line):^{width}}")
        elif isinstance(metric, dict):
            for key, value in metric.items():
                if key not in self.key_translate:
                    continue
                metric_str = self.fmt_str(f'    {self.fmt_key(key)}', value=False)
                metric_str = f'{metric_str:<20}│'
                value_str  = self.fmt_str(str(value), value=True)
                value_str  = f'{value_str:<{width - 21}}'
                output = f"{metric_str}{self._reset} {value_str}"
                print(output)
                self._write_to_log(f"{metric_str} {value_str}")
        else:
            metric_str = self.fmt_str(self.fmt_key(str(metric)), value=None)
            color = self.styles['metric']['color']
            output = f"{color}{metric_str:^{width}}{self._reset}"
            print(output)
            self._write_to_log(f"{metric_str:^{width}}")

    def p_info(self, info: str) -> None:
        color    = self.styles['info']['color']
        width    = self.styles['info']['width']
        char     = self.styles['info']['char']
        title    = f'{"INFO":<20}'
        info_str = self.fmt_str(info, value=True)
        string   = f'{info_str:<{width - 21}}'
        print(f"{color}{title}{self._reset} {string}")
        self._write_to_log(f"{title} {string}")

    def p_warning(self, warning: str) -> None:
        color       = self.styles['warning']['color']
        width       = self.styles['warning']['width']
        char        = self.styles['warning']['char']
        title       = f'{"WARNING":<20}'
        warning_str = self.fmt_str(warning, value=True)
        string      = f'{warning_str:<{width - 21}}'
        print(f"{color}{title}{self._reset} {string}")
        self._write_to_log(f"{title} {string}")

    def p_error(self, error: str) -> None:
        color     = self.styles['error']['color']
        width     = self.styles['error']['width']
        char      = self.styles['error']['char']
        title     = f'{"ERROR":<20}'
        error_str = self.fmt_str(error, value=True)
        string    = f'{error_str:<{width - 21}}'
        print(f"{color}{title}{self._reset} {string}")
        self._write_to_log(f"{title} {string}")

    def p_debug(self, debug: str) -> None:
        color     = self.styles['debug']['color']
        width     = self.styles['debug']['width']
        char      = self.styles['debug']['char']
        title     = f'{"DEBUG":<20}'
        debug_str = self.fmt_str(debug, value=True)
        string    = f'{debug_str:<{width - 21}}'
        print(f"{color}{title}{self._reset} {string}")
        self._write_to_log(f"{title} {string}")

    def p_success(self, success: str) -> None:
        color       = self.styles['success']['color']
        width       = self.styles['success']['width']
        char        = self.styles['success']['char']
        title       = f'{"SUCCESS":<20}'
        success_str = self.fmt_str(success, value=True)
        string      = f'{success_str:<{width - 21}}'
        print(f"{color}{title}{self._reset} {string}")
        self._write_to_log(f"{title} {string}")

    def p_progress(self, progress: str) -> None:
        color = self.styles['progress']['color']
        width = self.styles['progress']['width']
        print(f"{color}{progress:^{width}}{self._reset}", end='\r', flush=True)
        
        if self.log_handle:
            if not self._progress_started:
                clean_text = self._strip_ansi(progress)
                self.log_handle.write(f"PROGRESS START: {clean_text}\n")
                self.log_handle.flush()
                self._progress_started = True
            self._last_progress_msg = progress

    def set_nocolor(self, nocolor: bool) -> None:
        self.nocolor = nocolor
        if nocolor:
            self.hc = ''
            self.bc = ''
        self._init_styles()

    def p_final(self, files: dict) -> None:
        r = self._reset
        def draw_box(lines, hcolor):
            width  = 78
            top    = f'{hcolor}╭{"─" * width}╮{r}'
            bottom = f'{hcolor}╰{"─" * width}╯{r}'
            print(top)
            for line in lines:
                print(f'{hcolor}│ {r}{line:<{width-1}}{hcolor}│{r}')
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
        
        if self.log_handle:
            width = 78
            self._write_to_log(f"+{'-' * width}+")
            for line in out_list:
                clean_line = self._strip_ansi(line)
                self._write_to_log(f"| {clean_line:<{width-1}}|")
            self._write_to_log(f"+{'-' * width}+")

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
