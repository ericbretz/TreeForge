import time
import functools
import os
from datetime              import datetime
from pathlib               import Path
from typing                import Set, Tuple, Any, Callable, Dict, List, Optional
from pathlib               import Path as PathType
from core.utils.printout   import PrintOut
from core.utils.formatters import format_bytes, format_number, get_file_size
from core.utils.constants  import ITERATIVE_STEPS, STEP_METRICS_CONFIG, FUNCTION_PARAMS_MAP


def _truncate_for_terminal(message: str, width: int = 80) -> str:
    if len(message) <= width:
        return message
    return '...' + message[-(width - 3):]

def get_filesystem_state(path: Path, target_dirs: Optional[List[Path]] = None) -> Tuple[Set[Path], Set[Path]]:
    files = set()
    dirs  = set()
    
    if target_dirs is None:
        for root, directories, filenames in os.walk(path):
            root_path = Path(root)
            dirs.update(root_path / d for d in directories)
            files.update(root_path / f for f in filenames)
    else:
        for target_dir in target_dirs:
            if not target_dir.exists():
                continue
            for root, directories, filenames in os.walk(target_dir):
                root_path = Path(root)
                dirs.update(root_path / d for d in directories)
                files.update(root_path / f for f in filenames)
    
    return files, dirs

def get_stage_target_dirs(datahub, stage_name: str) -> List[Path]:
    stage_lower = stage_name.lower()
    target_dirs = []
    
    if stage_lower == 'blast':
        if hasattr(datahub, 'master_dict') and 'blast' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['blast']['dir'])
        if hasattr(datahub, 'hcluster_enabled') and datahub.hcluster_enabled:
            if 'hcluster_blast' in datahub.master_dict:
                target_dirs.append(datahub.master_dict['hcluster_blast']['dir'])
    
    elif stage_lower == 'mcl':
        if 'mcl' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['mcl']['dir'])
        iter_key = f'iter_{datahub.current_iter}'
        if iter_key in datahub.master_dict and 'mafft' in datahub.master_dict[iter_key]:
            target_dirs.append(datahub.master_dict[iter_key]['mafft']['dir'])
    
    elif stage_lower == 'mafft':
        iter_key = f'iter_{datahub.current_iter}'
        if iter_key in datahub.master_dict and 'mafft' in datahub.master_dict[iter_key]:
            target_dirs.append(datahub.master_dict[iter_key]['mafft']['dir'])
    
    elif stage_lower == 'tree':
        iter_key = f'iter_{datahub.current_iter}'
        if iter_key in datahub.master_dict:
            if 'tree' in datahub.master_dict[iter_key]:
                target_dirs.append(datahub.master_dict[iter_key]['tree']['dir'])
                if 'trimmed' in datahub.master_dict[iter_key]['tree']:
                    target_dirs.append(datahub.master_dict[iter_key]['tree']['trimmed'])
        if 'prune' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['prune']['dir'])
    
    elif stage_lower == 'prune':
        if 'prune' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['prune']['dir'])
            if 'ortho1to1' in datahub.master_dict['prune']:
                target_dirs.append(datahub.master_dict['prune']['ortho1to1'])
    
    elif stage_lower == 'prank':
        if 'prank' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['prank']['dir'])
    
    elif stage_lower == 'astral':
        if 'super' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['super']['dir'])
        if 'results' in datahub.master_dict:
            if 'species_trees' in datahub.master_dict['results']:
                target_dirs.append(datahub.master_dict['results']['species_trees']['dir'])
            if 'gene_trees' in datahub.master_dict['results']:
                target_dirs.append(datahub.master_dict['results']['gene_trees']['dir'])
    
    elif stage_lower in ['busco_extract', 'genecluster_for_hcluster', 'hcluster_guide_tree', 'hcluster']:
        if 'hcluster_enabled' in datahub.master_dict:
            target_dirs.append(datahub.master_dict['hcluster_enabled']['dir'])
    
    if hasattr(datahub, 'dir_logs'):
        target_dirs.append(datahub.dir_logs)
    
    return [d for d in target_dirs if d is not None]

def convert_paths(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: convert_paths(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_paths(item) for item in obj]
    return obj

class BaseDecorator:
    """Handles common stuff that both decorators need to do."""
    
    def __init__(self, skip_metrics: bool = False):
        self.skip_metrics = skip_metrics
    
    def update_dict(self, datahub, step_name: str, step_dict: Any, execution_time: float) -> None:
        """Add step results to the master dictionary."""
        step_name_lower = step_name.lower()
        
        if step_name_lower in ITERATIVE_STEPS:
            iter_key = f'iter_{datahub.current_iter}'
            if iter_key not in datahub.master_dict:
                datahub.master_dict[iter_key] = {}
            if step_name_lower not in datahub.master_dict[iter_key]:
                datahub.master_dict[iter_key][step_name_lower] = {}
            
            target_dict = datahub.master_dict[iter_key][step_name_lower]
        else:
            if step_name_lower not in datahub.master_dict:
                datahub.master_dict[step_name_lower] = {}
            target_dict = datahub.master_dict[step_name_lower]
        
        if 'metrics' not in target_dict:
            target_dict['metrics'] = {}
        target_dict['metrics']['time'] = f'{execution_time:.2f} s'
        
        if step_dict is not None:
            if isinstance(step_dict, dict):
                for key, value in step_dict.items():
                    existing_value = target_dict.get(key)
                    if (existing_value is not None and 
                        isinstance(existing_value, dict) and 
                        isinstance(value, dict)):
                        target_dict[key].update(value)
                    else:
                        target_dict[key] = value
            else:
                target_dict['data'] = step_dict
    
    def print_nested_dict(self, datahub, d: Dict, exclude_keys: Optional[Set] = None) -> List[Tuple]:
        """Print nested dict as sorted list, exclude certain keys."""
        if exclude_keys is None:
            exclude_keys = set()
        
        result = []
        for key, value in sorted(d.items()):
            if key in exclude_keys:
                continue
                
            if isinstance(value, dict):
                for nested_key, nested_value in sorted(value.items()):
                    if isinstance(nested_value, (str, int, float, bool)):
                        result.append((nested_key, nested_value))
                    elif isinstance(nested_value, list):
                        result.append((nested_key, len(nested_value)))
            elif isinstance(value, (str, int, float, bool)):
                if key not in datahub.printClass.key_translate:
                    continue
                result.append((key, value))
            elif isinstance(value, list):
                if key not in datahub.printClass.key_translate:
                    continue
                result.append((key, len(value)))
        return result
    
    def print_metrics(self, datahub, step_name: str) -> None:
        """Shows the metrics for this step, if its supposed to."""
        if self.skip_metrics:
            return
            
        step_name_lower = step_name.lower()
        
        if step_name_lower in ITERATIVE_STEPS:
            iter_key = f'iter_{datahub.current_iter}'
            metrics  = datahub.master_dict[iter_key][step_name_lower]
            datahub.printout('metric', metrics)
        else:
            metrics = datahub.master_dict[step_name_lower]
            config = STEP_METRICS_CONFIG.get(step_name_lower, {})
            
            if step_name_lower == 'prune' and 'prune' in metrics:
                metrics = metrics['prune']
            
            filtered_metrics = self.print_nested_dict(
                datahub, metrics, 
                exclude_keys=config.get('exclude_keys', set())
            )
            
            if filtered_metrics:
                datahub.printout('metric', filtered_metrics)
    
    def _get_function_params(self, datahub, func_name: str) -> Dict[str, Any]:
        """Extract relevant parameters for a function from datahub."""
        params = {}
        param_map = FUNCTION_PARAMS_MAP.get(func_name, [])
        
        for attr_name, display_name in param_map:
            if hasattr(datahub, attr_name):
                value = getattr(datahub, attr_name)
                params[display_name] = value
        
        return params
    
    def _get_function_description(self, func: Callable) -> str:
        """Extract first line of docstring as description."""
        if func.__doc__:
            lines = func.__doc__.strip().split('\n')
            return lines[0].strip()
        return "No description available"
    
    def _determine_status(self, result: Any, exception: Optional[Exception]) -> str:
        """Determine execution status based on result and exceptions."""
        if exception:
            return "FAILED"
        
        if result is None:
            return "SUCCESS"
        
        if isinstance(result, dict):
            if result.get('status') == 'skipped':
                return "SKIPPED"
            if result.get('status') == 'failed':
                return "FAILED"
            if 'error' in result or 'errors' in result:
                return "WARNING"
        
        return "SUCCESS"
    
    def _format_metrics(self, func_name: str, result: Any) -> List[str]:
        """Format metrics from function return dictionary."""
        lines = []
        
        if not isinstance(result, dict):
            return lines
        
        func_lower = func_name.lower()
        
        # BLAST
        if func_lower == 'blast':
            if 'contig_count' in result:
                lines.append(f"  Total sequences concatenated: {format_number(result['contig_count'])} sequences")
        
        # MCL
        elif func_lower == 'mcl':
            if 'mcl_count' in result:
                lines.append(f"  Clusters generated: {format_number(result['mcl_count'])} clusters")
        
        # MAFFT
        elif func_lower == 'mafft':
            if 'mafft' in result:
                mafft_data = result['mafft']
                if 'files_processed' in mafft_data:
                    lines.append(f"  Files processed: {format_number(mafft_data['files_processed'])} files")
                if 'aln_files' in mafft_data and isinstance(mafft_data['aln_files'], list):
                    lines.append(f"  Alignments created: {format_number(len(mafft_data['aln_files']))} files")
                if 'cln_files' in mafft_data and isinstance(mafft_data['cln_files'], list):
                    lines.append(f"  Cleaned alignments: {format_number(len(mafft_data['cln_files']))} files")
                if 'iqtree_files' in mafft_data and isinstance(mafft_data['iqtree_files'], list):
                    lines.append(f"  Trees generated: {format_number(len(mafft_data['iqtree_files']))} files")
        
        # Tree
        elif func_lower == 'tree':
            if 'trim_tips' in result:
                trim_data = result['trim_tips']
                if 'trees_processed' in trim_data:
                    lines.append(f"  Trees processed (trim): {format_number(trim_data['trees_processed'])} trees")
                if 'tips_removed' in trim_data:
                    lines.append(f"  Tips trimmed: {format_number(trim_data['tips_removed'])} tips")
            
            if 'mask_tips' in result:
                mask_data = result['mask_tips']
                if 'trees_processed' in mask_data:
                    lines.append(f"  Trees processed (mask): {format_number(mask_data['trees_processed'])} trees")
                if 'tips_masked' in mask_data:
                    lines.append(f"  Tips masked: {format_number(mask_data['tips_masked'])} tips")
            
            if 'cut_branches' in result:
                cut_data = result['cut_branches']
                if 'trees_processed' in cut_data:
                    lines.append(f"  Trees processed (cut): {format_number(cut_data['trees_processed'])} trees")
            
            if 'write_tree' in result:
                write_data = result['write_tree']
                if 'final_files_count' in write_data:
                    lines.append(f"  Final tree files: {format_number(write_data['final_files_count'])} files")
        
        # Prune
        elif func_lower == 'prune':
            if 'prune' in result:
                prune_data = result['prune']
                if 'ortho1to1_count' in prune_data:
                    lines.append(f"  1-to-1 orthologs: {format_number(prune_data['ortho1to1_count'])} genes")
                if 'ortho_mi_count' in prune_data:
                    lines.append(f"  Many-to-many orthologs: {format_number(prune_data['ortho_mi_count'])} genes")
        
        # PRANK
        elif func_lower == 'prank':
            if 'prank' in result:
                prank_data = result['prank']
                if 'files_processed' in prank_data:
                    lines.append(f"  Files processed: {format_number(prank_data['files_processed'])} files")
        
        # ASTRAL
        elif func_lower == 'astral':
            if 'num_trees' in result:
                lines.append(f"  Gene trees used: {format_number(result['num_trees'])} trees")
        
        # BUSCO
        elif func_lower == 'busco_extract':
            if 'busco_count' in result:
                lines.append(f"  BUSCO genes extracted: {format_number(result['busco_count'])} genes")
        
        # HCluster
        elif func_lower == 'hcluster_guide_tree':
            if 'cluster_count' in result:
                lines.append(f"  Clusters created: {format_number(result['cluster_count'])} clusters")
            if 'clustered_sequences' in result:
                lines.append(f"  Sequences clustered: {format_number(result['clustered_sequences'])} sequences")
        
        # GeneCluster
        elif func_lower == 'genecluster_for_hcluster':
            if 'gene_count' in result:
                lines.append(f"  Gene clusters created: {format_number(result['gene_count'])} clusters")
        
        return lines
    
    def create_log_entry(self, datahub, original_func, execution_time: float, 
                        new_files: Set[Path], new_dirs: Set[Path],
                        status: str = "SUCCESS",
                        result: Any = None,
                        exception: Optional[Exception] = None) -> List[str]:
        base_dir = datahub.output_dir if datahub.output_dir != datahub.dir_base else datahub.dir_base
        func_name = original_func.__name__
        
        lines = []
        
        lines.append("=" * 80)
        lines.append(f"TREEFORGE PIPELINE LOG - {func_name.upper()} STAGE")
        lines.append("=" * 80)
        lines.append(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Status: {status}")
        lines.append(f"Duration: {execution_time:.2f} seconds")
        
        if func_name.lower() in ITERATIVE_STEPS:
            lines.append(f"Iteration: {datahub.current_iter}")
        
        lines.append("")
        
        description = self._get_function_description(original_func)
        lines.append("DESCRIPTION:")
        lines.append(f"{description}")
        lines.append("")
        
        params = self._get_function_params(datahub, func_name)
        if params:
            lines.append("CONFIGURATION:")
            for param_name, param_value in params.items():
                if isinstance(param_value, float):
                    if param_value < 0.01:
                        formatted_value = f"{param_value:.2e}"
                    else:
                        formatted_value = f"{param_value}"
                elif isinstance(param_value, bool):
                    formatted_value = "Yes" if param_value else "No"
                else:
                    formatted_value = str(param_value)
                
                lines.append(f"  {param_name}: {formatted_value}")
            lines.append("")
        
        metric_lines = self._format_metrics(func_name, result)
        if metric_lines:
            lines.append("PROCESSING METRICS:")
            lines.extend(metric_lines)
            lines.append("")
        
        lines.append("FILE OPERATIONS:")
        
        if new_dirs:
            lines.append(f"  New directories created ({len(new_dirs)}):")
            for d in sorted(new_dirs)[:20]:
                try:
                    rel_path = d.relative_to(base_dir)
                    lines.append(f"    - {rel_path}")
                except ValueError:
                    lines.append(f"    - {d}")
            if len(new_dirs) > 20:
                lines.append(f"    ... and {len(new_dirs) - 20} more directories")
        else:
            lines.append("  New directories created: (none)")
        
        lines.append("")
        
        if new_files:
            lines.append(f"  New files created ({len(new_files)}):")
            
            files_with_sizes = []
            for f in new_files:
                size = get_file_size(f)
                files_with_sizes.append((f, size))
            
            files_with_sizes.sort(key=lambda x: str(x[0]))
            
            for f, size in files_with_sizes[:30]:
                try:
                    rel_path = f.relative_to(base_dir)
                    size_str = format_bytes(size)
                    lines.append(f"    - {rel_path} ({size_str})")
                except ValueError:
                    size_str = format_bytes(size)
                    lines.append(f"    - {f} ({size_str})")
            
            if len(new_files) > 30:
                remaining = len(new_files) - 30
                total_size = sum(size for _, size in files_with_sizes[30:])
                lines.append(f"    ... and {remaining} more files ({format_bytes(total_size)})")
        else:
            lines.append("  New files created: (none)")
        
        lines.append("")
        
        if exception:
            lines.append("ERRORS:")
            lines.append(f"  {type(exception).__name__}: {str(exception)}")
            lines.append("")
        elif status == "WARNING":
            lines.append("WARNINGS:")
            if isinstance(result, dict) and 'warnings' in result:
                warnings = result['warnings']
                if isinstance(warnings, list):
                    for warning in warnings[:10]:
                        lines.append(f"  - {warning}")
                    if len(warnings) > 10:
                        lines.append(f"  ... and {len(warnings) - 10} more warnings")
                else:
                    lines.append(f"  {warnings}")
            else:
                lines.append("  See output for details")
            lines.append("")
        elif status == "SUCCESS":
            lines.append("WARNINGS: (none)")
            lines.append("")
        
        lines.append("SUMMARY:")
        summary = self._generate_summary(func_name, status, result, execution_time, new_files)
        lines.append(f"{summary}")
        lines.append("")
        
        lines.append("=" * 80)
        
        return lines
    
    def _generate_summary(self, func_name: str, status: str, result: Any, 
                         execution_time: float, new_files: Set[Path]) -> str:
        func_lower = func_name.lower()
        
        if status == "FAILED":
            return f"Failed to complete {func_name} stage. Check errors above for details."
        
        if status == "SKIPPED":
            return f"Skipped {func_name} stage as requested."
        
        if func_lower == 'blast':
            if isinstance(result, dict) and 'contig_count' in result:
                count = result['contig_count']
                return (f"Successfully completed all-by-all BLAST search for {format_number(count)} "
                       f"sequences in {execution_time:.1f} seconds. Ready for MCL clustering.")
            return f"Successfully completed BLAST stage in {execution_time:.1f} seconds."
        
        elif func_lower == 'mcl':
            if isinstance(result, dict) and 'mcl_count' in result:
                count = result['mcl_count']
                return (f"Successfully generated {format_number(count)} ortholog clusters using MCL. "
                       f"Ready for alignment and tree building.")
            return f"Successfully completed MCL clustering in {execution_time:.1f} seconds."
        
        elif func_lower == 'mafft':
            if isinstance(result, dict) and 'mafft' in result:
                mafft_data = result['mafft']
                if 'files_processed' in mafft_data:
                    count = mafft_data['files_processed']
                    return (f"Successfully aligned {format_number(count)} clusters and generated "
                           f"preliminary phylogenetic trees in {execution_time:.1f} seconds.")
            return f"Successfully completed MAFFT alignment stage in {execution_time:.1f} seconds."
        
        elif func_lower == 'tree':
            tree_count = len([f for f in new_files if f.suffix == '.tt'])
            if tree_count > 0:
                return (f"Successfully trimmed and processed {format_number(tree_count)} phylogenetic "
                       f"trees in {execution_time:.1f} seconds. Ready for pruning.")
            return f"Successfully completed tree processing stage in {execution_time:.1f} seconds."
        
        elif func_lower == 'prune':
            if isinstance(result, dict) and 'prune' in result:
                prune_data = result['prune']
                ortho1to1  = prune_data.get('ortho1to1_count', 0)
                ortho_mi   = prune_data.get('ortho_mi_count', 0)
                total      = ortho1to1 + ortho_mi
                return (f"Successfully identified {format_number(ortho1to1)} 1-to-1 orthologs and "
                       f"{format_number(ortho_mi)} many-to-many orthologs ({format_number(total)} total) "
                       f"in {execution_time:.1f} seconds.")
            return f"Successfully completed pruning stage in {execution_time:.1f} seconds."
        
        elif func_lower == 'prank':
            return f"Successfully completed PRANK refinement in {execution_time:.1f} seconds. Ready for final tree."
        
        elif func_lower == 'astral':
            if isinstance(result, dict) and 'num_trees' in result:
                count = result['num_trees']
                return (f"Successfully generated species tree from {format_number(count)} gene trees using "
                       f"ASTRAL in {execution_time:.1f} seconds. Pipeline complete!")
            return f"Successfully completed ASTRAL species tree generation in {execution_time:.1f} seconds."
        
        else:
            return f"Successfully completed {func_name} stage in {execution_time:.1f} seconds."
    
    def write_log(self, datahub, original_func, log_entry: List[str]) -> None:
        log_dir = datahub.dir_logs / datahub.run_timestamp
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / f"{original_func.__name__}.log"
        
        try:
            with open(log_file, "a") as f:
                f.write("\n".join(log_entry) + "\n")
        except Exception as e:
            msg = f"Warning: Failed to write log file {log_file}: {e}"
            print(_truncate_for_terminal(msg))
    
    
    def __call__(self, func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            original_func = getattr(func, '_original_func', func)
            
            datahub = args[0] if args else None
            if not datahub or not hasattr(datahub, '__class__') or datahub.__class__.__name__ != 'DataHub':
                return func(*args, **kwargs)
            
            stage_name                  = original_func.__name__
            target_dirs                 = get_stage_target_dirs(datahub, stage_name) #So TF only scans directories that are relevant, instead of the whole directory every stage
            initial_files, initial_dirs = get_filesystem_state(datahub.dir_treeforge, target_dirs)
            
            result     = None
            exception  = None
            start_time = time.time()
            
            try:
                result = func(*args, **kwargs)
            except Exception as e:
                exception = e
                result    = {'status': 'failed', 'error': str(e)}
            
            end_time = time.time()
            execution_time = end_time - start_time
            
            self.remove_phyx_logs(datahub)
            
            final_files, final_dirs = get_filesystem_state(datahub.dir_treeforge, target_dirs)
            new_files               = final_files - initial_files
            new_dirs                = final_dirs - initial_dirs

            status                  = self._determine_status(result, exception)

            log_entry               = self.create_log_entry(
                datahub        = datahub,
                original_func  = original_func,
                execution_time = execution_time,
                new_files      = new_files,
                new_dirs       = new_dirs,
                status         = status,
                result         = result,
                exception      = exception
            )
            
            self.write_log(datahub, original_func, log_entry)
            
            step_name = original_func.__name__
            
            self.update_dict(datahub, step_name, result, execution_time)
            self.print_metrics(datahub, step_name)
            
            if exception:
                raise exception
            
            return result
        
        wrapper._original_func = func
        return wrapper

    def remove_phyx_logs(self, datahub):
        phyx_logs = datahub.dir_base / 'phyx.logfile'
        if phyx_logs.exists():
            phyx_logs.unlink()

def timefunc(func) -> Callable[..., Any]:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        printClass     = PrintOut('info', '', '')
        printout       = printClass.printout
        start_time     = time.time()
        result         = func(*args, **kwargs)
        end_time       = time.time()
        execution_time = end_time - start_time
        printout('metric', {'Time': f'{execution_time:.2f} s'})
        return result
    return wrapper

def bilge_crew() -> Callable[..., Any]:
    return BaseDecorator(skip_metrics=False)

def tree_crew() -> Callable[..., Any]:
    return BaseDecorator(skip_metrics=True)
