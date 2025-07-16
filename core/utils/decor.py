import time
import json
import functools
import os
from datetime            import datetime
from pathlib             import Path
from typing              import Set, Tuple, Any, Callable, Dict, List, Optional
from core.utils.printout import PrintOut

ITERATIVE_STEPS     = {'mafft', 'tree'}
STEP_METRICS_CONFIG = {
    'blast' : {'exclude_keys': {'files'}},
    'mcl'   : {'exclude_keys': set()},
    'prune' : {'nested_key': 'prune'},
    'prank' : {'exclude_keys': set()},
    'astral': {'exclude_keys': set()}
}

def get_filesystem_state(path: Path) -> Tuple[Set[Path], Set[Path]]:
    """Scans a directory and returns all files and subdirectories found."""
    files = set()
    dirs  = set()
    
    for root, directories, filenames in os.walk(path):
        root_path = Path(root)
        dirs.update(root_path / d for d in directories)
        files.update(root_path / f for f in filenames)
    
    return files, dirs

def convert_paths(obj: Any) -> Any:
    """Turns Path objects into strings so they can be saved as JSON."""
    if isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: convert_paths(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_paths(item) for item in obj]
    return obj

class BaseDecorator:
    """Handles the common stuff that both decorators need to do."""
    
    def __init__(self, save_flag: Optional[bool] = None, skip_metrics: bool = False):
        self.save_flag    = save_flag
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
        """Print nested dict as sorted list, excluding certain keys."""
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
    
    def create_log_entry(self, datahub, original_func, execution_time: float, 
                        new_files: Set[Path], new_dirs: Set[Path]) -> List[str]:
        """Builds the log message that gets written to the file."""
        base_dir = datahub.output_dir if datahub.output_dir != datahub.dir_base else datahub.dir_base
        
        return [
            f"Date: {datetime.now().strftime('%Y-%m-%d')}",
            f"Time: {time.strftime('%H:%M:%S', time.localtime())}",
            f"\nFunction: {original_func.__name__}",
            f"Time taken: {execution_time:.2f} seconds",
            f"\nNew directories created:",
            *[f"  - {d.relative_to(base_dir)}" for d in sorted(new_dirs)],
            f"\nNew files created:",
            *[f"  - {f.relative_to(base_dir)}" for f in sorted(new_files)],
            "-" * 50
        ]
    
    def write_log(self, datahub, original_func, log_entry: List[str]) -> None:
        """Writes the log entry to a file."""
        log_dir = datahub.dir_logs / datahub.run_timestamp
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / f"{original_func.__name__}.log"
        
        try:
            with open(log_file, "a") as f:
                f.write("\n".join(log_entry) + "\n")
        except Exception as e:
            print(f"Warning: Failed to write log file {log_file}: {e}")
    
    def save_master_dict(self, datahub) -> None:
        """Saves the master dictionary to a JSON file."""
        try:
            datahub.dir_logs.mkdir(parents=True, exist_ok=True)
            serializable_dict = convert_paths(datahub.master_dict)
            with open(os.path.join(datahub.dir_logs, "previous_run.json"), "w") as f:
                json.dump(serializable_dict, f, indent=4)
        except Exception as e:
            print(f"Warning: Failed to save master dict: {e}")
    
    def __call__(self, func: Callable) -> Callable:
        """The main decorator that wraps the function."""
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            original_func = getattr(func, '_original_func', func)
            
            datahub = args[0] if args else None
            if not datahub or not hasattr(datahub, '__class__') or datahub.__class__.__name__ != 'DataHub':
                return func(*args, **kwargs)
            
            initial_files, initial_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            start_time     = time.time()
            result         = func(*args, **kwargs)
            end_time       = time.time()
            execution_time = end_time - start_time
            
            self.remove_phyx_logs(datahub)
            
            final_files, final_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            new_files = final_files - initial_files
            new_dirs  = final_dirs  - initial_dirs
            
            log_entry = self.create_log_entry(datahub, original_func, execution_time, new_files, new_dirs)
            self.write_log(datahub, original_func, log_entry)
            
            should_save = self.save_flag if self.save_flag is not None else getattr(datahub, 'save', True)
            step_name = original_func.__name__
            
            self.update_dict(datahub, step_name, result, execution_time)
            self.print_metrics(datahub, step_name)
            
            if should_save:
                self.save_master_dict(datahub)
            
            return result
        
        wrapper._original_func = func
        return wrapper

    def remove_phyx_logs(self, datahub):
        """Removes the phyx logs."""
        phyx_logs = os.path.join(datahub.dir_base, 'phyx.logfile')
        if os.path.exists(phyx_logs):
            os.remove(phyx_logs)

def timefunc(func) -> Callable[..., Any]:
    """Simple timing decorator."""
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

def bilge_crew(save_flag=None) -> Callable[..., Any]:
    """Main decorator that times functions, tracks file changes, saves results to logs and dict."""
    return BaseDecorator(save_flag=save_flag, skip_metrics=False)

def tree_crew(save_flag=None) -> Callable[..., Any]:
    """Same as bilge_crew but doesn't print metrics."""
    return BaseDecorator(save_flag=save_flag, skip_metrics=True)
