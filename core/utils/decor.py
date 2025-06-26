import time
from datetime import datetime
import os
from pathlib import Path
from typing import Set, Tuple, Any, Callable
import json
from core.utils.printout import PrintOut

def get_filesystem_state(path: Path) -> Tuple[Set[Path], Set[Path]]:
    """Get files and directories in the given path."""
    files = set()
    dirs  = set()
    
    for root, directories, filenames in os.walk(path):
        root_path = Path(root)
        dirs.update(root_path / d for d in directories)
        files.update(root_path / f for f in filenames)
    
    return files, dirs

def timefunc(func) -> Callable[..., Any]:
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
    """
    Main decorator that times functions, tracks file changes, saves results to logs and dict.
    """
    def decorator(func) -> Callable[..., Any]:
        def convert_paths(obj):
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            return obj

        def update_dict(datahub, step_name, step_dict, execution_time):
            """Add step results to the master dictionary."""
            # Handle iteration steps (mafft and tree)
            if step_name.lower() in ['mafft', 'tree']:
                iter_key = f'iter_{datahub.current_iter}'
                if iter_key not in datahub.master_dict:
                    datahub.master_dict[iter_key] = {}
                if step_name.lower() not in datahub.master_dict[iter_key]:
                    datahub.master_dict[iter_key][step_name.lower()] = {}
                
                # Add timing info
                if 'metrics' not in datahub.master_dict[iter_key][step_name.lower()]:
                    datahub.master_dict[iter_key][step_name.lower()]['metrics'] = {}
                datahub.master_dict[iter_key][step_name.lower()]['metrics']['time'] = f'{execution_time:.2f} s'
                
                # Update with step results
                if isinstance(step_dict, dict):
                    for key, value in step_dict.items():
                        existing_value = datahub.master_dict[iter_key][step_name.lower()].get(key)
                        if existing_value is not None and isinstance(existing_value, dict) and isinstance(value, dict):
                            datahub.master_dict[iter_key][step_name.lower()][key].update(value)
                        else:
                            datahub.master_dict[iter_key][step_name.lower()][key] = value
                else:
                    datahub.master_dict[iter_key][step_name.lower()]['data'] = step_dict
            else:
                # Handle regular steps
                if step_name.lower() not in datahub.master_dict:
                    datahub.master_dict[step_name.lower()] = {}
                
                # Add timing info
                if 'metrics' not in datahub.master_dict[step_name.lower()]:
                    datahub.master_dict[step_name.lower()]['metrics'] = {}
                datahub.master_dict[step_name.lower()]['metrics']['time'] = f'{execution_time:.2f} s'
                
                # Update with step results
                if isinstance(step_dict, dict):
                    for key, value in step_dict.items():
                        existing_value = datahub.master_dict[step_name.lower()].get(key)
                        if existing_value is not None and isinstance(existing_value, dict) and isinstance(value, dict):
                            datahub.master_dict[step_name.lower()][key].update(value)
                        else:
                            datahub.master_dict[step_name.lower()][key] = value
                else:
                    datahub.master_dict[step_name.lower()]['data'] = step_dict

        def wrapper(*args, **kwargs):
            # Get the original function
            original_func = getattr(func, '_original_func', func)
            
            # Check if first arg is DataHub
            datahub = args[0] if args else None
            if not datahub or not hasattr(datahub, '__class__') or datahub.__class__.__name__ != 'DataHub':
                return func(*args, **kwargs)
            
            def print_nested_dict(d, exclude_keys=None):
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
            
            initial_files, initial_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time
            
            final_files, final_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            new_files = final_files - initial_files
            new_dirs  = final_dirs - initial_dirs
            
            log_entry = [
                f"Date: {datetime.now().strftime('%Y-%m-%d')}",
                f"Time: {time.strftime('%H:%M:%S', time.localtime())}",
                f"\nFunction: {original_func.__name__}",
                f"Time taken: {execution_time:.2f} seconds",
                f"\nNew directories created:",
                *[f"  - {d.relative_to(datahub.dir_base)}" for d in sorted(new_dirs)],
                f"\nNew files created:",
                *[f"  - {f.relative_to(datahub.dir_base)}" for f in sorted(new_files)],
                "-" * 50
            ]
            
            # Create log directory with timestamp
            log_dir = datahub.dir_logs / datahub.run_timestamp
            log_dir.mkdir(parents=True, exist_ok=True)
            
            log_file = log_dir / f"{original_func.__name__}.log"
            
            try:
                with open(log_file, "a") as f:
                    f.write("\n".join(log_entry) + "\n")
            except Exception as e:
                pass
            
            # Save to master dictionary
            should_save = save_flag if save_flag is not None else getattr(datahub, 'save', True)
            
            # Always update the master dictionary and print metrics, regardless of save flag
            step_name = original_func.__name__
            update_dict(datahub, step_name, result, execution_time)
            
            # Handle iteration steps
            if step_name.lower() in ['mafft', 'tree']:
                iter_key = f'iter_{datahub.current_iter}'
                datahub.printout('metric', datahub.master_dict[iter_key][step_name.lower()])
            else:
                # Use the same filtering logic as _prev_print
                metrics = datahub.master_dict[step_name.lower()]
                
                if step_name.lower() == 'blast':
                    blast_metrics = print_nested_dict(metrics, exclude_keys={'files'})
                    if blast_metrics:
                        datahub.printout('metric', blast_metrics)
                elif step_name.lower() == 'mcl':
                    mcl_metrics = print_nested_dict(metrics)
                    if mcl_metrics:
                        datahub.printout('metric', mcl_metrics)
                elif step_name.lower() == 'prune':
                    if 'prune' in metrics:
                        prune_metrics = print_nested_dict(metrics['prune'])
                        if prune_metrics:
                            datahub.printout('metric', prune_metrics)
                elif step_name.lower() == 'prank':
                    prank_metrics = print_nested_dict(metrics)
                    if prank_metrics:
                        datahub.printout('metric', prank_metrics)
                elif step_name.lower() == 'astral':
                    astral_metrics = print_nested_dict(metrics)
                    if astral_metrics:
                        datahub.printout('metric', astral_metrics)
                else:
                    # Default case for other steps
                    filtered_metrics = print_nested_dict(metrics)
                    if filtered_metrics:
                        datahub.printout('metric', filtered_metrics)
            
            # Only save to JSON if should_save is True
            if should_save:
                # Save to JSON
                datahub.dir_logs.mkdir(parents=True, exist_ok=True)
                serializable_dict = convert_paths(datahub.master_dict)
                with open(datahub.dir_logs / "previous_run.json", "w") as f:
                    json.dump(serializable_dict, f, indent=4)
            
            return result
        
        # Store original function
        wrapper._original_func = func
        
        # Preserve metadata
        import functools
        return functools.update_wrapper(wrapper, func)
    return decorator

def tree_crew(save_flag=None) -> Callable[..., Any]:
    """
    Special decorator for Tree function - same as bilge_crew but skips printing metrics.
    """
    def decorator(func) -> Callable[..., Any]:
        def convert_paths(obj):
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            return obj

        def update_dict(datahub, step_name, step_dict, execution_time):
            """Add step results to the master dictionary."""
            # Handle iteration steps (mafft and tree)
            if step_name.lower() in ['mafft', 'tree']:
                iter_key = f'iter_{datahub.current_iter}'
                if iter_key not in datahub.master_dict:
                    datahub.master_dict[iter_key] = {}
                if step_name.lower() not in datahub.master_dict[iter_key]:
                    datahub.master_dict[iter_key][step_name.lower()] = {}
                
                # Add timing info
                if 'metrics' not in datahub.master_dict[iter_key][step_name.lower()]:
                    datahub.master_dict[iter_key][step_name.lower()]['metrics'] = {}
                datahub.master_dict[iter_key][step_name.lower()]['metrics']['time'] = f'{execution_time:.2f} s'
                
                # Update with step results (only if not None)
                if step_dict is not None and isinstance(step_dict, dict):
                    for key, value in step_dict.items():
                        existing_value = datahub.master_dict[iter_key][step_name.lower()].get(key)
                        if existing_value is not None and isinstance(existing_value, dict) and isinstance(value, dict):
                            datahub.master_dict[iter_key][step_name.lower()][key].update(value)
                        else:
                            datahub.master_dict[iter_key][step_name.lower()][key] = value
                elif step_dict is not None:
                    datahub.master_dict[iter_key][step_name.lower()]['data'] = step_dict

        def wrapper(*args, **kwargs):
            # Get original function
            original_func = getattr(func, '_original_func', func)
            
            # Check if first arg is DataHub
            datahub = args[0] if args else None
            if not datahub or not hasattr(datahub, '__class__') or datahub.__class__.__name__ != 'DataHub':
                return func(*args, **kwargs)
            
            initial_files, initial_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time
            
            final_files, final_dirs = get_filesystem_state(datahub.dir_treeforge)
            
            new_files = final_files - initial_files
            new_dirs  = final_dirs - initial_dirs
            
            log_entry = [
                f"Date: {datetime.now().strftime('%Y-%m-%d')}",
                f"Time: {time.strftime('%H:%M:%S', time.localtime())}",
                f"\nFunction: {original_func.__name__}",
                f"Time taken: {execution_time:.2f} seconds",
                f"\nNew directories created:",
                *[f"  - {d.relative_to(datahub.dir_base)}" for d in sorted(new_dirs)],
                f"\nNew files created:",
                *[f"  - {f.relative_to(datahub.dir_base)}" for f in sorted(new_files)],
                "-" * 50
            ]
            
            # Create log directory with timestamp
            log_dir = datahub.dir_logs / datahub.run_timestamp
            log_dir.mkdir(parents=True, exist_ok=True)
            
            log_file = log_dir / f"{original_func.__name__}.log"
            
            try:
                with open(log_file, "a") as f:
                    f.write("\n".join(log_entry) + "\n")
            except Exception as e:
                pass
            
            # Save to master dictionary
            should_save = save_flag if save_flag is not None else getattr(datahub, 'save', True)
            if should_save:
                step_name = original_func.__name__
                update_dict(datahub, step_name, result, execution_time)
                
                # Save to JSON (no metrics printing for tree)
                datahub.dir_logs.mkdir(parents=True, exist_ok=True)
                serializable_dict = convert_paths(datahub.master_dict)
                with open(datahub.dir_logs / "previous_run.json", "w") as f:
                    json.dump(serializable_dict, f, indent=4)
            
            return result
        
        # Store original function
        wrapper._original_func = func
        
        # Preserve metadata
        import functools
        return functools.update_wrapper(wrapper, func)
    return decorator
