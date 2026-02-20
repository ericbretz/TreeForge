import sys
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Optional, Union, List, Callable


def run_logged_subprocess(
    cmd           : Union[str, list],
    subprocess_dir: Path,
    log_name      : str,
    shell         : bool = True,
    check         : bool = True,
    **kwargs
) -> subprocess.CompletedProcess:

    subprocess_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp     = datetime.now().strftime('%H%M%S')
    stdout_file   = subprocess_dir / f"{log_name}_{timestamp}.stdout.txt"
    stderr_file   = subprocess_dir / f"{log_name}_{timestamp}.stderr.txt"
    
    result = subprocess.run(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        shell  = shell,
        check  = False,
        **kwargs
    )
    
    if result.stdout:
        with open(stdout_file, 'wb') as f:
            f.write(result.stdout)
    
    if result.stderr:
        with open(stderr_file, 'wb') as f:
            f.write(result.stderr)
    
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, 
            cmd, 
            output = result.stdout,
            stderr = result.stderr
        )
    
    return result

def get_subprocess_dir(dir_logs: Path, run_timestamp: str) -> Path:
    return dir_logs / run_timestamp / 'subprocesses'

def run_stage_subprocess(
    cmd                : str,
    stage_name         : str,
    subprocess_dir     : Optional[Path],
    printout_func      : Callable,
    check              : bool = True,
    allowed_returncodes: List[int] = None
) -> subprocess.CompletedProcess:
    if allowed_returncodes is None:
        allowed_returncodes = [0]
    
    try:
        if subprocess_dir:
            result = run_logged_subprocess(cmd, subprocess_dir, stage_name.replace(' ', '_'), 
                                          shell=True, check=False)
        else:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                   shell=True, check=False)
        
        if check and result.returncode not in allowed_returncodes:
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            
        return result
        
    except subprocess.CalledProcessError as e:
        printout_func('error', f'{stage_name} failed')
        if e.stderr:
            stderr_text = e.stderr.decode('utf-8') if isinstance(e.stderr, bytes) else str(e.stderr)
            printout_func('error', stderr_text)
        sys.exit(1)
