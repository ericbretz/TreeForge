from pathlib import Path
from typing import Union


def format_bytes(bytes_num: int) -> str:
    units      = ['B', 'KB', 'MB', 'GB', 'TB']
    value      = float(bytes_num)
    unit_index = 0
    
    while value >= 1024 and unit_index < len(units) - 1:
        value /= 1024
        unit_index += 1
    
    if unit_index == 0:
        return f"{int(value)} {units[unit_index]}"
    else:
        return f"{value:.1f} {units[unit_index]}"

def format_number(num: int) -> str:
    return f"{num:,}"

def get_file_size(file_path: Union[str, Path]) -> int:
    try:
        path = Path(file_path) if isinstance(file_path, str) else file_path
        if path.exists() and path.is_file():
            return path.stat().st_size
        return 0
    except:
        return 0

def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    
    minutes           = int(seconds // 60)
    remaining_seconds = seconds % 60
    
    if minutes < 60:
        if remaining_seconds > 0:
            return f"{minutes}m {remaining_seconds:.1f}s"
        return f"{minutes}m"
    
    hours             = int(minutes // 60)
    remaining_minutes = minutes % 60
    
    if remaining_seconds > 0:
        return f"{hours}h {remaining_minutes}m {remaining_seconds:.1f}s"
    elif remaining_minutes > 0:
        return f"{hours}h {remaining_minutes}m"
    return f"{hours}h"
