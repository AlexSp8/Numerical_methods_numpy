"""
Utilities for file I/O and data processing.
"""

from pathlib import Path
from typing import List


def read_file_float_data(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            row = [float(v) for v in line.split()]
            data.append(row)
    return data

def write_file(file_path: str | Path, content: str, mode: str = 'w') -> None:
    """
    Write content to a file.

    Args:
        file_path: Path to the file to write
        content: Content to write
        mode: File mode ('w' for write, 'a' for append)
    """
    with open(file_path, mode) as f:
        f.write(content + '\n')
