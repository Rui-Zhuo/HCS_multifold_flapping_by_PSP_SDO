"""Small, dependency-light helpers shared by command-line scripts."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterator


def ensure_directory(path: str | Path) -> Path:
    """Create an output directory and return its normalized path."""
    directory = Path(path).expanduser()
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def require_file(path: str | Path, description: str = "input file") -> Path:
    """Return an existing file or raise an actionable error."""
    candidate = Path(path).expanduser()
    if not candidate.is_file():
        raise FileNotFoundError(f"{description} does not exist: {candidate}")
    return candidate


def require_directory(path: str | Path, description: str = "input directory") -> Path:
    """Return an existing directory or raise an actionable error."""
    candidate = Path(path).expanduser()
    if not candidate.is_dir():
        raise FileNotFoundError(f"{description} does not exist: {candidate}")
    return candidate


def iter_cadence(start: datetime, end: datetime, minutes: int) -> Iterator[datetime]:
    """Yield inclusive timestamps at a fixed cadence."""
    if minutes <= 0:
        raise ValueError("Cadence must be a positive number of minutes")
    if end < start:
        raise ValueError("End time must not precede start time")
    current = start
    step = timedelta(minutes=minutes)
    while current <= end:
        yield current
        current += step
