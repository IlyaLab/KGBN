"""Helpers for accessing packaged KGBN data files."""

from importlib import resources
from pathlib import Path


DATA_PACKAGE = "KGBN.data"


def data_path(filename: str) -> Path:
    """Return a filesystem path for a packaged data file."""
    return Path(resources.files(DATA_PACKAGE).joinpath(filename))


def data_files() -> list[str]:
    """Return the names of packaged data files."""
    return sorted(
        item.name
        for item in resources.files(DATA_PACKAGE).iterdir()
        if item.is_file() and item.name != "__init__.py"
    )
