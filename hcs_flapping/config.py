"""TOML-backed project configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

DEFAULT_CONFIG_PATH = Path("config.toml")


@dataclass(frozen=True)
class ProjectConfig:
    """Read-only access to nested project configuration."""

    values: dict[str, Any]
    source: Path | None = None

    def get(self, section: str, key: str, default: Any = None) -> Any:
        """Return one configuration value with a fallback."""
        return self.values.get(section, {}).get(key, default)

    def path(self, key: str, default: str | Path | None = None) -> Path:
        """Return a normalized path from the paths section."""
        value = self.get("paths", key, default)
        if value is None:
            raise KeyError(f"Missing required configuration key: paths.{key}")
        return Path(value).expanduser()


def load_config(path: str | Path | None = None) -> ProjectConfig:
    """Load a TOML file, or return empty configuration when it is absent."""
    config_path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    if not config_path.exists():
        return ProjectConfig(values={}, source=None)
    with config_path.open("rb") as stream:
        values = tomllib.load(stream)
    return ProjectConfig(values=values, source=config_path.resolve())
