"""Plot a radial-field slice from an FDIPS IDL output file."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hcs_flapping.config import load_config
from hcs_flapping.utils import ensure_directory, require_file

DEFAULT_INPUT_DIR = Path("data/FDIPS")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--input", type=Path, help="FDIPS field output file")
    parser.add_argument("--radial-index", type=int, default=-2)
    parser.add_argument("--limit", type=float, default=0.001, help="Symmetric field limit")
    parser.add_argument("--save", type=Path)
    return parser.parse_args()


def plot_fdips_slice(input_path: Path, radial_index: int, field_limit: float):
    """Build the FDIPS radial-field slice figure."""
    from spacepy.pybats import IdlFile

    fdips_field = IdlFile(str(require_file(input_path, "FDIPS field output")))
    longitude = np.asarray(fdips_field["Longitude"])
    latitude = np.asarray(fdips_field["Latitude"])
    radial_field = np.asarray(fdips_field["Br"])
    radial_slice = radial_field[radial_index, :, :]

    figure, axes = plt.subplots(figsize=(8, 4))
    mesh = axes.pcolormesh(longitude, latitude, radial_slice.T, cmap="bwr",
                           vmin=-field_limit, vmax=field_limit, shading="auto")
    figure.colorbar(mesh, ax=axes, label=r"$B_r$")
    axes.plot([70, 100], [-2, -4], "k", label="PSP trajectory")
    axes.set_xlabel("Carrington longitude (deg)")
    axes.set_ylabel("Carrington latitude (deg)")
    axes.set_aspect("equal")
    axes.legend()
    figure.tight_layout()
    return figure


def main() -> None:
    """Plot the selected FDIPS output."""
    args = parse_args()
    config = load_config(args.config)
    default_file = config.path("fdips_output", DEFAULT_INPUT_DIR) / "fdips_field_360_720.out"
    figure = plot_fdips_slice(args.input or default_file, args.radial_index, args.limit)
    if args.save:
        output = ensure_directory(args.save.parent) / args.save.name
        figure.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved: {output}")
        plt.close(figure)
    else:
        plt.show()


if __name__ == "__main__":
    main()
