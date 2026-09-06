"""Plot SDO/HMI SHARP CEA magnetograms for one active region."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt

from hcs_flapping.config import load_config
from hcs_flapping.utils import ensure_directory, require_directory, require_file

DEFAULT_INPUT = Path("data/SDO/HMI/SHARP/AR2796")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sharp", type=int)
    parser.add_argument("--date", default="20210117", help="Record date as YYYYMMDD")
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Plot matching SHARP magnetograms."""
    args = parse_args()
    import sunpy.map

    config = load_config(args.config)
    input_dir = require_directory(args.input or config.path("hmi_sharp", DEFAULT_INPUT),
                                  "SHARP directory")
    output_dir = ensure_directory(args.output or config.path("output", "outputs") / "sharp")
    sharp = args.sharp or int(config.get("jsoc", "sharp_number", 7532))
    pattern = re.compile(
        rf"hmi\.sharp_cea_720s\.{sharp}\.(\d{{8}})_(\d{{6}})_TAI\.Br\.fits$",
        re.IGNORECASE,
    )

    count = 0
    for br_path in sorted(input_dir.glob("*.fits")):
        match = pattern.match(br_path.name)
        if not match or match.group(1) != args.date:
            continue
        record_date, record_time = match.groups()
        magnetogram_path = require_file(
            input_dir
            / f"hmi.sharp_cea_720s.{sharp}.{record_date}_{record_time}_TAI.magnetogram.fits",
            "paired SHARP magnetogram",
        )
        magnetogram = sunpy.map.Map(magnetogram_path)
        figure, axes = plt.subplots(figsize=(8, 3))
        image = axes.pcolormesh(magnetogram.data, cmap="gray", shading="auto")
        figure.colorbar(image, ax=axes, label="magnetogram (G)")
        axes.set_aspect("equal")
        axes.set_title(f"{record_date}_{record_time}", fontsize=14)
        figure.tight_layout()
        output_path = output_dir / f"{record_date}_{record_time}_Br_mag.png"
        if args.show:
            plt.show()
        else:
            figure.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close(figure)
        count += 1
        print(f"Processed: {output_path}")
    if count == 0:
        raise FileNotFoundError(f"No SHARP Br records for {args.date} found in {input_dir}")


if __name__ == "__main__":
    main()
