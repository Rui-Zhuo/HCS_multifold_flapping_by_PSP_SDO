"""Download and plot one GOES/SUVI FITS image."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import requests

from hcs_flapping.utils import ensure_directory

DEFAULT_URL = (
    "https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/"
    "goes/goes18/l2/data/suvi-l2-ci195/2023/03/26/"
    "dr_suvi-l2-ci195_g18_s20230326T000000Z_e20230326T000400Z_v1-0-2.fits"
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--data-dir", type=Path, default=Path("outputs/goes"))
    parser.add_argument("--save", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Download the FITS file if needed and render the SUVI map."""
    args = parse_args()
    import sunpy.map

    data_dir = ensure_directory(args.data_dir)
    destination = data_dir / args.url.rsplit("/", 1)[-1]
    if args.overwrite or not destination.exists():
        response = requests.get(args.url, timeout=60)
        response.raise_for_status()
        destination.write_bytes(response.content)
        print(f"Downloaded: {destination}")
    goes_map = sunpy.map.Map(destination)
    figure = plt.figure(figsize=(8, 8))
    axes = figure.add_subplot(projection=goes_map)
    goes_map.plot(axes=axes)
    figure.colorbar(axes.images[0], ax=axes)
    if args.save:
        save_path = ensure_directory(args.save.parent) / args.save.name
        figure.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
