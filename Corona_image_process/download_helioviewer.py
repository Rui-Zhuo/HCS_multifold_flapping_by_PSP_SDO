"""Compare near-simultaneous SOHO/LASCO images from Helioviewer and VSO."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt

from hcs_flapping.utils import ensure_directory


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time", default="2011-01-14T12:00:00")
    parser.add_argument("--output", type=Path, default=Path("outputs/helioviewer"))
    parser.add_argument("--detector", default="C2", choices=("C2", "C3"))
    parser.add_argument("--save", type=Path)
    return parser.parse_args()


def main() -> None:
    """Download, compare, and plot both LASCO products."""
    args = parse_args()
    from sunpy.io import read_file
    from sunpy.map import Map
    from sunpy.net import Fido, attrs as a, helioviewer

    observation_time = datetime.fromisoformat(args.time)
    output_dir = ensure_directory(args.output)
    client = helioviewer.HelioviewerClient()
    hv_file = client.download_jp2(
        observation_time,
        observatory="SOHO",
        instrument="LASCO",
        detector=args.detector,
        directory=str(output_dir),
    )
    hv_data, hv_header = read_file(hv_file)[0]
    hv_map = Map(hv_data, hv_header)

    result = Fido.search(
        a.Time(observation_time, observation_time + timedelta(minutes=10)),
        a.Instrument.lasco,
        a.Detector(args.detector.lower()),
    )
    if len(result) == 0:
        raise RuntimeError("VSO returned no LASCO records")
    downloaded = Fido.fetch(result[0, 0], path=str(output_dir))
    vso_data, vso_header = read_file(downloaded[0])[0]
    vso_map = Map(vso_data, vso_header)

    print(f"Helioviewer CROTA2: {hv_header.get('CROTA2')}")
    print(f"VSO CROTA2: {vso_header.get('CROTA2')}")
    figure = plt.figure(figsize=(12, 5))
    hv_axes = figure.add_subplot(1, 2, 1, projection=hv_map)
    vso_axes = figure.add_subplot(1, 2, 2, projection=vso_map)
    hv_map.plot(axes=hv_axes)
    vso_map.plot(axes=vso_axes)
    hv_axes.set_title("Helioviewer JP2")
    vso_axes.set_title("VSO FITS")
    figure.tight_layout()
    if args.save:
        save_path = ensure_directory(args.save.parent) / args.save.name
        figure.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
