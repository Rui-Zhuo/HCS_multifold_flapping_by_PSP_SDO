"""Download pre-exported SDO/HMI SHARP CEA FITS records."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import requests

from hcs_flapping.config import load_config
from hcs_flapping.constants import HMI_CADENCE_MINUTES, HMI_COMPONENTS
from hcs_flapping.utils import ensure_directory, iter_cadence

DEFAULT_OUTPUT = Path("data/SDO/HMI/SHARP/AR2796")
DEFAULT_URL_TEMPLATE = (
    "https://jsoc1.stanford.edu/SUM38/D1937912030/S00000/"
    "hmi.sharp_cea_720s.{sharp}.{timestamp}_TAI.{component}.fits"
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--start", default="2021-01-18T00:00:00")
    parser.add_argument("--end", default="2021-01-19T23:59:59")
    parser.add_argument("--cadence", type=int, default=HMI_CADENCE_MINUTES)
    parser.add_argument("--sharp", type=int)
    parser.add_argument("--url-template", default=DEFAULT_URL_TEMPLATE)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def download_sharp(output_dir: Path, sharp_number: int, start: datetime, end: datetime,
                   cadence_minutes: int, url_template: str,
                   overwrite: bool = False) -> None:
    """Download SHARP components from a previously submitted JSOC export."""
    output_dir = ensure_directory(output_dir)
    with requests.Session() as session:
        for timestamp in iter_cadence(start, end, cadence_minutes):
            timestamp_text = timestamp.strftime("%Y%m%d_%H%M%S")
            for component in HMI_COMPONENTS:
                url = url_template.format(sharp=sharp_number, timestamp=timestamp_text,
                                          component=component)
                destination = output_dir / url.rsplit("/", 1)[-1]
                if destination.exists() and not overwrite:
                    print(f"Skipping existing file: {destination.name}")
                    continue
                print(f"Downloading {timestamp_text} {component}")
                response = session.get(url, timeout=60)
                response.raise_for_status()
                destination.write_bytes(response.content)


def main() -> None:
    """Run the SHARP downloader."""
    args = parse_args()
    config = load_config(args.config)
    output = args.output or config.path("hmi_sharp", DEFAULT_OUTPUT)
    sharp = args.sharp or int(config.get("jsoc", "sharp_number", 7532))
    download_sharp(output, sharp, datetime.fromisoformat(args.start),
                   datetime.fromisoformat(args.end), args.cadence,
                   args.url_template, args.overwrite)


if __name__ == "__main__":
    main()
