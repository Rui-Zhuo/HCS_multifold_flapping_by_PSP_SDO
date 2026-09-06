"""Download SDO/HMI 720 s full-disk line-of-sight magnetograms from JSOC."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from hcs_flapping.config import load_config
from hcs_flapping.constants import HMI_CADENCE_MINUTES
from hcs_flapping.utils import ensure_directory, iter_cadence

DEFAULT_OUTPUT = Path("data/SDO/HMI/full_disk")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--output", type=Path, help="HMI FITS output directory")
    parser.add_argument("--email", help="JSOC notification email")
    parser.add_argument("--start", help="UTC start time in ISO format")
    parser.add_argument("--end", help="UTC end time in ISO format")
    parser.add_argument("--cadence", type=int, help="Sampling cadence in minutes")
    return parser.parse_args()


def download_hmi(output_dir: Path, email: str, start: datetime, end: datetime,
                 cadence_minutes: int) -> None:
    """Download the first HMI record found around each requested timestamp."""
    import astropy.units as u
    from astropy.time import Time
    from sunpy.net import Fido, attrs as a

    output_dir = ensure_directory(output_dir)
    for timestamp in iter_cadence(start, end, cadence_minutes):
        search_time = Time(timestamp, scale="utc")
        print(f"Searching HMI record near {search_time.isot}")
        query = Fido.search(
            a.Time(search_time - 1 * u.s, search_time + 1 * u.s),
            a.Sample(cadence_minutes * u.min),
            a.jsoc.Series.hmi_m_720s,
            a.jsoc.Notify(email),
        )
        if len(query) == 0 or len(query[0]) == 0:
            print("No data found")
            continue
        record = query[0][0]
        record_time = str(record["T_REC"])
        compact_time = (
            f"{record_time[0:4]}{record_time[5:7]}{record_time[8:10]}_"
            f"{record_time[11:13]}{record_time[14:16]}{record_time[17:19]}"
        )
        expected = output_dir / f"hmi.m_720s.{compact_time}_TAI.1.magnetogram.fits"
        if expected.exists():
            print(f"Skipping existing file: {expected.name}")
            continue
        files = Fido.fetch(record, path=str(output_dir), overwrite=False)
        print(f"Saved: {files}")


def main() -> None:
    """Run the HMI downloader from CLI and TOML configuration."""
    args = parse_args()
    config = load_config(args.config)
    email = args.email or config.get("jsoc", "email")
    if not email or email == "your-email@example.com":
        raise ValueError("Set jsoc.email in config.toml or pass --email")
    output = args.output or config.path("hmi_full_disk", DEFAULT_OUTPUT)
    start_text = args.start or config.get("jsoc", "start", "2021-01-16T00:00:00")
    end_text = args.end or config.get("jsoc", "end", "2021-01-17T23:59:59")
    cadence = args.cadence or config.get("jsoc", "cadence_minutes", HMI_CADENCE_MINUTES)
    download_hmi(output, email, datetime.fromisoformat(start_text),
                 datetime.fromisoformat(end_text), int(cadence))


if __name__ == "__main__":
    main()
