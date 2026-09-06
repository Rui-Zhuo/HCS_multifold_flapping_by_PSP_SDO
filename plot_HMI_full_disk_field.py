"""Plot HMI full-disk magnetograms around a Carrington footpoint."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.colors as colors
import matplotlib.pyplot as plt

from hcs_flapping.config import load_config
from hcs_flapping.utils import ensure_directory, require_directory

DEFAULT_INPUT = Path("data/SDO/HMI/full_disk")
HMI_PATTERN = re.compile(
    r"hmi\.m_720s\.(\d{8})_(\d{6})_TAI\.\d+\.magnetogram\.fits$",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--longitude", type=float)
    parser.add_argument("--latitude", type=float)
    parser.add_argument("--fov", type=float, help="Field of view in arcsec")
    parser.add_argument("--all-cadences", action="store_true")
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def plot_magnetogram(file_path: Path, output_dir: Path, longitude_deg: float,
                     latitude_deg: float, fov_arcsec: float, show: bool) -> Path:
    """Plot one HMI map and mark the specified Carrington coordinate."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.map
    from sunpy.coordinates import frames
    from sunpy.coordinates.ephemeris import get_earth

    hmi_map = sunpy.map.Map(file_path)
    observer = get_earth(hmi_map.date)
    carrington_coord = SkyCoord(
        lon=longitude_deg * u.deg,
        lat=latitude_deg * u.deg,
        frame=frames.HeliographicCarrington,
        obstime=hmi_map.date,
        observer=observer,
    )
    x_pixel, y_pixel = hmi_map.world_to_pixel(
        carrington_coord.transform_to(hmi_map.coordinate_frame)
    )

    figure = plt.figure(figsize=(10, 8))
    axes = figure.add_subplot(111, projection=hmi_map)
    norm = colors.SymLogNorm(linthresh=5, vmin=-500, vmax=500, base=10)
    image = axes.imshow(hmi_map.data, cmap="bwr", norm=norm)
    figure.colorbar(image, ax=axes, shrink=0.9, pad=0.05, label=r"$B_{LOS}$ (G)",
                    extend="both")
    axes.scatter(
        x_pixel,
        y_pixel,
        color="lime",
        edgecolors="black",
        s=100,
        linewidths=2,
        label=f"Carrington footpoint\n({longitude_deg:.1f} deg, {latitude_deg:.1f} deg)",
    )
    axes.legend(loc="upper right", fontsize=10)
    half_fov_pixels = ((fov_arcsec * u.arcsec) / hmi_map.scale.axis1).value / 2
    axes.set_xlim(x_pixel.value - half_fov_pixels, x_pixel.value + half_fov_pixels)
    axes.set_ylim(y_pixel.value - half_fov_pixels, y_pixel.value + half_fov_pixels)
    axes.set_aspect("equal")
    axes.set_title(f"{hmi_map.date.isot} HMI magnetogram")
    axes.invert_xaxis()
    axes.invert_yaxis()
    figure.tight_layout()

    match = HMI_PATTERN.match(file_path.name)
    record_id = "_".join(match.groups()) if match else file_path.stem
    output_path = output_dir / f"{record_id}_Br.png"
    if show:
        plt.show()
    else:
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(figure)
    return output_path


def main() -> None:
    """Plot all selected magnetograms in the configured directory."""
    args = parse_args()
    config = load_config(args.config)
    input_dir = require_directory(args.input or config.path("hmi_full_disk", DEFAULT_INPUT),
                                  "HMI full-disk directory")
    output_dir = ensure_directory(args.output or config.path("output", "outputs") / "hmi")
    longitude = args.longitude if args.longitude is not None else float(
        config.get("hmi_plot", "footpoint_longitude_deg", 101.4)
    )
    latitude = args.latitude if args.latitude is not None else float(
        config.get("hmi_plot", "footpoint_latitude_deg", -27.2)
    )
    fov = args.fov or float(config.get("hmi_plot", "fov_arcsec", 200.0))

    selected = []
    for path in sorted(input_dir.glob("*.fits")):
        match = HMI_PATTERN.match(path.name)
        if match and (args.all_cadences or match.group(2)[2:] == "0000"):
            selected.append(path)
    if not selected:
        raise FileNotFoundError(f"No matching HMI FITS files found in {input_dir}")
    for path in selected:
        result = plot_magnetogram(path, output_dir, longitude, latitude, fov, args.show)
        print(f"Processed: {result}")


if __name__ == "__main__":
    main()
