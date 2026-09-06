"""Compute and visualize a PFSS model for the 2021-01-17 PSP interval."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.lines import Line2D
import numpy as np

from hcs_flapping.config import load_config
from hcs_flapping.constants import PSP_CARRINGTON_TRACK_DEG
from hcs_flapping.utils import ensure_directory, require_file

DEFAULT_HMI = Path("data/SDO/HMI/hmi.Synoptic_Mr.2239.fits")
DEFAULT_GONG = Path("data/GONG/202101171314.fits")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config.toml"))
    parser.add_argument("--magnetogram", type=Path)
    parser.add_argument("--nrho", type=int)
    parser.add_argument("--rss", type=float, help="Source-surface radius in solar radii")
    parser.add_argument("--field-limit", type=float, help="Input field clipping limit in gauss")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def prepare_magnetogram(path: Path, field_limit: float):
    """Remove monopole offset, replace NaNs, and clip the magnetogram."""
    import sunpy.map

    original = sunpy.map.Map(require_file(path, "synoptic magnetogram"))
    data = np.nan_to_num(original.data.astype(float), nan=0.0)
    data -= np.mean(data)
    data = np.clip(data, -field_limit, field_limit)
    norm = SymLogNorm(linthresh=0.01, vmin=-field_limit, vmax=field_limit)
    return sunpy.map.Map(data, original.meta, plot_settings={"norm": norm})


def plot_source_surface(model_output, rss: float):
    """Plot source-surface radial field, PIL, and PSP Carrington track."""
    source_br = model_output.source_surface_br
    figure = plt.figure()
    axes = figure.add_subplot(projection=source_br, label="Neutral Line")
    source_br.plot(axes=axes)
    if len(model_output.source_surface_pils) > 0:
        axes.plot_coord(model_output.source_surface_pils[0])
    longitudes, latitudes = zip(*PSP_CARRINGTON_TRACK_DEG)
    axes.add_line(Line2D(longitudes, latitudes, linewidth=2, color="black",
                         label="PSP trajectory"))
    figure.colorbar(axes.images[0], ax=axes, orientation="horizontal")
    axes.set_title(f"Source-surface magnetic field at {rss:g} solar radii")
    axes.set_xlim(0, 360)
    axes.set_ylim(-90, 90)
    axes.legend(loc="upper right", frameon=True)
    return figure


def trace_meridional_field(model_input, model_output):
    """Trace field lines seeded in one meridional plane."""
    import astropy.constants as const
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import pfsspy

    figure, axes = plt.subplots()
    axes.set_aspect("equal")
    radius = 1.01 * const.R_sun
    longitude = np.pi / 2 * u.rad
    latitude = np.linspace(-np.pi / 2, np.pi / 2, 20) * u.rad
    seeds = SkyCoord(longitude, latitude, radius, frame=model_output.coordinate_frame)
    tracer = pfsspy.tracing.PythonTracer()
    for field_line in tracer.trace(seeds, model_output):
        coordinates = field_line.coords
        coordinates.representation_type = "cartesian"
        color = {0: "black", -1: "tab:red", 1: "tab:blue"}[field_line.polarity]
        axes.plot(coordinates.y / const.R_sun, coordinates.z / const.R_sun, color=color)
    axes.add_patch(patches.Circle((0, 0), 1, color="black", fill=False))
    axes.add_patch(patches.Circle((0, 0), model_input.grid.rss, color="black",
                                  linestyle="--", fill=False))
    axes.set_xlabel(r"$y/R_\odot$")
    axes.set_ylabel(r"$z/R_\odot$")
    axes.set_title("PFSS meridional-plane solution")
    return figure


def trace_field_3d(model_output):
    """Trace a sparse global set of PFSS field lines."""
    import astropy.constants as const
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import pfsspy

    figure = plt.figure()
    axes = figure.add_subplot(111, projection="3d")
    radius = 1.2 * const.R_sun
    latitude = np.linspace(-np.pi / 2, np.pi / 2, 10, endpoint=False)
    longitude = np.linspace(0, 2 * np.pi, 10, endpoint=False)
    latitude, longitude = np.meshgrid(latitude, longitude, indexing="ij")
    seeds = SkyCoord(longitude.ravel() * u.rad, latitude.ravel() * u.rad, radius,
                     frame=model_output.coordinate_frame)
    tracer = pfsspy.tracing.PythonTracer()
    for field_line in tracer.trace(seeds, model_output):
        coordinates = field_line.coords
        coordinates.representation_type = "cartesian"
        color = {0: "black", -1: "tab:red", 1: "tab:blue"}[field_line.polarity]
        axes.plot(coordinates.x / const.R_sun, coordinates.y / const.R_sun,
                  coordinates.z / const.R_sun, color=color, linewidth=1)
    axes.set_title("PFSS global solution")
    return figure


def main() -> None:
    """Run the configured PFSS calculation and create three figures."""
    args = parse_args()
    import pfsspy

    config = load_config(args.config)
    if args.magnetogram:
        magnetogram_path = args.magnetogram
    elif config.get("pfss", "magnetogram", "hmi").lower() == "gong":
        magnetogram_path = config.path("gong_synoptic", DEFAULT_GONG)
    else:
        magnetogram_path = config.path("hmi_synoptic", DEFAULT_HMI)
    nrho = args.nrho or int(config.get("pfss", "nrho", 10))
    rss = args.rss or float(config.get("pfss", "source_surface_radius", 2.5))
    field_limit = args.field_limit or float(config.get("pfss", "field_limit_gauss", 20.0))
    magnetogram = prepare_magnetogram(magnetogram_path, field_limit)
    model_input = pfsspy.Input(magnetogram, nrho, rss)
    model_output = pfsspy.pfss(model_input)
    figures = [
        ("source_surface", plot_source_surface(model_output, rss)),
        ("meridional", trace_meridional_field(model_input, model_output)),
        ("global_3d", trace_field_3d(model_output)),
    ]
    if args.output:
        output_dir = ensure_directory(args.output)
        for name, figure in figures:
            figure.savefig(output_dir / f"pfss_{name}.png", dpi=300, bbox_inches="tight")
            plt.close(figure)
    if args.show or not args.output:
        plt.show()


if __name__ == "__main__":
    main()
