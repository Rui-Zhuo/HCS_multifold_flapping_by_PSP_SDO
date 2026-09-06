"""Apply multiscale Gaussian normalization (MGN) to a solar FITS image."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hcs_flapping.constants import MGN_PARAMETERS
from hcs_flapping.utils import ensure_directory, require_file


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input FITS file")
    parser.add_argument("--instrument", choices=tuple(MGN_PARAMETERS), default="lasco_c2")
    parser.add_argument("--save", type=Path)
    return parser.parse_args()


def apply_mgn(image_map, instrument: str):
    """Apply the project MGN parameter set and preserve map metadata."""
    import sunkit_image.enhance as enhance
    import sunpy.map

    parameters = MGN_PARAMETERS[instrument]
    normalized = enhance.mgn(
        image_map.data.astype(float, copy=False),
        sigma=parameters["sigma"],
        weights=parameters["weights"],
        k=parameters["k"],
        gamma=parameters["gamma"],
        h=parameters["h"],
    )
    return sunpy.map.Map(normalized, image_map.meta)


def mask_lasco_c2(mgn_map):
    """Mask the LASCO/C2 image outside 2.4--6 observed solar radii."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import sunpy.map
    from sunpy.map.maputils import all_coordinates_from_map

    pixel_coords = all_coordinates_from_map(mgn_map)
    solar_center = SkyCoord(0 * u.deg, 0 * u.deg, frame=mgn_map.coordinate_frame)
    pixel_radii = np.hypot(pixel_coords.Tx - solar_center.Tx,
                           pixel_coords.Ty - solar_center.Ty)
    mask = (pixel_radii < mgn_map.rsun_obs * 2.4) | (pixel_radii > mgn_map.rsun_obs * 6)
    return sunpy.map.Map(mgn_map.data, mgn_map.meta, mask=mask)


def main() -> None:
    """Render the raw and normalized solar images."""
    args = parse_args()
    import sunpy.map

    image_map = sunpy.map.Map(require_file(args.input, "solar FITS image"))
    mgn_map = apply_mgn(image_map, args.instrument)
    maps = [("Raw", image_map), ("MGN", mgn_map)]
    if args.instrument == "lasco_c2":
        maps.append(("Masked MGN", mask_lasco_c2(mgn_map)))
    figure = plt.figure(figsize=(6 * len(maps), 5))
    for index, (title, current_map) in enumerate(maps, start=1):
        axes = figure.add_subplot(1, len(maps), index, projection=current_map)
        current_map.plot(axes=axes)
        axes.set_title(title)
        if title == "Masked MGN":
            current_map.draw_limb(axes=axes)
    figure.tight_layout()
    if args.save:
        save_path = ensure_directory(args.save.parent) / args.save.name
        figure.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
