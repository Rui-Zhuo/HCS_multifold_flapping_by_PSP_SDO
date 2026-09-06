"""Print the Earth position in Heliocentric Inertial coordinates."""

from __future__ import annotations

import argparse

def main() -> None:
    """Calculate and print the Earth position at one UTC time."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("time", nargs="?", default="2021-01-17 12:00:00")
    args = parser.parse_args()
    from sunpy.coordinates import frames
    from sunpy.coordinates.ephemeris import get_earth

    earth_position = get_earth(args.time)
    earth_hci = earth_position.transform_to(frames.HeliocentricInertial)
    print("Earth in HCI coordinates:", earth_hci)


if __name__ == "__main__":
    main()
