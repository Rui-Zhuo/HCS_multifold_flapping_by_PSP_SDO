"""Project-wide numerical and plotting constants."""

from __future__ import annotations

RANDOM_SEED = 20210117

HMI_CADENCE_MINUTES = 12
HMI_SAMPLE_MINUTES = (0, 12, 24, 36, 48)
HMI_COMPONENTS = ("Bp", "Bt", "Br", "magnetogram")

PSP_CARRINGTON_TRACK_DEG = ((67.923843, -2.4845026), (97.512878, -3.8445756))

MGN_PARAMETERS = {
    "lasco_c2": {
        "sigma": (1.25, 2.5, 5.0, 10.0, 20.0, 40.0),
        "weights": (0.907, 0.976, 1.0, 1.0, 1.0),
        "k": 0.8,
        "gamma": 1.0,
        "h": 0.9,
        "vmin": -0.1,
        "vmax": 1.1,
    },
    "eui_fsi": {
        "sigma": (1.25, 2.5, 5.0, 10.0),
        "weights": (0.907, 0.976, 1.0, 1.0),
        "k": 0.7,
        "gamma": 2.0,
        "h": 0.9,
        "vmin": 0.07,
        "vmax": 1.0,
    },
    "aia": {
        "sigma": (1.25, 2.5, 5.0, 10.0, 20.0, 40.0),
        "weights": (0.907, 0.976, 1.0, 1.0, 1.0),
        "k": 0.7,
        "gamma": 3.2,
        "h": 0.7,
        "vmin": -0.01,
        "vmax": 1.0,
    },
}
