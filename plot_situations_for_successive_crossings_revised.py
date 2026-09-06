"""Plot four schematic HCS-crossing and flapping configurations."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from hcs_flapping.utils import ensure_directory

LINE_WIDTH = 3
FONT_SIZE = 16
QUIVER_LENGTH = 2
HEAD_WIDTH = 5
POINT_SIZE = 20
X_LIMITS = (-10, 10)
Y_LIMITS = (-6, 6)


def hcs_normal(x_cross: float, wavelength: float) -> tuple[float, float]:
    """Return the unit normal of a cosine-shaped current sheet."""
    slope_term = 6 * np.pi / wavelength * np.sin(2 * np.pi / wavelength * x_cross)
    magnitude = np.hypot(slope_term, 1)
    return slope_term / magnitude, 1 / magnitude


def plot_hcs(axes, x: np.ndarray, y: np.ndarray, alpha: float = 1.0) -> None:
    """Draw an HCS curve with a position-dependent colour."""
    colour_norm = Normalize(vmin=x.min(), vmax=x.max())
    axes.scatter(x, y, c=plt.cm.winter(colour_norm(x)), s=POINT_SIZE,
                 alpha=alpha, edgecolors="none")
    axes.set_xlim(*X_LIMITS)
    axes.set_ylim(*Y_LIMITS)
    axes.set_aspect("equal")


def plot_hcs_normals(axes, x_cross: np.ndarray, y_cross: np.ndarray, wavelength: float):
    """Draw the local HCS normal at every crossing."""
    quiver = None
    for x_value, y_value in zip(x_cross, y_cross):
        x_normal, y_normal = hcs_normal(x_value, wavelength)
        quiver = axes.quiver(
            x_value, y_value, QUIVER_LENGTH * x_normal, QUIVER_LENGTH * y_normal,
            color="red", linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH,
            angles="xy", scale_units="xy", scale=1,
        )
    return quiver


def plot_psp_orbit(axes, x_psp: np.ndarray):
    """Draw the PSP path and propagation direction."""
    axes.plot(x_psp, np.zeros_like(x_psp), "k", linewidth=LINE_WIDTH)
    return axes.quiver(
        1, 0, -2, 0, color="black", linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH,
        angles="xy", scale_units="xy", scale=1,
    )


def plot_steady_motion(axes, x: float, y_up: float, y_down: float):
    """Draw opposing steady-flapping directions."""
    upper = axes.quiver(
        x, y_up, 0, QUIVER_LENGTH, edgecolor="orange", facecolor="none",
        linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH * 2,
        angles="xy", scale_units="xy", scale=1,
    )
    axes.quiver(
        x, y_down, 0, -QUIVER_LENGTH, edgecolor="orange", facecolor="none",
        linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH * 2,
        angles="xy", scale_units="xy", scale=1,
    )
    return upper


def plot_kink_motion(axes, x_up: float, y_up: float, x_down: float, y_down: float):
    """Draw the propagation of a kink-like HCS perturbation."""
    upper = axes.quiver(
        x_up, y_up, QUIVER_LENGTH, 0, edgecolor="orange", facecolor="none",
        linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH * 2,
        angles="xy", scale_units="xy", scale=1,
    )
    axes.quiver(
        x_down, y_down, QUIVER_LENGTH, 0, edgecolor="orange", facecolor="none",
        linewidth=LINE_WIDTH, headwidth=HEAD_WIDTH * 2,
        angles="xy", scale_units="xy", scale=1,
    )
    return upper


def label_times(axes, labels: list[tuple[float, float, str]]) -> None:
    """Add schematic crossing-time labels."""
    for x_value, y_value, text in labels:
        axes.text(x_value, y_value, text, color="darkblue", fontsize=FONT_SIZE,
                  fontweight="bold")


def build_figure():
    """Construct the four-panel HCS scenario figure."""
    x_hcs = np.linspace(-10, 10, 2001)
    x_psp = np.linspace(-10, 10, 2001)
    wavelength = 10.0
    folded_hcs = 3 * np.cos(2 * np.pi / wavelength * x_hcs)
    flat_hcs = np.zeros_like(x_hcs) - 0.3
    crossings = np.array([-1.5, -0.5, 0.5, 1.5]) / 2 * wavelength
    zero_crossings = np.zeros_like(crossings)

    figure, axes = plt.subplots(2, 2, figsize=(10, 6))
    figure.subplots_adjust(hspace=0.3, wspace=0.1, bottom=0.15)

    axes[0, 0].set_title("(a) Static folded CS", fontsize=FONT_SIZE)
    plot_hcs(axes[0, 0], x_hcs, folded_hcs)
    plot_hcs_normals(axes[0, 0], crossings, zero_crossings, wavelength)
    plot_psp_orbit(axes[0, 0], x_psp)

    axes[0, 1].set_title("(b) Steady flapping of unfolded CS", fontsize=FONT_SIZE)
    plot_hcs(axes[0, 1], x_hcs, flat_hcs)
    plot_hcs(axes[0, 1], x_hcs, flat_hcs + 3, alpha=1e-2)
    plot_hcs(axes[0, 1], x_hcs, flat_hcs - 3, alpha=5e-3)
    plot_psp_orbit(axes[0, 1], x_psp)
    plot_hcs_normals(axes[0, 1], crossings, np.full_like(crossings, -0.3), 1e10)
    plot_steady_motion(axes[0, 1], x=-8, y_up=1, y_down=-1)
    label_times(axes[0, 1], [(-9.8, 0.4, r"$t_1$"), (-9.8, 3.3, r"$t_2$"),
                              (-9.8, -2.7, r"$t_3$")])

    axes[1, 0].set_title("(c) Kink-like flapping of folded CS", fontsize=FONT_SIZE)
    plot_hcs(axes[1, 0], x_hcs, folded_hcs)
    plot_hcs(axes[1, 0], x_hcs + 3, folded_hcs, alpha=1e-2)
    plot_hcs(axes[1, 0], x_hcs + 6, folded_hcs, alpha=5e-3)
    plot_hcs_normals(axes[1, 0], crossings, zero_crossings, wavelength)
    plot_psp_orbit(axes[1, 0], x_psp)
    plot_kink_motion(axes[1, 0], x_up=-7.5, y_up=2, x_down=-9.5, y_down=-2)
    label_times(axes[1, 0], [(-9.8, 3.3, r"$t_1$"), (-6.8, 3.3, r"$t_2$"),
                              (-3.8, 3.3, r"$t_3$")])

    axes[1, 1].set_title("(d) Steady flapping of folded CS", fontsize=FONT_SIZE)
    plot_hcs(axes[1, 1], x_hcs, folded_hcs)
    plot_hcs(axes[1, 1], x_hcs, folded_hcs + 1.8, alpha=1e-2)
    plot_hcs(axes[1, 1], x_hcs, folded_hcs - 1.8, alpha=5e-3)
    normal_key = plot_hcs_normals(axes[1, 1], crossings, zero_crossings, wavelength)
    psp_key = plot_psp_orbit(axes[1, 1], x_psp)
    motion_key = plot_steady_motion(axes[1, 1], x=-8, y_up=2.2, y_down=0.2)
    label_times(axes[1, 1], [(-9.7, 3.1, r"$t_1$"), (-9.7, 4.9, r"$t_2$"),
                              (-9.7, 1.3, r"$t_3$")])

    for panel in axes.flat:
        panel.axis("off")
    axes[1, 1].quiverkey(psp_key, -0.5, -0.1, 2, label="PSP orbit", labelpos="E",
                         fontproperties={"size": FONT_SIZE - 5})
    axes[1, 1].quiverkey(normal_key, 0, -0.1, 2, label="CS normal", labelpos="E")
    axes[1, 1].quiverkey(motion_key, 0.5, -0.1, 2, label="CS motion", labelpos="E")
    figure.text(0.99, 0.01, "plotted by plot_situations_for_successive_crossings_revised.py",
                ha="right", fontsize=FONT_SIZE - 4)
    return figure


def main() -> None:
    """Render or save the HCS scenario figure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", type=Path)
    args = parser.parse_args()
    figure = build_figure()
    if args.save:
        output = ensure_directory(args.save.parent) / args.save.name
        figure.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved: {output}")
        plt.close(figure)
    else:
        plt.show()


if __name__ == "__main__":
    main()
