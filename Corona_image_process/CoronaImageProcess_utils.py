"""Numerical utilities for coronagraph image enhancement and slit extraction."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
from scipy import signal
from scipy.interpolate import interp2d

def build_atrous_coef(j: int, method: str = 'B_spline') -> np.ndarray:
    """Build the separable 2-D a-trous smoothing kernel for scale ``j``."""
    if not isinstance(j, (int, np.integer)) or j < 0:
        raise ValueError("j must be a non-negative integer")
    if method == 'B_spline':
        if j == 0:
            atrous_coef = np.array([1/16.0, 1/4.0, 3/8.0, 1/4.0, 1/16.0])
            return np.outer(atrous_coef, atrous_coef)
        elif j >= 1:
            len_atrous_eff = 5 + 2**(j+1)
            atrous_coef = np.zeros(len_atrous_eff)
            atrous_coef[0] = 1/16.0
            atrous_coef[int(2**(j-1)+1)] = 1/4.0
            atrous_coef[int((len_atrous_eff-1)/2)] = 3/8.0
            atrous_coef[int(3*2**(j-1)+3)] = 1/4.0
            atrous_coef[-1] = 1/16.0
            return np.outer(atrous_coef, atrous_coef)
    elif method == 'linear':
        if j == 0:
            atrous_coef = np.array([0.25, 0.5, 0.25])
            return np.outer(atrous_coef, atrous_coef)
        elif j >= 1:
            len_atrous_eff = 3 + 2**j
            atrous_coef = np.zeros(len_atrous_eff)
            atrous_coef[0] = 0.25
            atrous_coef[int((len_atrous_eff-1)/2)] = 0.5
            atrous_coef[-1] = 0.25
            return np.outer(atrous_coef, atrous_coef)
    else:
        raise ValueError("method must be 'B_spline' or 'linear'")
        
def a_trous_wavelet_2D(input_data: np.ndarray, level_num: int,
                       method: str = 'B_spline') -> np.ndarray:
    """Return a 2-D a-trous wavelet decomposition without downsampling."""
    if not isinstance(level_num, (int, np.integer)) or level_num < 1:
        raise ValueError("level_num must be a positive integer")
    input_data = np.asarray(input_data)
    if input_data.ndim != 2:
        raise ValueError("input_data must be a two-dimensional array")
    if method == 'B_spline':
        if level_num == 1:
            len_atrous_coef = 5
        elif level_num >= 2:
            len_atrous_coef = 5 + 2**(level_num)
    elif method == 'linear':
        if level_num == 1:
            len_atrous_coef = 3
        elif level_num >= 2:
            len_atrous_coef = 3 + 2**(level_num-1)
    else:
        raise ValueError("method must be 'B_spline' or 'linear'")
        
    data_shape = np.shape(input_data)

    if np.min(data_shape) > len_atrous_coef:
        output_c = np.empty(shape=(data_shape[0], data_shape[1], level_num))
        wavelet_coef = np.empty(shape=(data_shape[0], data_shape[1], level_num))
        output_c[:,:,0] = np.copy(input_data)
        for i in range(1,level_num):
            output_c[:, :, i] = signal.convolve(
                input_data, build_atrous_coef(i, method), mode="same"
            )
        wavelet_coef[:,:,:-1] = -np.diff(output_c, axis=2)
        wavelet_coef[:,:,-1] = np.copy(output_c[:,:,-1])
        return wavelet_coef
    raise ValueError(
        f"level_num={level_num} requires a kernel smaller than both image dimensions {data_shape}"
    )

def interp_to_slit(image: np.ndarray, x_slit: np.ndarray,
                   y_slit: np.ndarray) -> np.ndarray:
    """Linearly interpolate an image along paired pixel coordinates."""
    image = np.asarray(image)
    x_slit = np.asarray(x_slit)
    y_slit = np.asarray(y_slit)
    if image.ndim != 2:
        raise ValueError("image must be a two-dimensional array")
    if x_slit.shape != y_slit.shape:
        raise ValueError("x_slit and y_slit must have the same shape")
    # Construct Cartesian Coordinate
    x = np.linspace(0, image.shape[0] - 1, image.shape[0])
    y = np.linspace(0, image.shape[1] - 1, image.shape[1])
    # Interpolate with interp2d
    interpfun = interp2d(x, y, image, kind='linear')
    interp_matrix = interpfun(x_slit, y_slit)
    interp_array = np.diagonal(interp_matrix)
    slit_pixels = interp_array.reshape(x_slit.shape)
    
    return slit_pixels

def radial_slit(beg_point: Sequence[float], end_point: Sequence[float],
                min_x: float, min_y: float, step: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """Sample a ray from ``beg_point`` toward ``end_point`` to a lower boundary."""
    beg_x, beg_y = beg_point
    end_x, end_y = end_point
    direct = np.array([end_x - beg_x, end_y - beg_y])
    distance = np.linalg.norm(direct)
    if distance == 0:
        raise ValueError("beg_point and end_point must differ")
    if step <= 0:
        raise ValueError("step must be positive")
    if direct[0] >= 0 and direct[1] >= 0:
        raise ValueError("ray must approach min_x or min_y to terminate")
    e_direct = direct / distance
    x_slit, y_slit = [beg_x], [beg_y]
    while True:
        new_x = x_slit[-1] + step * e_direct[0]
        new_y = y_slit[-1] + step * e_direct[1]
        if new_x < min_x or new_y < min_y:
            break
        x_slit.append(new_x)
        y_slit.append(new_y)
    return np.array(x_slit), np.array(y_slit)

def insert_nan_columns(image: np.ndarray,
                       insertions: Iterable[tuple[int, int]]) -> np.ndarray:
    """Insert missing-cadence columns represented by NaN values."""
    image = np.asarray(image, dtype=float)
    for col, num_insertions in insertions:
        for _ in range(num_insertions):
            image = np.insert(image, col, np.nan, axis=1)
            col += 1
    return image
