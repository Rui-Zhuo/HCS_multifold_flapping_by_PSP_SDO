# Coronal image processing

Event-oriented scripts for enhancing SOHO/LASCO and related EUV/coronagraph
images, extracting radial or lateral slits, and constructing time-distance maps.

Main components:

- `multiscale_gaussian_normalization_for_fits.py`: reusable FITS MGN CLI;
- `a_trous_wavelet_for_jp2.py`: a-trous and MGN comparison for JP2 sequences;
- `slit_observation.py`: slit extraction and spectral filtering;
- `CoronaImageProcess_utils.py`: shared wavelet, interpolation, and gap helpers;
- `MGN/`: vendored legacy MGN/tomography implementation, kept unchanged.

Run `python -m Corona_image_process.multiscale_gaussian_normalization_for_fits --help`
from the repository root for the portable entry point. Older event scripts retain their
historical absolute path defaults; edit them or migrate the path to the root
`config.toml` before use.

The MGN parameter sets are centralized in `hcs_flapping/constants.py`. They
follow the existing project settings and call `sunkit_image.enhance.mgn`.
Consult [Morgan & Druckmuller (2014)](https://link.springer.com/article/10.1007/s11207-014-0523-9),
*Solar Physics*, 289, 2945--2955, and the
sunkit-image documentation when changing Gaussian widths or weights. Verify the
licence and attribution of files under `MGN/` before redistribution.
