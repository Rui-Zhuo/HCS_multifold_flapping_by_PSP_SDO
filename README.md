# Multi-fold HCS flapping observed by PSP and SDO

Research code for analysing successive heliospheric current sheet (HCS)
crossings observed by Parker Solar Probe (PSP), together with photospheric and
coronal magnetic-field context from SDO/HMI and potential-field models.

## Repository layout

```text
.
|-- hcs_flapping/                 Shared Python configuration and utilities
|-- tests/                        Fast tests independent of mission data
|-- Corona_image_process/         Coronagraph enhancement and slit analysis
|   `-- MGN/                      Vendored legacy MGN/tomography code
|-- FDIPS/                        External FDIPS source (not modified here)
|-- FDIPS_SWMF/                   Local FDIPS/SWMF products (not modified here)
|-- PSI-master/                   External PSI source (ignored by Git)
|-- *.m                           PSP/HCS analysis and MATLAB figures
`-- *.py                          Download, field-model and plotting entry points
```

The scripts retain their historical filenames so existing workflows continue to
work. Machine-specific paths and event parameters belong in `config.toml`; copy
the documented template before running data-dependent scripts:

```powershell
Copy-Item config.example.toml config.toml
python -m pip install -e ".[dev]"
python -m unittest discover -s tests -v
```

Python 3.10 is the reference interpreter. Dependency versions in
`pyproject.toml` reproduce the environment in which this cleanup was validated.
No stochastic algorithm is currently used; `random_seed = 20210117` is recorded
in the configuration for future reproducible extensions.

## Typical workflows

```powershell
# Inspect available options without accessing data
python download_HMI_full_disk_fits.py --help
python plot_HMI_full_disk_field.py --help
python pfss.py --help

# Use the default local configuration
python get_Earth_position.py
python plot_situations_for_successive_crossings_revised.py --save outputs/hcs_scenarios.png
```

All data downloaders and plotting scripts now create output directories when
needed, validate their inputs, and avoid performing work when imported.

## Data and reproducibility

Mission data are intentionally not committed. The main external data products
are PSP in-situ observations, SDO/HMI magnetograms, GONG synoptic maps,
Helioviewer/SOHO-LASCO imagery, and FDIPS outputs. Update `[paths]` in
`config.toml` for the local archive. JSOC downloads additionally require a valid
email in `[jsoc]`; the template contains no personal address.

The large MATLAB scripts preserve the original event-selection and physical
analysis logic. Paths are centralized in `hcs_config.m`. For a specific machine,
copy `hcs_config_local_template.m` to the ignored `hcs_config_local.m`, rename
the function accordingly, and update the paths. The same values may instead be
set through the `HCS_*` environment variables documented in `hcs_config.m`.

## External code

`FDIPS`, `FDIPS_SWMF`, and `PSI-master` are excluded from the cleanup scope.
`Corona_image_process/MGN` contains legacy code attributed in-file to Huw Morgan
and includes its own notices. Check the redistribution terms of every external
component before making the GitHub repository public.

## Related publication and citation

This repository contains the analysis code associated with:

> Zhuo, R., He, J., Duan, D., Wu, Z., & Hou, C. (2026).
> “Flapping of Multifold Heliospheric Current Sheet Observed by Parker Solar
> Probe.” *The Astrophysical Journal*, **1000**, 243.
> [https://doi.org/10.3847/1538-4357/ae40ff](https://iopscience.iop.org/article/10.3847/1538-4357/ae40ff/meta)

BibTeX:

```bibtex
@article{Zhuo2026HCSFlapping,
  author  = {Zhuo, Rui and He, Jiansen and Duan, Die and Wu, Ziqi and Hou, Chuanpeng},
  title   = {Flapping of Multifold Heliospheric Current Sheet Observed by Parker Solar Probe},
  journal = {The Astrophysical Journal},
  year    = {2026},
  volume  = {1000},
  number  = {2},
  eid     = {243},
  doi     = {10.3847/1538-4357/ae40ff},
  url     = {https://doi.org/10.3847/1538-4357/ae40ff}
}
```

GitHub-compatible citation metadata are provided in `CITATION.cff`. No
project-wide open-source licence is asserted by this cleanup. Add a licence only
after confirming that it is compatible with the external code and data products
included in the repository.
