# GitHub release checklist

## Scientific and technical checks

- Run `python -m unittest discover -s tests -v`.
- Run the representative event with local PSP, HMI, and FDIPS products and
  compare key crossing times, LMN directions, and figures with archived results.
- Start MATLAB successfully and run `checkcode` on modified `.m` files.
- Confirm `config.toml` and `hcs_config_local.m` remain ignored.
- Confirm no FITS/CDF/HDF mission data, credentials, or personal paths are staged.

## Provenance and licensing

- Decide whether `FDIPS`, `FDIPS_SWMF`, and `PSI-master` should be Git submodules,
  release assets, or omitted from the public repository.
- Review the redistribution notices under `FDIPS` and
  `Corona_image_process/MGN`.
- Select a project-wide licence only after confirming compatibility with the
  bundled external code.
- Confirm the author and repository fields in `CITATION.cff`.
- Confirm that GitHub renders the preferred article citation from
  `CITATION.cff`.

## GitHub settings

- Protect `main` and require the CI workflow before merging.
- Enable issue/PR templates only if external contributions will be accepted.
- Create a `v0.1.0` release after tagging the exact validated commit.
