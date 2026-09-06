# Contributing

1. Create a focused branch and keep mission data outside Git.
2. Copy `config.example.toml` to the ignored `config.toml` for local paths.
3. Run `python -m unittest discover -s tests -v` and `python -m compileall` before committing.
4. Keep physical units explicit in names, docstrings, labels, and configuration.
5. Do not modify vendored code in `FDIPS`, `FDIPS_SWMF`, `PSI-master`, or
   `Corona_image_process/MGN` without documenting its provenance and licence.

Changes to event selection or plasma/field calculations should include a short
scientific justification and a comparison with a known event product.
