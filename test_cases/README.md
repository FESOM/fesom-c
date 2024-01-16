## Model Test Cases

This directory contains tests for the model.

### LE_medium:
- Test for Lock Exchange on a medium grid, approximately 4500 nodes.

#### MESH - Grid for the Experiment:
- `nod2d.out` - Grid nodes.
- `elem2d.out` - Elements.
- `depth.out` - Depths.
- `aux3d.out` - Required for partitioning.

Other files and directories are created using `partit`. Directories named `*_XX` represent grid partitions for a specified number of CPUs.

- `fv_run.dat` - Parameter file.
- `namelist.ice` - Namelist for the ice model (not used in the LE experiment).
- `outputfilelist.dat` - File to control output into NetCDF files.

#### postproc
- energy_tocheck.dat  - file with total energy, to control
- interp_lines_z_sigma.py  - python script to regrid output on transect aslon X (NetCDF output)
- lines_regular_0_T25.nc - example of output from interp_lines_z_sigma.py
