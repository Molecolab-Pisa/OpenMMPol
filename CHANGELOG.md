# Changelog: OpenMMPol 1.2.0 (August 2026)

## 🚀 New Feature: Density Fitting Module

Complete implementation of Density Fitting Charges (DFC) approximation to represent the QM electrostatic potential via fitted charges.

**What changed:**
- Added a new DFC module with grid generation, SVD charge fitting, energy and gradient computation
- Support for atom-centered grid and rigid multi-points grids (Fibonacci grid and cubic grid)
- Exposed DC through C, Fortran, and Python interfaces
- Added JSON configuration for DFC grid types and parameters
- Integrated DFC into the system lifecycle (init, update, terminate)
- Extended memory allocation to support 4D arrays needed by DFC

### Detailed Changes

- **`src/mod_density_fit.F90`** — Full DFC module: grid generation, SVD charge fitting, X/Xinv matrices, gradient computation
- **`include/openmmpol.h`** — DFC C API declarations
- **`src/mod_c_interface.F90`** — C bindings 
- **`src/mod_interface.F90`** — Fortran interface: `ommp_set_fit_potential`, `ommp_df_compute_induced_dipoles`, `ommp_df_geomgrad`, `ommp_init_density_fit`
- **`src/mod_mmpol.F90`** — `df` field added to system type, coordinate update → `df_update`, lifecycle management
- **`pythonapi/pyopenmmpol/pyommp_interface.cpp`** — pybind11 bindings (`get_df_*`, `df_geomgrad`, all nabla matrices, `df_get_nabla`)
- **`src/smartinput.c`** — Full `density_fit` JSON section parsing
- **`config/openmmpol_const.h.in`** — preprocessor defined parameter
- **`src/mod_constants.F90`** — DFC Fortran parameters (`ommp_df_*`)
- **`src/CMakeLists.txt`** — Added `mod_density_fit.F90`
- **`src/mod_memory.F90`** — Added `r_alloc4`, `i_alloc4`, `r_free4`, `i_free4` for 4D arrays

---

## 🛠️  Improvements

Several minor improvements mainly in the interface layer.

**What changed:**
- Improved polarization LS configuration the convergence threshold and the use of the previous solution as initial guess can now be decided from C, Fortran and Python interfaces or by the user using JSON keys `polarization_ls_conv_thr` and `polarization_ls_use_guess`
- `ommp-ip` 
  - supports some new arguments to extend the control of the output JSON file  (`--never-use-fmm`, `--use-fmm`, and `--fail-on-duplicated-prm` ).
  - precision of QM coordinates is truncated to 5 digits for consistency with the `.xyz` QM file written by MDAnalysis.

- `ommp_xyz2mmp` now automatically ignores duplicated parameters in Tinker `.prm` forcefield definition, as they are very common.



### Detailed Changes

- **`pythonapi/pyopenmmpol/ommp_ip/ommp_ip.py`** — CLI flag `--fail-on-duplicated-prm`, `--never-use-fmm` and `--use-fmm`
- **`app/ommp_xyz2mmp.c`** — Added `ignore_duplicated_prm()` calls
- **`src/mod_solvers.F90`** — Exported `OMMP_DEFAULT_SOLVER_TOL`
- **`src/mod_electrostatics.F90`** — Added `def_conv_thr` and `def_use_guess` to electrostatics type
- **`src/mod_polarization.F90`** — Added optional `arg_tol` / `arg_use_guess` to `polarization()`; passed tolerance to iterative solvers
- **`src/smartinput.c`** — Added parsing of `polarization_ls_conv_thr` and `polarization_ls_use_guess` JSON keys
- Interfaces across C / Fortran / Python for `ommp_set_polarization_conv_thr` and `ommp_set_polarization_use_guess`

---

## 🐛 Bug Fixes

**What changed:**
- Fixed buffer overflow and string comparison bugs in smartinput version parsing
- Ensured fatal error messages are always printed regardless of verbosity setting
- Fixed issues for partially polarizable systems (related to screening lists).

### Detailed Changes

- **`src/smartinput.c`** — Fixed buffer overflow: `commit[8]` → `commit[32]`; replaced `strcpy` with `sprintf` for commit hash; fixed `strcmp` → `strncmp` (8-char) in version comparison
- **`src/mod_io.F90`** — Fatal error messages now always printed (verbosity `-1` overrides filtering)
- **`src/mod_electrostatics.F90`** — Fixed `polar_mm` → `mm_polar` variable name in FMM; fixed `remove_null_pol` signature (`rebuild_list` arg); fixed FMM tree cleanup ordering; added `free_screening_list`
- **`src/mod_link_atom.F90`** — Updated `remove_null_pol` call with new `rebuild_list` arg

---

## 🧹 Cleanup

**What changed:**
- Removed `input-preparation` sub-project and its Learning Kit data. Now this has its own life as [GraphAssign](https://github.com/mattiabondanza/GraphAssign)
- Removed large test case and model data files (~120k lines total)

---

## 📝 Documentation

**What changed:**
- Added a comprehensive new user guide covering installation, JSON interface, `ommp-ip`, MMP format etc.
- Rewrote the README: moved installation details to user guide, added project structure diagram
