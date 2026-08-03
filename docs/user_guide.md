# OpenMMPol User Guide

## Polarizable QM/MM with OpenMMPol

**Version:** 1.1.4
**License:** LGPL-3.0

---

## Table of Contents

- [1. What Is OpenMMPol?](#1-what-is-openmmpol)
- [2. Dependencies](#2-dependencies)
- [3. Building from Source](#3-building-from-source)
- [4. Conda](#4-conda)
- [5. The Smart Input JSON Interface](#5-the-smart-input-json-interface)
- [6. Input Preparation with OMMP-IP](#6-input-preparation-with-ommp-ip)
- [7. The MMP File Format](#7-the-mmp-file-format)
- [8. Running a Calculation](#8-running-a-calculation)
- [9. Troubleshooting](#9-troubleshooting)
- [Further Reading](#further-reading)

---

## 1. What Is OpenMMPol?

OpenMMPol is a __library__ that computes QM/MM energies and forces of a molecular system modeled with polarizable molecular mechanics in induced-point-dipoles formalism (eg. MMPol, AMOEBA). It is built specifically to allow an easy implementation of QM/MMPol and QM/AMOEBA multiscale models within QM software. 
The MMPol region can be described by a Tinker-compatible force field (tested: AMOEBA, AMBER).

**Key features:**
- QM/MM boundary with link atoms
- Fast Multipole Method (FMM) for large systems
- Analytical geometry gradients
- QM-DFC model to simplify and speedup QM/MMPol

- HDF5 checkpoint/restart (**not complete**)

**Tested Force fields:** AMOEBA (`amoeba09.prm`, `amoebabio18.prm`), AMBER (`amber99sb.prm`). See the [Tinker `.prm` format](https://dasher.wustl.edu/tinker/).

---

## 2. Dependencies

### Required
Fortran 2008 compiler (GCC ≥7.5, Intel ≥2024, NVIDIA HPC ≥24), CMake ≥3.20, LAPACK, OpenSSL, cJSON, OpenMP.

### Optional
HDF5 (HDF5 I/O, `ommp_pp`), pybind11 + NumPy (Python bindings).

### ommp-ip (input tool)
MDAnalysis, NumPy, SciPy — installed automatically with Python support.

### Install Dependencies (OpenSuse)
```bash
zypper in cJSON-devel gcc gcc-c++ gcc-fortran make cmake \
  python lapack-devel liblapack3 hdf5 hdf5-devel zlib-devel \
  openssl-devel python3-MDAnalysis python3-numpy python3-pybind11
```

---

## 3. Building from Source

```bash
git clone https://github.com/Molecolab-Pisa/OpenMMPol.git
cd OpenMMPol
cp config.cmake custom.cmake   # customize options
cmake -C custom.cmake -B build
cmake --build build -j
cmake --install build
```


### CMake Options

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Debug, Release, RelWithDebInfo | Release | Compiler flags |
| `CMAKE_INSTALL_PREFIX` | any path | — | Install directory |
| `WITH_HDF5` | ON/OFF | ON | HDF5 support |
| `WITH_PYTHON` | ON/OFF | ON | Python bindings + OMMP-IP |
| `WITH_TESTS` | ON/OFF | ON | Compile tests |
| `TESTLANG` | C, F03 | C | Test language |
| `CMAKE_Fortran_COMPILER` | path | — | Fortran compiler |
| `CMAKE_C_COMPILER` | path | — | C compiler |
| `CMAKE_CXX_COMPILER` | path | — | C++ compiler |

### Install Components
```bash
cmake --install build --component Library    # lib + headers
cmake --install build --component Python     # pyopenmmpol
cmake --install build --component UtilityApp # ommp_xyz2mmp, ommp_pp
```

---

## 4. Conda
Install basic dependecies in e new environment:
```bash
conda create -n ommp python=3.xx numpy>2 cmake gcc gxx gfortran \
                     make lapack liblapack hdf5 zlib openssl \
                     libcurl pybind11 openmp libgomp
conda activate ommp
pip install build
```
Install cJSON library:
```bash
git clone https://github.com/DaveGamble/cJSON.git
cd cJSON
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make -j
make install
```

```bash
git clone https://github.com/Molecolab-Pisa/OpenMMPol.git
cd OpenMMPol
cp config.cmake custom.cmake   # customize options
cmake -C custom.cmake -B build \
      -DCMAKE_Fortran_COMPILER=`which gfortran`\
      -DCMAKE_C_COMPILER=`which gcc` \
      -DCMAKE_CXX_COMPILER=`which g++` \
      -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

cmake --build build -j
cmake --install build
```
Finaly install python package:
```bash
cd build
make PythonPackage
pip install pythonapi/dist/pyopenmmpol*.whl
```

## 5. The Smart Input JSON Interface

The **Smart Input (SI) JSON** file is the primary and prefereed mechanism to setup a calculation using OpenMMPol. It specifies input files, simulation parameters, QM region, and optional features — all in one human-readable file and computer-writable file.

### Complete JSON Schema

```json
{
    "name": "string",                  // Optional: calculation name
    "description": "string",           // Optional: description
    "version": "string",               // Optional: required OMMP version (semver)
    "verbosity": "string",             // [none, low, high, debug]
    "solver": "string",                // [default, cg, conjugate gradient, inversion, diis]
    "matrix_vector": "string",         // [default, direct, incore]
    "vdw_cutoff": number,              // Number, in Å (default: 14.0)
    "polarization_ls_conv_thr": number,  // Number, positive (default: 1d-8)
    "polarization_ls_use_guess": "true", // ["true", "false"] (default: true)

    // === INPUT/OUTPUT FILES ===
    // Choose ONE of: xyz_file+prm_file, mmpol_file, or hdf5_file

    "xyz_file": {
        "path": "string",              // Required: Tinker XYZ/ARC file
        "mode": "read",                // Always "read" for input files
        "md5sum": "string"             // Optional: MD5 checksum for validation
    },

    "prm_file": {
        "path": "string",              // Required if using xyz_file
        "mode": "read",
        "md5sum": "string"
    },

    "mmpol_file": {
        "path": "string",              // Required: MMP file (self-contained)
        "mode": "read",
        "md5sum": "string"
    },

    "hdf5_file": {
        "path": "string",              // HDF5 checkpoint file (WITH_HDF5)
        "mode": "read",
        "md5sum": "string"
    },

    "output_file": {
        "path": "string",              // Output log file path
        "mode": "write"                // Always "write" for output
    },

    // === QM REGION ===
    // EITHER specify qm_atoms + qm_coords OR leave absent for pure-MM

    "qm": {
        "qm_atoms": [                  // Array of element strings: "C", "N", "O", ...
            "C", "C", "N", "O", "H", "H"
        ],
        "qm_coords": [                 // Array of [x, y, z] in Å (3 × natoms)
            [1.234, 2.345, 3.456],
            [4.567, 5.678, 6.789]
        ],
        "qm_atom_types": [             // Optional: Tinker atom type integers
            133, 127, 129, 128
        ],
        "prm_file": {                  // Optional: .prm for VdW (required if qm_atom_types set)
            "path": "path/to/forcefield.prm",
            "mode": "read"
        },
        "qm_frozen_atoms": [           // Optional: subset of QM atoms to freeze (0-indexed)
            0, 1
        ]
    },

    // === GLOBAL SETTINGS ===
    "frozen_atoms": [                  // Array: atom indices (1-based) to freeze
        1, 2, 3, 10, 11, 12
    ],

    "remove_pol": [                    // Array: atom indices (1-based) to turn off polarizability
        5, 6
    ],

    // === LINK ATOMS ===
    "link_atoms": [                    // Array of link atom objects (requires qm section)
        {
            "MM_id": 400,            // MM atom index (1-based, bonded to QM)
            "QM_id": 15,             // QM atom index (1-based, at boundary)
            "LA_id": 16,             // Link atom index (1-based, new atom)
            "bond_length": 1.0,      // Optional: link atom distance from QM in Å (default: 1.0)
            "eel_remove": 0          // Optional: order of neighbour that will have electrostatic terms removed (default: 2)
        }
    ],

    // === FMM (Fast Multipole Method) ===
    "use_fmm": "string",               // [true, false]
    "fmm_max_l": 3,                    // Optional: FMM maximum multipole order (default: 8)
    "fmm_pol_max_l": 2,                // Optional: FMM for polarization max order (default: 8)
    "fmm_min_cell_radius": 14.0,       // Optional: min cell radius in Å (default: 7.0)
    "fmm_distance_thr": 14.0,          // Optional: FMM transition distance in Å (default: 5.0)

    // === DENSITY FITTING ===
    "density_fit": {
        "charge_points": {
            "type": "string",          // "atoms", "fibonacci", or "cubic"
            "n_pts_per_atom": 1,       // Optional: points per atom (fibonacci only)
            "radius": 0.0,             // Optional: radius (Å, cubic and fibonacci only)
            "source": "string"         // Optional: "qm" or "mm" (default: "qm")
        },
        "fit_points": {
            "type": "string",          // "atoms" or "cubic"
            "n_pts_per_atom": 1,       // Optional: points per atom
            "radius": 0.0,             // Optional: radius (Å, cubic and fibonacci only)
            "source": "string"         // Optional: "qm" or "mm" (default: "mm")
        }
    },

    // === PARAMETER LOOSENESS ===
    "ignore_duplicated_prm": "string", // [true, false]: ignore dup angle+opb
    "ignore_duplicated_angle_prm": "string", // [true, false]: ignore dup angles
    "ignore_duplicated_opb_prm": "string"    // [true, false]: ignore dup OOP
}
```

### Field Details

#### Top-Level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | No | Human-readable name, printed to log |
| `description` | string | No | Description, printed to log |
| `version` | string | No | Required OMMP version in semver format (e.g. `1.1.4`). If current build version is older, initialization fails |
| `verbosity` | string | No | Log verbosity: `none` (no output), `low` (minimal), `high` (default), `debug` (verbose) |
| `solver` | string | No | Induced dipole solver: `default` (CG), `cg`/`conjugate gradient`, `inversion` (direct solve), `diis` (DIIS extrapolation) |
| `matrix_vector` | string | No | MatVec strategy: `default` (auto), `direct` (compute on fly), `incore` (store matrix) |
| `vdw_cutoff` | number | No | Lennard-Jones cutoff in Å (default: no-cutoff). Stored internally in a.u. |
| `polarization_ls_conv_thr` | number | No | Convergence threshold for the polarization iterative solvers (CG, DIIS). Positive number (default: 1d-8) |
| `polarization_ls_use_guess` | string | No | Use the previous induced dipole as initial guess: `"true"` or `"false"` (default: true) |
| `output_file` | object | No | Log output file path with `mode: "write"` |

#### Input File Fields

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `xyz_file` | object | See below | Requires `prm_file`. Tinker XYZ/ARC format. |
| `prm_file` | object | Conditional | Required when using `xyz_file`. Tinker force-field format. |
| `mmpol_file` | object | See below | Self-contained system state. Alternative to xyz+prm. |
| `hdf5_file` | object | See below | HDF5 checkpoint. Requires `WITH_HDF5`. Alternative to mmp. |

**Mutual exclusivity:** Exactly one of `xyz_file+prm_file`, `mmpol_file`, or `hdf5_file` must be specified.

#### QM Section

The `qm` block defines the QM region. **All fields within `qm` must be present together** — `qm_atoms` with `qm_coords`, optionally `qm_atom_types`, `prm_file` (for VdW), and `qm_frozen_atoms`.

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `qm_atoms` | string[] | Yes (with coords) | Element symbols: `"C"`, `"N"`, `"O"`, `"H"`, etc. |
| `qm_coords` | number[][] | Yes (with atoms) | 3×N matrix `[x, y, z]` in Ångström |
| `qm_atom_types` | int[] | No | Tinker atom type integers. Requires `prm_file` for VdW. |
| `prm_file` | object | Conditional | Required if `qm_atom_types` is set; provides VdW parameters |
| `qm_frozen_atoms` | int[] | No | Indices (0-based) of QM atoms to freeze (non-polarizable) |

#### Frozen Atoms vs Remove Pol

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `frozen_atoms` | int[] | No | Atoms whose coordinates are fixed during gradients. 1-based indices |
| `remove_pol` | int[] | No | Atoms with polarizability turned off (become rigid). 1-based indices |

#### Link Atoms

Each entry in `link_atoms` defines one link atom at a QM/MM boundary:

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `MM_id` | int | Yes | 1-based index of MM atom (donor, bonded to QM boundary atom) |
| `QM_id` | int | Yes | 1-based index of QM boundary atom |
| `LA_id` | int | Yes | 1-based index of the link atom to be inserted |
| `bond_length` | float | No | Distance QM→LA in Å (default: 1.0 Å) |
| `eel_remove` | int | No | Electrostatic interaction terms to remove (default: 2) |

#### FMM Settings

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `use_fmm` | string | No | `"true"` or `"false"`. Only needed to override default |
| `fmm_max_l` | int | No | Maximum multipole expansion order for total electrostatics (default: 8) |
| `fmm_pol_max_l` | int | No | Max order for polarization (default: 8) |
| `fmm_min_cell_radius` | float | No | Minimum cell radius in Å (default: 7.0) |
| `fmm_distance_thr` | float | No | Transition distance between direct/FMM in Å (default: 5.0) |

#### Density Fitting

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `charge_points.type` | string | Yes | Grid type for electrostatic charge points: `"atoms"` (atom centers), `"fibonacci"` (spherical grid), `"cubic"` (cubic grid) |
| `charge_points.n_pts_per_atom` | int | Conditional | Points per atom for fibonacci (default: 1) |
| `charge_points.radius` | float | Conditional | Grid radius in Å for cubic (default: 0.0 = auto) |
| `charge_points.source` | string | No | Source topology: `"qm"` (default) or `"mm"` |
| `fit_points.type` | string | Yes | Grid type for fitting: `"atoms"` or `"cubic"` |
| `fit_points.n_pts_per_atom` | int | Conditional | Points per atom |
| `fit_points.radius` | float | Conditional | Radius in Å for cubic |
| `fit_points.source` | string | No | Source topology: `"mm"` (default) or `"qm"` |

#### Parameter Looseness

| Field | Type | Required | Notes |
|-------|------|----------|-------|
| `ignore_duplicated_prm` | string | No | `"true"` ignores both angle and OOP duplicates |
| `ignore_duplicated_angle_prm` | string | No | `"true"` ignores duplicated angle parameters |
| `ignore_duplicated_opb_prm` | string | No | `"true"` ignores duplicated out-of-plane parameters |

---

### JSON Examples

#### Minimal .mmp loading 
```json
{
    "mmpol_file": {
        "path": "system.mmp",
        "md5sum": "abc123..."
    },
    "verbosity": "high"
}
```

#### QM/MM with QM from JSON Arrays
```json
{
    "xyz_file": {"path": "mm.arc"},
    "prm_file": {"path": "amoeba09.prm"},
    "qm": {
        "qm_atoms": ["C", "N", "O"],
        "qm_coords": [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [2.2, 0.0, 0.0]],
        "qm_atom_types": [133, 133, 129],
        "prm_file": {"path": "amoeba09.prm"}
    },
    "verbosity": "low"
}
```

---

## 6. Input Preparation with OMMP-IP

**OMMP-IP** (`ommp-ip`) is the recommended tool for preparing OpenMMPol inputs. It takes Tinker XYZ (+ a PDB of the exact same system just to simplify selection with residue information etc), separates QM from MM, detects link atoms, and produces all files needed by the library.

### Invocation
```bash
ommp-ip -i system.xyz [-p system.pdb] [-d forcefield.prm] \
        [-q "QM_selection"] [-f "frozen_selection"] [-o out]
```

| Flag | Required | Description |
|------|----------|-------------|
| `-i <xyz>` | Yes | Combined Tinker XYZ/ARC (QM + MM) |
| `-p <pdb>` | No | PDB for residue/atom name metadata |
| `-d <prm>` | No | Tinker force-field file |
| `-q <sel>` | No | MDAnalysis selection for QM atoms (default: `not all`) |
| `-f <sel>` | No | MDAnalysis selection for frozen atoms (default: `not all`) |
| `-o <base>` | No | Output basename (default: `out`) |

### What It Produces

| File | Content |
|------|---------|
| `<base>_mm.arc` | MM atoms in Tinker XYZ (ARC format) |
| `<base>_qm.xyz` | QM atoms (+ link atoms) in Tinker XYZ |
| `<base>_si.json` | Smart Input JSON — ready to run |
| `<base>_mm.pdb` | MM atoms (PDB format) |
| `<base>_frozen.pdb` | Frozen atoms (PDB format) |

### Examples

**Ligand as QM, protein backbone frozen:**
```bash
ommp-ip -i complex.xyz -p complex.pdb \
  -d amoebabio18.prm -q "resname LIG" \
  -f "(backbone or name O*) and not resname LIG" \
  -o ligand_calc
```

**Simple small molecule in water:**
```bash
ommp-ip -i mol_water.xyz -d amoeba09.prm \
  -q "resname MOL" -o mol_calc
```

### Link Atom Detection

OMMP-IP scans QM-MM bonds and places a hydrogen link atom 1.0 Å from the QM atom along the bond direction. It guesses the Tinker atom type from other H atoms on the QM atom; if none found, type is set to 0 (requires manual correction).

### MDAnalysis Selection Strings

```
"resname HIS"                          # by residue name
"resname HIS and not name H*"          # exclude atoms
"name CA"                              # by atom name
"resid 45 to 50"                       # by residue number
"element C"                            # by element
"(resname ALA or resname GLY) and backbone"  # compound
"not resname SOL"                      # complement
```
For anything more complex refer to the [official documentation](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).

---

## 7. The MMP File Format

> **Note:** The **MMP** file is an old legacy format for QM/MMPol and QM/AMOEBA calculation. It's supperted by OpenMMPol with several limitation (using this input format neither bonded terms nor VdW are computed, linka toms are not available etc.) **it should be used only for legacy calculations!**  

`.mmp` is a text format for storing the electrostatic part of the forcefield for a system: coordinates, charges/multipoles, polarizabilities, connectivity, and solver settings. It is a compact input — everything needed for a polarizable electrostatics calculation is in one file and it is compatible with the internal `g16-molecolab` Gaussian version. It can be generated using using some internal tools or from a Tinker `xyz` and `prm` using the utility 



> **Note:** The MMP file stores only the electrostatic system (no bonded terms, no QM data). Bonded parameters are loaded separately from `.prm` files when building the full system. For QM/MM workflows, use OMMP-IP to generate MMP + JSON together.

### Version
The MMP file starts with a **revision number**: `2` or `3`. Version 3 includes additional fields (PCM parameters). Both are fully supported.

### Format Layout

#### Header Block (system-wide parameters)

| # | Field | Type | Used in OpenMMPol? | Description |
|---|:-----:|:----:|:-----------:|:----------|
| 1 | `version` | int | ✅ | MMP revision (2 or 3) |
| 2 | `job_type` | int | ❌ | 0=QM/MM, 1=Tinker QM/MM, 2=EET |
| 3 | `verbosity` | int | ❌ | Verbosity flag |
| 4 | `ff_type` | int | ✅ | 0=AMBER-like, 1=AMOEBA-like |
| 5 | `ff_subtype` | int | ✅ | AMBER: 0=Wang-AL, 1=Wang-DL, 2=Thole; AMOEBA: 0 |
| 6 | `disabled` | int | ❌ | Disabled flag |
| 7 | `damping` | real | ❌ | Dipole-dipole damping parameter (Å) |
| 8 | `solver` | int | ❌ | 0=default, 1=matrix inversion, 2=Jacobi/DIIS, 3=CG |
| 9 | `matvec` | int | ❌ | 0=default, 1=incore, 2=O(N²), 3=FMM large, 4=FMM all |
| 10 | `conv_thresh` | int | ❌ | Convergence threshold 10^-N for iterative solvers |
| 11 | `fmm_accuracy` | int | ❌ | FMM accuracy (<0=10⁻⁵, 0=10⁻⁶, 1=10⁻⁷) |
| 12 | `fmm_box_mm` | real | ❌ | FMM box size for MM (Å) |
| 13 | `fmm_box_pcm` | real | ❌ | FMM box size for MM/PCM (Å, v3 only) |
| 14 | `solv_model` | int | ❌ | 0=none, 1=ddCOSMO, 2=ddPCM |
| 15 | `pcm_lmax` | int | ❌ | Max angular momentum for PCM |
| 16 | `pcm_leb_points` | int | ❌ | Number of Lebedev integration points |
| 17 | `pcm_conv` | int | ❌ | PCM convergence threshold 10^-N |
| 18 | `eps_solvent` | real | ❌ | Dielectric constant of solvent |
| 19 | `eps_optical` | real | ❌ | Optical dielectric constant (v3 only) |
| 20 | `switch_region` | real | ❌ | Switching region size for PCM |
| 21 | `cavity_type` | int | ❌ | 0=default, 1=SAS/UFF, 2=from-file, 3=VdW/Bondi, 4=SAS/Bondi |
| 22 | `probe_radius` | real | ❌ | Probe radius for SAS cavity (Å) |
| 23 | `n_spheres` | int | ❌ | Number of cavity spheres (0 if auto) |
| 24 | `natoms` | int | ✅ | Number of MM atoms |

#### Per-Atom Blocks

| # | Block | Count | Format | Used? | Description |
|---|:-----:|:-----:|:------:|:-----:|:----------|
| 25 | `atomic_numbers` | N | int[N] | ❌ | Atomic numbers (read for counting only, not stored) |
| 26 | `coordinates` | N | real[N][3] | ✅ | Coordinates in **Å** (converted to Bohr internally) |
| 27 | `residue_numbers` | N | int[N] | ❌ | Residue numbers (read for structure, not stored) |
| 28 | `charges/multipoles` | N | real[N] or real[N][10] | ✅ | Charges (AMBER) or multipoles (AMOEBA) |
| 29 | `polarizabilities` | N | real[N] | ✅ | Polarizabilities in Bohr³ (0.0 for non-polarizable) |
| 30 | `connectivity` | N | int[N][8] | ✅ | 1-2 neighbor connectivity (max 8 neighbors) |
| 31 | `pol_groups` | N | int[N][120] | ✅ (ff_type=1) | Polarization group members (read only when `ff_type=1`, i.e. AMOEBA) |
| 32 | `rot_frames` | N | int[N][4] | ✅ (ff_type=1) | Multipole rotation frame: mol_frame, iz, ix, iy (read only when `ff_type=1`, i.e. AMOEBA) |

### AMOEBA Multipoles (Block 28)
When `ff_type=1` (AMOEBA), each atom stores 10 values:
`[q, μx, μy, μz, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz]` — monopole, dipole (3), quadrupole (5), all in atomic units.

### Key Conventions
- Coordinates stored in **Å** in the file (converted to Bohr internally)
- Polarizabilities in **Å³** (converted to A.U. internally )
- Atom indices are **1-based** in connectivity blocks
- Non-polarizable atoms have polarizability = 0.0 (**not omitted**)
- Zero-filled arrays pad unused neighbors/group members

### Creating an MMP File

**From XYZ + force-field:**: `ommp_xyz2mmp input.xyz forcefield.prm output.mmp`

**From within a software using OpenMMPol**: use `ommp_save_mmp(sys, filename, version)`.

---

## 8. Running a Calculation

Running a calculation on OpenMMPol alone is possible but does make sense just for test pourposes. To this end there are several simple programs in C and Fortran in `tests/test_programs`. 

Since it is a library meant to be interfaced with QM software, it is normally used from within a **host QM software**. In QM software interface, the path to the appropriate json file is provided and the QM/MMPol calculation is run from there.

---

## 9. Troubleshooting

### "File does not correspond to md5sum"
The file content doesn't match the `md5sum` in JSON. Either update the checksum or remove it.

### "No input for MM system found"
Specify one input source: `xyz_file+prm_file`, `mmpol_file`, or `hdf5_file`.

### Induced dipoles not converging
- Check for overlapping atoms
- Verify polarizability values (especially AMOEBA multipoles)
- Switch to `diis` solver
- Run with `verbosity: "debug"` for detailed output

### Link atom type is 0
OMMP-IP couldn't guess the Tinker atom type. Edit the JSON manually: set `qm_atom_types` to the correct values, or re-run OMMP-IP with a force-field `.prm` file that has known atom types.

### MMP file errors
- Verify atom count matches between blocks
- Check that all atom indices in connectivity blocks are valid (1 to NATOMS)
- Ensure the file ends after the last expected block

---

## Further Reading

- **FORD API docs:** https://molecolab-pisa.github.io/OpenMMPol
- **Tinker format:** https://dasher.wustl.edu/tinker/
