# Ariadne Set

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/set/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/set/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/set/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/set/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/set/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/set/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/set/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/set/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/set/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/set)

Ariadne Set is the standalone C++20 set-representation layer of Ariadne. It implements constrained and validated image sets, affine constrained image sets, curves and interpolated curves, and list sets. Grid pavings and paver algorithms belong to the **downstream Paving layer**, not Set.

## Layout

- Public headers: `include/set/`, with umbrella header `include/ariadne-set.hpp`.
- C++ implementations: `src/`.
- C++ tests: `tests/`.
- Python binding, tests and tutorial: `python/bindings/`, `python/tests/` and `python/tutorials/`.
- Standalone installed C++ example: `tutorials/tutorial_set/`.

The source tree retains the historical `polyhedron.cpp`, `polytope.cpp`, and `zonotope.cpp` implementation files, but they were not compiled by the original Kernel Set object target and are **not enabled** in this extraction.

## Dependencies

Set depends directly on:

- [ariadne-cps/solving](https://github.com/ariadne-cps/solving) (Git submodule `ss/solving`).
- [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration) (Git submodule `ss/configuration`).

Solving provides Function, Geometry and their lower dependencies transitively. Set does not depend on Kernel or Paving. The architecture is `geometry <- function <- solving <- set <- paving`.

The private LP state header under `include/set/detail/` allows the affine separation algorithm to build without requiring an unpublished header from Kernel.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/set.git
cd set
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

For installation and Python bindings:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-set bindings
cmake --install build-release --prefix /path/to/install
```

## Python bindings

With Python 3.11 or newer and a matching development environment, the build contributes the Set bindings to the `pyariadne` extension via the lower Solving, Function and Geometry modules. The extension installs under the chosen prefix.

## Tutorials

- C++: `tutorials/tutorial_set/tutorial_set.cpp`, built separately against an installed `ariadne-set`.
- Python: `python/tutorials/tutorial_set.py`.

## Contribution guidelines

Please contact the developer: Luca Geretti <luca.geretti@univr.it>.

## License

Ariadne Set is released under the GNU General Public License v3.0 (or later where individual source files specify it).
