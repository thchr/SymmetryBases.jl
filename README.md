# SymmetryBases.jl

This package provides access to computation of Hilbert bases associated with the space of band structure symmetries (see `compatibility_bases`, `nontopological_bases`, and `split_fragiletrivial`), using Normaliz (which must be installed separately).

The package additionally provides a number of utilities to easily check the topology of a symmetry vector, using elementary band representations accessed from [Crystalline.jl](https://github.com/thchr/Crystalline.jl) (see `calc_detailed_topology`, `calc_topology`, `isbandstruct`, and `indicators`).

## License

The [Normaliz](https://github.com/Normaliz/Normaliz) library as well as its Python bindings, [PyNormaliz](https://github.com/Normaliz/PyNormaliz), are required to run this software and must be installed separately.
Note that Normaliz and PyNormaliz are licensed under the GPLv3 (see their license files), but the bindings to the library in this package, SymmetryBases.jl, as well as any other source code in this package, are licensed under the MIT License.
This means that code using the Normaliz library via the SymmetryBases.jl bindings is subject to Normaliz's licensing terms. If you distribute a derived or combined work, i.e. a program that links to and is distributed with the Normaliz library, then that distribution falls under the terms of the GPLv3. If you just distribute source code that links to SymmetryBases.jl, and users have to download Normaliz and PyNormaliz separately to provide a backend, then the GPL probably doesn't have much effect on you.