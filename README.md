# FastScape.jl
[![CI](https://github.com/boriskaus/FastScape.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/boriskaus/FastScape.jl/actions/workflows/CI.yml)

Julia interface to [FastScape](https://fastscape.org) - a flexible and modular landscape evolution model developed at GFZ Potsdam.

This translates nearly all routines of the fortran version of the library. As of now, you have to compile this library yourself and place it in the main directory of repository (but a BinaryBuilder version of FastScape seems faeasible).

Examples of how to use it are in the `/test` directory.
