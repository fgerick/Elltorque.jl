# Elltorque

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fgerick.github.io/Elltorque.jl/dev) -->
[![Build Status](https://travis-ci.com/fgerick/Elltorque.jl.svg?token=NJNkFC9qALxxCxMBhjwi&branch=master)](https://travis-ci.com/fgerick/Elltorque.jl)

Reproducing results of [Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle](https://doi.org/10.1093/gji/ggaa166).

## Prerequisits

Installed working texlive, python/python3 with matplotlib >v2.1 for support of latest colormaps. A working Julia >v1.0 without linked MKL library for Arpack.jl support.


## Install

<!-- Go to Julia REPL and run `]add https://github.com/fgerick/Elltorque.jl.git`. -->
Open the repository directory and run `julia install_local.jl` to install the package.

## Run

To run the calculations and the plots you simply run

```julia
using Elltorque
Elltorque.run(true)
```

To run without calculating use

```julia
Elltorque.run(false)
```
