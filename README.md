# Elltorque

Results and figures used in [Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle](https://doi.org/10.1093/gji/ggaa166).

## Prerequisits

Installed working texlive, python/python3 with matplotlib >v2.1 for support of latest colormaps. A working Julia >v1.3.


## Install

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

## Citation
If you use this software please cite the article [Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle](https://doi.org/10.1093/gji/ggaa166).
