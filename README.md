# Elltorque

Results and figures used in [Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle](https://doi.org/10.1093/gji/ggaa166).

## Prerequisites

Installed texlive, python/python3 with matplotlib >v2.1 for support of latest colormaps. A working Julia >v1.3.


## Install

Open the repository directory and run `julia install_local.jl` to install the package.

## Run

To run the calculations and the plots you simply run

```julia
using Elltorque
Elltorque.run(true)
```
from within the Julia REPL (takes around 2-3h). To run without calculating the data use

```julia
Elltorque.run(false)
```

## Citation
If you use this data and/or software please cite the article [Pressure torque of torsional Alfvén modes acting on an ellipsoidal mantle](https://doi.org/10.1093/gji/ggaa166).
