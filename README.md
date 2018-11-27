# P. vivax model
This is the C++ source code for the _Plasmodium vivax_ transmission model developed by Dr. Michael WHITE (<mwhite@pasteur.fr>).

## What's in this repo ?
The actual model code is found in the pv\_mod folder, in the [Source.cpp](../blob/master/Pv_mod/Source.cpp) file. Numerical optimisation routines are found in the same directory.
The root of this repository contains example of model parameter files, as well as some R code to visualize model outputs.

## Build the model
This section provides details about how to build the model on different operating systems.
### On your laptop
Under UNIX/Linux, use GCC 7+ (I know for sure it'll fail with GCC 4 because of syntaxical errors with '\<' and '\>' signs).

````bash
g++ -O3 -o Pv_model Source.cpp com.cpp linpack.cpp randlib.cpp
```

Under Windows... No idea. To be detailed.

### On Pasteur TARS cluster
You need to load GCC into your environment as it's not available by default.
```bash
module load gcc/7.2.0
g++ -O3 -o Pv_model Source.cpp com.cpp linpack.cpp randlib.cpp
```
## Use the model

TODOs: make a much more elaborated README file
