% Some notes on using the code
% Jānis Erdmanis
% November 18 2016

# About

  Magnetic droplet surface evolution boundary integral algorithms (BIE or BEM) according to magnetosatics and Stokes equations in general 3D case are presented in this repository.  

# Requirements

  The code is written in `julia 0.4` and should work fine with latter versions. The code depends on `JLD.jl` which stores and loads meshes and simulations; `SurfaceGeometry.jl` which deals with mesh generation, stabilisation and others which are used for processing simulation results. 

Assuming you are `Ubuntu` user install required dependencies `sudo apt-get install -y libcgal-dev liblapack-dev libblas-dev`, `sudo apt-get install csh hdf5-tools python-matplolib cmake` and if you want to use Distmesh mesh you also need Matlab. In julia REPL execute following lines
```
Pkg.init()
Pkg.add("JLD")
Pkg.add("ArgParse")
Pkg.add("FastGaussQuadrature")
Pkg.clone("https://github.com/akels/SurfaceGeometry.jl")
Pkg.build("SurfaceGeometry")
```
Since packages might break features in the future I provide backup in `libs.zip`. Then instalation is just unziping this folder in julia package directory. 

Optional dependencies which are needed to view and process simulation data are installed with `sudo apt-get install python-matplotlib` and in julia REPL
```
Pkg.add("ThreeJS")
Pkg.add("Escher")
Pkg.add("PyPlot")
Pkg.add("Elliptic")
```

# Running simulation

  For running simulation execute in MDrop directory `julia -L .juliarc.jl mdrop.jl --simulation=FastField --Bm=20 --mu=10` which will run `FastField.jl` and sore results in `~/SimulationData`. Instead of passing command line arguments you can run simulation also as `julia -L .juliarc.jl -L config.jl Fastfield.jl`.
  
# For citing this work refer to article

```
Erdmanis, J. & Kitenbergs, G. & Perzynski, R. & Cebers, A. (2017)
Magnetic droplet in rotating field: numerical simulation and comparison with experiment
```  
