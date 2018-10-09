% Some notes on using the code
% Jānis Erdmanis
% November 18 2016

# About

  Magnetic droplet surface evolution boundary integral algorithms (BIE or BEM) according to magnetostatics and Stokes equations in general 3D case are presented in this repository.  

# SetUP

The code originaly was written with `julia 0.4` but now essentials are ported for `julia 0.7` and will run on `julia 1.0` when warnings will be adressed. It is now required that LINUX is being used for running simulation which is due to a binary dependency ElTopo.

To run the code first repoitory needs to be cloned whch will be in `MDrop` folder. Then in `julia 0.7` one executes a following commands in the REPL (assuming that `.` is in repository folder)
```
]activate .
resolve
```
which will install all the necessary dependencies for the code. To proceed then with a calculation one then executes a following commands:
```
include(".juliarc.jl")
include("config.jl")
include("FastField.jl")
```

# Requirements OLD

  The code is written in `julia 0.4` and tested extensively  on `Ubuntu 16.04`. Main dependencies of the code are `JLD.jl` which stores and loads meshes and simulations; `SurfaceGeometry.jl` which deals with mesh generation, stabilization and others which are used for processing simulation results; 

Assuming you are `Ubuntu` install required dependencies with `sudo apt-get install -y libcgal-dev liblapack-dev libblas-dev`, `sudo apt-get install csh hdf5-tools cmake` and if you want to use Distmesh for mesh generation you also need Matlab. In julia REPL execute following lines
```
Pkg.init()
Pkg.add("JLD")
Pkg.add("ArgParse")
Pkg.add("FastGaussQuadrature")
Pkg.clone("https://github.com/akels/SurfaceGeometry.jl")
Pkg.build("SurfaceGeometry")
```
Since packages might break features in the future I provide backup in `libs.zip`. Then installation is just unzipping this folder in julia package directory. 

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
