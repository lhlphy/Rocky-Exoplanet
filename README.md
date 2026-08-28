This project can calculate phase curve of exoplanet in different surface conditions, such as specular and lambert surface.

Phase curve can be calculated by simulation or a complex integration or a simple optical principle, but only the simulation and original integration are accurate.

./lava_lib : Lava albedo data used for Fig.1

./lib : main numerical code

./temp
    Lambert_0.1 : diffuse model with As=An=0.1
    Specular_0.1 : specular model with An=0.1
    Specular_0.2 : specular model with An=0.2
    Specular_As0.1 : specular model with As=0.1
    Specular_As0.2 : specular model with As=0.2
    Specular_As0.4 : specular model with As=0.4
    NonFresnel_flat : specular model with non-Fresnel reflection, validation of the flat phase curve

params.txt : configuration file for the simulation
PS.csv : planetary system parameters from NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/)
run.sh : run script
