# Summer Project Code Repository
Skimmer flow simulation code for summer project with smf group
MATLAB script for interfacing with Bird DS2V software and custom DSMC code in MATLAB and MEX C++

# prepare_data.m 
Transfers data from Bird DS2V software output DS2FF.DAT to MATLAB array

# plot_data.m
Plots centerline flux through disk of skimmer radius from skimmer tip (set as x = 0 in the simulation) to the exit

# Bird.m
Prepares DS2VD.DAT file for simulation in Bird DS2V software with focus on the skimmer tip 

# simulation.m and .cpp files
Custom DSMC simulation, used along with the stlread file

# stlread.m and license.txt
stlread function written by Eric Johnson and used to import external mesh of skimmer surface in custom DSMC simulation
