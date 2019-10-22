# Description
Analysis of seismic data from the Piton de la Fournaise in the Réunion Island. The data was given by the Piton de la Fournaise Volcano Observatory and the project was directed by Alexandrine Gesret.

# Explanation
The purpose of this project is to find the epicenter and time of seisms based on the measures made by 20 stations.

For each station, the Eikonal uses a basic speed model in the ground to compute the time taken to go to the station from any point of the map.

The best epicenter and time is found by a complete optimisation: a L2 and L1 minimisation of the difference between the observations and the predicted arrivals of S and P waves to the stations.

# Usage
- Modify the main.f90 file to update the solver.
- Run `gfortran main.f90 -o main.exe` to compile it
- Run `./main.exe` to run it

`Données sismiques.ipynb` contains the analysis
