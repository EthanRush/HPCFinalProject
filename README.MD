Scene options:

1/Default: Shirley Orbs

2: Cornell Box

3: Shierly Box


To run the MPI Implementation:
>cd src

>module load openmpi

>mpiCC -std=c++11 MPI.cpp -o mpi -lm

>sbatch MPI_sceneX.script //X represents the scene you want rendered

To run the OMP Implementation

>cd src

>g++ -std=c++11 -fopenmp OMP.cpp -o omp 

>sbatch OMP_sceneX.script //X represents the scene you want rendered

To run the Single Threaded Implementation

>cd src

>g++ -std=c++11 singleThread.cpp -o sing 

>sbatch SING_sceneX.script //X represents the scene you want rendered
