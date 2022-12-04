To run the MPI Implementation:

module load openmpi
mpiCC -std=c++11 MPI.cpp -o mpi -lm
sbatch MPI.script

To run the OMP Implementation
g++ -fopenmp OMP.cpp -o omp 
sbatch OMP.script