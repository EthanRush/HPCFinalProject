To run the MPI Implementation:
cd src
module load openmpi
mpiCC -std=c++11 MPI.cpp -o mpi -lm
sbatch MPI.script

To run the OMP Implementation

cd src
srun -p --partition=short -N 1 -n 1 -cpuspertask=12 --pty --time=02:00:00 /bin/bash 
g++ -std=c++11 -fopenmp OMP.cpp -o omp 
./omp

To run the Single Threaded Implementation

cd src
srun -p --partition=short -N 1 -n 1 -cpuspertask=12 --pty --time=02:00:00 /bin/bash 
g++ -std=c++11 singleThread.cpp -o sing 
./sing