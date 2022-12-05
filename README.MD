To run the MPI Implementation:
cd src
module load openmpi
mpiCC -std=c++11 MPI.cpp -o mpi -lm
sbatch MPI_scene<X>.script //X represents the scene you want rendered

To run the OMP Implementation

cd src
srun --partition=express -N 1 -n 1 --cpus-per-task=12 --pty --time=01:00:00 /bin/bash 
g++ -std=c++11 -fopenmp OMP.cpp -o omp 
sbatch OMP_scene<X>.script //X represents the scene you want rendered

To run the Single Threaded Implementation

cd src
srun --partition=short -N 1 -n 1 --pty --time=02:00:00 /bin/bash 
g++ -std=c++11 singleThread.cpp -o sing 
./sing