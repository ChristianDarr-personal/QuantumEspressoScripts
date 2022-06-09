#! /bin/bash
#SBATCH -J FH # name of the job
#SBATCH -o FH.%J.out #output file %J=job number
#SBATCH -e FH.%J.err #output error file; %J=job number
#SBATCH -n 8
#SBATCH -p main --qos main
#SBATCH --mem-per-cpu 8G


# Load the modules that we need.
#module load fftw3/gcc5.4.0
module load hdf5/gcc/1.10.4
module load mpi/openmpi/gcc/3.1.1
module load physical/espresso/6.1

# Needed for intel mpi to work.
#export I_MPI_PMI_LIBRARY="/lib64/libpmi.so"

export ESPRESSO_PSEUDO='/stuperm/sbudhathoki/QE/Heusler'

srun sumpdos.x *\(Co\)* > Co_tot.dat
srun sumpdos.x *\(Ti\)* > Ti_tot.dat
srun sumpdos.x *\(Sn\)* > Sn_tot.dat
srun sumpdos.x *\(Co\)*\(s\) > Co_s.dat
srun sumpdos.x *\(Co\)*\(p\) > Co_p.dat
srun sumpdos.x *\(Co\)*\(d\) > Co_d.dat
srun sumpdos.x *\(Ti\)*\(s\) > Ti_s.dat
srun sumpdos.x *\(Ti\)*\(p\) > Ti_p.dat
srun sumpdos.x *\(Ti\)*\(d\) > Ti_d.dat
srun sumpdos.x *\(Sn\)*\(s\) > Sn_s.dat
srun sumpdos.x *\(Sn\)*\(p\) > Sn_p.dat
srun sumpdos.x *\(Sn\)*\(d\) > Sn_d.dat


wait $!
echo "========================="
echo "Job complete!"
echo "$(date)"