#! /bin/bash
#SBATCH -J FH # name of the job
#SBATCH -o FH.%J.out #output file %J=job number
#SBATCH -e FH.%J.err #output error file; %J=job number
#SBATCH -n 8
#SBATCH -p main --qos main
#SBATCH --mem-per-cpu 8G

cat > FH.pdos.in << EOF

&PROJWFC
    prefix        = 'CTS',
    outdir        = './outdir',
    degauss       =  0.001,
    deltae        =  0.001,
    emax          =  20,
    emin          = -15,
    filpdos        = 'FH.pdos.dat',
/

EOF

# Load the modules that we need.
#module load fftw3/gcc5.4.0
module load hdf5/gcc/1.10.4
module load mpi/openmpi/gcc/3.1.1
module load physical/espresso/6.1

# Needed for intel mpi to work.
#export I_MPI_PMI_LIBRARY="/lib64/libpmi.so"

export ESPRESSO_PSEUDO='/stuperm/sbudhathoki/QE/Heusler'

srun projwfc.x < FH.pdos.in > FH.pdos.out

wait $!
echo "========================="
echo "Job complete!"
echo "$(date)"