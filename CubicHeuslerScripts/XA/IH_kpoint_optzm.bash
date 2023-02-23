#! /bin/bash
#SBATCH -J IH # name of the job
#SBATCH -o IH.%J.out #output file %J=job number
#SBATCH -e IH.%J.err #output error file; %J=job number
#SBATCH -n 8
#SBATCH -p main --qos main
#SBATCH --mem-per-cpu 8G

rm -f IH_scf.out Inverse.fcc
touch Inverse.fcc
for k in 6 8 10 12 14 16 18 20
do
cat > IH_scf.in << EOF

&CONTROL
  calculation   = 'scf',
  prefix        = 'X2YZ',               !**Please use name of your choice**!
  restart_mode  = 'from_scratch',
  verbosity     = 'high',
  nstep         = 120, 
  outdir        = './outdir',
  forc_conv_thr = 1.00000e-05,
  etot_conv_thr = 1.00000e-05, 
/
&SYSTEM
  ibrav       = 2,                           !**FCC crystal**!
  celldm(1)   = 10.56,                       !use sensible lattice parameter!
  nat         = 4,                           !**NUMBER OF ATOMS**!
  ntyp        = 3,                           !**NUMBER OF TYPES OF ATOMS**!
  ecutwfc     = 80,                          !**First test convergence**!
  ecutrho     = 800,                         !**Default = ecutrho ~ 4*ecutwfc**!
  occupations = 'smearing',
  smearing    = 'gaussian',
  degauss     = 0.002,
  nspin       = 2,
  starting_magnetization(1)= 0.01,           !**check for magnetization of each atom**!
  starting_magnetization(2)= 0.01,           !**check for magnetization of each atom**!
  starting_magnetization(3)= 0.01,           !**check for magnetization of each atom**!
/
&ELECTRONS
  electron_maxstep  = 200,
  mixing_mode       = 'plain',               !**charge density Broyden mixing**!
  conv_thr          = 1.00000e-07,
  mixing_beta       = 7.00000e-01,
  diagonalization   = 'cg',
/
&IONS
  ion_dynamics = 'bfgs'                      !**use BFGS quasi-newton algorithm**!
/
&CELL
  cell_dynamics='bfgs'
/
ATOMIC_SPECIES
X   51.9961 Cr.pbe-spn-kjpaw_psl.1.0.0.UPF    !**Symbol, atomic weight, and potential for X-element**!
Y   55.845  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF    !**Symbol, atomic weight, and potential for Y-element**!
Z   28.085  Si.pbe-n-kjpaw_psl.1.0.0.UPF      !**Symbol, atomic weight, and potential for Z-element**!

ATOMIC_POSITIONS {alat}
X     0.000000000         0.000000000         0.000000000      !X_1  A and please use symbol for X-element!
X     0.250000000         0.250000000         0.250000000      !X_2  B and please use symbol for X-element!
Y     0.500000000         0.500000000         0.500000000      !Y_1  C and please use symbol for Y-element!
Z     0.750000000         0.750000000         0.750000000      !Z_1  D and please use symbol for Z-element!

K_POINTS {automatic}
$k $k $k  0 0 0

EOF

# Load the modules that we need.
module load fftw3/gcc5.4.0
module load hdf5/gcc/1.10.4
module load mpi/openmpi/gcc/3.1.1
module load physical/espresso/6.1

# Needed for intel mpi to work.
export I_MPI_PMI_LIBRARY="/lib64/libpmi.so"

#path to direcotry where you have potential files
export ESPRESSO_PSEUDO='/stuperm/sbudhathoki/QE/X2YZ'   

srun pw.x < IH_scf.in > IH_scf.out

# extract Etot from output
etot=`grep -e ! IH_scf.out | awk '{print $(NF-1)}'`
echo $k $etot >> Inverse.fcc
done
cat Inverse.fcc

wait $!
echo "========================="
echo "Job complete!"
echo "$(date)"

