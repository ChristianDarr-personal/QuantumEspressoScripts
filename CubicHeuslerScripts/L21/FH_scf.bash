#! /bin/bash
#SBATCH -J FH # name of the job
#SBATCH -o FH.%J.out #output file %J=job number
#SBATCH -e FH.%J.err #output error file; %J=job number
#SBATCH -n 8
#SBATCH -p main --qos main
#SBATCH --mem-per-cpu 8G

rm -f FH_scf.out Full.fcc
touch Full.fcc
for k in 11.508432099
do
cat > FH_scf.in << EOF

&CONTROL
  calculation   = 'scf',
  prefix        = 'CTS',               
  restart_mode  = 'from_scratch',
  verbosity     = 'high',
  nstep         = 120, 
  outdir        = './outdir',
  forc_conv_thr = 1.00000e-04,
  etot_conv_thr = 1.00000e-04, 
/
&SYSTEM
  ibrav       = 2,                           
  celldm(1)   = $k,                       
  nat         = 4,                           
  ntyp        = 3,                           
  ecutwfc     = 200,                          
  ecutrho     = 1000, 
  nbnd        = 52,                        
  occupations = 'smearing',
  smearing    = 'gaussian',
  degauss     = 0.002,
  nspin       = 2,
  starting_magnetization(1)= 0.1,           
  starting_magnetization(2)= 0.1,           
  starting_magnetization(3)= 0.1,           
/
&ELECTRONS
  electron_maxstep  = 200,
  mixing_mode       = 'plain',               
  conv_thr          = 1.00000e-06,
  mixing_beta       = 7.00000e-01,
  diagonalization   = 'cg',
/
ATOMIC_SPECIES
Co   58.93319    Co.pbe-spn-kjpaw_psl.1.0.0.UPF    
Ti   47.86700    Ti.pbe-spn-kjpaw_psl.1.0.0.UPF    
Sn   118.7100    Sn.pbe-dn-kjpaw_psl.1.0.0.UPF      

ATOMIC_POSITIONS {alat}
Co     0.000000000         0.000000000         0.000000000      !X_1  A and please use symbol for X-element!
Co     0.500000000         0.500000000         0.500000000      !X_2  C and please use symbol for X-element!
Ti     0.250000000         0.250000000         0.250000000      !Y_1  B and please use symbol for Y-element!
Sn     0.750000000         0.750000000         0.750000000      !Z_1  D and please use symbol for Z-element!

K_POINTS {automatic}
20 20 20  0 0 0

EOF

# Load the modules that we need.
module load fftw3/gcc5.4.0
module load hdf5/gcc/1.10.4
module load mpi/openmpi/gcc/3.1.1
module load physical/espresso/6.1

# Needed for intel mpi to work.
export I_MPI_PMI_LIBRARY="/lib64/libpmi.so"

#path to direcotry where you have potential files
export ESPRESSO_PSEUDO='/stuperm/sbudhathoki/QE/Heusler'   

srun pw.x < FH_scf.in > FH_scf.out

# extract Etot from output
etot=`grep -e ! FH_scf.out | awk '{print $(NF-1)}'`
mtot=`grep 'total magnetization' FH_scf.out | tail -1 | awk '{print $(NF-2)}'`
echo $k $etot $mtot >> Full.fcc
done
cat Full.fcc

wait $!
echo "========================="
echo "Job complete!"
echo "$(date)"