# Uncomment the next 3 lines if not using modules
module load modules/0-traditional
module load gcc/9.3.0
module load matlab
IK=2
NU=10
NP=5
matlab -nodesktop -nosplash -nodisplay -r "run test1(${IK},${NU},${NP})"
