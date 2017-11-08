# Pas réussi à intégrer le .f90 dans le makefile

# Compilation en (lapack) parallèle
#echo 'ifort -openmp -I/opt/intel/composer_xe_2013.2.146/mkl/include -c diagonaliser.f90'
#ifort -openmp -I/opt/intel/composer_xe_2013.2.146/mkl/include -c diagonaliser.f90

# Compilation en (lapack) séquentiel
echo 'ifort -I$MKLROOT/include -c diagonaliser.f90'
ifort -I$MKLROOT/include -c diagonaliser.f90

make

echo '###############################################################################'
echo 'Attention à copier le bin résultant dans le bon dossier --->   cp ../bin/bd ../../bin/'
echo '###############################################################################'