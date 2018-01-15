#!/bin/sh
# @ job_name = cipsi
# @ output = $(step_name).$(jobid).stdout
# @ error = $(step_name).$(jobid).stdout
 
# @ step_name = get_file
# @ job_type = serial
# @ class = archive
# @ queue
 
# @ step_name = exec
# @ dependency = (get_file == 0)
# @ job_type = serial
# @ parallel_threads = 8
# @ as_limit=56gb
# @ wall_clock_limit=90:00:00
# @ queue
 
# @ step_name = put_file
# @ dependency = (exec >= 0)
# @ job_type = serial
# @ class = archive
# @ notify_user = romain.vexiau@u-psud.fr
# @ notification= complete    
# @ queue


case ${LOADL_STEP_NAME} in
  get_file )
    set -ex
    cd $TMPDIR
    mfget input/* cip_run/    
    cp $WORKDIR/cipsi/bin/Autocip13.exe cip_run/
  ;;
  exec )
    set -x
    module load idrmem
    cd $TMPDIR/cip_run
    date > idris.out
    ./Autocip13.exe<auto.in >> idris.out 2>>idris.out 
    date >> idris.out
    cat idris.out
  ;;
  put_file )
    set -x
    cd $TMPDIR/cip_run/
    find . -name "*_mat*" -type f -delete
    find . -name "*fort*" -type f -delete
    find . -name "*BD_SCRATCH*" -type f -delete
    find . -name "COMB_VCPP" -type f -delete
    find . -name "*_ijkl*" -type f -delete
    find . -name "*_f50*" -type f -delete
    find . -name "*_scr*" -type f -delete
    find . -name "*_pqrs*" -type f -delete   
    find . -name "*_f25*" -type f -delete   
    find . -name "*_det*" -type f -delete
    find . -name "*_bda*" -type f -delete
    find . -name "F10_*" -type f -delete
    find . -name "*_dcla*" -type f -delete
    find . -name "*_dspina*" -type f -delete       
    result=`sed -n 4p auto.in`
    tar -zcf idris.tar.gz *
    mfput ./idris.tar.gz output/$result/

  ;;
 
esac

