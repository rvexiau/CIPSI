#!/bin/bash
#modifier pour faire des raccourcis claviers

Remote=ada

echo "Command to execute : ?"
echo "m = molpro // f = fetch results // c = cipsi // r = recompile "
read prog

if [ "$prog" == "m" ] ; then
  name=molpro.ll
  rsync -az --progress Molpro/data/ ergon:Molpro/
  rsync -az --progress Molpro/molpro.ll adapp:Molpro/molpro.ll
  ssh -t -t $Remote << EOL1
  cd Molpro
  llsubmit $name
  llq -u rxbl004
  exit
EOL1
fi

if [ "$prog" == "r" ] ; then
   echo "Machine : ?"
   echo "l = Local // i = Idris // u = lumat"
   read Task
  
  if [ "$Task" == "l" ] ; then
   echo "use lumat to compile"
  fi
  
  if [ "$Task" == "u" ] ; then
   rsync -az $CIPSI_ROOT/compile.sh lumat:cipsi/     
   rsync -az $CIPSI_ROOT/Autocip_13/* lumat:cipsi/Autocip_13/   
   rsync -az $CIPSI_ROOT/CIPSI_sources/* lumat:cipsi/CIPSI_sources/   
   ssh -t -t lumat << EOL2
   cd cipsi
   ./compile.sh
   exit
EOL2
  fi
  
  if [ "$Task" == "i" ] ; then
   rm -r $CIPSI_ROOT/idris
   rsync -az $CIPSI_ROOT/Autocip_13 $CIPSI_ROOT/idris/
   rsync -az $CIPSI_ROOT/CIPSI_sources/* $CIPSI_ROOT/idris/CIPSI_sources/
   rsync $CIPSI_ROOT/compile.sh $CIPSI_ROOT/idris/ 
   cd $CIPSI_ROOT/idris
   rm cipsi.tar.gz
   rm -f ergon:cipsi.tar.gz
   tar -zcf cipsi.tar.gz *
   rsync $CIPSI_ROOT/idris/cipsi.tar.gz ergon:
   rm -r $CIPSI_ROOT/idris   
   ssh -t -t ada << EOL3
     mfget cipsi.tar.gz "\$WORKDIR"
     cd "\$WORKDIR"
     tar -zxf cipsi.tar.gz -C cipsi/      
     cd cipsi
     ./compile.sh
     exit
EOL3
  fi
fi

if [ "$prog" == "f" ] ; then
   echo "Machine : ?"
   echo "i = Idris // u = lumat"
   read Task
   
   if [ "$Task" == "u" ] ; then
     result=`sed -n 4p input/auto.in`
     ssh -t -t lumat << EOL4
       cd /data/vexiau/
       tar -zcf lumat.tar.gz $result
       exit   
EOL4

#     rsync -az lumat:results/Cipsi/$result/ ./$result  
     rsync -az --progress lumat:/data/vexiau/lumat.tar.gz ./  
     tar -zxf lumat.tar.gz
     rm lumat.tar.gz
     ssh -t -t lumat << EOL7
       cd /data/vexiau/
       rm -r $result
       rm lumat.tar.gz
       exit   
EOL7
  fi
  if [ "$Task" == "i" ] ; then
     result=`sed -n 4p input/auto.in`
     rsync -az --progress ergon:output/$result/ ./$result
     cd $result
     tar -zxf idris.tar.gz
     ssh -t -t ergon << EOL5
       cd output
       rm -r $result
       exit   
EOL5
  fi
fi


if [ "$prog" == "c" ] ; then
   echo "Machine : ?"
   echo "l = Local // i = Idris // u = lumat"
   read Task
  
  if [ "$Task" == "l" ] ; then
   result=`sed -n 4p input/auto.in`  
   rsync -az input/* $result
   rsync -az $CIPSI_ROOT/bin/Autocip13 $result
   rsync $CIPSI_ROOT/bin/data/def_grid.dat $result
   cd $result
   export OMP_NUM_THREADS=4
#   valgrind --leak-check=yes --track-origins=yes --suppressions=/home/romain/valgrind.supp --log-file="valgrind.out" ./Autocip12<auto.in 
   ./Autocip13<auto.in>cipsi.out   
  fi
  
  if [ "$Task" == "u" ] ; then
   result=`sed -n 4p input/auto.in`
   rm -rf lumat:results/Cipsi/$result
#   rsync input/* lumat:results/Cipsi/$result   
#   rsync $CIPSI_ROOT/bin/cipsi.sh lumat:results/Cipsi/$result
#   rsync $CIPSI_ROOT/bin/def_grid.dat lumat:results/Cipsi/$result  
   rsync input/* lumat:/data/vexiau/$result   
   rsync $CIPSI_ROOT/bin/cipsi.sh lumat:/data/vexiau/$result 
   rsync $CIPSI_ROOT/bin/data/def_grid.dat lumat:/data/vexiau/$result 
   ssh -t -t lumat << EOL6
   cd /data/vexiau/$result
   rsync ~/cipsi/bin/Autocip13 ./Autocip13
   sbatch cipsi.sh
   squeue
   exit
EOL6

  fi
  
  if [ "$Task" == "i" ] ; then
    rsync -az input/ ergon:input/
    rsync -az $CIPSI_ROOT/bin/data/def_grid.dat ergon:input/
    rsync -az $CIPSI_ROOT/bin/cipsi.ll ada:.
    ssh -t -t $Remote << EOL7
    llsubmit cipsi.ll
    llq -u rgdw003
    exit
EOL7
  fi
fi
