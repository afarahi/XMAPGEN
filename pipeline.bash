#!/bin/bash

echo "Start running pipeline"

FILES=/scratch/physdoeproject_fluxoe/aryaf/Chandra/halo_catalogs/*.fit
#REALIZATION_ID=0
realization_id=5

for fname in $FILES
do
  echo "Processing ${fname:59:100} file..."
  # take action on each file. $f store current file name
  #python main.py 1 ${fname:59:100} 
  #python main.py 2 ${fname:59:100} 
  python main.py 3 ${fname:59:100} 
  #rm Catalog/Output_File/${fname:59:100}
  exit 0
done

#for IIND in {35..36}
#do
#   #echo $IIND
#   FNAME=$(printf 'XCAT_Aardvark_1.0_%i.fit' $IIND)
#   JOBONE=`qsub -v fname=$FNAME,realization_id=$REALIZATION_ID PBSmap`
#done
#
#exit 0



