 #!/usr/bin/bash 
clear

# to clear current plot folder to get rid of any old plots make clear_flag = 1
#---------------------------------------------------------------------

clear_flag=1
if [ "$clear_flag" -eq 1 ]; then
rm ../plot/*.png
rm ../plot/*.eps
rm -r ../plot/
rm -r ../data/
rm ../output/out*
echo "old results are deleted"
fi

#---------------------------------------------------------------------
# create folders for data, plot, output text files, image files, results

mkdir -p ../data 
mkdir -p ../plot
mkdir -p ../output
mkdir -p ../im
mkdir -p ../result

echo "run matlab code"
#---------------------------------------------------------------------
# run analysis over each channel in a loop

# number of channels in im/ directory
ndata=$(<ndata)
echo 'number of tif data files ->' $ndata

ls ../im -1 | sed -e 's/\.tif$//' > filenames

#--------------------------- DATA LOOP -------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for j in $(seq 1 1 $ndata) 	
do
dataname=`sed "${j}q;d" filenames`
printf "%s" $dataname>'dataname'

echo 'running for data ->' ${dataname}

time matlab -nodesktop -nodisplay < SAM.m > ../output/out_${dataname}


# ----------------------------organize plots--------------------------
cd ../plot
mkdir -p "plots_${dataname}"
mv *_${dataname}*.* plots_${dataname}/
cd ../code


done
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#--------- save code in compressed form ------------

#tar -czvf compCode.tar.gz code *.sh *.m filenames
#mv compCode.tar.gz ../plot/compCode.tar.gz








