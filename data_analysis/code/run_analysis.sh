 #!/usr/bin/bash 
clear
# clear current plot folder------------------------------------------


# create folders fro data, plot etc..
 
mkdir -p ../../result
mkdir -p ../plot

cp ../../code/ndata ndata

echo "run analysis codes to get cell division statistics"
#---------------------------------------------------------------------
# run analysis over tif files in a loop

# dataname from filenames
ndata=$(<ndata)
#printf "%s" $ndata>'ndata'
echo 'no of tif data ->' $ndata


#--------------------------- DATA LOOP -------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# get the filenames for the structured cell data
cd ../../data/
ls cell_id* -1 | sed -e 's/\.txt$//' > datafilenames
mv datafilenames ../data_analysis/code/datafilenames
cd ../data_analysis/code/

for j in $(seq 1 1 $ndata) 	
do
dataname=`sed "${j}q;d" datafilenames`
printf "%s" $dataname>'dataname'

echo 'running for data ->' ${dataname}

time matlab -nodesktop -nodisplay < ft.m > ftout

done


#--------------------------- DATA LOOP -------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# dataname from filenames
ndata=$(<ndata)
#printf "%s" $ndata>'ndata'
echo 'no of tif data ->' $ndata

# get the filenames for the structured cell data
cd ../../data/
ls cellone* -1 | sed -e 's/\.txt$//' > datafilenames2
mv datafilenames2 ../data_analysis/code/datafilenames2
cd ../data_analysis/code/

for j in $(seq 1 1 $ndata) 	
do
dataname=`sed "${j}q;d" datafilenames2`
printf "%s" $dataname>'dataname'

echo 'running for data ->' ${dataname}

time matlab -nodesktop -nodisplay < growth_rate.m > grout




done



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#--------------------- PROCESS WHOLE DATA SET OF ALL CHANNELS------------------------
cd ../../result

cat alltau_*.txt > alltau.txt

cat mdd_tau_*.txt > allmdd_tau.txt
cat mdd_lb_*.txt > allmdd_lb.txt
cat mdd_ld_*.txt > allmdd_ld.txt

cat mdc_tau_*.txt > allmdc_tau.txt
cat mdc_lb_*.txt > allmdc_lb.txt
cat mdc_ld_*.txt > allmdc_ld.txt

cat allgrowth_*.txt > allgrowth.txt

cd ../data_analysis/code


time matlab -nodesktop -nodisplay < ftalldata.m > allout


time matlab -nodesktop -nodisplay < gen_correlation.m > genout









