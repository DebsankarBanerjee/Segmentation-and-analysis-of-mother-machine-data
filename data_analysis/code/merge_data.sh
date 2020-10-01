

a=49
for i in ../data/set2/cell_id_ch_*.txt; do
  echo $i
  printf -v j "%03d" $a			 #03 pad to length of 3
  mv ${i} ../data/set2/cell_id_ch_${j}.txt
  let a=a+1
done


a=49
for i in ../data/set2/cellone_time_ch_*.txt; do
  echo $i
  printf -v j "%03d" $a			 #03 pad to length of 3
  mv ${i} ../data/set2/cellone_time_ch_${j}.txt
  let a=a+1
done

