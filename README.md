# Segmentation and analysis of mother machine data (SAM)
Bacterial cell segmentation in mother machine data using MATLAB
Console version of SAM
# To run in MATLAB console:
	make four folders named "code", "data", "plot" and "im"
	put all files except the .tif ones inside /code
	Open MATLAB in /code folder path
1. all .tif files should be in /im folder

2. write the filenames in the file code/dataname
   example: if the files in "im" folder are xy01_01.tif and xy01_02.tif then the dataname file should be

xy01_01

xy01_02

3. run SAM.m in the MATLAB console

4. the folders "data" and "plot" should exist for writing data and plotting frames.

5. the data file cell_id_* will contain the main data in the format:
cell_id_array : current_index / time@Birth / time@Division / length@Birth / length@Division / parent ID / index at division / position at division
    
