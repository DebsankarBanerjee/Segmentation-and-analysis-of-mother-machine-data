# Segmentation and analysis of mother machine data (SAM)
Bacterial cell segmentation in mother machine data using MATLAB\
Console version of SAM
# To run in MATLAB console:
	Make four folders named "code", "data", "plot" and "im"
	Put all files except the .tif ones inside /code
	Put the .tif files inside /im
	Open MATLAB in /code folder path
	Run the code >SAM
# For new data files:
	All .tif files should be in /im folder
	Write the filenames in the file code/dataname for example: if the files in "im" folder are xy01_01.tif and xy01_02.tif then the dataname file should be
	
	xy01_01
	xy01_02
	
	run SAM.m in the MATLAB console
# Results:
	The folder "plot" contains frames with segmentation
	The folder "data" contains the data files cell_id_*
	Data structure in cell_id_* 
	array : current_index // time@Birth // time@Division // length@Birth // length@Division // parent ID // index at division // position at division
