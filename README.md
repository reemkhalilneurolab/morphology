Overview
========
This code accompanies the paper "TOPOLOGICAL SHOLL DESCRIPTORS FOR NEURONAL CLUSTERING AND CLASSIFICATION"
Data is available https://doi.org/10.5281/zenodo.6118983


Requierments
============
1. MATLAB release R2020a
2. Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Installation
============

1. Extract the ZIP file or clone the git repository in a directory of your choice on your local machine. 
2. Place your swc files in the morphology/data folder 

        ├── morphology/data
        │   ├── class1 (folder with swc files for class1)
        │   ├── class2 
            .
            .
        │   └── classn
     
   To run the swc files included in the paper download the data folder from https://doi.org/10.5281/zenodo.6118983
   navigate to data/set(n)/swc_files and copy all the folders into the morphology/data folder
   
3. Run the main.m script morphology/main.m file from the editor tab or by excuting run('main') from your command window.
   The morphology/config.m file should add and set all the necessary paths and the script should start executing. 
  
4. To manually add the `morphology/` folder to your path in MATLAB :
    - use the "Set Path" dialog in MATLAB, or 
    - by running the addpath(pwd) function from your command window or `startup` script.

5- All the output matrices and the images will be generated in morphology/save


