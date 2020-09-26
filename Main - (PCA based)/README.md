# HRTF synthesis based on anthropometric data 

These are the main files used to build process the SOFA HRTFs, generate the PCA sub-space, train the regression model and rebuild HRTFs to their raw form.

* *The main routines are ordered in the recommended workflow.*


## Anthropometry_Datasets.m

Select the morphological data used to train the regression model and remove subjects with incomplete or nonexistent measurements.


## Preprocess_Datasets.m

Process all the SOFA HRTF to have the same signal properties, spatial resolution and remove the subjecs without anthropometric data across all the different dataset. 


## PCA_DTF.m

Principal Component Analysis over the DTFs, to reduce the number of parameters necessary to the regression model. 


## NeuralNet_Datasets.m

Neural Network regression model to determine the relationship between the anthropometry and the low-dimensional DTFs.


## Rebuild_HRTF.m 

Example of the reconstruction procedure for a subject removed from neural net training processing and some interesting graphs.


### Error_Analysis.m

Spectral and time domain error evaluated for the exemple case defined at the previous function. 


## Simulation_Erro.m

Generates individualized HRTFs for the HUTUBs dataset and performs analysis of the spectral and time domain features in comparisson with the simulated (BEM) data in the dataset. The same error evaluation was done between generic HRTFs and the database. 


## Model_test.m

Virtuallly the same as Rebuild_HRTF.m, but you can write your own anthropometry. Can be used to test a new generated model. *Becareful to input the correct data at the right place - the use of the Individualization App is highly recomended instead of this routine.* 
