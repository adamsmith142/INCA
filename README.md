# INCA
INCA (INverse method for ConcAvity index) Code to calculate the concavity index using inverse methodology.

This code calculates the concavity index for a region based on the methodology of Smith and Fox (In review). The methodology uses an inversion to calculate channel steepness index (see Smith et al. 2022, https://doi.org/10.1016/j.earscirev.2022.103970) at different values of the concavity index. Calculating the misfit between the predicted channel, and the actual channel, allows assessment of which concavity index is correct. The DEM of the clearwater catchment can be downloaded from the following repository - Smith, Adam (2024), “Clearwater DEM”, Mendeley Data, V1, doi: 10.17632/pbppkd7gwn.1. The concavity index of this catchment is given by the minimum of the curve of misfit vs concavity index. In this code, the concavity index is measured as 0.52, however as the spacing between concavity index values is 0.02 it could be between 0.51 and 0.53. 

This code takes approximately 45 minutes to run in its present form. The code can be made to run faster by using a smaller river network. 
