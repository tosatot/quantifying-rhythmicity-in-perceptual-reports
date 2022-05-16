In this repository there is the code relative to the paper:

"Quantifying rhythmicity in perceptual reports" 
by Tommaso Tosato, Gustavo Rohenkohl, Jarrod Robert Dowdall, Pascal Fries


The repository contains the following functions:

simulation.m:  guides the user through the generation, analysis, statistics, and visualization of a simulated dataset. 
genData.m:     generates data
atcDFT.m:      apply the DFT to the mean accuracy time course
stLSS.m:       apply the LSS to the mean accuracy time course
sinFIT.m:      fit a sine (and a dampened oscillation) to the mean accuracy time course
randEff.m:     perform a random effect statistical test
pairedTTEst.m: obtain the t-values of a paired t test 
calcFDR.m:     perform an FDR multiple comparison correction



Dependencies:

'Matlab R2018b' 
'Signal Processing Toolbox'
'Statistics and Machine Learning Toolbox'
'Curve Fitting Toolbox'
'Circular Statistics Toolbox'

the Circular Statistics Toolbox can be freely downloaded here:
https://se.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics




Author: Tommaso Tosato
tommaso.tosato@esi-frankfurt.de
