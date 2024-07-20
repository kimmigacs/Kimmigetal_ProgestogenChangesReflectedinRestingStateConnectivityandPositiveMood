### **Analysis data and code for Kimmig et al.: "To start or to discontinue the pill – changes in progestogens are reflected by resting-state connectivity and positive mood"**

<br>**Boxplots**: R code and anonymized source data provided. Running the scripts create and save the boxplot figures provided in the manuscript.

*The boxplot figures in the manuscript were created using R (version 4.3.2) with the following libraries: ggplot2 (version 3.4.4.), ggpubr (version 0.6.0), ggprism (version 1.0.5), glue (version 1.7.0), exactRankTests (version  0.8-35, for indications of significance).* 

<br>**IS-RSAs**: Three Jupyter Notebook files are provided for the different types of intersubject RSAs presented in the manuscript, including all output.txt files: 
1) hormonal change by parcel-wise resting state functional connectivity
2) hormonal change by behavioral change
3) parcel-wise resting state functional connectivity by behavioral change

Note: Step-wise use of the Jupyter Notebooks is required to select the specific dataframes (e.g. progestogens, estrogens, positive mood, negative mood, parcelwise connectivity etc.) for the respective IS-RSA. 

**Example 1: IS-RSA between changes in progestogen levels and parcel-wise** (RSA_withDeltaBrainandHormones_ManuscriptFinalAnonymized.ipynb)
- after importing packages and loading 'Parcel-wise Resting State Connectivity dataframe', skip 'Estrogen Dataframe' and directly proceed with running 'Progestogen Dataframe'
- Run all subsequent cell until 'Save results of Estrogen-RSFC RSA', skip this cell and proceed with 'Save results of Progestogen-RSFC RSA' instead
- Results will be saved in a textfile called: outputdelta_progestogens_seuclidean_twotailedFWEcorrected.txt

**Example 2: IS-RSA between changes in positive mood and parcel-wise** (RSA_withDeltaBrainandBehavior_ManuscriptFinalAnonymized.ipynb)
-  after importing packages and loading 'Parcel-wise Resting State Connectivity dataframe', run 'Load Itemwise positive mood data' and skip all other cells loading itemwise behavioral dataframes until
-  Run all subsequent cell starting from 'Drop participants with missing behavioral data from both dataframes' until 'Save results of Positive Mood-RSFC RSA', skip all remaining cells for saving data of the other behavioral measures
-  Results will be saved in a textfile called: outputdelta_posmood_seuclidean_twotailedFWEcorrected.txt

*IS-RSA scripts are dependent on the following python packages, versions used in our work are denoted in the brackets: pandas (version 1.2.4), numpy (version 1.20.1), math, sklearn (version 0.24.1) for metrics and pairwise distances, and scipy (version 1.8.0) for correlation of distance matrices and permutation testing. Representational Distance Matrices were plotted using seaborn (version 0.11.1) and numpy (version 1.20.1). The scripts were ran in Jupyter Notebook (version 6.3.0) with Python (version 3.8.8).*
