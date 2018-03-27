# GMNB

  GMNB (gamma Markov negative binomial) is proposed for differential expression analysis of temporal RNA-seq count data across different phenotypes or treatment conditions. The inputs of GMNB are RNA-seq count data taken at multiple time points under two conditions. It does not need for ad-hoc normalization step. The output provides Bayes Factor for each gene in the study to rank the significance of dynamic gene differential expression across conditions. 

```
Main function:
      GMNB.m

Depends:
      CRT.m
      get_NB_loglikelihood.m
```

## Running the tests

All tables and figures from the paper can be reproduced based on the code detailed as follows.

### Simulation study:

The presented results can be reproduced by the code in Sim Directory:

* **simNB.m:** To generate the synthetic RNA-seq count data for 1000 genes under two conditions according to the GMNB model with the gamma-Markov temporal dependencies between dispersion parameters.

* **simGP.m:** To generate the synthetic RNA-seq count data for 1000 genes under two conditions according to the DyNB model assumptions.

* **simAR.m:** To generate the synthetic RNA-seq count data for 1000 genes under two conditions according to the NB-AR(1) model assumptions.

* **Demo.m:** This demo shows 20 runs of GMNB for each generative model and plots both receiver operating characteristic (ROC) and precision-recall (PR) curves.

* **rocm.m:** Calculate the mean and variance of ROC and PR for 20 runs.

### Real-world RNA-seq Data:

Directory contains the code to generate the plots and tables for real-world RNA-seq data:

* **Demo.m:** This demo shows how GMNB can be implemented to analyze two real-world temporal RNA-seq data.

#### Data Directory:

* [GSE52260_rc_human_timeSeriesRNAseq.mat](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52260): Human Th17 cell induction dataset used in the paper. 

* [DyNB_example_data.mat](http://research.cs.aalto.fi/csb/software/dynb/DyNB_example_data.mat): RNA-seq data in Aijo et al. \[2014\]( Human-activated T and Th17 cells). This data has been used to compare the order of validated genes between GMNB and DyNB methods.

The figures were plotted by the code in Plot Directory:

* **plot_GMNB.m:** This file can be used to plot the normalized gene expression profile of different genes over time estimated by the proposed GMNB model. The input “k” for plotting the reported genes in the paper has been commented in this file.

* **plot_GMNB_rawCount.m:** This file can be used to plot the gene expression profile of different genes over time estimated by GMNB model based on the observed read count. The input “k” for plotting the reported genes in the paper has been commented in this file.

In addition, the following directories contain the implementation code for other existing methods for performance comparison:

**DyNB Directory:**
    
* DyNB.m: A MATLAB implementation of DyNB.
* DyNB_caller.m: An example script demonstrating the use of the DyNB.

**deseq2 Directory:**
    
* main_AR.R: An example R script demonstrating the use of DESeq2.
