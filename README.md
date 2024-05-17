# Pan-cancer proteogenomics signatures associated with HRD, MSI, APOBEC, and smoking.

This is a repo with the scripts for CPTAC pan-cancer proteogenomics signatures project. Also see here for other CPTAC pan-cancer papers: 

  * https://doi.org/10.1016/j.cell.2023.07.014

  * https://doi.org/10.1016/j.cell.2023.07.013

  * https://doi.org/10.1016/j.ccell.2023.06.009


# Analysis scripts:

1. ```Define_High_Low_groups``` -- this directory contains scripts that are used to define "High" and "Low" groups of samples for each mutational signature based on that signature mutational count and fraction. 


2. ```Expression_markers``` -- this directory contains scripts to run differential expression markers between "High" and "Low" samples groups using cancer cohort as a covariate.

   + The analysis was done using pan-cancer harmonised dataset, that was also used in the 1st and 2nd CPTAC papers listed above.

   + In the folder, there are three scripts for markers based on: gene, protein and phosphosite expression.


3. ```Proteogenomic_markers_Figs_5_6``` -- this directory contains analysis and plotting scripts for proteogenomic markers associated wih tht HRD, MSI, APOBEC, and smoking cancer phenotypes.

   + ```Proteogenomic_markers_Figs_5_6/Volcano_plots``` -- scripts for plotting protein and phosphosite expression markers. Top markers are highlighted in those plots.

   + ```Proteogenomic_markers_Figs_5_6/Heatmap_plots``` -- scripts for plotting top gene, protein and phosphosite expression, that are also used for subtyping using ```ConsensusClusterPlus``` R package. For each of the phenotype/cancer pair, there are two steps/scripts. For example, for HRD-BRCA pair the scripts should be used in the following order:

     * ```Proteogenomic_markers_Figs_5_6/Heatmap_plots/Generate_matrix.BRCA.20240105.R``` -- script to select top markers, perform subtyping based on those, and also creating the matrix and annotations needed for plotting.


     * ```Proteogenomic_markers_Figs_5_6/Heatmap_plots/Make_heatmap_BRCA_HRD.20240105.R``` -- script to make a heatmap plot with ```ComplexHeatmap``` R package, using the results from the previous step.