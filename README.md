
This repository provides all scripts and a detailed walk-through to repeat analyses shown in Figures 8 and 10 in "Loss of CTRP10 results in female obesity with preserved metabolic health" by Chen, Wong et al 2024


Analyses are broken down into 2 sections:

CTRP10 Differential expression analysis and Integration (Fig 8)

CTRP10 Sex-specific connectivities (Fig 10)

Each analysis should operate independently with the provided links to required files

## Link to tpm kallisto counts, for CTRP10 Differential expression analysis and Integration (Fig 8)
## https://drive.google.com/file/d/1RmYstO-Vo5ZvGd2rjgPZZ7p3MGeW3PmQ/view?usp=drive_link

Publicly available data used were acquired from (and pre processed for further analysis):
GTEx V8: https://gtexportal.org/home/datasets (the raw data were filtered where individuals were required to show counts > 0 in 1.2e6 gene_tissue combinations across all data. Given that our goal was to look across tissues at enrichments, this was done to limit spurious influence of genes only expressed in specific tissues. Post-filtering consists of 310 individuals and 1.8e7 gene_tissue combinations). This results in 52% NA values in the final matrix used.
## Link to pre processed GTEx data used in this study
## https://drive.google.com/file/d/16Wh-6nagVIBXwrFINCoFKBB5zuUqGFeL/view?usp=drive_link
### Link to metadata including sex datatable
## https://drive.google.com/file/d/13Yg1a0uXagMffX0fTimLUhGhSgw-04I1/view?usp=drive_link

MGI list of mouse-human orthologs: http://www.informatics.jax.org/homology.shtml
## Link to human/mouse orthologs table used in this study
## https://drive.google.com/file/d/1-6R13c546yWgHMlbunW5snDlkflZGicp/view?usp=drive_link

Additional Considerations
While the tissues shown in this study capture notable metabolic processes, all available tissues are included in the filtered dataset to enable exapnsion of the analysis. Specifically, these can be identified as columns in 'working_dataset'. The scripts can be easily repurposed to incorporate other tissues (for example, coronary artery)

We apologize for bulky and recursive portions of the scripts. Realizing that there and more efficent mechanisms of executing these analyses, we felt that breaking down step-by-step is helpful to facilitate learning.

All feedback is welcome, please email: lmvelez@uci.edu & mseldin@uci.edu
