# Statistical Analysis of the association between miRNAs in breast milk, perinatal probiotic supplement and the development of atopic dermatitis


This is a Github repository for the code used in Lene Tillerli Omdal's master thesis at NTNU. Below follws a breif description of the code files. The "testing" folder was mostly used for exploring functions and methods, thus we do not provide a description of these files. The "Code" folder contains all code used to perform the anlalysis and the "Data" folder contains data stored along the way of the analysis. The miRNA breast milk data that was analyzed is not provided here.


## Code file explanations

1. *load-count-data.R*: This R file contains code for reading the data, check for missing values and some cleaning of the data.
2. *normalize-and-filter.R*: This R file contains code for performinf pre-processing of the data. This includes TMM and cpm normalization, filtering, and joining of identical sequences.
3. *heatmaps.R*: This R file containes code for performing hierachical clustering and creating heatmaps
4. *voom-analysis.R*: This R file containes code for perfomring the voom analysis.
5. *heatmap-top-mirna.R*: This R filse contaies code for performing hierarcical clutering of the 20 top ranked miRNAs by the voom analysis
6. *cut-trees.R*: This R file includes code for exploring how samples and miRNAs were clustered.
7. *enet-mod-functions.R*: This R file containes the code for the nested CV and repeated CV functions for determining the model parameters of an elastic net model. Additionally, it containes funtions for performing bootstraping of both of these functions.
8. *elasticnet-model.R*: This R file containes the code for running nested CV and repeated CV on the original data and for performing the bootstraping.
9. *bca.R* This R file containes the code for calculationg the modified bias correction factor, jackknife estimtes and bias-corrected accelerated confidence intervals of the coefficients of the elasticnet model.

## Project thesis abstract

The Probiotics in the Prevention of Allergy among Children in Trondheim (ProPACT) study showed a $40\%$ reduction in the risk of developing atopic dermatitis in children whose mothers received a probiotic supplement before and whilst breastfeeding, compared to a placebo alternative. The biological explanation for this risk reduction has not yet been fully understood. In this thesis, we analyze if microRNAs in breast milk, $10$ days postpartum, are possible contributors to the risk reduction.
We perform two analyses; one in which we examine the effect of the probiotic supplement on microRNAs and one in which we examine whether any microRNAs are associated with development of atopic dermatitis by $2$ years. In addition, we perform a clustering analysis to explore patterns and groupings in the data. 

The analysis in this thesis is based on data from $60$ mother-child pairs which were semi-randomly selected from the ProPACT study.
We use differential expression analysis to investigate if the probiotic supplement has an effect on the expression values of individual microRNAs. 
To perform this analysis we employ the statistical method \textit{voom}.
To investigate whether any microRNAs are associated with the development of atopic dermatitis, we perform variable selection using an elastic net model.
To estimate confidence intervals of the coefficients, we employ bootstrapping and the bias-corrected accelerated method.
In the exploratory clustering analysis we use hierarchical clustering with euclidean and correlation based dissimilarity measures.


The probiotic supplement is associated with differential expression of one microRNA, $miR\text{-}577$, when taking into account the multiplicity of tests by controlling the false discovery rate at $10\%$ using the Benjamini and hochberg method.
In total, $47$ microRNAs have a raw $p$-value below $0.05$ but except for $miR\text{-}577$ none of them have an acceptable false discovery rate adjusted $p$-value. 
The five microRNAs, $miR$-$342\text{-}3p, miR\text{-}3605\text{-}3p,\\miR\text{-}500a/b\text{-}5p, miR\text{-}625\text{-}3p$ and $miR$-$6515\text{-}5p$, are associated with development of atopic dermatitis. 
The microRNA $miR$-$3605\text{-}3p$ is also one of the microRNAs with a raw $p$-value below $0.05$ in the analysis of the effect of probiotics.
The cluster analysis groups $miR$-$3605\text{-}3p$, $miR$-$6515\text{-}5p$ and $miR$-$577$ together, indicating that they are highly correlated. 
Thus, the probiotic supplement affects one microRNA, and five microRNAs are found to be associated with atopic dermatitis in breast milk at $10$ days postpartum.
However, we found no conclusive evidence that probiotics affect the same microRNAs that are associated with atopic dermatitis. 
Some were grouped together in the cluster analysis and may be highly correlated, and further studies may consider focusing on the microRNAs that were most promising.

