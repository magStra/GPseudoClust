# GPseudoClust
Code for the GPseudoClust method, including pre- and postprocessing steps. 

Magdalena E Strauß, Paul DW Kirk, John E Reid, Lorenz Wernisch (2019); 
GPseudoClust: deconvolution of shared pseudo-trajectories at single-cell resolution.

Author of the code: ME Strauss

A tutorial is provided in .mlx and .pdf format in the files GPseudoClustByExample.mlx and GPseudoClustByExample.pdf (
inside GPseudoClust folder).

The clustering method can be run without downloading additional software. 

However, for postprocessing the following is required:

https://github.com/magStra/nonparametricSummaryPSM

For plotting we use the subaxis function by Aslak Grinsted:

https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot


The folder lmkk_summaryMatrixRepresentation contains additional methods for postprocessing, which use the following software, which requires 
separate download:

1) Code implementing the localised kernel k-means method available at https://github.com/mehmetgonen/lmkkmeans,

Gönen, M. and Margolin, A.A. (2014). Localized data fusion for kernel k-means clustering with application to cancer biology. 
In Advances in Neural Information Processing Systems 27, pages 1305-1313.

2) The Mosek optimisation software (https://www.mosek.com/).

3) The SIMLR software, where one of the functions is used for the estimation of the optimal number of clusters for the summary 
clustering in the post-processing. 

Wang, B. et al. (2017). Visualization and analysis of single-cell RNA-seq data by
kernel-based similarity learning. Nat Meth, 14, 414-416.

https://github.com/BatzoglouLabSU/SIMLR


