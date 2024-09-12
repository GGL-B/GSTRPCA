
GSTRPCA

Title: Irregular tensor singular value decomposition for single-cell multi-omics data clustering
Description: GSTRPCA is written in the MATLAB programming language. To use, please download the GSTRPCA folder and follow the instructions provided in the README.


To operate:
	* for the example of GSTRPCA:
		we have provided some pre-data for GSTRPCA, if you want to run GSTRPCA framwork, please run the script "GSTRPCA_Demo" directly.      

   *  Some of the data in the paper is included in the "Data" file.

Files：

GSTRPCA_Demo.m - The main function.

prox_l1.m - Soft threshold operator for computing irregular sparse tensors.

Imsvd1.m - t-SVD calculation of regular low-rank tensors.(fill with zeros)

t_SVT2.m - GSVD calculation of regular low-rank tensors.(fill with zeros)

t_SVT1.m - GSVD calculation of irregular low-rank tensors.

TSN1.m - The iterative process of GSTRPCA algorithm.

SS.m - Downstream clustering experiments.

10X_inhouse.mat - A real single-cell multi-omics data used in the cell type clustering example. The dataset has been preprocessed.

compute_NMI.m - A program designed to compute the Normalized Mutual Information (NMI) between two sets of clusterings.

AMI.m - A program designed to compute the Adjusted Mutual Information (AMI) between two sets of clusterings.

ARI.m - A program designed to compute the Adjusted Rand Index (ARI) between two sets of clusterings.


