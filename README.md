# IRS_normalization
An exploration of internal reference scaling (IRS) normalization in isobaric tagging proteomics experiments. Also, examples of how IRS-normalized data affects statistical testing, and how to avoid using ratios in the analyses. 

The IRS method was first described in this publication:
> Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

The analysis is of a mouse lens development time course (6 points 3 days apart from E15 to P9) where three replicates of the time points were done in 3 separate TMT labelings. The lens is a unique system that has been studied for many years and the prior knowledge can be used to guide some analysis steps. The data is from this publication:
> Khan, S.Y., Ali, M., Kabir, F., Renuse, S., Na, C.H., Talbot, C.C., Hackett, S.F. and Riazuddin, S.A., 2018. Proteome Profiling of Developing Murine Lens Through Mass Spectrometry. Investigative Ophthalmology & Visual Science, 59(1), pp.100-107.

## Contents:
Four [jupyter notebook](http://jupyter.org) files (R kernel). If you click on the notebook files (*.ipynb extensions), they will render and display in your bowser. Please be patient as they can take a minute to render.
* understanding_IRS.ipynb is Part 1 (normalizations)
* statistical_testing.ipynb is Part 2 (edgeR testing)
* statistical_testing_ratios.ipynb is Part 3 (taking ratios and using limma)
* statistical_testing_take2.ipynb is Part 4 (testing P0 vs P3)

Data from Kahn, et al.:
* iovs_58-13-55-s01.csv

Sample information for design matrix:
* design.csv

Saved results from the statisticl testing:
* final_part3.csv (and final_part3.xlsx)

Added HTML renderings of the notebooks for those who just want to see the analysis steps and figures (these may load faster):
* [understanding_IRS.html](https://pwilmart.github.io/IRS_normalization/understanding_IRS.html) is Part 1 (needs the png file)
* irs_diagram.png is an image for the understanding_IRS.html script
* [statistical_testing.html](https://pwilmart.github.io/IRS_normalization/statistical_testing.html) is Part 2
* [statistical_testing_ratios.html](https://pwilmart.github.io/IRS_normalization/statistical_testing_ratios.html) is Part 3
* [statistical_testing_take2.html](https://pwilmart.github.io/IRS_normalization/statistical_testing_take2.html) is Part 4

Added R scripts extracted from the notebooks. These can be used in RStudio or modified for your own analyses.

## Other repositories that may be of interest:
* [TMT data analysis examples](https://github.com/pwilmart/TMT_analysis_examples.git)
