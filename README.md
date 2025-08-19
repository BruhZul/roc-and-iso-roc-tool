# ROC Tool

ROC Tool is a tool that allows the user to visualize ROC curves and easily perform more advanced analysis on them. Specifically, the tool allows users to plot multiple ROC curves and their Regions of Interest, compute AUC and RRA, and plot multiple iso-PM curves.

For more information on how Regions of Interest, RRA, and iso-PM curves work, refer to the following paper:

Luigi Lavazza, Sandro Morasca, Gabriele Rotoloni. Software Defect Prediction evaluation: New metrics based on the ROC curve.Information and Software Technology, 2025, https://doi.org/10.1016/j.infsof.2025.107865

If you use this tool for your research, please reference that paper.

You can find more details on how to use it and how to install it in the Documentation.pdf file.


The Tool uses the following libraries:
- Shiny (https://cran.r-project.org/package=shiny)
- ggplot2 (https://cran.r-project.org/package=ggplot2)
- shinyjs (https://cran.r-project.org/package=shinyjs)
- RColorBrewer (https://cran.r-project.org/package=RColorBrewer)
- xml2 (https://cran.r-project.org/package=xml2)
- iRRA (Morasca, Sandro, and Luigi Lavazza. "On the assessment of software defect prediction models via ROC curves." Empirical Software Engineering 25.5 (2020): 3977-4019. DOI: https://doi.org/10.1007/s10664-020-09861-4)

This project is licensed under the GNU General Public License 3.0. See the LICENSE file for more details. This license applies to all past and present versions and commits of the repository.
