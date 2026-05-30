# A Statistical Framework for Spatial Boundary Estimation and Change Detection: Application to the Sahel Sahara Climate Transition
### Authors: Stephen Tivenan, Indranil Sahoo, Yanjun Qian

In this paper, we introduce a novel, unified methodology which models boundaries and detects differences in boundary lines. Our method uses heteroskedastic Gaussian process (GP) regression and the scaled Maximum Absolute Difference (MAD) Global Envelope Test (GET) to estimate spatial boundaries and distinguish differences between boundary lines. The framework enables local boundary modeling, while establishing a procedure to detects changes or natural fluctuations in a line. The proposed procedure was applied to arid climates' transitional boundaries in the Sahel and the Sahara between  1960-1989 and in a simulation study.


The Code for the analysis and pictures shown in the paper can be found in the Boundary-Exploratory-Code folder. The Decade_Boundary_Line_Primary_Secondary.R file applies the MAD-GET for two decadal lines (assumed to be independent from each other) for both the primary and secondary boundary lines. While the Specific_Year_GET_Boundary_code.R and Variance_Plots_Primary_Secondary.R files contains functions that generate the GET for a training set and applies it to a specific and calculates the variance of the boundary line over a period of time. Each of the functions can be applied in on the data found in the Boundary_Data_Set folder inside the Boundary-Exploratory-Code.   The last piece of code is shown as a Rmarkdown as a pdf. Here, the decade boundary of the 1960s, 1970s, and 1980s are displayed with an introduction of the 1960s data analysis section of the paper.


Beyond highlighting the main application of the defined methodology, a simulation study was  done to validate GET's size under the null and to demonstrated the test's power in identifying boundary shifts. he Boundary-MAD-GET repository highlights simulation study procedure in "Simulation Example"  page and in an R Shiny "Interactive Simulation" page (https://tivenansm.github.io/Boundary-Scaled-MAD-GET/). The interactive simulation generates two random observations of two defined separate boundaries fromB-splines and Gaussian Process (GP) noises. Then the scaled MAD GET is impplemented to access if there is a significant difference between the mean lines. Simulation settings for the constructing the line can be found in the coefficient slider, the number of grid points and the noise settings. The scaled MAD GET settings include both the Monte Carlo Simulation and the size of the distribution for each simulation.


To run this application locally make sure you have R installed and download the code from  myapp folder with the following packages:

```R
install.packages(c("shiny", "MASS", "splines2", "plotly", "dplyr"))
shiny::runGitHub("YOUR_USERNAME/Global-Envelope-Simulation")
