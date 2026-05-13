# A Statistical Framework for Spatial Boundary Estimation and Change Detection: Application to the Sahel Sahara Climate Transition
### Authors: Stephen Tivenan, Indranil Sahoo, Yanjun Qian

In this paper, we introduce a novel, unified methodology which models boundaries and detects differences in boundary lines. Our method uses heteroskedastic Gaussian process (GP) regression and the scaled Maximum Absolute Difference (MAD) Global Envelope Test (GET) to estimate spatial boundaries and distinguish differences between boundary lines. The framework enables local boundary modeling, while establishing a procedure to detects changes or natural fluctuations in a line. The proposed procedure was applied to arid climates' transitional boundaries in the Sahel and the Sahara between  1960-1989 and in a simulation study. Primarily the simulation was conducted to validate GET's size under the null and to demonstrated the test's power in identifying boundary shifts. 

The Boundary-MAD-GET repository highlights simulation study procedure in "Simulation Example"  page and in an R Shiny "Interactive Simulation" page. The interactive simulation generates two random observations of boundary lines from two separate B-splines and Gaussian Process (GP) noises, then evaluates the scaled MAD GET to access if there is a significant difference between the mean lines. Simulation settings for the constructing the line can be found in the coefficient slider, the number of grid points and the noise setting. The scaled MAD GET settings include both the Monte Carlo Simulation and the size of the distribution for each simulation.


To run this application locally, ensure you have R installed and run:

```R
install.packages(c("shiny", "MASS", "splines2", "plotly", "dplyr"))
shiny::runGitHub("YOUR_USERNAME/Global-Envelope-Simulation")
