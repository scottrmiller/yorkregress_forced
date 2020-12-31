# yorkregress_forced
**Yorkregress_forced** performs a York regression, a weighted bivariate regression in which both variables have errors (Model II or errors-in-variables regression), and forces the solution through the origin.  This code uses the maximum likelihood estimation (MLE) algorithm of York et al. (2004).  The output of this function includes the estimate of slope, standard error of slope, correlation coefficient, t-statistic, p-value of the slope, and the number of iterations in the MLE routine.  To assess goodness of fit, the function also outputs the mean squared weighted deviation (MSWD) or reduced weighted chi-square, the standard deviation of MSWD, and the p-value of chi-square (Wendt and Carl, 1991).  This code adopts the method of Thirumalai et al. (2011) and Trappitsch et al. (2018) for forcing the solution through a point.  York regression has a long history of use in the geosciences but has wide applicability in empirical studies broadly where both variables commonly contain errors (Wehr and Saleska, 2017).
 
**References:**

Thirumalai, K., A. Singh, and R. Ramesh (2011), A MATLAB code to perform weighted linear regression with (correlated or uncorrelated) errors in bivariate data, Journal of the Geological Society of India, 77(4), 377-380, https://doi.org/10.1007/s12594-011-0044-1 (for code, see https://github.com/planktic/RegressBivariate)

Trappitsch, R., Boehnke, P., Stephan, T., Telus, M., Savina, M. R., Pardo, O., Davis, A.M., Dauphas, N., Pellin, M.J., & Huss, G. R. (2018). New constraints on the abundance of 60Fe in the early solar system. The Astrophysical Journal Letters, 857(2), L15, http://doi.org/10.3847/2041-8213/aabba9

Wehr, R., & Saleska, S. R. (2017). The long-solved problem of the best-fit straight line: application to isotopic mixing lines. *Biogeosciences*, 14, 17-29, https://doi.org/10.5194/bg-14-17-2017

Wendt, I., & Carl, C. (1991). The statistical distribution of the mean squared weighted deviation. Chemical Geology: Isotope Geoscience Section, 86(4), 275-285, https://doi.org/10.1016/0168-9622(91)90010-T

York, D., Evensen, N. M., Martınez, M. L., & De Basabe Delgado, J. (2004). Unified equations for the slope, intercept, and standard errors of the best straight line. *American Journal of Physics*, 72(3), 367-375, https://doi.org/10.1119/1.1632486
