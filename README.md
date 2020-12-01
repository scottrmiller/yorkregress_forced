# york_forcefit
**York_forcefit** performs a weighted bivariate regression in which both variables have errors (Model II or errors-in-variables regression), and is specifically written to force the solution through the origin.  This code uses the maximum likelihood estimation algorithm of York et al. (2004).  The output of this function includes the estimate of slope, standard error of slope, correlation coefficient, t-statistic, and p-value.  York regression has a long history of use in the geosciences but has wide applicability in empirical studies where both variables commonly contain errors (Wehr and Saleska, 2017). 
 
**References:**

Wehr, R., & Saleska, S. R. (2017). The long-solved problem of the best-fit straight line: application to isotopic mixing lines. *Biogeosciences*, 14, 17-29, https://doi.org/10.5194/bg-14-17-2017

York, D., Evensen, N. M., Martınez, M. L., & De Basabe Delgado, J. (2004). Unified equations for the slope, intercept, and standard errors of the best straight line. *American Journal of Physics*, 72(3), 367-375, https://doi.org/10.1119/1.1632486
