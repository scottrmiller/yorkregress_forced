function [b,sigma_b,r,p,t,MSWD,sigma_MSWD,p_chi,i] = yorkregress_forced(X,Y,sigma_X,sigma_Y,ri)

% YORKREGRESS_FORCED York regression forced through origin
%
% Description 
%     Regression forced through the origin where both x and y variables
%     have errors (errors-in-variables regression), based on the maximum
%     likelihood estimation algorithm of York et al. (2004).  This function
%     follows the method of Thirumalai et al. (2011) and Trappitsch et al.
%     (2018) for forcing the solution through a point, in this case by
%     adding a point corresponding to the origin (0,0) at the end of the
%     vectors for x and y and assigning each very small errors (1e-15), but
%     corrects the estimates of the standard error of slope, correlation
%     coefficient, and the t-statistic.  Goodness-of-fit parameter from
%     York et al. (2004) and standard deviation of chi-square from Wendt
%     and Carl (1991).  For a Matlab code of the general York regression
%     not forced through the origin or forced through one of the
%     datapoints, see, for example, Thirumalai et al. (2011).
%
% Input arguments
%     X             Observed x variable (column vector)
%     Y             Observed y variable (column vector)
%     sigma_X       Error in x (column vector)
%     sigma_Y       Error in y (column vector)
%     ri            [Optional] correlation coefficient of errors. Equals 
%                   zero if not assigned.
%
% Output arguments
%     b             Slope of best fit line
%     sigma_b       Standard error of slope
%     wr            Correlation coefficient
%     p             p value, for H0 that slope of line = 0. Values <
%                   significance level (say 0.05) indicate slope is
%                   significant and non-zero.
%     t             t statistic
%     MSWD          Mean weighted squared deviation (MSWD), a measure of
%                   goodness-of-fit.  Also known as reduced weighted chi
%                   square
%     sigma_MSWD    Standard deviation of MSWD
%     p_chi         p value for chi-square, for H0 that the observed data 
%                   are well represented by a linear relationship.  Values
%                   greater than the significance level support linear
%                   relationship.
%     i             Number of iterations
%
% References:
%
%     Mahon, K. I. (1996). The New “York” regression: Application of an
%     improved statistical method to geochemistry. International Geology
%     Review, 38(4), 293-303, https://doi.org/10.1080/00206819709465336
%
%     Thirumalai, K., A. Singh, and R. Ramesh (2011), A MATLAB code to
%     perform weighted linear regression with (correlated or uncorrelated)
%     errors in bivariate data, Journal of the Geological Society of India,
%     77(4), 377-380, https://doi.org/10.1007/s12594-011-0044-1
%     https://github.com/planktic/RegressBivariate
%
%     Trappitsch, R., Boehnke, P., Stephan, T., Telus, M., Savina, M. R., %
%     Pardo, O., Davis, A.M., Dauphas, N., Pellin, M.J., & Huss, G. R.
%     (2018). New constraints on the abundance of 60Fe in the early solar
%     system. The Astrophysical Journal Letters, 857(2), L15,
%     http://doi.org/10.3847/2041-8213/aabba9
%
%     Wehr, R., & Saleska, S. R. (2017). The long-solved problem of the
%     best-fit straight line: application to isotopic mixing lines.
%     Biogeosciences, 14, 17-29, https://doi.org/10.5194/bg-14-17-2017
%
%     Wendt, I., & Carl, C. (1991). The statistical distribution of the
%     mean squared weighted deviation. Chemical Geology: Isotope Geoscience
%     Section, 86(4), 275-285, https://doi.org/10.1016/0168-9622(91)90010-T
% 
%     York, D., Evensen, N. M., Martinez, M. L., & De Basabe Delgado, J.
%     (2004). Unified equations for the slope, intercept, and standard
%     errors of the best straight line. American Journal of Physics, 72(3),
%     367-375, https://doi.org/10.1119/1.1632486
%
% Author: Scott R. Miller, Department of Geology and Geophysics,
% University of Utah, 
% Date: 23 December 2020



%% Tolerance (tol) for determining number of iterations in maximum likelihood 
% estimate (MLE).  Iteration procedure converges when successive
% estimations of slope (b) agree within this tolerance.  Usually converges
% in about 10 iterations.  However, occasionally it can go on so a maximum
% number of iterations (maxiter) can be set.
tol = 1e-15;
maxiter = 100;


%% If input for ri is not provided, make it zero.
if nargin < 5
    ri = 0;
end


%% Modify input data to force through origin.  
% Add coordinates for origin (0,0) to X and Y vectors and add very small
% values (1e-15) to sigma vectors.
Xi = [X; 0];
Yi = [Y; 0];
sigma_Xi = [sigma_X;1e-15];
sigma_Yi = [sigma_Y;1e-15];


%% Ordinary least squares regression, for approximating initial value of slope (b1)
tmp = polyfit(Xi,Yi,1);
b1 = tmp(1);


%% Maximum likelihood estimation using weights

% Error weights
omega_X = 1./(sigma_Xi.^2);                     
omega_Y = 1./(sigma_Yi.^2);
alpha = sqrt(omega_X.*omega_Y);

b = b1;                             
delta = tol;                            
i = 0; 

% MLE loop.  Iterations continue while difference between successive
% estimates of slope (b) are greater than or equal to specified tolerance,
% AND while greater than 100
while delta >= tol && i <= maxiter
    i = i+1;                            % Count iterations
    b_previous = b;                     % Save previous estimate of slope
    Wi = omega_X.*omega_Y./(omega_X + b^2*omega_Y - 2*b*alpha*ri);
    X_bar = sum(Wi.*Xi)/sum(Wi);   % Weighted mean of observations of X
    Y_bar =  sum(Wi.*Yi)/sum(Wi);  % Weighted mean of observations of Y
    U = Xi - X_bar;
    V = Yi - Y_bar;
    beta = Wi.*((U./omega_Y)+((b.*V)./omega_X) - (b.*U + V).*(ri./alpha));
    b = sum(Wi.*beta.*V)/sum(Wi.*beta.*U);    % slope of best fit line, y = bx
    delta = abs(b - b_previous);
end

% To correctly estimate the standard error of slope, correlation
% coefficient, the p-value, and t-statistic, we re-calculate some of the
% above variables after using original X and Y (i.e., without point at
% origin) and resizing Wi and beta to sizes of X and Y.
n = numel(U)-1;
Wi = Wi(1:n);
beta = beta(1:n);
X_bar = sum(Wi.*X)/sum(Wi);
Y_bar = sum(Wi.*Y)/sum(Wi);
U = X - X_bar;
V = Y - Y_bar;
xi = X_bar + beta;
x_bar = sum(Wi.*xi)/sum(Wi);
ui = xi - x_bar;
variance_b = 1./(sum(Wi.*(ui.*ui)));
sigma_b = sqrt(variance_b);    % standard error of slope b


%% Calculate t statistic, p-value, and correlation coefficient (r)
B = 0;
t = (b-B)/sigma_b;
p = 2*(1-tcdf(abs(t),n-2));
r = sum(U.*V)./sqrt(sum(U.^2)*sum(V.^2)); 


%% Calculate goodness-of-fit parameter (MSWD)  
% The mean squared weighted deviation (MSWD), or reduced weighted
% chi-square, is calculated following York et al. (2004).  In the case here
% where slope is the only free parameter, there are n-1 degrees of freedom
% (Mahon, 1996). If the data are well represented by a straight-line model
% through the origin, the reduced weighted chi-square, or MSWD, is expected
% to have a value about 1 +/-2 standard deviations (Wendt and Carl, 1991;
% Wehr and Saleska, 2017).  Deviations significantly different from one
% indicate that the variables x and y are not well related by a straight
% line and/or the assigned errors are incorrect (possibly not normally
% distributed).
S = sum((Wi).*(Y - b.*X).^2);   % Weighted sum of deviations from best-fit 
                                % line, or chi square.  Equivalent to
                                % chi-square formulation in Wendt and Carl
                                % (1991) (summation in equation 1).
MSWD = S/(n-1);                 % Reduced weighted chi-square, or MSWD 
sigma_MSWD = sqrt(2/(n-1));     % Standard deviation in the weighted 
                                % reduced chi squared statistic (Wendt and
                                % Carl, 1991; Wehr and Saleska, 2017).
                                

%% Calculate p value for chi-square statistic
% This is done by calculating the area under the curve (probability) to the
% left and right of the chi-square value.  That is, these are the
% probabilities of realize more extreme, lower or higher, chi-square
% values.  The smaller of these is taken as the final value.
p_chi_left = chi2cdf(S,n-1);    
p_chi_right = chi2cdf(S,n-1,'upper'); 
p_chi = min(p_chi_left,p_chi_right);


