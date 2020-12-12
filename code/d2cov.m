function cov = d2cov(d,covpars,covtype)
arguments
    d double
    covpars double
    covtype char = 'squaredexponential'
end
% d2cov  distance array to covariance based on covtype and covpars. covpars
% follows documentation in 'KernelParameters' section of documentation for fitrgp
%
% for a tutorial/description on a few covariance functions, see
% https://www.mathworks.com/help/stats/kernel-covariance-function-options.html
% or https://www.cs.toronto.edu/~duvenaud/cookbook/
switch covtype
    case 'squaredexponential'
        L = covpars(2);
        sigma = covpars(1);
        cov = sigma.^2*exp(-0.5*d.^2)./L.^2;
end