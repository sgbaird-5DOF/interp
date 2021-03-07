function y = tmvnpdf(X,mu,Sigma,l,u)
arguments
    X double
    mu double
    Sigma double
    l double
    u double
end
%TMVNPDF  truncated multi-variate normal probability density function
y = mvnpdf(X,mu,Sigma);
y(y < l) = 0;
y(y > u) = 0;