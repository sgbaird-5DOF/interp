function y = tmvn(mu,Sigma,l,u,n,method,zerofloorQ)
arguments
    mu(:,1) double
    Sigma(:,:) double
    l(:,1) double = 0.75*mu;
    u(:,1) double = 1.5*mu;
    n(1,1) double = 1;
    method char {mustBeMember(method,{'mvrandn','slicesample'})} = 'slicesample'
    zerofloorQ(1,1) logical = true
end
% TMVN  sample from truncated multivariate normal distribution using mean and covariance matrix
%--------------------------------------------------------------------------
% Inputs:
%  mu - mean of distribution
%  Sigma - covariance matrix (or vector of diagonals along covariance matrix)
%  l - lower bound(s)
%  u - upper bound(s)
%  n - number of posterior samples
%  zerofloorQ - whether to constrain the sampling to non-negative values or not
%
% Outputs:
%  y - new sampled values
%
% Usage:
%  y = tmvn(mu,Sigma,l,u,n,zerofloorQ);
%
% Dependencies:
%  Truncated Normal and Student's t-distribution toolbox (FEX) (1)
%
% see also MVRANDN, MVNRND, TRANDN
%
% Author(s): Sterling Baird
%
% Date: 2020-01-21
%
% References:
%  (1) Botev, Z. I. The Normal Law under Linear Restrictions: Simulation
%  and Estimation via Minimax Tilting. J. R. Stat. Soc. B 2017, 79 (1),
%  125???148. https://doi.org/10.1111/rssb.12162.
%--------------------------------------------------------------------------
sz = size(mu);

switch method
    case 'mvrandn'
        addpathdir('mvrandn.m')
        
        % allow for scalar inputs to l and u
        if isscalar(l)
            l = repelem(l,sz(1),sz(2));
        end
        if isscalar(u)
            u = repelem(u,sz(1),sz(2));
        end
        
        if zerofloorQ
            % output will never be lower than 0
            l = max([-mu,l-mu],[],2);
        end
        sig = mvrandn(l,u,Sigma,n); % sample from posterior distribution with mu == 0
        y = mu+sig; %add mean (mu) back in
    case 'slicesample'
        y = slicesample(mu,n,'pdf',@(X) tmvnpdf(X,mu,Sigma,l,u));
end
end