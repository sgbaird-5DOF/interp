function y = tmvn(mu,Sigma,l,u,n,zerofloorQ)
arguments
   mu(:,1) double
   Sigma(:,:) double
   l(:,1) double = 0.75*mu;
   u(:,1) double = 1.5*pred;
   n(1,1) double = 1;
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
%  125–148. https://doi.org/10.1111/rssb.12162.
%--------------------------------------------------------------------------
sz = size(mu);

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

end