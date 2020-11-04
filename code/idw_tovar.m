function Fint = idw_tovar(X0,F0,Xint,p,rad,L)
% % IDW_TOVAR  Original inverse distance weight function by Andrew Tovar
% function Fint = idw(X0,F0,Xint,p,rad,L)
%
% Fint = idw(X0,F0,Xint) uses input coordinates X0 and input values F0
% where X0 is a N by M input matrix of N samples and M number of variables.
% F0 is vector of N responses. Xint is a Q by M matrix of coordinates to be
% interpolated. Fint is the vector of Q interpolated values.
%
% Fint = idw(X0,F0,Xint,p,rad) uses the power p (default p = 2) and radius
% rad (default rad = inf).
%
% Fint = idw(X0,F0,Xint,p,rad,L) uses L-distance. By defaults L=2
% (Euclidean norm).
%
% Example:
%
% X1 = [800;2250;3250;2250;900;500];
% X2 = [3700;4200;5000;5700;5100;4900];
% F =  [13.84;12.15;12.87;12.68;14.41;14.59];
% Q = 100;
% [X1int,X2int] = meshgrid(0:4000/(Q-1):4000, 3200:(5700-3200)/(Q-1):5700);
% Fint = idw([X1,X2],F,[X1int(:),X2int(:)]);
% contourf(X1int, X2int, reshape(Fint,Q,Q), 20)
%
% Contact info:
%
% Andres Tovar
% tovara@iupui.edu
% Indiana University-Purdue University Indianapolis
%
% Code developed for the course Design of Complex Mechanical Systems (ME
% 597) offered for the first time in Spring 2014
%
% Default input parameters
%
% Modified by Sterling Baird (2020-10-17)

if nargin < 6
    L = 2;
    if nargin < 5
        rad = inf;
        if nargin < 4
            p = 2;
        end
    end
end

% Basic dimensions
N = size(X0,1); % Number of samples
M = size(X0,2); % Number of variables
Q = size(Xint,1); % Number of interpolation points

% Inverse distance weight output
Fint = zeros(Q,1);

pd = pdist2(X0,Xint);
pd(pd == 0) = eps;
pd(pd > r) = Inf;
W = 1./(pd.^p);
Fint = sum(W.*F0)./sum(W)

for ipos = 1:Q    
    % Distance matrix
    DeltaX = X0 - repmat(Xint(ipos,:),N,1);
    DabsL = zeros(size(DeltaX,1),1);
    for ncol = 1:M
        DabsL = DabsL + abs(DeltaX(:,ncol)).^L;
    end
    Dmat = DabsL.^(1/L);
    
    Dmat(Dmat==0) = eps;
    Dmat(Dmat>rad) = inf;
    
    % Weights
    W = 1./(Dmat.^p);
    
    % Interpolation
    Fint(ipos) = sum(W.*F0)/sum(W);
end

end