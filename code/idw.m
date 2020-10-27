function [yq,W,r,nints,numnonints,int_fraction] = idw(X,Xq,y,r,L)
arguments
    X double %rows of input points
    Xq double %rows of query points
    y(:,1) double %property values
    r double = [] %radius
    L(1,1) double = 2 %Euclidean or 2-norm
end
% IDW inverse-distance weighting interpolation
%referenced Andres Tovar's FEX implementation (idw.m)

%default radius or user-supplied radius
if isempty(r)
    [~,~,mu] = get_knn(X,'norm',1);
    r = mu*sqrt(2); %note: since alen = ~0.5*omega, mu in octonion is ~2x this (2*mu*sqrt(2))
end

pd = pdist2(X,Xq); %pairwise distance matrix

%parse pd based on radius
pd(pd == 0) = eps;
pd(pd > r) = Inf;

W = 1./(pd.^L); %weight matrix

%store as sparse matrix if more than half of elements are zero 
numnonzero = nnz(W);
if numnonzero < 0.5 * numel(pd)
    W = sparse(W);
end

%compute interpolated values
yq = sum(W.*y)./sum(W);

%assign NN value if no GBs fall into radius for a particular property
nnIDs = find(isnan(yq));
nnList = dsearchn(X,Xq(nnIDs,:));
%assign NN property values
yq(nnIDs) = y(nnList);

if issparse(yq)
    yq = full(yq);
end
%make into column vector
yq = yq.';

%idw parameters
numnonints = length(nnIDs);
npredpts = size(Xq,1);
nints = npredpts-numnonints;
int_fraction = nints/(nints+numnonints);
disp(['numnonints = ' int2str(numnonints)])

end

%% CODE GRAVEYARD
%{
% yavg = mean(y);
% yq(isnan(yq)) = yavg;
% yq(yq <= eps) = yavg;

if isempty(r)
    pdX = pdist(X);
%     sig = std(pdX);
    mu = mean(pdX);
%     if sig > 0
%         r = mu-2*sig;
%     else
        r = mu;
%     end
    clear pdX
end

    r double = [] %radius, if rand(size(X,1)^2) can't fit in memory, be sure to specify

%     pdX = pdist(X);
%     mu = mean(pdX)/sqrt(2);
%     r = mu;
%     clear pdX

%}