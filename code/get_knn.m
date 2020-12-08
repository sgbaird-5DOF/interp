function [nnX,D,mu,sigma,idxtmp] = get_knn(X,dtype,K,NV)
arguments
    X double
    dtype char {mustBeMember(dtype,{'omega','norm','alen'})} = 'omega'
    K(1,1) double {mustBePositive,mustBeInteger} = 1
    NV.dispQ(1,1) logical = false
    NV.ID = []
    NV.Y = []
end
% GET_KNN  k-nearest neighbor points, distances, mean, and std
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Inputs:
%  pts - rows of points (euclidean)
%
%  dtype - distance type
%
% Outputs:
%  creates a histogram figure
%
% Usage:
%	pts = sqrt2norm(normr(rand(1000,8))); %octonion case (omega), input must have norm == sqrt(2)
%  nnhist(pts)
%
%  pts = nnhist(rand(1000,5)); % euclidean distance
%  nnhist(pts)
%
% Dependencies:
%  get_omega.m
%
% Notes:
%  *
%--------------------------------------------------------------------------
%get nearest neighbor IDs and euclidean distances
ID = NV.ID;
Y = NV.Y;
assert(isempty(ID) || isempty(Y),'ID and Y should not be supplied simultaneously')
if ~isempty(ID)
    Y = X(ID,:);
elseif isempty(Y)
    Y = X; %check NNs within set of points (otherwise check NN against specific pts in NV.Y)
end
    
[idxtmp,Dtmp] = knnsearch(X,Y,'K',K+1);

%remove "self" column
idxtmp = idxtmp(:,2:end);
Dtmp = Dtmp(:,2:end);

%initialize
[mu,sigma] = deal(zeros(K,1));
D = cell(K,1);

for k = 1:K
    idx = idxtmp(:,k);
    
    %nearest neighbor pts
    nnX = X(idx,:);
    
    %get distances for plotting
    switch dtype
        case 'omega'
            Drad = get_omega(nnX,Y);
            D{k} = rad2deg(Drad);
        case 'norm'
            D{k} = Dtmp(:,k);
        case 'alen' %general arc length formula
            assert(norm(X(1,:))-1 < 1e-6,'norm(pts,2) must be 1 (within tol)')
            Drad = real(acos(dot(nnX,Y,2)));
            D{k} = rad2deg(Drad);
    end
    
    mu(k) = mean(D{k});
    sigma(k) = std(D{k});
    
    if NV.dispQ
        disp(['nn: ' int2str(k) ', mu = ' num2str(mu(k)) ', sigma = ' num2str(sigma(k))])
    end
end