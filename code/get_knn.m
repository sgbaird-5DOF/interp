function [nnpts,D,mu,sigma,idx] = get_knn(pts,dtype,K,NV)
arguments
    pts
    dtype char {mustBeMember(dtype,{'omega','norm','alen'})} = 'omega'
    K(1,1) double {mustBePositive,mustBeInteger} = 1
    NV.dispQ(1,1) logical = false
end
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Description: nearest neighbor distance histogram
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
[idxtmp,Dtmp] = knnsearch(pts,pts,'K',K+1);

%remove "self" column
idxtmp = idxtmp(:,2:end);
Dtmp = Dtmp(:,2:end);

%initialize
[mu,sigma] = deal(zeros(K,1));
D = cell(K,1);

for k = 1:K
    idx = idxtmp(:,k);
    
    %nearest neighbor pts
    nnpts = pts(idx,:);
    
    %get distances for plotting
    switch dtype
        case 'omega'
            Drad = get_omega(pts,nnpts);
            D{k} = rad2deg(Drad);
        case 'norm'
            D{k} = Dtmp(:,k);
        case 'alen' %general arc length formula
            assert(norm(pts(1,:))-1 < 1e-6,'norm(pts,2) must be 1 (within tol)')
            Drad = real(acos(dot(pts,nnpts,2)));
            D{k} = rad2deg(Drad);
    end
    
    mu(k) = mean(D{k});
    sigma(k) = std(D{k});
    
    if NV.dispQ
        disp(['nn: ' int2str(k) ', mu = ' num2str(mu(k)) ', sigma = ' num2str(sigma(k))])
    end
end