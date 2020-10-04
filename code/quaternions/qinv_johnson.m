function qi = qinv_johnson(q)
arguments
   q(:,4)
end
% QINV takes the inverse of quaternions.
%-------------------------------------------------------------------------%
%
% Inputs:
%   q - An npts-by-4 array containing quaternion components in 4-vector 
%       format.
%
% Outputs:
%   qi - An npts-by-1 vector containing the inverse of each of the npts
%        quaternions in q, i.e., qi(i) = q(i,:)#-1, where #-1 denotes the
%        operation of quaternion inversion.
%
%Filename:  qinv.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%-------------------------------------------------------------------------%
%---perform quaternion inversion---%
a = qnorm(q).^-2;
a = a(:,ones(4,1));
qi = a.*qconj(q);

%---correct for numerical errors---%
qi(abs(qi) < eps) = 0;

end

%--------------------------CODE GRAVEYARD----------------------------------
%{

%---check inputs---%
assert((size(q,2) == 4),'q must be an npts-by-4 array.')

%}