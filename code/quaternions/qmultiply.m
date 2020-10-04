function qc = qmultiply(qa,qb)
arguments
    qa(:,4)
    qb(:,4)
end
% QMULTIPLY performs quaternion multiplication.
%-------------------------------------------------------------------------%
%
% Inputs:
%   qa - An npts-by-4 array containing quaternion components in 4-vector 
%        format.
%   qb - An npts-by-4 array containing quaternion components in 4-vector 
%        format.
%
% Outputs:
%   qc - An npts-by-4 array containing the quaternion components resulting
%        from multiplication of quaternions in qa with quaternions in qb by
%        row, i.e., qc(i,:) = qa(i,:)#qb(i,:), where # is here used to
%        denote quaternion multiplication.
%
%Filename:  qmultiply.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%-------------------------------------------------------------------------%

%---check inputs---%
assert((size(qa,1) == size(qb,1)),'qa and qb must have the same number of points.')

%---perform quaternion multiplication---%
qc(:,1) = qb(:,1).*qa(:,1)-qb(:,2).*qa(:,2)-qb(:,3).*qa(:,3)-qb(:,4).*qa(:,4);
qc(:,2) = qb(:,2).*qa(:,1)+qb(:,1).*qa(:,2)+qb(:,4).*qa(:,3)-qb(:,3).*qa(:,4);
qc(:,3) = qb(:,3).*qa(:,1)-qb(:,4).*qa(:,2)+qb(:,1).*qa(:,3)+qb(:,2).*qa(:,4);
qc(:,4) = qb(:,4).*qa(:,1)+qb(:,3).*qa(:,2)-qb(:,2).*qa(:,3)+qb(:,1).*qa(:,4);

%-----------------------------CODE GRAVEYARD-------------------------------
%{
assert((size(qa,2) == 4),'qa must be an npts-by-4 array.')
assert((size(qb,2) == 4),'qb must be an npts-by-4 array.')
%}