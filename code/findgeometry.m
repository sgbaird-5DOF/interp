function geometry = findgeometry(qlist,precision)
arguments
	qlist(:,4) double {mustBeReal,mustBeFinite}
	precision(1,1) double {mustBeInteger} = 3
end
%-------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-30
%
% Description: Output the misorientation FZ geometry of a point given a
%	quaternion. (Assumes FCC 2020-06-30). The basic idea is check if in
%	surface. If in surface, check if in line. If in line, check if on point.
%	When nothing is found in the next geometry type (surface, line, point),
%	the geometry name from the next highest set is taken. If it's not on a
%	surface, the geometry is assumed to be 'interior'.
%
% Inputs:
%		qlist	=== rows of quaternion (according to convention in [1])
%
% Outputs:
%		geometry === geometry corresponding to the point in the
%			misorientation FZ (interior or the specific surface, line, or
%			point). Possible outputs to geometry for FCC: 'OAB', %'OBCE',
%			'OADE', 'CDE', 'OB', 'CE', 'ED', 'OE', 'OA', 'AC', 'B', %'E', 'A',
%			'C', 'O'. 3 or 4 letters == surface, 2 letters == line, 1 letter
%			== point
%
% References:
%		Misorientation FZ geometry equations obtained from [1] S. Patala,
%		C.A. Schuh, Symmetries in the representation of grain boundary-plane
%		distributions, Philos. Mag. 93 (2013) 524–573.
%		https://doi.org/10.1080/14786435.2012.722700.
%--------------------------------------------------------------------------

nq = size(qlist,1);
geometry = cell(1,nq);
for qnum = 1:nq
	q = qlist(qnum,:);
	%assign dummy variables for each of the quaternion components
	qcell = num2cell(q);
	[q0,q1,q2,q3] = qcell{:};
	
	crystal = 'FCC';
	
	%create helper function to compare values
% 	precision = 3;
	r = @(n1,n2) round(n1-n2,precision);
	
	%find geometry
	switch crystal
		case {'FCC','BCC'}
			
			d = q2rod(q);
			misFZ = InCubicFZ(d);
			
			if misFZ
				
				%define points
				k = sqrt(2)-1;
				
				qA = [cos(pi/8),sin(pi/8),0,0];
				qB = [1/sqrt(1+2*k^2),k/sqrt(1+2*k^2),k/sqrt(1+2*k^2),0];
				qC = 1/(2*sqrt(2))*[1/k,1,1,k];
				
				qD = [(7+(-4).*2.^(1/2)).^(-1/2),...
					((1/17).*(5+(-2).*2.^(1/2))).^(1/2),...
					((1/34).*(5+(-2).*2.^(1/2))).^(1/2),...
					((1/34).*(5+(-2).*2.^(1/2))).^(1/2)];
				
				qE = [sqrt(3)/2,1/(2*sqrt(3)),1/(2*sqrt(3)),1/(2*sqrt(3))];
				
				checknames = {...
					'OAB','OBCE','OADE','CDE',... surfaces
					'OB','CE','ED','OE','OA','AC',... lines
					'B','E','A','C','O','D'}; %points
				checklist = {... surfaces
					[r(q3,0)==0,r(q1,q2)~=0,r(q2,0)~=0], ... %OAB
					[r(q1,q2)==0,r(q3,0)~=0,r(q3,q1)~=0,r(2*q1+q3,q0)~=0], ... %OBCE
					[r(q2,q3)==0,r(q2,0)~=0,r(q1,q2)~=0,r(q1+2*q2,q0)~=0], ... %OADE
					[r(q1+q2+q3,q0)==0,r(q1,q2)~=0,r(q2,q3)~=0],... %CDE
					... lines
					[r(q1,q2)==0,r(q3,0)==0,r(q2-1/sqrt(2)*sin(pi/4),0)~=0,r(q0,1)~=0],... %OB
					[r(q1,q2)==0,r(2*q1+q3,q0)==0,r(q3,1/(2*sqrt(3)))~=0],... %CE
					[r(q2,q3)==0,r(q1+2*q2,q0)==0,r(q3,1/(2*sqrt(3)))~=0],... %ED
					[r(q1,q2)==0,r(q2,q3)==0,r(q3,1/(2*sqrt(3)))~=0,r(q0,1)~=0],... %OE
					[r(q2,0)==0,r(q3,0)==0,r(q1,sin(pi/8))~=0,r(q0,1)~=0],... %OA
					[r(q1,(sqrt(2)-1)*q0)==0,r(q3,(sqrt(2)-1)*q2)==0,r(q1,q2)~=0,r(q1,0)~=0],... %AC
					... points
					[r(q,qB)==0],... %B
					[r(q,qE)==0],... %E
					[r(q,qA)==0],... %A
					[r(q,qC)==0],... %C
					[r(q0,1)==0],...
					[r(q,qD)==0]}; %O
				
				for i = 1:length(checklist)
					check(i) = all(checklist{i});
				end
				
				geomID = find(check);
				
				if ~isempty(geomID)
					%take last geometry (assumes list is in order of surface, line,
					%and point, and that it doesn't register as true on two of the
					%same feature type, e.g. a true for 2 separate lines. If a point
					%and a line register as true, the point will be taken since it is
					%the finer feature.)
					geometry{qnum} = checknames{geomID(end)};
				else
					geometry{qnum} = 'interior';
				end
				
			else
				geometry{qnum} = 'exterior';
			end
			
	end
	
end

end