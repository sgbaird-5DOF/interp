function [A,R] = symaxis(q,geometry)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-30
%
% Description: Return the symmetry axes (if multiple) and the geometry
% given a quaternion (assumes FCC or BCC crystal structure 2020-06-30).
%
% Input:
%		q			===	quaternion, according to convention in [1], where
%			dummy variables q0, q1, q2, q3 correspond to q(1),q(2),q(3),q(4),
%			respectively.
%
%		geometry	===	which feature on the misorientation FZ the point is
%			on. Possible inputs to geometry for FCC: 'OAB', %'OBCE', 'OADE',
%			'CDE', 'OB', 'CE', 'ED', 'OE', 'OA', 'AC', 'B', %'E', 'A', 'C',
%			'O'. 3 or 4 letters == surface, 2 letters == line, 1 letter ==
%			point
%
% Output:
%		A			===	[ax,ay,az]: the three possible symmetry axes that
%			correspond to the boundary plane FZ stereograms. If there is no
%			symmetry axis for one or two of these, the output will be [0,0,0].
%			Symmetry axis equations obtained from [1]
%
%		R			===	rotation matrix to go from az to [0,0,1]
%
% References:
%		[1] S. Patala, C.A. Schuh, Symmetries in the representation of grain
%		boundary-plane distributions, Philos. Mag. 93 (2013) 524–573.
%		https://doi.org/10.1080/14786435.2012.722700.
%
%--------------------------------------------------------------------------

%assign dummy variables for each of the quaternion components
qcell = num2cell(q);
[q0,q1,q2,q3] = qcell{:};

crystal = 'FCC';

%%
%assign symmetry axes based on crystal and geometry types
switch crystal
	case {'FCC','Bcc'}
		switch geometry
			case 'OAB'
				%sigma 15, 27b
				az = [q2,-q1,q0];
				
			case 'OBCE'
				%sigma 17b, 25b
				az = [(q0+q3)/sqrt(2),(q3-q0)/sqrt(2),-sqrt(2)*q1];
				
			case 'OADE'
				%sigma 21b, 23
				az = [-sqrt(2)*q2,(q0+q1)/sqrt(2),(q1-q0)/sqrt(2)];
				
			case 'CDE'
				az = [q0-q3,q0-q1,q0-q2];
				
			case 'OB'
				%sigma 9, 11, 19a, 27a
				ax = [q1,-q1,q0];
				%ay = [q0/sqrt(2),q0/sqrt(2),-sqrt(2)*q1]; from Patala2013, but seems to be incorrect (issues with signs)
				%by taking symbolic cross product of ax & az (these are perpendicular to each other), we get:
				ay = [-q0/sqrt(2),q0/sqrt(2),sqrt(2)*q1];
				az = [1/sqrt(2),1/sqrt(2),0];
				
			case 'CE'
				ax = [2*q1,q1+q3,q1+q3];
				ay = sqrt(2)*[-(q1+q3),q1,q1];
				az = [0,-1/sqrt(2),1/sqrt(2)];
				
			case 'ED'
				ax = [q1+q2,2*q2,q1+q2];
				ay = sqrt(2)*[q2,-(q1+q2),q2];
				az = [1/sqrt(2),0,-1/sqrt(2)];
				
			case 'OE'
				%sigma 7, 13b, 19b, 21a
				ax = 1/sqrt(2)*[q0+q1,q1-q0,-2*q1];
				ay = 1/sqrt(6)*[-3*q1+q0,3*q1+q0,-2*q0];
				az = [1/sqrt(3),1/sqrt(3),1/sqrt(3)];
				
			case 'OA'
				%sigma 5, 13a, 17a, 25a
				ax = [0,q0,q1];
				ay = [0,-q1,q0];
				az = [1,0,0];
				
			case 'AC'
				az = [0,1/sqrt(2),1/sqrt(2)];
				
			case {'B','D'} % (D == B)??
				ax = [q1,-q1,q0];
				ay = [-q0/sqrt(2),q0/sqrt(2),sqrt(2)*q1]; %sign error corrected from Patala 2013 paper. Same error as 'OB'
				az = [1/sqrt(2),1/sqrt(2),0];
				
			case 'E'
				%sigma 3
				ax = 1/sqrt(2)*[2/sqrt(3),-1/sqrt(3),-1/sqrt(3)];
				ay = [0,1/sqrt(2),-1/sqrt(2)];
				az = [1/sqrt(3),1/sqrt(3),1/sqrt(3)];
				
			case 'A'
				ax = [0,cos(pi/8),sin(pi/8)];
				ay = [0,-sin(pi/8),cos(pi/8)];
				az = [1,0,0];
				
			case 'C'
				ax = [1,0,0];
				ay = [0,1/sqrt(2),1/sqrt(2)];
				az = [0,-1/sqrt(2),1/sqrt(2)];
				
			case 'O'
				%sigma 1
				ax = [1,0,0];
				ay = [0,1,0];
				az = [0,0,1];				
			
			case 'interior' %is this right? Are there are no symmetry axes (so anywhere goes for meshing), or is it a hemisphere so I choose arbitrary axes?
				ax = [1,0,0];
				ay = [0,1,0];
				az = [0,0,1];
				
			case 'twosphere'
				ax = [1,0,0];
				ay = [0,1,0];
				az = [0,0,1];
		end
end

% assign [0,0,0] for the axes that are empty
if exist('ax','var') == 0 && exist('ay','var') == 0
	nsymaxes = 1;
	ax = [0,0,0];
	ay = [0,0,0];
else
	nsymaxes = 3;
end


%coordinate system rotation matrix
%tranpose into column vectors
ax = ax.';
ay = ay.';
az = az.';

%combine column vectors
A = [ax,ay,az];

switch nsymaxes
	case 1
		zaxis = [0,0,1].';
		R = vecpair2rmat(zaxis,az); %R*zaxis == az, only valid if single axis specified
	case 3
		R = A;
end





end

%----------------------------HELPER FUNCTIONS------------------------------
function R = vecpair2rmat(v1,v2)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Description: Compute rotation matrix. If v1 == v2 or v1 == -v2 within the
% given numerical precision, then the identity matrix or -1*identity matrix
% is given, respectively.
%
% Inputs:
%		v1, v2	===	two vectors to find the rotation matrix to go from
%							v1->v2
%
% Outputs:
%		R			=== Rotation matrix such that R*v1 == v2, and R\v2 == v1
%
% References
%	https://math.stackexchange.com/a/897677
%	https://math.stackexchange.com/a/476311/76513

precision = 12;
r = @(a,b) round(a-b,precision);

isEqual = all(r(v1,v2)==0); %pointing in same direction
isOpposite = all(r(v1,-v2)==0); %pointing in opposite direction
isParallel = isEqual || isOpposite; %either of previous two cases

if ~isParallel
	ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
	R = eye(3) + ssc(cross(v1,v2)) + ssc(cross(v1,v2))^2*(1-dot(v1,v2))/(norm(cross(v1,v2))^2);
	
elseif isEqual
	% same vector
	R = eye(3);
	
elseif isOpposite
	% vectors pointing in opposite directions
	R = -eye(3);
end

end
%----------------------------vecpair2rmat----------------------------------

%Notes:
%	Consider moving geometry calculation to separate function.


%--------------------------CODE GRAVEYARD----------------------------------
%{
if exist('az','var')==0
	az = [0,0,1];
end

if exist('ay','var')==0
	
end
%}