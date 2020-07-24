%{
Possible inputs to geometry for FCC: 'OAB', 'OBCE', 'OADE', 'CDE', 'OB',
'CE', 'ED', 'OE', 'OA', 'AC', 'B', %'E', 'A', 'C', 'O'
%}
geometry = 'CE';

switch geometry
	case 'ED'
		q = [3 1 2 2];
		q = q/norm(q); %normalize quaternion
end

[A,R] = symaxis(q,geometry);

ax = A(:,1);
ay = A(:,2);
az = A(:,3);

check = all(R\az == [0,0,1].');

if check
	disp('inverting rotation matrix resulted in +z vector')
else
	disp('inverting rotation matrix did not result in +z vector')
end
