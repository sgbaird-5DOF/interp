% sphere stereographic inverse test

% pts = eye(3);
pts = rand(5,3);
newpts = sphere_stereograph(pts);

newpts2 = sphere_stereograph_inverse(newpts);

t=num2cell(newpts2,1);
plot3(t{:},'*')

if all((pts-newpts) < 1e-6,'all')
	disp('projection mapped back properly')
else
	disp('projection did not map back properly')
end