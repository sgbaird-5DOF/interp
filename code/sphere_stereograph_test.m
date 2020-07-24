%sphere stereograph test

% pts = eye(3);
pts = rand(5,3);
newpts = sphere_stereograph(pts);

t=num2cell(newpts,1);
plot3(t{:},'*')

proj_down(newpts,1e-6)