d = 3;

switch d
	case 3
		X = [...
			1 1 0
			1 0 1
			0 1 1]; %rows of vectors
		
		xA = [1 1 1];
	case 4
		X = [...
			1 1 0 1
			1 1 1 0
			0 1 1 1
			1 0 1 1 ];
		
		xA = [1 1 1 1];
end

lambda = numStabBary(X,xA)

if d == 3
	%scatter3(X(1,:),X(2,:),X(3,:))
end

