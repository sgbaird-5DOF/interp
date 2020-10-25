%avgrepeats test
seed = 11;
rng(seed)

meshList = [ ...
	1 1 1
	1 1 1
	0 1 0
	1 1 1
	0 1 0
	1 0 1
	2 2 2];

npts = size(meshList,1);
propList = randi(100,npts,1);

[meshListNew,propListNew] = avgrepeats(meshList,propList)

t=n2c(meshList);
table(t{:},propList,'VariableNames',{'x','y','z','property'})

t=n2c(meshListNew);
table(t{:},propListNew,'VariableNames',{'x','y','z','property'})