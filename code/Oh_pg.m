% OH_PG  Oh point group load/testing function
pgnum = 32 %cubic Oh
symnames = load('PGnames.mat'); %need to add crystal_symmetry_ops to path in order for this to work
symops = load('PGsymops.mat');
qpt = symops.Q{pgnum}; %choose point group symmetry

pgname = symnames.PG_names{pgnum};
disp(['loading point group: ' pgname])

t = num2cell(qpt,2);

Rmats = cellfun(@(q) qu2om(q),t,'UniformOutput',false);

Rcombs = allcomb(Rmats,Rmats);

axz = [0 0 1];
for i = 1:length(Rcombs)
	nA(i,:) = (Rcombs{i,1}*Rcombs{i,2}*axz.').';
end

t=num2cell(nA,1);
plot3(t{:},'*')