clear; close all
load('misFZfeatures.mat','qlist','dlist')

figure
plotFZrodriguez();

scl = 0.1;
x = [1,0,0]*scl;
y = [0,1,0]*scl;
z = [0,0,1]*scl;

w = [0,0,0];
hold on
quiver3(w,w,w,x,y,z,1,'linewidth',1)
text(x(1),x(2)+0.01,x(3),'d_1','FontWeight','bold')
text(y(1),y(2),y(3)+0.01,'d_2','FontWeight','bold')
text(z(1)+0.01,z(2),z(3),'d_3','FontWeight','bold')

featurelist = {'A','B','C','D','E','O'};

for i = 1:length(featurelist)
	feature = featurelist{i};
	tmp = num2cell(dlist.(feature));
	hold on
	plot3(tmp{:},'ko')
	text(tmp{1:2},tmp{3}+0.02,feature)
	%quiver3(tmp{1:2},0,0,0,tmp{3},'Color','k','ShowArrowhead','off')
end
disp(' ')

%%
q = qlist.interior;
nA = normr([1 1 1]);

fnames ={'GBfive2oct.m','ax2qu.m'};
for i = 1:length(fnames)
	fname = fnames{i};
	file = dir(fullfile('**',fname));
	addpath(file.folder)
end

o = 1/sqrt(2)*GBfive2oct(q,nA);

five = GBoct2five(o);
q2 = five.q;
nA2 = five.nA;

q
q2
nA
nA2

if all(abs(q2-q) < 1e-6)
	disp('q matches')
else
	disp('q unequal')
	
end
if all(abs(nA2-nA) <1e-6)
	disp('nA matches')
else
	disp('nA unequal')
end

%%
addpathdir('PGnames.mat')
o
five
osyms = osymsets(o,32);

otemp1 = osyms{1}(1,:)
fivetemp = GBoct2five(otemp1)
otemp2 = GBfive2oct(fivetemp.q,fivetemp.nA)
fivetemp2 = GBoct2five(otemp2)