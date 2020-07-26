function hAx = plotFZrodriguez_vtx()

[~,hAx] = plotFZrodriguez();

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

addpathdir({'misFZfeatures.mat','q2rod.m'})
load('misFZfeatures.mat','dlist')

featurelist = {'A','B','C','D','E','O'};
for i = 1:length(featurelist)
	feature = featurelist{i};
	tmp = num2cell(dlist.(feature));
	hold on
	plot3(tmp{:},'ko')
	text(tmp{1:2},tmp{3}+0.02,feature)
	%quiver3(tmp{1:2},0,0,0,tmp{3},'Color','k','ShowArrowhead','off')
end

ac=[dlist.A;dlist.C];
plot3(ac(:,1),ac(:,2),ac(:,3),'k:')
% hold off

hAx = gca;


end