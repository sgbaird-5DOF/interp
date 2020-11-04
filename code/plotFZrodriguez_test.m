%PLOTFZRODRIGUEZ_TEST  
clear; close all;
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

%%
k = sqrt(2)-1;

% 'B'
qlist.B = [1/sqrt(1+2*k^2),k/sqrt(1+2*k^2),k/sqrt(1+2*k^2),0];
% 'E'
qlist.E = [sqrt(3)/2,1/(2*sqrt(3)),1/(2*sqrt(3)),1/(2*sqrt(3))];
% 'A'
qlist.A = [cos(pi/8),sin(pi/8),0,0];
% 'C'
qlist.C = 1/(2*sqrt(2))*[1/k,1,1,k];
% 'D'
qlist.D = [2.^(-1/2),1+(-1).*2.^(-1/2),(-1/2)+2.^(-1/2),(-1/2)+2.^(-1/2)];
% 'O'
qlist.O = [1,0,0,0];


qnames = fields(qlist);
for i = 1:length(qnames)
	qname = qnames{i};
	q = qlist.(qname);
	q = q/norm(q); %normalize quaternion
	qlist.(qname) = q;
	%q = disorientation(q,'cubic');
	
	d = q(2:4)./q(1);
	dlist.(qname) = d;
	
	hold on
	plot3(d(1),d(2),d(3),'ro')
	text(d(1),d(2),d(3)+0.01,qname)
end

%%
lnames = {'AB','OE','CE','BC','OA','AC','OB','ED'}; %line names

for i = 1:length(lnames)
	lname = lnames{i};
	dlist.(lname) = 0.5*(dlist.(lname(2)) + dlist.(lname(1)));
	d = dlist.(lname);
	qlist.(lname) = rod2q(d);
	hold on
% 	plot3(d(1),d(2),d(3),'ro')
% 	text(d(1),d(2),d(3)+0.01,lname)
end

snames = {'OAB','OBCE','OADE','CDE'}; %surface names
for i = 1:length(snames)
	sname = snames{i};
	dlist.(sname) = 1/3*(dlist.(sname(3)) + dlist.(sname(2)) + dlist.(sname(1)));
	d = dlist.(sname);
	qlist.(sname) = rod2q(d);
	hold on
	%plot3(d(1),d(2),d(3),'ro')
	%text(d(1),d(2),d(3)+0.01,sname)
end

cnames = {'CDE','O','A','B'};
d = [0,0,0];
for i = 1:length(cnames)
	cname = cnames{i};
	d = [d+dlist.(cname)];
end
d = d/length(cnames);
dlist.interior = d;
qlist.interior = rod2q(d);

plot3(d(1),d(2),d(3),'ko')
text(d(1),d(2),d(3)+0.01,'interior')

dlist.twosphere = d;
qlist.twosphere = rod2q(d);



ac=[dlist.A;dlist.C];
plot3(ac(:,1),ac(:,2),ac(:,3),':')

save('misFZfeatures.mat','qlist','dlist','qnames','lnames','snames')
% 
% %%
% qD = [sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))];
% 
% d = qD;
% q(1) = sqrt(1/(1+sum(d.^2)));
% q(2) = q(1)*d(1);
% q(3) = q(1)*d(2);
% q(4) = q(1)*d(3);
