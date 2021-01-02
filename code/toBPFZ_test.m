function toBPFZ_test(Kminus1)
arguments
    Kminus1(1,1) double = 1
end
K = Kminus1 + 1;

npts = 500;
qlist = get_cubo(npts);
nAlist = zeros(npts,3);
rng(12)
nAlist(1,:)= normr(rand(3,1).');
nAlist(2:npts,:) = normr(rand(3,npts-1).'-0.5);
posIDs = find(all(nAlist(2:npts,:) > 0,2)); %positive IDs
exID = posIDs(2);
nAex = nAlist(exID,:); %example
addpathdir('q2rod.m')

%get symmetrically equivalent ones
[nAout,nArot] = toBPFZ(qlist,nAlist);

%example
% nAoutex = nAout(posID1,:);
nArotex = nArot{exID+1};

%reference normal
nAref = nArot{1}(18,:);


%% plotting
paperfigure(1,2);
tnum = 1;
nexttile(tnum) % voronoi cell point cloud
nAsub = nArot(2:end);
v = num2cell(randi(24,npts-1,1));
nAdata = cellfun(@(x,i) x(i,:),nAsub,v,'UniformOutput',false);
nAdata = vertcat(nAdata{:});
t = n2c(nAdata);
plot3(t{:},'r.')

ax = gca;
ax.View = [140 35];
hold on

[k,D] = cellfun(@(x) knnsearch(x,nAref,'K',K,'IncludeTies',true),nAsub);
klen = cellfun(@length,k);
D1 = cellfun(@(x) x(1),D);
D2 = cellfun(@(x) x(2),D);

[pts,t] = cellvertcat(nAsub,k,1);
% kex = k{posID1}(1); %example
% nAsubex = nAsub{posID1};
ptex = pts(exID,:);

if any(D1 == D2)
   warning('more than one minimum found for at least one point')
end

% nnID = 1;
% pts = cellfun(@(x,i) x(i(nnID),:),nAsub,k,'UniformOutput',false);
% pts = vertcat(pts{:});

% clist = {'b','g','y','c','m','w'};
clist = {'b','','','','','',''};
for i = 1:K-1
    cellvertcatplot(nAsub,k,i,clist{i});
end

%reference unit normal
t = n2c(nAref*1.02);
plot3(t{:},'ko','MarkerFaceColor','w')
% quivplot(1.5)
axis equal tight off
test_voronoisphere(nArot{1}.')
if K == 2
    lgdlbltmp = {'$p_i^*$ (symmetrized)'};
    lgdlbl = ['$p_{i,1}$ (starting)',lgdlbltmp,'$p_\mathrm{ref}$ (reference)'];
    legend(lgdlbl,'Location','north','Interpreter','latex')
%     legend('starting',lgdlbltmp{:},'reference','Location','north')
else
    lgdlbltmp = strcat('NN',{' '},num2cell(num2str((1:K-1).')));
    legend('starting',lgdlbltmp{:},'reference','Location','southeast')
end

papertext(tnum);

tnum = 2;
nexttile(tnum) %voronoi cell, single intersection example
hold on
t = n2c(setdiff(nArotex,ptex,'rows'));
plot3(t{:},'m.')
t = n2c(ptex);
plot3(t{:},'b.')
t = n2c(nAref*1.02);
plot3(t{:},'ko','MarkerFaceColor','w')

papertext(tnum);
ax = gca;
ax.View = [140 35];
test_voronoisphere(nArot{1}.')

axis equal tight off
% lgdlbl = {'p_i','p_i^*','p_\rmmath{ref}'};
% lgdlbl = {'starting','symmetrized','reference'};
lgdlbltmp = {'$p_1^*$ (symmetrized)'};
lgdlbl = ['$p_{1,j}$ (starting)',lgdlbltmp,'$p_\mathrm{ref}$ (reference)'];
legend(lgdlbl,'Location','north','Interpreter','latex')

end
function cellvertcatplot(x,k,nnID,c)
    [~,t] = cellvertcat(x,k,nnID);
    plot3(t{:},[c '.'])
end

function [pts,t] = cellvertcat(x,k,nnID)
pts = cellfun(@(x,i) x(i(nnID),:),x,k,'UniformOutput',false);
pts = vertcat(pts{:});
t = n2c(pts);
end

%% CODE GRAVEYARD
%{
% t=n2c(pts);
% plot3(t{:},'.')
% 
% nnID = 2;
% pts2 = cellfun(@(x,i) x(i(nnID),:),nAsub,k,'UniformOutput',false);
% pts2 = vertcat(pts2{:});
% 
% t=n2c(pts2);
% plot3(t{:},'.')

% if ~all(klen == K)
%     warning('more than one minimum found for at least one point')
% end
%}
