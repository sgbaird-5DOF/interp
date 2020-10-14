%% result parity plot
methodlist = {'sphgpr','gpr','sphbary','pbary','nn','avg'};
cellfun(@(type) median(mdlparstbl(strcmp(mdlparstbl.method,type),:).rmse),methodlist)
cellfun(@(type) median(mdlparstbl(strcmp(mdlparstbl.method,type),:).mae),methodlist)
uuid = 'ad8908c4';
tbltmp = mdltbl(strcmp(mdltbl.uuid,uuid),:);
barypars = tbltmp.barypars{1};
yactual = padcat(tbltmp.data.props(barypars.ids),tbltmp.data.props(barypars.ilist));
ypred = padcat(tbltmp.propOut{1}(barypars.ids),tbltmp.propOut{1}(barypars.ilist));
parityplot(yactual,ypred,'title','sphbary','legend',{'interp','nn'})

%% distance histogram
pd1 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts).';
pd2 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_omega).';
pd3 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_alen).';
parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','')

%% distance histogram
pts = normr(rand(388,3)-0.5);
pts2 = normr(rand(388,3)-0.5);
pd1 = pdist(pts);
pd3 = pdist(pts,@get_alen);
parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','','title','388 3D points')
figure
t=n2c(pts);
scatter3(t{:})
axis equal