function [tpredlist,tsdlist,propList,methodlist,A,B] = ...
    tunnelplot_test(testnum,ninputpts,n,tpredlist,tsdlist,propList,methodlist,A,B)
arguments
    testnum(1,1) double = 2
    ninputpts(1,1) double = 50000
    n(1,1) double = 300
    tpredlist = []
    tsdlist = []
    propList = []
    methodlist = []
    A = []
    B = []
end
switch testnum
    case 1
        %%
        mdlnum = 1;
        mdl = mdlcat(mdlnum);
        npts = mdl.mesh.npts;
        ids = randi(npts,18,1);
        paperfigure(3,3);
        cstAQ = true;
        for i = 1:2:length(ids)
            nexttile
            if cstAQ
                %         A = mdl.mesh.pts(ids(1),:);
                A = mdl.oref;
            else
                A = mdl.mesh.pts(ids(i),:);
            end
            B = mdl.mesh.pts(ids(i+1),:);
            tunnelplot(mdl,A,B);
        end
        
    case 2
        %%
%         ninputpts = 1000;
        switch ninputpts
            case 50000
                fnames = {...
                    'gpr50000_gitID-6ede824_puuID-8ecab3fe_paper-data5.mat',...
                    'pbary50000_gitID-6ede824_puuID-06272075_paper-data5.mat',...
                    'nn50000_gitID-6ede824_puuID-1f06700b_paper-data5.mat',...
                    'idw50000_gitID-6ede824_puuID-ba2b5202_paper-data5.mat'};
            case 1000
                fnames = {...
                    'pbary1000_gitID-6ede824_puuID-bdd7211c_paper-data5.mat',...
                    'gpr1000_gitID-6ede824_puuID-e6366570_paper-data5.mat',...
                    'nn1000_gitID-6ede824_puuID-73c9a5b3_paper-data5.mat',...
                    'idw1000_gitID-6ede824_puuID-ae3ba253_paper-data5.mat'};
        end
        nnames = length(fnames);
        mdllist = cell(nnames,1);
        for i = 1:nnames
            fname = fnames{i};
            load(fname,'mdl');
            mdllist{i} = mdl;
        end
        %%
        paperfigure();
        [tpredlist,tsdlist,propList,methodlist,A,B] = ...
            tunnelplot(mdllist,A,B,n,'nnQ',false,'nnQ2',false,...
            'tpredlist',tpredlist,'tsdlist',tsdlist,'propList',propList,...
            'methodlist',methodlist);
end

%% CODE GRAVEYARD
%{
%         if all(cellfun(@isempty,{tpredlist,tsdlist,propList,methodlist}))
%         [tpredlist,tsdlist,propList,methodlist,A,B] = tunnelplot(mdllist,A,B,10,'nnQ',false,'nnQ2',false);
%         else
            [tpredlist,tsdlist,propList,methodlist,A,B] = ...
                tunnelplot(mdllist,A,B,300,'nnQ',false,'nnQ2',false,...
                'tpredlist',tpredlist,'tsdlist',tsdlist,'propList',propList,...
                'methodlist',methodlist);
%         end
%}