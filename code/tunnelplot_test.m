function [tpredlist,tsdlist,propList,methodlist,A,B] = ...
    tunnelplot_test(testnum,ninputpts,n,tpredlist,tsdlist,propList,methodlist,A,B,nv)
arguments
    testnum(1,1) double = 2
    ninputpts(1,1) double = 1000
    n(1,1) double = 100
    tpredlist = []
    tsdlist = []
    propList = []
    methodlist = []
    A = []
    B = []
    nv.extend = 0.1
end
extend = nv.extend;

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
        if isempty(tpredlist)
            switch ninputpts
                case 50000
%                     fnames = {...
%                         'gpr50000_gitID-6ede824_puuID-8ecab3fe_paper-data5.mat',...
%                         'pbary50000_gitID-6ede824_puuID-06272075_paper-data5.mat',...
%                         'nn50000_gitID-6ede824_puuID-1f06700b_paper-data5.mat',...
%                         'idw50000_gitID-6ede824_puuID-ba2b5202_paper-data5.mat'};
                    fnames = {...
                        'gpr50000_gitID-dd6fa64_puuID-88aa16b9_paper-data9.mat',...
                        'pbary50000_gitID-dd6fa64_puuID-53fd75d6_paper-data9.mat',...
                        'nn50000_gitID-dd6fa64_puuID-5b81c95f_paper-data9.mat',...
                        'idw50000_gitID-dd6fa64_puuID-26fa8dde_paper-data9.mat'};
                case 1000
%                     fnames = {...
%                         'pbary1000_gitID-6ede824_puuID-bdd7211c_paper-data5.mat',...
%                         'gpr1000_gitID-6ede824_puuID-e6366570_paper-data5.mat',...
%                         'nn1000_gitID-6ede824_puuID-73c9a5b3_paper-data5.mat',...
%                         'idw1000_gitID-6ede824_puuID-ae3ba253_paper-data5.mat'};
                    fnames = {...
                        'pbary1000_gitID-9ebfa46_puuID-32707b2d_paper-data10.mat',...
                        'gpr1000_gitID-9ebfa46_puuID-4bc2b881_paper-data10.mat',...
                        'nn1000_gitID-9ebfa46_puuID-01e60698_paper-data10.mat',...
                        'idw1000_gitID-9ebfa46_puuID-e0b47f81_paper-data10.mat'};
            end
            nnames = length(fnames);
            mdllist = cell(nnames,1);
            for i = 1:nnames
                fname = fnames{i};
                load(fname,'mdl');
                mdllist{i} = mdl;
            end
        else
            mdllist = [];
        end
        %%
        [tpredlist,tsdlist,propList,methodlist,A,B] = ...
            tunnelplot(mdllist,A,B,n,'nnQ',false,'nnQ2',false,...
            'tpredlist',tpredlist,'tsdlist',tsdlist,'propList',propList,...
            'methodlist',methodlist,'extend',extend,'lgdloc','south');
    case 3
        npts = 100;
        qm = get_cubo(npts); nA = normr(rand(npts,3)); %random (qm,nA) pairs
        propList = 1:npts; %property values
        qm2 = get_cubo(npts); nA2 = normr(rand(npts,3)); %random (qm,nA) pairs
        method = 'gpr'; %interpolation method
        [propOut,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,propList,qm2,nA2,method);
    case 4
        load("SymmetrizedKernelTest_TMS.mat")
        tunnelplot(mdllist2, A, B)
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