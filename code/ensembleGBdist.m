function dmin = ensembleGBdist(o,o2,K,dtype,nv)
arguments
   o(:,8) double = get_ocubo(100)
   o2(:,8) double = []
   K(1,1) double = 10
   dtype char = 'euclidean'
   nv.pdtype char {mustBeMember(nv.pdtype,{'dist','pdist','pdist2'})} = 'pdist'
   nv.orefs = get_orefs(K)
   nv.disp {mustBeLogical} = true
   nv.pgnum(1,1) double {mustBeInteger} = 32
end
%% ensemble grain boundary distance (take minimum)
pdtype = nv.pdtype;
orefs = nv.orefs;
dispQ = nv.disp;
pgnum = nv.pgnum;
% d = cell(1,K);

% if strcmp(pdtype,'dist')
switch dtype
    case 'euclidean'
        fn = @(a,b) norm(a-b);
    case 'arclength'
        fn = @get_alen;
    case 'omega'
        fn = @get_omega;
end
% end
for i = 1:K
    oref = orefs(i,:);
    
    if dispQ
        disp(['K: ' int2str(i)])
    end
    osym = get_octpairs(o,'oref',oref,'dispQ',false); %symmetrize w.r.t. to oref
    osym = normr(osym);
    
    if ~isempty(o2)
        osym2 = get_octpairs(o2,'oref',oref,'dispQ',false,'pgnum',pgnum);
        osym2 = normr(osym2);
    end
    
    switch pdtype
        case 'dist'
            d = fn(osym,osym2);
        case 'pdist'
            d = squareform(pdist(osym,fn));
        case 'pdist2'
            d = squareform(pdist2(osym,osym2,fn));
    end
    if i == 1
        dmin = d;
    end
    dmin = min(d,dmin);
end

end

%% CODE GRAVEYARD
%{
%     if i == 1
%         oref(i,:) = get_ocubo(1,'random',[],10); %just because I use this frequently
%     else
%         oref(i,:) = get_ocubo(); %random reference octonion
%     end

%}