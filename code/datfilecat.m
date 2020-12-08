function datfilecat(fpaths,fpathout)
arguments
    fpaths cell = strcat('../../TJ2GBE/TJdata/',{...
        'triples_Ni_0131_21520.dat',...
        'triples_Ni_0319_25444.dat',...
        'triples_Ni_0331_13980.dat',...
        'triples_Ni_0707_13274.dat',...
        'triples_Ni_0717_13874.dat'})
    fpathout char = '../../TJ2GBE/TJdata/triples_Ni_cat_88092.dat'
end

nfiles = length(fpaths);
[TJs,e1,e2,e3,m1,m2,m3,nTJs] = deal(cell(nfiles,1));
for i = 1:nfiles
    fpath = fpaths{i};
    [TJs{i},e1{i},e2{i},e3{i},m1{i},m2{i},m3{i},nTJs{i}] = datfile2em(fpath,0,false,false,true);
end

%concatenate
n2cvcat = @(x) num2cell(vertcat(x{:}));
TJscat = n2cvcat(TJs);
e1cat = n2cvcat(e1);
e2cat = n2cvcat(e2);
e3cat = n2cvcat(e3);
m1cat = n2cvcat(m1);
m2cat = n2cvcat(m2);
m3cat = n2cvcat(m3);

nTJstot = sum([nTJs{:}]);

fid = fopen(fpathout,'w');
fsp = '% .4g\t% .4g\t % .4g\n';
myfprintf = @(x) fprintf(fid,fsp,x{:});
for i = 1:nTJstot
    %unpack
    TJtmp = TJscat(i,:);
    e1tmp = e1cat(i,:);
    m1tmp = m1cat(i,:);
    e2tmp = e2cat(i,:);
    m2tmp = m2cat(i,:);
    e3tmp = e3cat(i,:);
    m3tmp = m3cat(i,:);
    
    %print
    fprintf(fid,'%i\n',i);
    myfprintf(TJtmp);
    myfprintf(e1tmp);
    myfprintf(m1tmp);
    myfprintf(e2tmp);
    myfprintf(m2tmp);
    myfprintf(e3tmp);
    myfprintf(m3tmp);
end

%% CODE GRAVEYARD
%{
% fpath = '../../TJ2GBE/TJdata/triples_Ni_0131_21520.dat';
% fid = fopen(fpath);
% 
% 
% while ~feof(fid)
%     A{i} = fgetl
% end

    fprintf(fid,fsp,TJtmp);
    fprintf(fid,fsp,e1tmp);
    fprintf(fid,fsp,m1tmp);
    fprintf(fid,fsp,e2tmp);
    fprintf(fid,fsp,m2tmp);
    fprintf(fid,fsp,e3tmp);
    fprintf(fid,fsp,m3tmp);
%}