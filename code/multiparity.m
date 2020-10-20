function multiparity(parity,ID,NV)
arguments
    parity(:,1) cell
    ID
    NV.charlblQ(1,1) logical = true
    NV.titleQ(1,1) logical = true
end

fig=figure;
fig.Position = [418.6 194.6 644.8 555.2];

tiledlayout('flow','TileSpacing','compact','Padding','compact')
nIDs = length(ID);
for i = 1:nIDs
    nexttile
    ptmp = parity{i};
    if NV.charlblQ
        alphabet = ('a':'z').';
        chars = num2cell(alphabet(1:nIDs));
        chars = chars.';
        charlbl = strcat('(',chars,')');
    else
        charlbl = repelem({''},length(ID));
    end
    if NV.titleQ
        t = ID;
    else
        t = repelem({''},length(ID));
    end
    plottype = 'hex';
    parityplot(ptmp.ytrue,ptmp.ypred,plottype,'charlbl',charlbl{i},'title',t{i},'scatterOpts',struct('MarkerEdgeAlpha',0.1))
end