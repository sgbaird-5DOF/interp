%----colorbar------
z_label = 'D_eff (m^2/s)';
clims = [1e-15 1e-13];

cb = colorbar(ax1);
cb.Label.String=z_label;
cb.Position = [0.814,0.402,0.0288,0.461];
caxis(ax1,clims);
cb.Limits=clims;
colormap(ax1,'jet');