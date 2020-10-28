function sphplot(n,scl,NV)
arguments
    n = 80;
    scl = 0.8;
    NV.axview = [15.314243969189249,29.158682228660489];
end
% SPHPLOT  plot a simple sphere for visualization purposes
[x,y,z]=sphere(n);
x = scl*x;
y = scl*y;
z = scl*z;
surf(x,y,z,'EdgeColor','none');
ax=gca;
ax.View = NV.axview;
shading interp
axis equal

end