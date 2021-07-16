function spinGIF(folder,fname,nv)
arguments
    folder = '.'
    fname = 'test.gif'
    nv.ntheta = 400
    nv.dt = 0.05
    nv.m = 1
    nv.n = 1
end

ntheta = nv.ntheta;
dt = nv.dt;
m = nv.m;
n = nv.n;
nsubplots = m*n;
fpath = fullfile(folder,[fname '.gif']);

fig = gcf;

for i = 1:ntheta
    for j = 1:nsubplots
        if nsubplots > 1
            subplot(m,n,j)
        end
        camorbit(360/ntheta,0)
    end
    drawnow
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,fpath,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,fpath,'gif','WriteMode','append','DelayTime', dt);
    end
end