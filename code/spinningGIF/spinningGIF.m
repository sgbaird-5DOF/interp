% spinningGIF(fname): makes a spinning GIF of the current plot and saves it
% Usage: make your 3D plot (using plot3(...) or scatter3(...) etc.) and
% then call SpinningGIF with the file name that you want
function spinningGIF(fname,axisOptions,dt,nv)
arguments
    fname
    axisOptions = {'on','equal','vis3d'}
    dt = 0.03
    nv.ntiles = 1
end
ntiles = nv.ntiles;
%     view(0,10)

fig = gcf;
[center,pos] = deal(cell(1,ntiles));
radius = zeros(1,ntiles);

for i = 1:ntiles
    if ntiles > 1
        nexttile(i)
    end
    
    axis(axisOptions{:})
    
    ax = gca;
    center{i} = get(ax, 'CameraTarget');
    pos{i} = get(ax, 'CameraPosition');
    radius(i) = norm(center{i}(1:2) - pos{i}(1:2));
    angles = 0:0.005*pi:2*pi;
%     set(ax,'CameraViewAngleMode','manual');
end

for ii=1:length(angles)
    angle = angles(ii);
    
    for j = 1:ntiles
        if ntiles > 1
            ax = nexttile(j);
        end
        v = [center{j}(1) + radius(j) * cos(angle),...
            center{j}(2) + radius(j) * sin(angle),...
            pos{j}(3)];
        view(v)
        ax.CameraTarget = [0,0,0];
    end
    
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1
        imwrite(imind,cm,fname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime', dt);
    end
end
end