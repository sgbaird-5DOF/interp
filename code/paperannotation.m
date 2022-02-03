function paperannotation(num,Interpreter,NV)
arguments
    num(1,1) double
    Interpreter = 'latex'
    NV.xypos(1,2) double = [0.05,0.85]
    NV.fontsize double = 12
    NV.ax = gca
end
% PAPERTEXT  do a text() command for the figure labels '(a)', '(b)', etc. for num == 1, 2, etc.
x = NV.xypos(1);
y = NV.xypos(2);
ax = NV.ax;
h = 0.1;
w = 0.1;
charlbl = get_charlbl();

if num > 1
    x = x + ax.Position(1) - w;
end

annotation('textbox',[x,y,w,h],'String',charlbl{num},'Units','normalized',...
    'FontSize',NV.fontsize,'Interpreter',Interpreter,'EdgeColor','none')
end