function papertext(num,Interpreter,NV)
arguments
    num(1,1) double
    Interpreter = 'latex'
    NV.xypos(1,2) double = [0.025,0.95]
    NV.fontsize double = 12
end
% PAPERTEXT  do a text() command for the figure labels '(a)', '(b)', etc. for num == 1, 2, etc.
x = NV.xypos(1);
y = NV.xypos(2);
charlbl = get_charlbl();

text(x,y,charlbl{num},'Units','normalized','FontSize',NV.fontsize,'Interpreter',Interpreter)
end