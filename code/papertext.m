function papertext(num,Interpreter,NV)
arguments
    num(1,1) double
    Interpreter = 'latex'
    NV.xypos(1,2) double = [0.025,0.95]
end
% PAPERTEXT  generate text() for the figure labels '(a)', '(b)', etc.
x = NV.xypos(1);
y = NV.xypos(2);
charlbl = get_charlbl();

text(x,y,charlbl{num},'Units','normalized','FontSize',12,'Interpreter',Interpreter)
end