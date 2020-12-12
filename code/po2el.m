function el = po2el(po,radQ)
arguments
    po(:,1) double
    radQ(1,1) logical = true
end
% PO2EL  convert polar angles to elevation angles (default is radians)
if ~radQ
    %degrees
    el = 90-po;
else
    %radians
    el = pi/2-po;
end
end