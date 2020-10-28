function el = po2el(po)
arguments
    po(:,1) double
end
% PO2EL  convert polar angles to elevation angles
el = 90-po;
end