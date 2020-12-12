function po = el2po(el,radQ)
arguments
   el(:,1) double 
   radQ(1,1) logical = true
end
% EL2PO  elevation angles to polar angles (default radians)
if ~radQ
    po = 90 - el;
else
    po = pi/2 - el;
end
end