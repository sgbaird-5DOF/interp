function po = el2po(el)
arguments
   el(:,1) double 
end
% EL2PO  elevation angles to polar angles
po = 90 - el;
end