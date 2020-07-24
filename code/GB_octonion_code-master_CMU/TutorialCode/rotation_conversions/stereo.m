function [Theta,R] = stereo(aa)
%Stereographic projection of rotation axes
% n x 3 matrix, columns are unit vectors [ax ay az]
n = length(aa(:,1));
% r = zeros(1,n);
R = zeros(1,n);
Theta = zeros(1,n);
%Phi = zeros(1,n);

for i = 1:n
    a = aa(i,1:3);
    r = norm(a);
    if a(1) < 0
        Theta(i) = atan(a(2)/a(1)) + pi;
    else 
        Theta(i) = atan(a(2)/a(1));
    end 
    phi = acos(a(3)/r);
    R(i) = sin(pi-phi)/(1-cos(pi-phi));
end 

R(isnan(R)) = 0;
Theta(isnan(Theta)) = 0;
% 
% figure
% pax = polaraxes;
% polarscatter(Theta,R)
% title('stereographic projection of rotation axes')
% pax.ThetaLim = [0 90];
% pax.RLim = [0 2];
% 
% figure
% polarhistogram(Theta)
% 
% figure
% histogram(R)
% figure
% histogram2(Theta,R)
% xlabel('Theta')
% ylabel('R')

end

