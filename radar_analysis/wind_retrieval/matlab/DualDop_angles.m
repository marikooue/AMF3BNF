function [alpha] = DualDop_angles(P1_x, P1_y, P2_x, P2_y, xdis, ydis)
% calculate angle between two lines at each gridpoint for dual base line
% P1_x, P1_y: (x,y) coordinate for point 1
% P2_x, P2_y: (x,y) coordinate for point 2
% xdis: 1D xdimension
% ydis: 1D y dimension

[xdis2,ydis2]=meshgrid(xdis,ydis);


a2=(P2_x-P1_x)^2+(P2_y-P1_y)^2;

b2=(xdis2-P2_x).^2+(ydis2-P2_y).^2;

c2=(xdis2-P1_x).^2+(ydis2-P1_y).^2;

b=sqrt(b2); c=sqrt(c2); 

cosa=(b2+c2-a2)./(2.0.*b.*c);
alpha=acosd(cosa);
