function [width] = Beamwidth_dis(P1_x, P1_y, beamwidth, xdis, ydis)
% calculate beam width (km) at each gridpoint 
% P1_x, P1_y: (x,y) coordinate for radar location
% beamwidth (degrees)
% xdis: 1D xdimension
% ydis: 1D y dimension

[xdis2,ydis2]=meshgrid(xdis,ydis);

a2=(xdis2-P1_x).^2+(ydis2-P1_y).^2; 

a=sqrt(a2);

width=a .* tand(beamwidth/2) * 2;




