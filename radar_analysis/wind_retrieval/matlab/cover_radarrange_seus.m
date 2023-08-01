%%%
% Plot a map of radar locations, distances, and multi-Doppler domain
clear all;
%- Finxed variables
Rearth = 6378.1; %km
%beamwidth
bw_nexrad=0.9;
bw_csapr=1.0;

RadarLocations;

Radar(1).lat=lat_UAH; Radar(1).lon=lon_UAH; Radar(1).name='UAH'; 
Radar(2).lat=lat_KHTX; Radar(2).lon=lon_KHTX; Radar(2).name='KHTX'; 
Radar(3).lat=lat_KGWX; Radar(3).lon=lon_KGWX; Radar(3).name='KGWX'; 
Radar(4).lat=lat_KBMX; Radar(4).lon=lon_KBMX; Radar(4).name='KBMX'; 
Radar(5).lat=lat_KMXX; Radar(5).lon=lon_KMXX; Radar(5).name='KMXX'; 
Radar(6).lat=lat_Huntsville; Radar(6).lon=lon_Huntsville; Radar(6).name='Huntsville'; 
Radar(7).lat=lat_Fayetteville; Radar(7).lon=lon_Fayetteville; Radar(7).name='Fayetteville'; 
Radar(8).lat=lat_CourlandAirport; Radar(8).lon=lon_CourlandAirport; Radar(8).name='CourlandAirport'; 
Radar(9).lat=lat_BlackWarriorWorkCenter; Radar(9).lon=lon_BlackWarriorWorkCenter; Radar(9).name='BlackWarriorWorkCenter'; 
Radar(10).lat=lat_Supplemental; Radar(10).lon=lon_Supplemental; Radar(10).name='Supplemental'; 
Radar(11).lat=lat_NEONMAYF; Radar(11).lon=lon_NEONMAYF; Radar(11).name='NEONMAYF'; 
Radar(12).lat=lat_Whitsitt; Radar(12).lon=lon_Whitsitt; Radar(12).name='Whitsitt'; 
Radar(13).lat=lat_RSA; Radar(13).lon=lon_RSA; Radar(13).name='RSA'; 
Radar(14).lat=lat_Cullman; Radar(14).lon=lon_Cullman; Radar(14).name='CullmanAirport'; %Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; 
Radar(15).lat=lat_CSAPR2; Radar(15).lon=lon_CSAPR2; Radar(15).name='CSAPR2'; Radar(15).bw=bw_csapr;%Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; Radar(13).bw=bw_csapr;
Radar(16).lat=lat_CSU; Radar(16).lon=lon_CSU; Radar(16).name='CSU C band'; Radar(16).bw=0.95;%Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; 
Radar(17).lat=lat_SEUSCSAPR2; Radar(17).lon=lon_SEUSCSAPR2; Radar(17).name='CSAPR2'; 
Radar(18).lat=lat_SACRsup; Radar(18).lon=lon_SACRsup; Radar(18).name='SACR candidate3'; 
%-- cities
lat_Nashville = 36.16266; lon_Nashville = -86.7816;
Cities=kml2struct('MajorCities.kml');
%-- map
[txlat,txlon]=borders('Texas');
[allat,allon]=borders('Alabama');
[tnlat,tnlon]=borders('Tennessee');

%%
load_radar_dist_data=1;
%--
lon_radar1=lon_KBMX;
lat_radar1=lat_KBMX;

lon_radar2=lon_KMXX;
lat_radar2=lat_KMXX;

lon_radar1=lon_UAH;
lat_radar1=lat_UAH;
lon_radar2=lon_KHTX;
lat_radar2=lat_KHTX;
lon_radar2=lon_SEUSCSAPR2;
lat_radar2=lat_SEUSCSAPR2;
%lon_radar2=lon_Cullman;
%lat_radar2=lat_Cullman;

% lon_radar1=lon_CSAPR2;
% lat_radar1=lat_CSAPR2;
% lon_radar2=lon_CSU;
% lat_radar2=lat_CSU;
maxrange1=122;
maxrange2=150;

bw_radar1=bw_csapr;
bw_radar2=0.95;

standard_range_km = 80;
retrieval_angle_deg = 10;

dlat=0.008993;
dlon=360.0 / (2.0 * pi * Rearth * cos(lat_radar1 /180.0 * pi));

ixlon=-80:1:80; dxkm=2.0; xlon = ixlon * dxkm * dlon + lon_radar1;
iylat=-80:1:80; dykm=2.0; ylat = iylat * dykm * dlat + lat_radar1;
izalt=0:1:15; dzkm=1.0; zalt = izalt * dzkm;
radarP1x=0;radarP1y=0;
radarP2x=(lon_radar2-lon_radar1)/(dxkm * dlon);
radarP2y=(lat_radar2-lat_radar1)/(dxkm * dlat);
%zalt=[0.5,5,10];
[xlon2,ylat2]=meshgrid(xlon,ylat);
[xlon3,ylat3,zalt3]=meshgrid(xlon,ylat,zalt);

%%
dist_from_rad1 = sqrt(((xlon2-lon_radar1)/dlon).^2.0 + ((ylat2-lat_radar1)/dlat).^2.0);
dist_from_rad2 = sqrt(((xlon2-lon_radar2)/dlon).^2.0 + ((ylat2-lat_radar2)/dlat).^2.0);

%-angle
dist_rad1rad2 = sqrt(((lon_radar2-lon_radar1)/dlon).^2.0 + ((lat_radar2-lat_radar1)/dlat).^2.0);
angle_rad1rad2 = asind((lat_radar2-lat_radar1)./dist_rad1rad2);
angle_rad2rad1 = asind((lat_radar1-lat_radar2)./dist_rad1rad2);
angle_from_rad1 = asind(((ylat2-lat_radar1)/dlat)./dist_from_rad1) + angle_rad1rad2;
angle_from_rad2 = asind(((ylat2-lat_radar2)/dlat)./dist_from_rad2) + angle_rad2rad1;
%-
cover_rad1rad2=zeros(size(dist_from_rad1));
cover_rad1rad2(dist_from_rad1<standard_range_km & dist_from_rad2<standard_range_km)=1;
cover_rad1rad2_all=cover_rad1rad2;
cover_rad1rad2(abs(angle_from_rad1)<retrieval_angle_deg & abs(angle_from_rad2)<retrieval_angle_deg)=0;
area_rad1rad2_radar=zeros(size(cover_rad1rad2));

for i=1:length(ixlon)
    for j=1:length(iylat)
        distkm=sqrt(((xlon2-xlon(i))/dlon).^2.0 + ((ylat2-ylat(j))/dlat).^2.0);
        cover1=zeros(size(distkm));
        cover1(distkm<standard_range_km)=1;
        cover3=cover1+cover_rad1rad2; num_grids=sum(sum(cover3==2));
        area_rad1rad2_radar(i,j) = num_grids * (dxkm*dykm);
    end
end
arearatio_rad1rad2_radar=area_rad1rad2_radar/(sum(sum(cover_rad1rad2))* (dxkm*dykm));

%% beta angles
beta=20;reso=1.0;%km
beamwidth1=bw_radar1;
beamwidth2=bw_radar2;
alpha_km=DualDop_angles(radarP1x*dxkm,radarP1y*dykm,radarP2x*dxkm,radarP2y*dykm,ixlon*dxkm,iylat*dykm); 
[width1] = Beamwidth_dis(radarP1x*dxkm,radarP1y*dykm, beamwidth1, ixlon*dxkm, iylat*dykm);
[width2] = Beamwidth_dis(radarP2x*dxkm,radarP2y*dykm, beamwidth2, ixlon*dxkm, iylat*dykm);
%area=(alpha_70km>beta & alpha_70km<(180-beta) & width1<1 & width2<1);

figure
pcolor(xlon,ylat,alpha_km); shading flat; colorbar;
hold on;
scatter([lon_radar1],[lat_radar1],20,'ok','filled');
scatter([lon_radar2],[lat_radar2],20,'ok','filled');
contour(xlon,ylat,alpha_km,[beta,180-beta],'-k','LineWidth',2,'DisplayName',['beta=',num2str(beta)]); 
contour(xlon,ylat,width1,[0,reso,reso+100],'-.k','LineWidth',1,'DisplayName',['bw=',num2str(reso)]); 
contour(xlon,ylat,width2,[0,reso,reso+100],'-.k','LineWidth',1); 
contour(xlon,ylat,dist_from_rad1,[maxrange1,maxrange1],'-.k','LineWidth',1,'DisplayName',['max range=',num2str(maxrange1)]); 
contour(xlon,ylat,dist_from_rad2,[maxrange2,maxrange2],'-.k','LineWidth',1,'DisplayName',['max range=',num2str(maxrange2)]); 
xlabel('Lon');ylabel('Lat');title(['Baseline=',num2str(dist_rad1rad2),' km ','Solid line: beta=',num2str(beta),'^o Dashed line: beamreso=',num2str(reso),'km ']);
plot(allon,allat,'-k','LineWidth',1);
plot(tnlon,tnlat,'-k','LineWidth',1);
plot(txlon,txlat,'-k','LineWidth',1);
for c=1:length(Cities)
    if(Cities(c).Lon>min(xlon) & Cities(c).Lon<max(xlon) & Cities(c).Lat>min(ylat) & Cities(c).Lat<max(ylat))
    scatter([Cities(c).Lon],[Cities(c).Lat],5,'ow','LineWidth',2,'DisplayName',Cities(c).Name)
    text([Cities(c).Lon],[Cities(c).Lat],Cities(c).Name,'FontSize',10)
    end
end
%-radar option
%for c=1:length(Radar)
for c=15:16
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    scatter([Radar(c).lon],[Radar(c).lat],20,'sr','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end
end
hold off

%% Save
save(['cover_radarrange_',standard_range_km,'km'], ...
    'xlon','ylat','arearatio_rad1rad2_radar','cover_rad1rad2','cover_rad1rad2_all',...
    'area_rad1rad2_radar')
%% radar corrdinate
nexrad_vcp_rangekm=2:0.25:180;
nexrad_vcp_azimuth=0:0.5:360;
nexrad_vcp_elevation=[0.48,0.88,1.32,1.8,2.42,3.12,4.0,5.1,6.42,8.0,10.02,12.48,15.6,19.51];
[nexrad_r3,nexrad_az3,nexrad_el3]=meshgrid(nexrad_vcp_rangekm,nexrad_vcp_azimuth,nexrad_vcp_elevation);
nexrad_vcp_z = nexrad_r3 .* sind(nexrad_el3);
nexrad_vcp_y = nexrad_r3 .* cosd(nexrad_el3).* sind(nexrad_az3);
nexrad_vcp_x = nexrad_r3 .* cosd(nexrad_el3).* cosd(nexrad_az3);
nexrad_vcp_lat = nexrad_vcp_y * dlat + lat_radar1;
nexrad_vcp_lon = nexrad_vcp_x * dlon + lon_radar1;
nexrad2_vcp_lat = nexrad_vcp_y * dlat + lat_radar2;
nexrad2_vcp_lon = nexrad_vcp_x * dlon + lon_radar2;

csapr_vcp_rangekm=0.06:0.12:150;
csapr_vcp_azimuth=0:1.0:360;
csapr_vcp_elevation=[0.75,1.2,1.9,2.6,3.5,4.4,5.3,6.4,7.8,9.6,11.7,14.3,17.5,21.4,26.1,33.,42.];
[csapr_r3,csapr_az3,csapr_el3]=meshgrid(csapr_vcp_rangekm,csapr_vcp_azimuth,csapr_vcp_elevation);
csapr_vcp_z = csapr_r3 .* sind(csapr_el3);
csapr_vcp_y = csapr_r3 .* cosd(csapr_el3).* sind(csapr_az3);
csapr_vcp_x = csapr_r3 .* cosd(csapr_el3).* cosd(csapr_az3);
csapr_vcp_lat1 = csapr_vcp_y * dlat + lat_radar1;
csapr_vcp_lon1 = csapr_vcp_x * dlon + lon_radar1;
csapr_vcp_lat2 = csapr_vcp_y * dlat + lat_radar2;
csapr_vcp_lon2 = csapr_vcp_x * dlon + lon_radar2;

csapr_rhi_rangekm=0.06:0.12:150;
csapr_rhi_azimuth=0:1.0:360;
csapr_rhi_elevation=0:1:30;
[rhi_r3,rhi_az3,rhi_el3]=meshgrid(csapr_rhi_rangekm,csapr_rhi_azimuth,csapr_rhi_elevation);
csapr_rhi_z = rhi_r3 .* sind(rhi_el3);
csapr_rhi_y = rhi_r3 .* cosd(rhi_el3).* sind(rhi_az3);
csapr_rhi_x = rhi_r3 .* cosd(rhi_el3).* cosd(rhi_az3);
csapr4_rhi_lat = csapr_rhi_y * dlat + lat_LaPorte;
csapr4_rhi_lon = csapr_rhi_x * dlon + lon_LaPorte;
csapr2_rhi_lat = csapr_rhi_y * dlat + lat_radar2;
csapr2_rhi_lon = csapr_rhi_x * dlon + lon_radar2;

%%
if(load_radar_dist_data)
    load('dist_radars_seus.mat');
else
dist_radarcor_csapr1=zeros(size(area_rad1rad2_radar));
dist_radarcor_csapr2=zeros(size(area_rad1rad2_radar));

for i=1:length(ixlon)
    for j=1:length(iylat)
    for k=1:length(zalt)
        distkm=sqrt(((csapr_vcp_lon1-xlon(i))/dlon).^2.0 + ((csapr_vcp_lat1-ylat(j))/dlat).^2.0 + (csapr_vcp_z-zalt(k)).^2.0);
        dist_radarcor_csapr1(i,j,k)=min(min(min(distkm)));
        distkm=sqrt(((csapr_vcp_lon2-xlon(i))/dlon).^2.0 + ((csapr_vcp_lat2-ylat(j))/dlat).^2.0 + (csapr_vcp_z-zalt(k)).^2.0);
        dist_radarcor_csapr2(i,j,k)=min(min(min(distkm)));
    end
    end
end

dist_sum_2radars=dist_radarcor_csapr1+dist_radarcor_csapr2;
dist_radarcor_NEONMAYF=dist_radarcor_csapr1;
dist_radarcor_Whitsitt=dist_radarcor_csapr2;
%Save
 save('dist_radars_seus', ...
     'xlon','ylat','dist_radarcor_NEONMAYF','dist_radarcor_Whitsitt',...
     'dist_sum_2radars')

end
%% plot
%%
figure
subplot(1,3,1)
%pcolor(xlon,ylat,area_rad1rad2_radar'); shading flat, colorbar;
pcolor(xlon,ylat,arearatio_rad1rad2_radar'*100); shading flat, colorbar;
title([num2str(standard_range_km),' km coverage % by 3rd radar'])
hold on;
%scatter([lon_radar1],[lat_radar1],20,'xk','LineWidth',4,'DisplayName','KHGX')
%scatter([lon_radar2],[lat_radar2],20,'og','LineWidth',4,'DisplayName','UAH')
%text(lon_radar1,lat_radar1,'KHGX','FontSize',12)
%text(lon_radar2,lat_radar2,'UAH','FontSize',12)
%scatter([lon_radar1],[lat_radar1],20,'xk','LineWidth',4,'DisplayName','KBMX')
%scatter([lon_radar2],[lat_radar2],20,'og','LineWidth',4,'DisplayName','KMXX')
text(lon_radar1,lat_radar1,'KBMX','FontSize',12);text(lon_radar2,lat_radar2,'KMXX','FontSize',12)
%text(lon_radar1,lat_radar1,'KHTX','FontSize',12);text(lon_radar2,lat_radar2,'UAH','FontSize',12)
%scatter([lon_CSAPR2],[lat_CSAPR2],20,'xr','LineWidth',4,'DisplayName','CSAPR2')
%scatter([lon_LaPorte],[lat_LaPorte],20,'xb','LineWidth',4,'DisplayName','LaPorte')
contour(xlon,ylat,cover_rad1rad2_all,[0.5,10000],'-k','LineWidth',5,'DisplayName',[num2str(standard_range_km),'km range']);
plot(allon,allat,'-k','LineWidth',1);
plot(tnlon,tnlat,'-k','LineWidth',1);
for c=1:length(Cities)
    if(Cities(c).Lon>min(xlon) & Cities(c).Lon<max(xlon) & Cities(c).Lat>min(ylat) & Cities(c).Lat<max(ylat))
    scatter([Cities(c).Lon],[Cities(c).Lat],5,'ow','LineWidth',2,'DisplayName',Cities(c).Name)
    text([Cities(c).Lon],[Cities(c).Lat],Cities(c).Name,'FontSize',10)
    end
end
%-radar option
for c=1:length(Radar)
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    scatter([Radar(c).lon],[Radar(c).lat],20,'sr','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end
end
hold off
%
subplot(1,3,2)
iz=10;
pcolor(xlon,ylat,squeeze(dist_sum_2radars(:,:,iz))'); shading flat, colorbar;
title(['z=',num2str(iz),' km nearest radar gate distance']); caxis([0 2])
hold on;
scatter([lon_radar1],[lat_radar1],20,'xk','LineWidth',4,'DisplayName','KHGX')
scatter([lon_radar2],[lat_radar2],20,'og','LineWidth',4,'DisplayName','AMF')
contour(xlon,ylat,cover_rad1rad2_all,[0.5,10000],'-k','LineWidth',5,'DisplayName',[num2str(standard_range_km),'km range']);
plot(allon,allat,'-k','LineWidth',1);
plot(tnlon,tnlat,'-k','LineWidth',1);
hold off

subplot(1,3,3)
iy=63;
pcolor(xlon,zalt,squeeze(dist_sum_2radars(:,iy,:))'); shading flat, colorbar;
title(['y=',num2str(iy),' km nearest radar gate distance'])
