%%%
% Plot a map of radar coverage and multi-Doppler domain
clear all;
%- Finxed variables
Rearth = 6378.1; %km
%beamwidth
bw_nexrad=0.9;
bw_csapr=1.0;

RadarLocations;

Radar(1).lat=lat_UAH; Radar(1).lon=lon_UAH; Radar(1).name='UAH ARMOR'; 
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
Radar(14).lat=lat_Cullman; Radar(14).lon=lon_Cullman; Radar(14).name='CullmanAirport'; %Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; Radar(13).bw=bw_csapr;
Radar(15).lat=lat_CSAPR2; Radar(15).lon=lon_CSAPR2; Radar(15).name='CSAPR2'; %Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; Radar(13).bw=bw_csapr;
Radar(16).lat=lat_CSU; Radar(16).lon=lon_CSU; Radar(16).name='CSU C band'; %Radar(14).maxrange=csapr_maxrange; Radar(13).vcp=vcp_csapr; Radar(13).bw=bw_csapr;
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
%load_radar_dist_data=1;
%--
lon_radar1=lon_SEUSCSAPR2;
lat_radar1=lat_SEUSCSAPR2;
lon_radar2=lon_UAH;
lat_radar2=lat_UAH;

lon_cradar=lon_SACRsup;
lat_cradar=lat_SACRsup;

standard_range_km = 100;
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
dist_from_crad = sqrt(((xlon2-lon_cradar)/dlon).^2.0 + ((ylat2-lat_cradar)/dlat).^2.0);

% sensitivity
dBZ_min=-50+(20.0.* log10(dist_from_crad));
dBZ_min(dist_from_crad>30)=nan;

%-angle
dist_rad1rad2 = sqrt(((lon_radar2-lon_radar1)/dlon).^2.0 + ((lat_radar2-lat_radar1)/dlat).^2.0);
angle_rad1rad2 = asind((lat_radar2-lat_radar1)./dist_rad1rad2);
angle_rad2rad1 = asind((lat_radar1-lat_radar2)./dist_rad1rad2);
angle_from_rad1 = asind(((ylat2-lat_radar1)/dlat)./dist_from_rad1) + angle_rad1rad2;
angle_from_rad2 = asind(((ylat2-lat_radar2)/dlat)./dist_from_rad2) + angle_rad2rad1;
%-

%% beta angles
beta=30;reso=1.0;%km
beamwidth=1.0;
alpha_km=DualDop_angles(radarP1x*dxkm,radarP1y*dykm,radarP2x*dxkm,radarP2y*dykm,ixlon*dxkm,iylat*dykm); 
[width1] = Beamwidth_dis(radarP1x*dxkm,radarP1y*dykm, beamwidth, ixlon*dxkm, iylat*dykm);
[width2] = Beamwidth_dis(radarP2x*dxkm,radarP2y*dykm, beamwidth, ixlon*dxkm, iylat*dykm);

figure
pcolor(xlon,ylat,dBZ_min); shading flat; colorbar;caxis([-50 -15])
hold on;
scatter([lon_radar1],[lat_radar1],40,'ok','filled');
scatter([lon_radar2],[lat_radar2],40,'ok','filled');
scatter([lon_cradar],[lat_cradar],50,'<r','LineWidth',2);
contour(xlon,ylat,alpha_km,[beta,180-beta],'-k','LineWidth',2,'DisplayName',['beta=',num2str(beta)]); 
contour(xlon,ylat,width1,[0,reso,reso+100],'-.k','LineWidth',1,'DisplayName',['bw=',num2str(reso)]); 
contour(xlon,ylat,width2,[0,reso,reso+100],'-.k','LineWidth',1); 
xlabel('Lon');ylabel('Lat');title(['Baseline=',num2str(dist_rad1rad2),' km ','Solid line: beta=',num2str(beta),'^o Dashed line: beamreso=',num2str(reso),'km ']);
plot(allon,allat,'-k','LineWidth',1);
plot(tnlon,tnlat,'-k','LineWidth',1);
plot(txlon,txlat,'-k','LineWidth',1);
%- location  
 c=1;
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    scatter([Radar(c).lon],[Radar(c).lat],30,'sr','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end
 c=17;
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    scatter([Radar(c).lon],[Radar(c).lat],30,'sr','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end
 c=9;
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    scatter([Radar(c).lon],[Radar(c).lat],40,'sk','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end
 c=18;
    if(Radar(c).lon>min(xlon) & Radar(c).lon<max(xlon) & Radar(c).lat>min(ylat) & Radar(c).lat<max(ylat))
    %scatter([Radar(c).lon],[Radar(c).lat],40,'sk','filled','DisplayName',Radar(c).name)
    text([Radar(c).lon],[Radar(c).lat],Radar(c).name,'FontSize',10)
    end

hold off

%% RHI

xdis = -30:1:30;
zdis = 0:0.2:6;

rangekm = 0:1:30;
elevs = 0:0.5:180;
azims = zeros(1,length(elevs));
rangekm2=repmat(rangekm,[length(elevs),1]);
rangekm2=rangekm2';
az2=repmat(azims,[length(rangekm),1]);
el2=repmat(elevs,[length(rangekm),1]);
az_rad = az2 * pi/180.0;
el_rad = el2 * pi/180.0;

dBZ_min=-50+(20.0.* log10(rangekm2));
dBZ_min(rangekm2>30)=nan;

[xkm, ykm, zkm] = EarthCurve(az_rad, el_rad, rangekm2);
dis = sqrt(xkm.^2+ykm.^2);dis(el_rad>pi/2)=dis(el_rad>pi/2)*-1;

el2(el2>90)=180-el2(el2>90);
figure
subplot(2,1,1)
pcolor(dis,zkm,dBZ_min); shading flat; colorbar; caxis([-50 -15]);
xlabel('Sistance from radar [km]');
ylabel('Height [km]');ylim([0 15])
title('Sensitivity [dB]')
subplot(2,1,2)
pcolor(dis,zkm,el2); shading flat; colorbar; caxis([0 90]);
hold on;
contour(dis,zkm,el2,[0,15],'-k','LineWidth',2); 
hold off
xlabel('Sistance from radar [km]');
ylabel('Height [km]');ylim([0 15])
title('Elevation [deg]')
