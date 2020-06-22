function projcoord = geog2projcoord(geogcoord,projsyst)
% projcoord = geog2projcoord(geogcoord,projsyst)
% Convert geographical (longitude,latitude) coordinates to projection (x,y)
% coordinates according to Lambert conformal conic projection
% Ref: https://fr.wikipedia.org/wiki/Projection_conique_conforme_de_Lambert
% Ref: https://www.linz.govt.nz/data/geodetic-system/coordinate-conversion/projection-conversions/lambert-conformal-conic-geographic
% Test values:
% projcoord = [840490 6522200; 840500 6522210; 840510 6522220; 840520 6522230; 840530 6522240; 840540 6522250; 840550 6522260; 840560 6522270; 840570 6522280; 840580 6522290; 840590 6522300; 840600 6522310; 840610 6522320; 840620 6522330; 840630 6522340];
% geogcoord = [4.8083731 45.7849386; 4.8085047 45.7850266; 4.8086364 45.7851145; 4.808768 45.7852025; 4.8088996 45.7852905; 4.8090312 45.7853784; 4.8091629 45.7854664; 4.8092945 45.7855543; 4.8094261 45.7856423; 4.8095577 45.7857303; 4.8096894 45.7858182; 4.809821 45.7859062; 4.8099526 45.7859942; 4.8100842 45.7860821; 4.8102159 45.7861701];
% projsyst = 'Lambert93';
%
% INPUTS
%---- geogcoord : n-by-2 matrix, [lambda phi] values (long,lat)
%---- projsyst  : string, projection name 'Lambert93', ...
%
% OUTPUTS
%---- projcoord : n-by-2 matrix, [x y] values (easting,northing)

if strcmp(projsyst,'Lambert93')
    a = 6378137.0; % equatorial radius [m]
    f = 1/298.257223563; % flattening
    lambda0 = 3; % reference longitude [deg]
    phi0 = 46.5; % reference latitude [deg]
    x0 = 700000; % reference easting (abscissa) [m]
    y0 = 6600000; % reference northing (ordinate) [m]
    phi1 = 44; % latitude of first standard parallel [deg]
    phi2 = 49; % latitude of second standard parallel [deg]
elseif strcmp(projsyst,'Lambert08')
    a = 6378137.0; % equatorial radius [m]
    f = 1/298.257223563; % flattening
    lambda0 = 4 + 21/60 + 33.177/3600;
    phi0 = 50 + 47/60 + 52.134/3600;
    x0 = 649328;
    y0 = 665262;
    phi1 = 49 + 50/60;
    phi2 = 51 + 10/60;
elseif strcmp(projsyst,'ETRS89-LCC')
    a = 6378137.0; % equatorial radius [m]
    f = 1/298.257223563; % flattening
    lambda0 = 10;
    phi0 = 52;
    x0 = 4000000;
    y0 = 2800000;
    phi1 = 35;
    phi2 = 65;
end

% First eccentricity
e2 = 2*f - f^2;
e = sqrt(e2);

m = @(phi) cos(phi)./sqrt(1 - e2.*sin(phi).^2);
t = @(phi) tan(pi/4 - phi/2)./((1 - e.*sin(phi))./(1 + e.*sin(phi))).^(e/2);

n = (log(m(pi/180*phi1)) - log(m(pi/180*phi2)))/(log(t(pi/180*phi1)) - log(t(pi/180*phi2)));
F = m(pi/180*phi1)/(n*t(pi/180*phi1)^n);

rho = @(phi) a.*F.*t(phi).^n;

% From geographic coordinates (lambda,phi) to projection coordinates (x,y)
Ncoord = size(geogcoord,1);
projcoord = zeros(Ncoord,2);
projcoord(:,1) = rho(pi/180*geogcoord(:,2)).*sin(n.*pi/180.*(geogcoord(:,1) - lambda0));
projcoord(:,2) = rho(pi/180*phi0) - rho(pi/180*geogcoord(:,2)).*cos(n.*pi/180.*(geogcoord(:,1) - lambda0));
projcoord = projcoord + [x0 y0];

end