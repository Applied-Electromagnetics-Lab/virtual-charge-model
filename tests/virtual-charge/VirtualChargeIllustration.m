clc; close all; clear;

% Colors
Red = [255,102,102]./255;
Blue = [102,102,255]./255;

LightingType = 'gouraud';
EdgeColor = 'none';
% Transparency value
alphavalue = 0.5;

%% Media properties and physical constants

qe = -1.60217662e-19; % Charge of electron [C]
kB = 1.38064852e-23; % Boltzman Constant [J/K]
eps0 = 8.85418782e-12; % Permitivity of free space [F/m]

epsr = 80;
sigma = 1;

eps = eps0*epsr;

wp = sigma / eps;
taup = 1/wp;

%% Trajectory

v = 20; % velocity [m/s]
tfinal = taup*3;
Nt = 100;

t = linspace(0,tfinal,Nt);
xp = 5e-9 .* sin(t./0.3e-9);
yp = 0 .* t;
zp = v .* t;

R = [xp;yp;zp].';

%% Virtual Charge Parameters

Nvc = 8;
Ntau = 3;
tauvc = Ntau * taup / Nvc; 

%% Locations and Sizes of Virtual Charges

I = (1:1:Nvc);
qvc = exp(-I.*tauvc/taup);
% The vector schould be scaled such that q + sum(qvc) = 0
scale = 1 / sum(qvc);
qvc = qvc.*scale;

tNOW = t(end);

for j = 1:Nvc
    % Calculate how long ago the virtual charge was created
    tvc = tNOW - j*tauvc;
    % If the virtual charge time is beyond the scope of the time
    % vector (i.e. it is from before t = 0) then we give it the
    % initial location. 
    if tvc < 0
        CV(j).Rvcj = R(1,:);
    else
        CV(j).Rvcj = interp1(t,R,tvc);
    end
    
    CV(j).q = qvc(j);
end

%% Plot Charges

% Radius of sphere representing tubulin monomer
SphereRadius = 5e-9;

% Generate Sphere
[X0,Y0,Z0] = sphere;

% Rezise sphere
X = X0.*SphereRadius;
Y = Y0.*SphereRadius;
Z = Z0.*SphereRadius;

% Plot Primary Charge
Rprimary = R(end,:);

% Translate Sphere
X_primary = X + Rprimary(1);
Y_primary = Y + Rprimary(2);
Z_primary = Z + Rprimary(3);
    
% Plot Spheres
hSurface=surf(X_primary,Y_primary,Z_primary);
 hold on
set(hSurface,'FaceColor',Red, ...
  'FaceAlpha',1,'FaceLighting',LightingType,'EdgeColor',EdgeColor)

for i = 1:Nvc
        
    % Rezise sphere
    X = X0.*SphereRadius*qvc(i)^(1/3);
    Y = Y0.*SphereRadius*qvc(i)^(1/3);
    Z = Z0.*SphereRadius*qvc(i)^(1/3);

    Xvc = X + CV(i).Rvcj(1);
    Yvc = Y + CV(i).Rvcj(2);
    Zvc = Z + CV(i).Rvcj(3);
    
    hSurface=surf(Xvc,Yvc,Zvc);
    hold on
    set(hSurface,'FaceColor',Blue, ...
    'FaceAlpha',alphavalue,'FaceLighting',LightingType,'EdgeColor',EdgeColor)
end

plot3(xp,yp,zp,'k--')

axis equal
set(gca,'visible','off')
camlight
view([0,-1,0])
  
