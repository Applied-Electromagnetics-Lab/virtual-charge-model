clc; close all; clear;

% Spatial Coordinates
Nx = 100;
Ny = 100;

DX = 10e-9;
DY = 10e-9;

X = linspace(-DX/2,DX/2,Nx);
Y = linspace(-DY/2,DY/2,Ny);
grid = meshgrid(X,Y);
Z0 = 1e-9;

% Material Properties
eps0 = 8.85418e-12;
epsr = 80;
sigma = 1;
eps = epsr*eps0;
wp = sigma / eps; % Plasma frequency [rad / s]
taup = 1/wp; % Plasma time constant [s]

% Path properties
f = 1e9;
T = 1/f;
w = 2*pi*f;
A = 2.5e-9;
q = 1;

% Time Vectors
Nt = 200;
t = linspace(0,2*T,Nt);
Pc = ([1;0;0] * A.*sin(w.*t)).';

% Generate a vector of virtual charges length Nvc. The charges should
% decay according to the plasma time constant. The decay per charge is
% determined by the virtual charge time constant.
Nvc = 20;
Ntau = 3;
tauvc = Ntau * taup / Nvc; 

i = (1:1:Nvc);
qvc = exp(-i.*tauvc/taup);
% The vector schould be scaled such that q + sum(qn) = 0
scale = -q / sum(qvc);
qvc = qvc.*scale;

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '1GHz.gif';

c1 = 'r';
c2 = 'b';
sz = 500;
alpha = 0.4;

for i = 1:Nt
    
    scatter(Pc(i,1),Pc(i,2),sz*q,c1,'filled','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
    xlim([-DX/2, DX/2]);
    ylim([-DY/2, DY/2]);
    
    hold on
    for n = 1:Nvc
        tvc = t(i) - tauvc*n;
        if tvc < 0
            Pvc = [0,0,0];
        else
            Pvc = interp1(t,Pc,tvc);
        end
        scatter(Pvc(1),Pvc(2),abs(sz*qvc(n)),c2,'filled','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha);
    end
    hold off
    
    set(gcf, 'color', 'white');
    set(gca,'visible','off')

    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif','DelayTime', 0.05, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime', 0.05, 'WriteMode','append'); 
    end 
end


