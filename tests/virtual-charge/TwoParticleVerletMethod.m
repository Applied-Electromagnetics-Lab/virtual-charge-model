clc; clear; close all

%% Add the src code
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end-1)-1);
folderadd = strcat(newdir,filesep,"src");
addpath(genpath(folderadd));

%% Physical Constants

qe = -1.60217662e-19; % Charge of electron [C]
kB = 1.38064852e-23; % Boltzman Constant [J/K]
eps0 = 8.85418782e-12; % Permitivity of free space [F/m]

T = 300; % Temperature [K]
epsr = 80; % Relative permitivity
eta = 8.9e-4; % Viscosity of water [Pa*s]

eps = eps0*epsr;
lambdaD = 1e-9;


%% External Field
E0 = 1e6; % V/m
f = 200e6;
w = 2*pi*f;

%% Time Frame
Tperiod = 1/f;
tend = 2*Tperiod;
Nt = 500;
t = linspace(0,tend,Nt);
h = t(2)-t(1);

%% Particle parameters

Np = 2;
a = 2.5e-9; % Particle radius
q = 10.*qe;
rho = 1e3; % kg/m^3
m = rho*4/3*pi*a^3;

%% Virtual Charge parameters
Nvc = 20;
Ntau = 4;
sigma = 1;
wp = sigma / eps; % Plasma frequency [rad / s]
taup = 1/wp; % Plasma time constant [s]

Nc = length(q);

%% Intialize Positions and variables
Dist0 = 10e-9;

Ex = E0.*cos(w.*t);
Ey = Ex.*0;
Ez = Ex.*0;

Eext = [Ey;Ex;Ez].';

for i = 1:Np
    P(i).X = zeros(Nt,3);
    P(i).F = zeros(Nt,3);
    P(i).V = zeros(Nt,3);
    P(i).U = zeros(Nt,1);
    P(i).q = q;
    P(i).a = a;
    P(i).m = m;
end

P(1).X(1,1) = -Dist0/2;
P(2).X(1,1) =  Dist0/2;

%% Verlet Integration (Time Loop)

for n = 1:(Nt-1)
    for i = 1:Np
        
        %%
        for j = 1:Np
            Pnow(j).X = P(j).X(1:n,:);
            Pnow(j).q = P(j).q;
        end
        Pnow = [Pnow(1:i-1),Pnow(i+1:end)];

        [Evc, ~, ~] = MyVirtualChargeEnsemble(Pnow, t(1:n), P(i).X(n,:), epsr, sigma, Nvc, Ntau);
        Fedynamic = Evc*P(i).q;

        %--------External Elecric Field--------
        Fext = P(i).q .* Eext(n,:);
        
        Force = Fedynamic + Fext;
        
        P(i).F(n,:) = Evc*P(i).q;
        
        P(i).X(n+1,:) = P(i).X(n,:) + h.*P(i).V(n,:) + h^2./(2*P(i).m).*Force;
        P(i).V(n+1,:) = P(i).V(n,:) + h/P(i).m.*Force;
    
    end
end

%% Plot Results

for n = 1:Nt
    r1 = P(1).X(n,:);
    r2 = P(2).X(n,:);
    IPD(n) = norm(r1 - r2);
end

figure
scatter( P(1).X(:,1), P(1).X(:,2));
hold on
scatter( P(2).X(:,1), P(2).X(:,2));

figure
plot(t*1e9,IPD.*1e9);
xlabel('Time [ns]')
ylabel('Interparticle Distance [nm]')

%% Compare with Theory

wp = sigma / eps;
Amag = abs( 1 ./ (1 - 1i*wp./w) );
epsbar = eps - 1i*sigma/w;

A1 = P(1).q*E0 / (P(1).m * w^2);
A2 = P(2).q*E0 / (P(2).m * w^2);

p1 = A1*P(1).q;
p2 = A2*P(2).q;

Ehorz = abs( p1 ./ (2*pi*epsbar*IPD.^3) );
Evert = abs( p2 ./ (4*pi*epsbar*IPD.^3) );

Fvert = abs( Evert*P(1).q );
Fhorz = abs( Ehorz*P(2).q );

F1avg = ( max(P(1).F(:,1)) - min(P(1).F(:,1)) )./2;
F2avg = ( max(P(2).F(:,1)) - min(P(2).F(:,1)) )./2;

figure
bar([Fvert(1),Fhorz(1),F1avg(1),F2avg(1)])
xticklabels({'Fv','Fh','F1','F2'})

% plot(t,Fhorz)
% hold on
% plot(t,Fvert)
% plot(t,F1avg.*ones(length(t)))






