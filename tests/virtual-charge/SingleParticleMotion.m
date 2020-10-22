
clc; close all; clear;

%% Moving Particle Test
% This script will move a single large particle and measure the E fields at
% different points in space

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

%% Medium properties
T = 300; % Temperature [K]
epsr = 80; % Relative permitivity of water
eta = 8.9e-4; % Viscosity of water [Pa*s]
sigma = 1; % Conductivity [S/m]
epsilon = eps0*epsr;

%% Particle Properties
R = 5e-9; % Particle Radius [m]
rho = 1000; % Particle Density [kg/m^3]
q = 10*qe; % Particle charge
v = 0.1; % Particle velocity [m/s]

%% Virtual charge parameters
Nvc = 10;
Ntau = 3;

%% Trajectory and initialization 
tstart = 5e-9; % Beginning of motion [s]
tend = 10e-9; % End of Motion [s]

dt = 100e-12; % Time step [s]
t = [0:dt:tend]; % Time vector
Nt = length(t);
P.X = zeros(Nt,3);
P.X(:,3) = v.*heaviside(t-tstart).*(t-tstart);
Rset = [20e-9,0,0];

%% Calculate Fields

for n = 1:Nt
    Pnow.X = P.X(1:n,:);
    Pnow.q = q;
    tnow = t(1:n);
    [Edynamic(n,:), ~, phiDynamic(n)] = MyVirtualChargeEnsemble(Pnow, tnow, Rset, epsr, sigma, Nvc, Ntau);
    phiStatic = 
    Emag(n) = norm(Edynamic(n,:));
end

%% Plot Results

figure
plot(t*1e9, Emag)
xlabel('Time [ns]')
ylabel('E_{mag} [V/m]')

