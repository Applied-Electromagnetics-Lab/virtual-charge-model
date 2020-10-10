
clc; close all; clear;

% This script is intended to test the normal brownian dynamics code by
% comparing data to that generated in that generated in the paper:
% A Combined Wormlike Chain and Bead Model for Dynamics Simulations of Long
% Linear DNA -- Jian 1997

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
epsr = 80; % Relative permitivity of water
eta = 8.9e-4; % Viscosity of water [Pa*s]
sigma = 1; % Conductivity [S/m]

eps = eps0*epsr;

%% Model Parameters

% DNA parameters
L0 = 5e-9; % Natural bond length [m]
h = 100*kB*T/L0^2; % Stretching force constant [J/m^2]
g = 9.82 * kB*T; % Bending rigidity constant [J]
nu = 0.243 * 1e10 * qe; % Linear charge density [C/m]
rD = 3.07e-9; % Debye length [m]
a = 1.59e-9; % Bead Radius [m]

% Time parameters
dt = 1000e-12; % Time step [s]
Nb = 41;
Nt = 10000;
t = [0:dt:dt*(Nt-1)];

% Initialize DNA location
Xinitial = zeros(Nb,3);
for i = 1:Nb
    Xinitial(i,3) = (i-1)*L0;
end

% External Field
% Because this is just comparing normal BD, the external field is set to
% zero. The virtual charge code will not activate if the external field is
% not present. 
E0 = [0,0,0]; % External Field amplitude [V/m]
E0 = 0;
f = 1e9; % External Field frequency [Hz]
w = 2*pi*f; % Angular frequency
Eext = sin(w.*t).' * E0;

% Virtual Charge Parameters
Ntau = 3;
Nvc = 10;

% Number of Trials
Ntrials = 1;

% Enter the Time rate for your simulation here. Since this will vary
% between machines and simulation parameters, it's suggested you run one
% quick trial and find how long it takes. Then enter the seconds per time
% steps in the variable below.
TimeRate = 62 / 1000; % Seconds / time step

%% Load all parameters into the Param Structure and initialize variables
Params.L0 = L0;
Params.h = h;
Params.g = g;
Params.nu = nu;
Params.rD = rD;
Params.a = a;
Params.dt = dt;
Params.Nb = Nb;
Params.Nt = Nt;
Params.epsr = epsr;
Params.sigma = sigma;
Params.T = T;
Params.eta = eta;
Params.Ntau = Ntau;
Params.Nvc = Nvc;
Params.Eext = Eext;
Params.Xinitial = Xinitial;

% Initialize Variables
D = zeros(Ntrials,Nt-1);
Re2e = zeros(Ntrials,1);
L = zeros(Ntrials,Nb-1);

% Compute how much time the simulation will take
TrialLength = TimeRate*Nt;
TrialH = floor(TrialLength/3600);
TrialM = floor(mod(TrialLength,3600)/60);
TrialS = mod(TrialLength,60);

TotalLength = TrialLength*Ntrials;
TotalH = floor(TotalLength/3600);
TotalM = floor(mod(TotalLength,3600)/60);
TotalS = mod(TotalLength,60);

% Print out simulation prediction
fprintf("The simulation will take approximately %i hours, %i minutes, %4.2f seconds per trial. \nIt will take approximately a total of %i hours, %i minutes, %4.2f seconds to compute all %i trials.\n",...
    TrialH, TrialM, TrialS, TotalH, TotalM, TotalS, Ntrials);

%% Run the simulation. Print out each trial
for Trial = 1:Ntrials
    tic
    [P, Re2e, L] = My_DNA_BD(Params);
    toc
    Trial
end

%% Plot Results and compare to paper

Lmean = mean(L);
Dmean = mean(D);

% save("June_23rd_N100_10us_Dij_2");
figure
histogram(L*1e9,30);
xlabel('Bond Length [nm]')

figure
histogram(Re2e*1e9,30);
xlabel('End to End Distance [nm]')

figure
plot(t*1e6,Re2e*1e9);
ylabel('End to End Distance [nm]')
xlabel('Time [\mu s]')

figure
plot(t(2:end)*1e6,Dmean.*1e12);
xlabel('Time \mu s')
ylabel('Diffusion Coefficient (\mu m)^2 / s')

