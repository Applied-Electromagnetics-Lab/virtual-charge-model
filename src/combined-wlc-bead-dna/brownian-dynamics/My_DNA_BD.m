function [P, Re2e, L, Dbulk] = My_DNA_BD(Params)
%% Description
% This function will run a brownian dynamics trial of DNA. The inputs are
% built into a structure called "Params". The fields of Params include:
% L0: Natural bond length [m]
% h: Stretching force constant [J/m^2]
% g: Bending rigidity constant [J]
% nu: Linear charge density [C/m]
% rD Debye length [m]
% a: Bead Radius [m]
% dt: Time step [s]
% Nb: Number of beads
% Nt: Number of time steps

% The outputs of the function include:
% P is the structure containing information on every particle in the
% simulation
% P(i).X is a Nt x 3 sized matrix giving the 3D location of particle i at
% every point in time [m]
% P(i).a is the radius of particle i [m]
% P(i).q is the charge of particle i [C]
% Dbulk is the translational diffusion coefficient of the DNA over time
% Re2e a vector length Nt giving the end to end distance of the DNA at the
% end of the simulation [m]
% L is a matrix size Nt x Nb giving the length of every bond of the
% simulation at every point in time

%% Model Parameters

L0 = Params.L0; % Natural bond length [m]
h = Params.h; % Stretching force constant [J/m^2]
g = Params.g; % Bending rigidity constant [J]
nu = Params.nu; % Linear charge density [C/m]
rD = Params.rD; % Debye length [m]
a = Params.a; % Bead Radius [m]

epsr = Params.epsr; % Relative Permitivity
sigma = Params.sigma; % Conductivity [S/m]
eta = Params.eta; % Viscosity [Pa*s]
T = Params.T; % Temperature [K]

dt = Params.dt; % Time step [s]
Nb = Params.Nb; % Number of beads
Nt = Params.Nt; % Number of time steps
Ntau = Params.Ntau; % Number of time constants for virtual charge calculation
Nvc = Params.Nvc; % Number of virtual charges
Xinitial = Params.Xinitial;
ComputeDynamicForce = Params.ComputeDynamicForce;

Eext = Params.Eext; % External Field

if Eext == 0
    NoDynamics = true;
else
    NoDynamics = false;
end

%% Physical Constants

qe = -1.60217662e-19; % Charge of electron [C]
kB = 1.38064852e-23; % Boltzman Constant [J/K]
eps0 = 8.85418782e-12; % Permitivity of free space [F/m]

eps = eps0*epsr;

%% Initialize Variables

t = [0:dt:dt*(Nt-1)];
D0 = kB*T / (6*pi*eta*a);
Re2e = zeros(1,Nt);
L = zeros(Nt,Nb-1);

for i = 1:Nb
    P(i).X = zeros(Nt,3);
    P(i).q = nu*L0;
    P(i).a = a;
end

for i = 1:Nb
    P(i).X(1,:) = Xinitial(i,:);
end
    
F = zeros(Nb,3);

%% Brownian Dynamics Code

for n = 1:(Nt-1)
    
    % Get the position of every particle in one matrix
    % Because we compute it often, we will also precompute the distance
    % bewteen each particle.
    Xn = zeros(Nb,3);
    rn = zeros(Nb,Nb,3);
    rmag = zeros(Nb,Nb);
    
    for i = 1:Nb
        Xn(i,:) = P(i).X(n,:);
        for j = 1:Nb
            rn(i,j,:) = P(i).X(n,:) - P(j).X(n,:);
            rmag(i,j) = norm(P(i).X(n,:) - P(j).X(n,:));
        end
    end    
    
    %% Get DNA Measurements
    [Re2e(n)] = MyEnd2EndDist(Xn); % End to end Distance
    [L(n,:)] = MyBondLengths(Xn); % Bond Lengths
    % Compute the translational diffusion coeffecient based on the initial
    % and current position of the DNA
    [Dbulk(n)] = MyTransDiffusionCoeff(Xn,Xinitial,t(n)); 
    
    % Calculate the Translational Diffusion tensor
    D = TransDiffusionTensor(P, eta, T, n);

    %% Calculate the force on each bead
    
    % Stretching Force
    [Fstretching] = StretchingForce(P, h, L0, n);
    
    % Bending Force
    [Fbending] = BendingForce(P, g, n);
    
    % Electrostatic Force
    [Felectro] = ElectrostaticForce(P, eps, rD, n);
    
    % Eletrodynamic Force (opt)
    if ComputeDynamicForce == true
        [Fedynamic] = ElectrodynamicForce(P, t(1:n), epsr, sigma, Nvc, Ntau, n);
    else
        Fedynamic = zeros(Nb,3);
    end
    
    % External Field Force
    Fext = zeros(Nb,3);
    for i = 1:Nb
        Fext(i,:) = P(i).q .* Eext(n,:);
    end
    
    % Sum forces
%     Fbending
%     Fstretching
%     Felectro
%     Fedynamic
%     Fext
    Ftotal = Fbending + Fstretching + Felectro + Fedynamic + Fext;
    
    %% Propagate Motion
    for i = 1:Nb
        DF = zeros(1,3);
        for j = 1:Nb
            DF = DF + Ftotal(j,:) * squeeze(D(i,j,:,:));
        end
        D0 = kB*T / (6*pi*eta*P(i).a);
        
        R = sqrt(D0*dt)*randn(1,3);
        
        P(i).X(n+1,:) = P(i).X(n,:) + (dt / (kB*T) * DF) + R;
    end
    
end

for i = 1:Nb
    Xfinal(i,:) = P(i).X(Nt,:);
end

Re2e(Nt) = MyEnd2EndDist(Xfinal);
L(Nt,:) = MyBondLengths(Xfinal);
[Dbulk(Nt)] = MyTransDiffusionCoeff(Xfinal,Xinitial,t(Nt));

end