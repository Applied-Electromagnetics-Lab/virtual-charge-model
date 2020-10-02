function [P, Re2e, L] = My_DNA_BD(Params)
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
% Re2e is the end to end distance of the DNA at the end of the simulation
% [m]
% L is a vector of length Nb giving the length of every bond at the end of
% the simulation [m]

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

for i = 1:Nb
    P(i).X = zeros(Nt,3);
    P(i).X(1,3) = (i-1)*L0;
    P(i).q = nu*L0;
    P(i).a = a;
end

for i = 1:Nb
    Xinitial(i,:) = P(i).X(1,:);
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
    
    %% Bulk Translation Diffusion Coeff
    % Compute the translational diffusion coeffecient based on the initial
    % and current position of the DNA
%     [Dbulk(n)] = MyTransDiffusionCoeff(Xn,Xinitial,t(n));
    
%% Calculate the Translational Diffusion tensor

    D = zeros(Nb,Nb,3,3);
    for i = 1:Nb
        for j = 1:Nb
            if i == j
                D(i,j,:,:) = kB*T / (6*pi*eta*P(i).a) * eye(3);
            else
                rij = reshape(rn(i,j,:),1,3);
                rijmag = rmag(i,j);
                D(i,j,:,:) = kB*T / (8*pi*eta*rijmag) * ...
                    ( (eye(3) + rij.'*rij/rijmag^2) + ...
                    2*P(i).a*P(j).a / rijmag^2 * (1/3 .* eye(3) - ...
                    rij.'*rij/rijmag^2) );
            end
        end
    end

%% Calculate the force on each bead
    for i = 1:Nb
        
        % --------Stretching Force--------
        % Calculate the stretching force based on every bond length
        % The first and last bead only experience a force from one bond.
        % Every other beand experiences a force from both bonds.
        if i == 1
            L1 = P(i).X(n,:) - P(i+1).X(n,:); % Bond lengths
            L1hat = L1./norm(L1); % Normalize vector
            Fstretching = -h/2 * ((norm(L1) - L0).*L1hat); % Compute force            
        elseif i == Nb
            L2 = P(i).X(n,:) - P(i-1).X(n,:);
            L2hat = L2./norm(L2); % Normalize vector
            Fstretching = -h/2 * ((norm(L2) - L0).*L2hat); % Compute force
        else
            L1 = P(i).X(n,:) - P(i+1).X(n,:); % Bond lengths
            L2 = P(i).X(n,:) - P(i-1).X(n,:);
            L1hat = L1./norm(L1); % Normalize vector
            L2hat = L2./norm(L2); % Normalize vector
            Fstretching = -h/2 * ((norm(L1) - L0).*L1hat + (norm(L2) - L0).*L2hat); % Compute force
        end

        % --------Calculate bending force--------
        % Calculate the bending force on each bead. This expression for
        % bending force is taken from "Multistep Brownian Dynamics: 
        % Application to Short Wormlike Chains" by S. A. ALLISON and J. A.
        % McCAMMON
        if i == 1
            Fbending = -g*( P(i).X(n,:) - 2*P(i+1).X(n,:) + P(i+2).X(n,:));
        elseif i == 2
            Fbending = -g*( -2*P(i-1).X(n,:) + 5*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
        elseif (i >= 3) && (i <= Nb - 2)
            Fbending = -g*( P(i-2).X(n,:) - 4*P(i-1).X(n,:) + 6*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
        elseif i == Nb - 1
            Fbending = -g*( -2*P(i+1).X(n,:) + 5*P(i).X(n,:) - 4*P(i-1).X(n,:) + P(i-2).X(n,:));
        elseif i == Nb
            Fbending = -g*( P(i).X(n,:) - 2*P(i-1).X(n,:) + P(i-2).X(n,:));
        end

        % --------Calculate Electrostatic force--------
        % Calculate the screened electrostatic force between charges. 
        Felectro = 0;
        for j = 1:Nb
            if i ~= j
                rij = reshape(rn(i,j,:),1,3);
                rijmag = rmag(i,j);
                rhat = rij./rijmag;
                Felectro = Felectro + ( P(i).q *P(j).q / (4*pi*eps*rijmag^2) * rhat * exp(-rijmag / rD) );
            end
        end
        
        % ---------Calculate Electrodynamic Force----------
        % Only bother doing this expensive computation if there is an
        % external field, otherwise skip it
        if NoDynamics == false
            
            for j = 1:Nb
                Pnow(j).X = P(j).X(1:n,:);
                Pnow(j).q = P(j).q;
            end
            Pnow = [Pnow(1:i-1),Pnow(i+1:end)];

            [Evc, ~, ~] = MyVirtualChargeEnsemble(Pnow, t(1:n), P(i).X(n,:), epsr, sigma, Nvc, Ntau);
            Fedynamic = Evc*P(i).q;
                    
            %--------External Elecric Field--------
            Fext = P(i).q .* Eext(n,:);
            
        else
            Fedynamic = [0,0,0];
            Fext = [0,0,0];
        end
        
        % Sum forces
        F(i,:) = Fbending + Fstretching + Felectro + Fedynamic + Fext;
        
    end
    
    %% Propagate Motion
    for i = 1:Nb
        DF = zeros(1,3);
        for j = 1:Nb
            DF = DF + F(j,:) * squeeze(D(i,j,:,:));
        end
        D0 = kB*T / (6*pi*eta*P(i).a);
        
        R = sqrt(D0*dt)*randn(1,3);
        
        P(i).X(n+1,:) = P(i).X(n,:) + (dt / (kB*T) * DF) + R;
    end
    
end

for i = 1:Nb
    Xfinal(i,:) = P(i).X(Nt,:);
end

[Re2e] = MyEnd2EndDist(Xfinal);
[L] = MyBondLengths(Xfinal);

end