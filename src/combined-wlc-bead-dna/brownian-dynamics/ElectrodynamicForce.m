function [Fedynamic] = ElectrodynamicForce(P, t, epsr, sigma, Nvc, Ntau, n)
% Computes the electrodynamic force between beads for the DNA model in
% "A Combined Wormlike-Chain and Bead Model for Dynamic Simulations of
% Long Linear DNA" by Hongmei Jian and Alexander V. Vologodskii and Tamar
% Schlick (1997)

% INPUTS
% P is the structure containing information on every particle in the
% simulation
% P(i).X is a Nt x 3 sized matrix giving the 3D location of particle i at
% every point in time [m]
% P(i).a is the radius of particle i [m]
% P(i).q is the charge of particle i [C]

% t is the time vector of the system [s]
% eps is the permitivity of the medium [F/m]
% sigma is the conductivity of the medium [S/m]
% Nvc is the number of virtual charges to approximate the wake
% Ntau is the number of time constants to approximate the wake
% n is the time step integer to use for P

% OUPUTS
% Fedynamic is a matrix size Nb x 3 giving the 3d dynamic force every
% bead is experiencing [N]

% Initialize variables
Nb = length(P);
Fedynamic = zeros(Nb,3);

% ---------Calculate Electrodynamic Force----------

% First we need to sort out the history of every particle
for j = 1:Nb
    Pnow(j).X = P(j).X(1:n,:);
    Pnow(j).q = P(j).q;
end
    
for i = 1:Nb
    % Exclude the particle we are computing the force on   
    Pnow2 = [Pnow(1:i-1),Pnow(i+1:end)];

    [Evc, ~, ~] = MyVirtualChargeEnsemble(Pnow2, t(1:n), P(i).X(n,:), epsr, sigma, Nvc, Ntau);
    Fedynamic(i,:) = Evc*P(i).q;
end

            
        
end

