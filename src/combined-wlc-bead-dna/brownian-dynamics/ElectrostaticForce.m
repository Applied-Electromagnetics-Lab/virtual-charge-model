function [Felectro] = ElectrostaticForce(P, eps, rD, n)
% Computes the electrostatic force between beads described in
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

% eps is the permitivity of the medium [F/m]
% rD is the debye length of the medium [m]
% n is the time step integer to use for P

% OUPUTS
% Felectro is a matrix size Nb x 3 giving the 3d electrostatic force every
% bead is experiencing [N]

% Initialize variables
Nb = length(P);
Felectro = zeros(Nb,3);

% --------Calculate Electrostatic force--------
% Calculate the screened electrostatic force between charges. 
for i = 1:Nb
    F = 0;
    for j = 1:Nb
        if i ~= j
            rij = P(i).X(n,:) - P(j).X(n,:);
            rijmag = norm(rij);
            rhat = rij./rijmag;
            F = F + ( P(i).q *P(j).q / (4*pi*eps*rijmag^2) * rhat * exp(-rijmag / rD) );
        end
    end
   Felectro(i,:) = F;
end

        
end

