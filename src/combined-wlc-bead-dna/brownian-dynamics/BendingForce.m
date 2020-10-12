function [Fbending] = BendingForce(P, g, n)
% Computes the bending force between beads described in
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

% g is a bending force constant [J]
% n is the time step integer to use for P

% OUPUTS
% Fbending is a matrix size Nb x 3 giving the 3d bending force every
% bead is experiencing [N]

% Initialize variables
Nb = length(P);
Fbending =  zeros(Nb,3);

% --------Calculate bending force--------
% Calculate the bending force on each bead. This expression for
% bending force is taken from "Multistep Brownian Dynamics: 
% Application to Short Wormlike Chains" by S. A. ALLISON and J. A.
% McCAMMON
for i = 1:Nb
    if i == 1
        Fbending(i,:) = -g*( P(i).X(n,:) - 2*P(i+1).X(n,:) + P(i+2).X(n,:));
    elseif i == 2
        Fbending(i,:) = -g*( -2*P(i-1).X(n,:) + 5*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
    elseif (i >= 3) && (i <= Nb - 2)
        Fbending(i,:) = -g*( P(i-2).X(n,:) - 4*P(i-1).X(n,:) + 6*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
    elseif i == Nb - 1
        Fbending(i,:) = -g*( -2*P(i+1).X(n,:) + 5*P(i).X(n,:) - 4*P(i-1).X(n,:) + P(i-2).X(n,:));
    elseif i == Nb
        Fbending(i,:) = -g*( P(i).X(n,:) - 2*P(i-1).X(n,:) + P(i-2).X(n,:));
    end    
end

end

