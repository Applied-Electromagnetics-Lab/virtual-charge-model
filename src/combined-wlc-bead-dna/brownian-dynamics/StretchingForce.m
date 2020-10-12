function [Fstretching] = StretchingForce(P,h,L0,n)
% Computes the stretching force between beads described in
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

% h is the spring constant [J/m^2]
% L0 is the resting length of the bonds [m]
% n is the time step integer to use for P

% OUPUTS
% Fstretching is a matrix size Nb x 3 giving the 3d stretching force every
% bead is experiencing [N]

% Initialize variables
Nb = length(P);
Fstretching = zeros(Nb,3);

% --------Stretching Force--------
% Calculate the stretching force based on every bond length
% The first and last bead only experience a force from one bond.
% Every other beand experiences a force from both bonds.
for i = 1:Nb
    if i == 1
        L1 = P(i).X(n,:) - P(i+1).X(n,:); % Bond lengths
        L1hat = L1./norm(L1); % Normalize vector
        Fstretching(i,:) = -h/2 * ((norm(L1) - L0).*L1hat); % Compute force            
    elseif i == Nb
        L2 = P(i).X(n,:) - P(i-1).X(n,:);
        L2hat = L2./norm(L2); % Normalize vector
        Fstretching(i,:) = -h/2 * ((norm(L2) - L0).*L2hat); % Compute force
    else
        L1 = P(i).X(n,:) - P(i+1).X(n,:); % Bond lengths
        L2 = P(i).X(n,:) - P(i-1).X(n,:);
        L1hat = L1./norm(L1); % Normalize vector
        L2hat = L2./norm(L2); % Normalize vector
        Fstretching(i,:) = -h/2 * ((norm(L1) - L0).*L1hat + (norm(L2) - L0).*L2hat); % Compute force
    end
end
        
end

