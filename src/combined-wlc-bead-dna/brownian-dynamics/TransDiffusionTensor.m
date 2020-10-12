function [D] = TransDiffusionTensor(P, eta, T, n)

% Computes the translational diffusion tensor based on the method described
% in "A Combined Wormlike-Chain and Bead Model for Dynamic Simulations of
% Long Linear DNA" by Hongmei Jian and Alexander V. Vologodskii and Tamar
% Schlick (1997)

% INPUTS
% P is the structure containing information on every particle in the
% simulation
% P(i).X is a Nt x 3 sized matrix giving the 3D location of particle i at
% every point in time [m]
% P(i).a is the radius of particle i [m]
% P(i).q is the charge of particle i [C]

% eta is the viscosity of the medium [Pa*s]

% OUPUTS
% D is the Translational Diffusion Tensor size Nb x Nb x 3 x 3 [m^2 / s]

kB = 1.38064852e-23; % Boltzman Constant [J/K]

Nb = length(P);

D = zeros(Nb,Nb,3,3);
for i = 1:Nb
    for j = 1:Nb
        if i == j
            D(i,j,:,:) = kB*T / (6*pi*eta*P(i).a) * eye(3);
        else
            rij = P(i).X(n,:) - P(j).X(n,:);
            rijmag = norm(rij);
            D(i,j,:,:) = kB*T / (8*pi*eta*rijmag) * ...
                ( (eye(3) + rij.'*rij/rijmag^2) + ...
                2*P(i).a*P(j).a / rijmag^2 * (1/3 .* eye(3) - ...
                rij.'*rij/rijmag^2) );
        end
    end
end
    
end

