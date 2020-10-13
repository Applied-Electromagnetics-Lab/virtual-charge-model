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
        r1 = P(i+1).X(n,:) - P(i).X(n,:);
        r2 = P(i+2).X(n,:) - P(i+1).X(n,:);
        u1 = r1./norm(r1);
        u2 = r2./norm(r2);
        
        Fbending(i,:) = -g*(u2 - u1)./norm(r1);
    elseif i == 2
        r1 = P(i).X(n,:) - P(i-1).X(n,:);
        r2 = P(i+1).X(n,:) - P(i).X(n,:);
        r3 = P(i+2).X(n,:) - P(i+1).X(n,:);
        u1 = r1./norm(r1);
        u2 = r2./norm(r2);
        u3 = r3./norm(r3);
        
        Fbending(i,:) = -g*(u2-u1) * (1/norm(r1) + 1/norm(r2)) - g*(u3-u2)/norm(r2);
    elseif (i >= 3) && (i <= Nb - 2)
        rim2 = P(i-1).X(n,:) - P(i-2).X(n,:);
        rim1 = P(i).X(n,:) - P(i-1).X(n,:);
        ri = P(i+1).X(n,:) - P(i).X(n,:);
        rip1 = P(i+2).X(n,:) - P(i+1).X(n,:);
        
        uim2 = rim2./norm(rim2);
        uim1 = rim1./norm(rim1);
        ui = ri./norm(ri);
        uip1 = rip1./norm(rip1);
        
        Fbending(i,:) = -g*(uim1 - uim2)./norm(rim1) + g*(ui - uim1)*(1/norm(rim1) + 1/norm(ri)) - g*(uip1 - ui)./norm(ri);
    elseif i == Nb - 1
        rim2 = P(i-1).X(n,:) - P(i-2).X(n,:);
        rim1 = P(i).X(n,:) - P(i-1).X(n,:);
        ri = P(i+1).X(n,:) - P(i).X(n,:);
        
        uim2 = rim2./norm(rim2);
        uim1 = rim1./norm(rim1);
        ui = ri./norm(ri);
        
        Fbending(i,:) = -g*(uim1 - uim2)./norm(rim1) + g*(ui - uim1)*(1/norm(rim1) + 1/norm(ri));
    elseif i == Nb
        rim2 = P(i-1).X(n,:) - P(i-2).X(n,:);
        rim1 = P(i).X(n,:) - P(i-1).X(n,:);
        
        uim2 = rim2./norm(rim2);
        uim1 = rim1./norm(rim1);
        
        Fbending(i,:) = -g*(uim1 - uim2)./norm(rim1);
    end    
end

% for i = 1:Nb
%     if i == 1
%         Fbending(i,:) = -g*( P(i).X(n,:) - 2*P(i+1).X(n,:) + P(i+2).X(n,:));
%     elseif i == 2
%         Fbending(i,:) = -g*( -2*P(i-1).X(n,:) + 5*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
%     elseif (i >= 3) && (i <= Nb - 2)
%         Fbending(i,:) = -g*( P(i-2).X(n,:) - 4*P(i-1).X(n,:) + 6*P(i).X(n,:) - 4*P(i+1).X(n,:) + P(i+2).X(n,:));
%     elseif i == Nb - 1
%         Fbending(i,:) = -g*( -2*P(i+1).X(n,:) + 5*P(i).X(n,:) - 4*P(i-1).X(n,:) + P(i-2).X(n,:));
%     elseif i == Nb
%         Fbending(i,:) = -g*( P(i).X(n,:) - 2*P(i-1).X(n,:) + P(i-2).X(n,:));
%     end    
% end

end

