function [E, B, phi] = MyVirtualChargeEnsemble(P, t, Rset, epsr, sigma, Nvc, Ntau)
% INPUTS
% P is a structure array of lenth Nc, where Nc is the number of charges
% It has the following fields:
% P(i).X(n,3) represents the 3D position of charge i at time n (size Nt x
% 3) [m]
% P(i).q represents the magnitude of charge i [C]
% t is a vector length Nt with the time stamp of each position labeled in
% P.X [s]
% Rset is a matrix size [Nr x 3] representing the 3D coordinates of every
% point where the fields and potentials will be computed.
% epsr is the relative permitivity of the medium
% sigma is the conductivity of the medium [S/m]
% Nvc is the number of virtual charges
% Ntau is the number of plasma time constants to consider before assuming
% the wake has settled to zero.

% OUTPUTS
% phi is a vector length Nr that gives the potential at every coordinate
% specified by Rset
% E is a matrix size [Nr,3] that gives the 3D E vector at each coordinate
% specified by Rset

% Material Properties
eps0 = 8.85418e-12;
eps = epsr*eps0;
wp = sigma / eps; % Plasma frequency [rad / s]
taup = 1/wp; % Plasma time constant [s]

% Compute the virtual charge time constant. This is the amount of time
% between each virtual charge.
tauvc = Ntau * taup / Nvc; 
tNOW = t(end);

% Get the Number of charges (Nc), number of time points (Nt) and number of
% locations to calculate the potential and E field (Nr)
Nc = length(P);
[Nt,~] = size(P(1).X);
[Nr,~] = size(Rset);

E = zeros(Nr,3);
B = zeros(Nr,3);
phi = zeros(Nr,1);

for i = 1:Nc
    
    % Generate a vector of virtual charges length Nvc. The charges should
    % decay according to the plasma time constant. The decay per charge is
    % determined by the virtual charge time constant.
    I = (1:1:Nvc);
    qvc = exp(-I.*tauvc/taup);
    % The vector schould be scaled such that q + sum(qvc) = 0
    scale = -P(i).q / sum(qvc);
    qvc = qvc.*scale;
    
    for k = 1:Nr
        % Calculate potentials from all virtual charges
        for j = 1:Nvc
            % Calculate how long ago the virtual charge was created
            tvc = tNOW - j*tauvc;
            % If the virtual charge time is beyond the scope of the time
            % vector (i.e. it is from before t = 0) then we give it the
            % initial location. 
            if tvc < 0
                Xvcj = P(i).X(1,:);
            else
                Xvcj = interp1(t,P(i).X,tvc);
            end

            % Calculate the fields from the virtual charge j from charge i
            [E_temp, B_temp, phi_temp] = ...
               MyMovingChargeField( Rset(k,:), Xvcj, [0,0,0], qvc(j), epsr, 1 ); 
           E(k,:) = E(k,:) + E_temp;
           B(k,:) = B(k,:) + B_temp;
           phi(k) = phi(k) + phi_temp;
        end
        % Now that every virtual charge contribution from charge i is
        % calculated, calculate potential from primary charge
        [E_temp, B_temp, phi_temp] = ...
               MyMovingChargeField( Rset(k,:), P(i).X(end,:), [0,0,0], P(i).q, epsr, 1 ); 
           
        % Sum all the fields together
        E(k,:) = E(k,:) + E_temp;
        B(k,:) = B(k,:) + B_temp;
        phi(k) = phi(k) + phi_temp;
    end
end

end
