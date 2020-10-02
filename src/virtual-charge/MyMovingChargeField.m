function [ E, B, phi ] = MyMovingChargeField( r, R, V, q, epsr, mur )

% Calculates the E and B field at point r from a point charge with q at
% point R moving with velocity V

% r is a 3 element vector describing the x,y,z coordinates of the point
% where the fields will be calculated [m]
% R is a 3 element vector describing the x,y,z coordinates of the position
% of the charge [m]
% V ia a 3 element vector describing the x,y,z directions of the velocity
% vector the dipole [m/s]
% q is a scalar representing the charge of the point [C]

% E = [Ex, Ey, Ez] [V/m]
% B = [Bx, By, Bz] [T]

% Define constants
eps0 = 8.85418e-12;
mu0 = 4*pi*10^(-7);

% Material Properties
eps = epsr*eps0;
mu = mu0*mur;

% R vectors
dR = norm(r - R);
dRhat = (r-R)./dR;

% E and B fields
E = q/(4*pi*eps) * dRhat/dR^2;
B = q*mu/(4*pi) * cross(V,dRhat)/dR^2;
phi = q/(4*pi*eps) /dR;

end

