function [L] = MyBondLengths(P)
% P is a matrix size N x 3
r = diff(P);
L = vecnorm(r,2,2);
end

