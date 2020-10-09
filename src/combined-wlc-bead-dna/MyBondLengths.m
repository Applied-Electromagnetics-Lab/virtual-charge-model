function [L] = MyBondLengths(P)
% P is a matrix size N x 3
[Nb,~] = size(P);
for i = 1:(Nb-1)
    r = P(i,:) - P(i+1,:);
    L(i) = norm(r);
end
% r = diff(P);
% L = vecnorm(r,2,2);
end

