function [R] = MyEnd2EndDist(P)
% P is a matrix size N x 3
r = P(end,:) - P(1,:);
R = norm(r);
end

