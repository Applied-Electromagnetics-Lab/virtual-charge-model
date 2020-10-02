function [E, phi] = MySingleChargeTimeSweep(X, t, q, R, Nvc, Ntau, epsr, sigma)

P(1).X = X;
P(1).q = q;

Nt = length(t);
phi = zeros(Nt,1);
E = zeros(Nt,3);

for i = 1:Nt

    P2(1).X = P(1).X(1:i,:);
    P2(1).q = P(1).q;
    tc = t(1:i);
    [E(i,:), ~, phi(i)] = MyVirtualChargeEnsemble(P2, tc, R, epsr, sigma, Nvc, Ntau);
    
end

end

