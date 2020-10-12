function [D] = MyTransDiffusionCoeff(Pt,P0,t)

[~,Np] = size(P0);
Nt = length(t);
L = 1;

G = 0;
for j = 1:Np
    for k = 1:Np
        G = G + norm( Pt(j,:) - P0(k,:) )^2;
    end
end
G = G / 1/(Np^2*L^2);
Xcm0 = sum(P0,2)/Np;
Xcmt = sum(Pt,2)/Np;
D = norm(Xcmt - Xcm0)^2 / (6*t);

end

