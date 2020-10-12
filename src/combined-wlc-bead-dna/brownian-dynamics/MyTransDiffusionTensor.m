function [D] = MyTransDiffusionTensor(P, eta)

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

