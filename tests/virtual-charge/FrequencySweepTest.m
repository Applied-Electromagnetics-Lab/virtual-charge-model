clc; close all; clear;

f = 100e6;
T = 1/f;
Nt = 500;
t = linspace(0,5*T,Nt);

w = 2*pi*f;
A = 1e-9;
P = sin(w.*t).'*[0,0,A];
Nvc = 40;
Ntau = 5;
epsr = 80;
sigma = 1;

eps0 = 8.854e-12; % permitivity of free space [F/m]
eps = eps0*epsr; % permitivity of medium [F/m]

R1 = [0,0,1000e-9];
R2 = [0,1000e-9,0];
Rset = [R1;R2];

q = 1;

[Nr,~] = size(Rset);

phi = zeros(Nt,Nr);
E = zeros(Nt,Nr,3);
Nc = 1;
Pc = zeros(Nt,3,Nc);

for i = 1:Nt
    
    Pc = P(1:i,:,:);
    tc = t(1:i);
    [E(i,:,:), ~, phi(i,:)] = MyVirtualChargeEnsemble(Pc, q, tc, Rset, epsr, sigma, Nvc, Ntau);
    
end

epsbar = eps - 1i*sigma/w;
phimax = abs(q*A/(4*pi*epsbar*norm(R1)^2));

figure
plot(t,phi(:,1))
hold on
plot(t,phimax.*ones(1,length(phi)),'k--')

% figure
% plot(t(Nt+1:2*Nt),Pc(1,:,3))

