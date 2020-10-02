clc; close all; clear;

f = 100e6;
T = 1/f;
Nt = 500;
t = linspace(0,10*T,Nt);
t0 = 5*T;
[~,t0I] = min(abs(t-t0));

w = 2*pi*f;
A = 1e-9;
Xp = sin(w.*t).'*[0,0,A];

q  =1;

P(1).X = Xp;
P(1).q = q;

Nvc = 40;
Ntau = 5;
epsr = 80;
sigma = 1;

eps0 = 8.854e-12; % permitivity of free space [F/m]
eps = eps0*epsr; % permitivity of medium [F/m]

R1 = [0,0,1000e-9];
R2 = [0,1000e-9,0];
Rset = [R1;R2];

[Nr,~] = size(Rset);

phi = zeros(Nt,Nr);
E = zeros(Nt,Nr,3);
Nc = 1;
Pc = zeros(Nt,3,Nc);

for i = 1:Nt

    P2(1).X = P(1).X(1:i,:);
    P2(1).q = P(1).q;
    tc = t(1:i);
    [E(i,:,:), ~, phi(i,:)] = MyVirtualChargeEnsemble(P2, tc, Rset, epsr, sigma, Nvc, Ntau);
    
end

epsbar = eps - 1i*sigma/w;
phimax = abs(q*A/(4*pi*epsbar*norm(R1)^2));

Emax1 = abs(2*q*A/(4*pi*epsbar*norm(R1)^3));
Emax2 = abs(  q*A/(4*pi*epsbar*norm(R2)^3));

figure
plot(t,phi(:,1))
hold on
plot(t,phimax.*ones(1,length(phi)),'k--')

figure
plot(t,E(:,1,3))
hold on
plot(t,Emax1.*ones(1,length(phi)),'k--')

figure
plot(t,E(:,2,3))
hold on
plot(t,Emax2.*ones(1,length(phi)),'k--')

% figure
% plot(t(Nt+1:2*Nt),Pc(1,:,3))

