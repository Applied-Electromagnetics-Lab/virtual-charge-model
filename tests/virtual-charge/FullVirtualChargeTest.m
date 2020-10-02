clc; close all; clear;

%% Material properties and constants

eps0 = 8.854e-12; % permitivity of free space [F/m]
qe = -1.60217662e-19; % Charge of electron [C]
kB = 1.38064852e-23; % Boltzman Constant [J/K]

epsr = 80;
eps = eps0*epsr; % permitivity of medium [F/m]
sigma = 1; % Conductivity of medium [S/m]
wp = sigma / eps;

%% Frequency Test
% Test the amplitude of the E field generated at a fixed distance around
% the dipole at both the Z and X axis as a function of frequency

Nfreq = 50;
freq = logspace(3,11,Nfreq);
w = 2*pi*freq;
epsbar = eps - 1i*sigma./w;
A = 1e-9;
q = 1*qe;
Nt = 100;

Nvc = 50;
Ntau = 5;
epsr = 80;

Dist  = 1e-6;
R1 = [0,0,Dist];
R2 = [Dist,0,0];

for i = 1:Nfreq
    T = 1/freq(i);
    t = linspace(0,5*T,Nt);
    X = sin(w(i).*t).' * [0,0,A];
    [E1, phi1] = MySingleChargeTimeSweep(X, t, q, R1, Nvc, Ntau, epsr, sigma);
    [E2, phi2] = MySingleChargeTimeSweep(X, t, q, R2, Nvc, Ntau, epsr, sigma);
    
    Emag1 = vecnorm(E1,2,2);
    Emax1(i) = max(Emag1);
    Edipole1(i) = abs( 2*q*A/(4*pi*Dist^3*epsbar(i)) );
    
    Emag2 = vecnorm(E2,2,2);
    Emax2(i) = max(Emag2);
    Edipole2(i) = abs( q*A/(4*pi*Dist^3*epsbar(i)) );
    
end

figure
loglog(freq,Emax1,'LineWidth',2);
hold on
loglog(freq,Edipole1,'k--','LineWidth',2);
xlabel('Frequency [Hz]')
ylabel('Emag [V/m]')
title('Off Z axis')
legend('Virtual Charge Fields','Dipole Approximation Fields')
xlim([min(freq),max(freq)])

figure
loglog(freq,Emax2,'LineWidth',2);
hold on
loglog(freq,Edipole2,'k--','LineWidth',2);
xlabel('Frequency [Hz]')
ylabel('Emag [V/m]')
title('Off X axis')
legend('Virtual Charge Fields','Dipole Approximation Fields')
xlim([min(freq),max(freq)])

figure
Error = abs( (Emax1 - Edipole1) ./ Emax1 ) * 100;
semilogx(freq,Error);
xlabel('Frequency [Hz]')
ylabel('Percent Error')
title('Off Z axis')
xlim([min(freq),max(freq)])

%% Spatial Sweep
% Test the amplitude of the electric field at a fixed frequency as a
% function of distancc for both Z and X directions

Nr = 50;
R = logspace(-9,-6,Nr);

freq = 1e9;
w = 2*pi*freq;
epsbar = eps - 1i*sigma./w;

A = 1e-9;
q = 1*qe;
Nt = 100;

Nvc = 40;
Ntau = 5;
epsr = 80;

T = 1/freq;
t = linspace(0,5*T,Nt);
X = sin(w.*t).' * [0,0,A];

for i = 1:Nr
    
    R1 = [0,0,R(i)];
    R2 = [R(i),0,0];
   
    [E1, phi1] = MySingleChargeTimeSweep(X, t, q, R1, Nvc, Ntau, epsr, sigma);
    [E2, phi2] = MySingleChargeTimeSweep(X, t, q, R2, Nvc, Ntau, epsr, sigma);
    
    Emag1 = vecnorm(E1,2,2);
    Emax1(i) = max(Emag1);
    Edipole1(i) = abs( 2*q*A/(4*pi*R(i)^3*epsbar) );
    
    Emag2 = vecnorm(E2,2,2);
    Emax2(i) = max(Emag2);
    Edipole2(i) = abs( q*A/(4*pi*R(i)^3*epsbar) );
    
end

figure
loglog(R,Emax1,'LineWidth',2);
hold on
loglog(R,Edipole1,'k--','LineWidth',2);
xlabel('Distance [m]')
ylabel('Emag [V/m]')
title('Off Z axis')
legend('Virtual Charge Fields','Dipole Approximation Fields')
xlim([min(R),max(R)])

figure
loglog(R,Emax2,'LineWidth',2);
hold on
loglog(R,Edipole2,'k--','LineWidth',2);
xlabel('Distance [m]')
ylabel('Emag [V/m]')
title('Off X axis')
legend('Virtual Charge Fields','Dipole Approximation Fields')
xlim([min(R),max(R)])

%% Angle Sweep
% Test the amplitude of the electric field at a fixed distance and
% an angle 360 degrees around the dipole

Ntheta = 50;
theta = linspace(0,2*pi,Ntheta);

freq = 1e9;
w = 2*pi*freq;
epsbar = eps - 1i*sigma./w;

A = 1e-9;
q = 1*qe;
Nt = 100;

Nvc = 40;
Ntau = 5;
epsr = 80;

Dist  = 1e-6;

T = 1/freq;
t = linspace(0,5*T,Nt);
X = sin(w.*t).' * [0,0,A];

for i = 1:Ntheta
    
    R = [sin(theta(i)),0,cos(theta(i))].*Dist;
   
    [E1, phi1] = MySingleChargeTimeSweep(X, t, q, R, Nvc, Ntau, epsr, sigma);
    
    Emag1 = vecnorm(E1,2,2);
    Emax1(i) = max(Emag1);
    
    EdipoleR = abs( 2*q*A/(4*pi*Dist^3*epsbar) * cos(theta(i)));
    EdipoleT = abs(  q*A/(4*pi*Dist^3*epsbar) * sin(theta(i)));
    Edipole(i) = sqrt(EdipoleR^2 + EdipoleT^2);
    
    Error(i) = abs( (Emax1(i) - Edipole(i)) / Emax1(i) ) * 100;
    
end

figure
polarplot(theta,Emax1,'LineWidth',2);
hold on
polarplot(theta,Edipole,'k--','LineWidth',2);

title('Theta Sweep')
legend('Virtual Charge Fields','Dipole Approximation Fields')

figure
plot(theta/pi*180,Error);
xlabel('Angle [degrees]')
ylabel('Percent Error')

%% Nvc and Ntau sweep
% Test the amplitude of the electric field at a fixed distance and
% frequency and orientation 

Nvc = [3:1:50];
Nn = length(Nvc);

freq = 1000e6;
w = 2*pi*freq;
epsbar = eps - 1i*sigma./w;

A = 1e-9;
q = 1*qe;
Nt = 100;

epsr = 80;

Dist  = 1e-6;
R = [0,0,Dist];

T = 1/freq;
t = linspace(0,5*T,Nt);
X = sin(w.*t).' * [0,0,A];

for i = 1:Nn
   
    [E0, phi0] = MySingleChargeTimeSweep(X, t, q, R, Nvc(i), 2, epsr, sigma);
    [E1, phi1] = MySingleChargeTimeSweep(X, t, q, R, Nvc(i), 3, epsr, sigma);
    [E2, phi2] = MySingleChargeTimeSweep(X, t, q, R, Nvc(i), 4, epsr, sigma);
    [E3, phi3] = MySingleChargeTimeSweep(X, t, q, R, Nvc(i), 5, epsr, sigma);
    [E4, phi4] = MySingleChargeTimeSweep(X, t, q, R, Nvc(i), 6, epsr, sigma);
    
    Emag0 = vecnorm(E0,2,2);
    Emax0(i) = max(Emag0);
    
    Emag1 = vecnorm(E1,2,2);
    Emax1(i) = max(Emag1);
    
    Emag2 = vecnorm(E2,2,2);
    Emax2(i) = max(Emag2);
    
    Emag3 = vecnorm(E3,2,2);
    Emax3(i) = max(Emag3);
    
    Emag4 = vecnorm(E4,2,2);
    Emax4(i) = max(Emag4);
    
    Edipole = abs( 2*q*A/(4*pi*Dist^3*epsbar) );
    
    Error0(i) = abs( (Emax0(i) - Edipole) / Emax0(i) ) * 100;
    Error1(i) = abs( (Emax1(i) - Edipole) / Emax1(i) ) * 100;
    Error2(i) = abs( (Emax2(i) - Edipole) / Emax2(i) ) * 100;
    Error3(i) = abs( (Emax3(i) - Edipole) / Emax3(i) ) * 100;
    Error4(i) = abs( (Emax4(i) - Edipole) / Emax4(i) ) * 100;
    
end

figure
plot(Nvc,Error0);
hold on
plot(Nvc,Error1);
plot(Nvc,Error2);
plot(Nvc,Error3);
plot(Nvc,Error4);
xlabel('Number of Virtual Charges')
ylabel('Percent Error')
legend('N\tau = 2', 'N\tau = 3','N\tau = 4','N\tau = 5','N\tau = 6')