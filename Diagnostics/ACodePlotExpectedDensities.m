%{
Written by: Daniel Duque
Last modified on 14 Mar 2020

This script plots the expected particle density
both for electrons and antiprotons
Malmberg Trap
%}

%-----------------------------
%Input
%-----------------------------

pathRead = [pwd, '\Data Files\'];
temperature = csvread([pathRead, 'Z Temperature.txt']);
densityE = csvread([pathRead, 'Z Expected Electron Density.csv']);
densityA = csvread([pathRead, 'Z Expected Antiproton Density.csv']);
parameters = csvread([pathRead, 'Z Trap Parameters.csv']); 

pathWrite = [pwd, '\Plots\'];

%-----------------------------
%End Input
%-----------------------------

chargeElectron = -1.602176634e-19;
KB = 1.380649e-23;
epsilon = 8.8541878128e-12;
massE = 9.1093837015e-31;

trapRadius = parameters(1,1);
trapLength = parameters(7,1);
Nz = parameters(5,1);
Nr = parameters(5,2);

dr = trapRadius/Nr;
dz = trapLength/Nz;
[X,Y] = meshgrid(0:dr:trapRadius-dr, 0:dz:trapLength);
   
densityE = reshape(densityE,Nz + 1, Nr) / chargeElectron;
densityA = reshape(densityA,Nz + 1, Nr) / chargeElectron;

n = max(densityE(:));
lambda = sqrt(epsilon * KB * temperature / (n * chargeElectron^2));
minGridsZ = trapLength / lambda
minGridsR = trapRadius / lambda

frequency = sqrt(n * chargeElectron^2 / (massE * epsilon)) / (2 *pi);
maxPeriod = 1 / frequency

figure
h=surf(X,Y,densityE);
set(h,'LineStyle','none')
title('Expected Electron Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('part / m^3');
savefig([pathWrite, 'Plot Expected Electron Density.fig'])

figure
h=surf(X,Y,densityA);
set(h,'LineStyle','none')
title('Expected Antiproton Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('part / m^3');
savefig([pathWrite, 'Plot Expected Antiproton Density.fig'])

electronProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        electronProfile(i) = electronProfile(i) + densityE(j,i);
    end
end

figure
h=plot(0:dr:trapRadius-dr,electronProfile);
title('Expected Electron Profile');
xlabel('r (m)');
ylabel('Line Density (part / m^2)');
savefig([pathWrite, 'Plot Expected Electron Profile.fig'])

antiprotonProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        antiprotonProfile(i) = antiprotonProfile(i) + densityA(j,i);
    end
end

figure
h=plot(0:dr:trapRadius-dr,antiprotonProfile);
title('Expected Antiproton Profile');
xlabel('r (m)');
ylabel('Line Density (part / m^2)');
savefig([pathWrite, 'Plot Expected Antiproton Profile.fig'])