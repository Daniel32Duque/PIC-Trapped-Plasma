%{
Written by: Daniel Duque
Last modified on 15 Mar 2020

This script plots the initial density of the plasmas
obtained from the macro-particle's initial distribution
And compares this with the expected equilibrium distribution.

This tells you how many electorns are in a central macro-particle
Recall that this number increases as 8 * indexR for macro-particles at
index different than 0

Also outputs the average percentaje error in the density from the
macroparticles compared to the expected density for the equilibrium

Run this, and increase/decrease the number of macro-particles until you are
happy with these results.
When to be happy is a total different question. Depends on what you want to
simulate.
%}

%-----------------------------
%Input
%-----------------------------

pathRead = [pwd, '\Data Files\'];
pathWrite = [pwd, '\Plots\'];

trapParameters = csvread([pathRead, 'Z Trap Parameters.csv']);
electronParameters = csvread([pathRead, 'Z Electron Parameters.csv']);
electronPositions = csvread([pathRead, 'Z PositionsElectrons.csv']);
densityE = csvread([pathRead, 'Z Expected Electron Density.csv']);

%-----------------------------
%End Input
%-----------------------------

trapRadius = trapParameters(1,1);
trapLength = trapParameters(7,1);
Nz = trapParameters(5,1);
Nr = trapParameters(5,2);

dr = trapRadius/Nr;
dz = trapLength/Nz;
[X,Y] = meshgrid(0:dr:trapRadius-dr, 0:dz:trapLength);

chargeElectron = -1.602176634e-19;
chargeMacroElectron = electronParameters(3);
macroDensityElectron = 4 * chargeMacroElectron / (pi * dz * dr * dr);
electronsInCentralMacroParticle = chargeMacroElectron / chargeElectron

macroElectrons = zeros(Nz + 1, Nr);
[rows, columns] = size(electronPositions);
for i = 1:rows
    indexR = electronPositions(i,1) + 1;
    indexZ = floor(electronPositions(i,2) / dz);
    z = electronPositions(i,2) - indexZ * dz;
    weightFactor = z / dz;
    indexZ = indexZ +1;
    macroElectrons(indexZ, indexR) = macroElectrons(indexZ, indexR) + (1 - weightFactor) * macroDensityElectron / chargeElectron;
    macroElectrons(indexZ + 1, indexR) = macroElectrons(indexZ + 1, indexR) + weightFactor * macroDensityElectron / chargeElectron;
end

figure
h=surf(X,Y,macroElectrons);
set(h,'LineStyle','none')
title('Macro-particle''s Electron Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('part / m^3');
savefig([pathWrite, 'Plot Macro-Electron Density.fig'])

densityE = reshape(densityE,Nz + 1, Nr) / chargeElectron;
percentElectrons = abs(densityE - macroElectrons);
percentElectrons = percentElectrons ./ densityE;
percentElectrons(isnan(percentElectrons)) = 0;
percentElectrons(percentElectrons == 1) = 0;
percentElectrons = percentElectrons * 100;

meanPercentErrorNonZeroRegionElectron = mean(nonzeros(percentElectrons(:,:)))

figure
h=surf(X,Y,percentElectrons);
set(h,'LineStyle','none')
title('Percentaje Error Electron Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Percentaje Error Electron Density.fig'])

percentRElectron = zeros(Nr,0);
for i = 1:Nr
    count = 0;
    total = 0;
    for j = 1:Nz + 1
        if percentElectrons(j,i) ~= 0
            count = count + 1;
            total = total + percentElectrons(j,i);
        end
    end
    percentRElectron(i) = total / count;
end
percentRElectron(isnan(percentRElectron)) = 0;

figure
h=plot(0:dr:trapRadius-dr,percentRElectron);
title('Average Percentaje Error Electron Density');
xlabel('r (m)');
ylabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Average Percentaje Error Electron Density.fig'])

macroProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        macroProfile(i) = macroProfile(i) + macroElectrons(j,i);
    end
end

figure
h=plot(0:dr:trapRadius-dr,macroProfile);
title('Macro-particle''s Electron Profile');
xlabel('r (m)');
ylabel('Line Density (part / m^2)');
savefig([pathWrite, 'Plot Macro Electron Profile.fig'])

electronProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        electronProfile(i) = electronProfile(i) + densityE(j,i);
    end
end

percentProfile = abs(electronProfile - macroProfile);
percentProfile = percentProfile ./ electronProfile;
percentProfile(isnan(percentProfile)) = 0;
percentProfile(percentProfile == 1) = 0;
percentProfile = percentProfile * 100;
figure
h=plot(0:dr:trapRadius-dr,percentProfile);
title('Percent Error Electron Profile');
xlabel('r (m)');
ylabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Percent Error Electron Profile.fig'])