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
antiprotonParameters = csvread([pathRead, 'Z Antiproton Parameters.csv']);
antiprotonPositions = csvread([pathRead, 'Z PositionsAntiprotons.csv']);
densityA = csvread([pathRead, 'Z Expected Antiproton Density.csv']);

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
chargeMacroAntiproton = antiprotonParameters(3);
macroDensityAntiproton = 4 * chargeMacroAntiproton / (pi * dz * dr * dr);
antiprotonsInCentralMacroParticle = chargeMacroAntiproton / chargeElectron

macroAntiprotons = zeros(Nz + 1, Nr);
[rows, columns] = size(antiprotonPositions);
for i = 1:rows
    indexR = antiprotonPositions(i,1) + 1;
    indexZ = floor(antiprotonPositions(i,2) / dz);
    z = antiprotonPositions(i,2) - indexZ * dz;
    weightFactor = z / dz;
    indexZ = indexZ +1;
    macroAntiprotons(indexZ, indexR) = macroAntiprotons(indexZ, indexR) + (1 - weightFactor) * macroDensityAntiproton / chargeElectron;
    macroAntiprotons(indexZ + 1, indexR) = macroAntiprotons(indexZ + 1, indexR) + weightFactor * macroDensityAntiproton / chargeElectron;
end

figure
h=surf(X,Y,macroAntiprotons);
set(h,'LineStyle','none')
title('Macro-particle''s Antiproton Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('part / m^3');
savefig([pathWrite, 'Plot Macro-Antiproton Density.fig'])

densityA = reshape(densityA,Nz + 1, Nr) / chargeElectron;
percentAntiprotons = abs(densityA - macroAntiprotons);
percentAntiprotons = percentAntiprotons ./ densityA;
percentAntiprotons(isnan(percentAntiprotons)) = 0;
percentAntiprotons(percentAntiprotons == 1) = 0;
percentAntiprotons = percentAntiprotons * 100;

meanPercentErrorNonZeroRegionAntiproton = mean(nonzeros(percentAntiprotons(:,:)))

figure
h=surf(X,Y,percentAntiprotons);
set(h,'LineStyle','none')
title('Percentaje Error Antiproton Density');
xlabel('r (m)');
ylabel('z (m)');
zlabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Percentaje Error Antiproton Density.fig'])

percentRAntiproton = zeros(Nr,0);
for i = 1:Nr
    count = 0;
    total = 0;
    for j = 1:Nz + 1
        if percentAntiprotons(j,i) ~= 0
            count = count + 1;
            total = total + percentAntiprotons(j,i);
        end
    end
    percentRAntiproton(i) = total / count;
end
percentRAntiproton(isnan(percentRAntiproton)) = 0;

figure
h=plot(0:dr:trapRadius-dr,percentRAntiproton);
title('Average Percentaje Error Antiproton Density');
xlabel('r (m)');
ylabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Average Percentaje Error Antiproton Density.fig'])

macroProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        macroProfile(i) = macroProfile(i) + macroAntiprotons(j,i);
    end
end

figure
h=plot(0:dr:trapRadius-dr,macroProfile);
title('Macro-particle''s Antiproton Profile');
xlabel('r (m)');
ylabel('Line Density (part / m^2)');
savefig([pathWrite, 'Plot Macro Antiproton Profile.fig'])

antiprotonProfile = zeros(Nr,1);
for i = 1:Nr
    for j = 1:Nz + 1
        antiprotonProfile(i) = antiprotonProfile(i) + densityA(j,i);
    end
end

percentProfile = abs(antiprotonProfile - macroProfile);
percentProfile = percentProfile ./ antiprotonProfile;
percentProfile(isnan(percentProfile)) = 0;
percentProfile(percentProfile == 1) = 0;
percentProfile = percentProfile * 100;
figure
h=plot(0:dr:trapRadius-dr,percentProfile);
title('Percent Error Antiproton Profile');
xlabel('r (m)');
ylabel('Percentaje Error (%)');
savefig([pathWrite, 'Plot Percent Error Antiproton Profile.fig'])