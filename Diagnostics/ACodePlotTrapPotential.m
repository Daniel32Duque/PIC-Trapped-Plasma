%{
Written by: Daniel Duque
Last modified on 14 Mar 2020

This script plots the initial electrostatic potential of the Trap
%}

pathRead = [pwd, '\Data Files\'];
pathWrite = [pwd, '\Plots\'];

potential = csvread([pathRead, 'Z Trap Potential.csv']);
parameters = csvread([pathRead, 'Z Trap Parameters.csv']); 
trapRadius = parameters(1,1);
trapLength = parameters(7,1);
Nz = parameters(5,1);
Nr = parameters(5,2);

dr = trapRadius/Nr;
dz = trapLength/Nz;
[X,Y] = meshgrid(0:dr:trapRadius-dr, 0:dz:trapLength);
   
potential = reshape(potential,Nz + 1, Nr);   

figure
h=surf(X,Y,potential);
set(h,'LineStyle','none')
title('Trap Potential');
xlabel('r (m)');
ylabel('z (m)');
zlabel('V');
savefig([pathWrite, 'Plot Trap Potential.fig'])