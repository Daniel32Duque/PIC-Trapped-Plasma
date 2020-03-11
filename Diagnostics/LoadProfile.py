# Written by Daniel Duque
# Last modified on 06/03/2020
# Diagnostic tools for a plasma loaded from an MCP profile

import csv
import matplotlib
import matplotlib.animation
import math
import copy
import numpy
import scipy.fftpack
from scipy import signal
from scipy.constants import Boltzmann
from mpl_toolkits import mplot3d

#----------------------------------------------
#INPUT:
#Name of the files where the data from the simulation was stored:
positionsFile = "PositionsElectrons.csv"
initialDensityFile = "Initial Density Estimate.csv"
trapParametersFile = "Trap Parameters.csv"
plasmaParametersFile = "ElectronsParameters.csv"

#Output parameters

#Do you want to plot the electrostatic potential of the trap
plotTrapPotential = True
saveFileNameInitialDensity = 'Initial Density 3.png'
#----------------------------------------------
#Data unpacking
#----------------------------------------------
#Unpack parameters of the trap. Look at the C++ output for the extract method to see how the parameters are formatted in the file
trapParametersMatrix = list(csv.reader(open(trapParametersFile), delimiter=','))
for row in trapParametersMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
#See in C++ extract what each number represents
trapRadius = trapParametersMatrix[0][0]
Nz, Nr = trapParametersMatrix[3]
hz, hr = trapParametersMatrix[4]
trapLength = trapParametersMatrix[5][0]
#Same as for trap parameters
plasmaParametersMatrix = list(csv.reader(open(plasmaParametersFile), delimiter=','))
for row in plasmaParametersMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
massSingle = plasmaParametersMatrix[0][0]
chargeSingle = plasmaParametersMatrix[1][0]
chargeMacro = plasmaParametersMatrix[2][0]
massMacro = massSingle * chargeMacro / chargeSingle;
macroChargeDensity = 4 * chargeMacro / (math.pi * hz * hr * hr)
#Read the positions, the first element corresponds to the radius index in the trap grid
#Remove the last element. We have one less speed from before, then make them all the same length as the speeds matrix
positionsMatrix = list(csv.reader(open(positionsFile), delimiter=','))
for row in positionsMatrix:
    if len(row) != 2:
        del row[len(row) - 1]
    for i in range(len(row)):
        row[i] = float(row[i])
#Read initial density estimate, remember it is embeded as a 1d array
transitionMatrix = list(csv.reader(open(initialDensityFile), delimiter=','))
for row in transitionMatrix:
    for i in range(len(row)):
        row[i] = float(row[i]) / chargeSingle
#Now arrange the 1D array as an appropriate 2D array as the trap grid
initialDensityEstimateMatrix = []
for i in range(len(transitionMatrix)):
    if i % (Nz + 1) == 0:
        aRow = []
        aRow.append(transitionMatrix[i][0])
        initialDensityEstimateMatrix.append(aRow)
    else:
        initialDensityEstimateMatrix[-1].append(transitionMatrix[i][0])
#------------------------------------------------
#Different things you can get:
#------------------------------------------------   
#------------------------------------------------   
#Initial Distribution
#------------------------------------------------
#Compare the estimated initial density with that obtained from the particle loading
#Plot initial estimate

x = numpy.linspace(0, trapLength, num = int(Nz + 1))
y = numpy.linspace(0, trapRadius, num = int(Nr), endpoint = False)
X, Y = numpy.meshgrid(x, y)
figInitialDensity = matplotlib.pyplot.figure()
ax1 = figInitialDensity.add_subplot(2, 1, 1, projection='3d')
ax1.plot_surface(X, Y, numpy.array(initialDensityEstimateMatrix), rstride=1, cstride=1,
                cmap='plasma', edgecolor='none')
ax1.set_title('Estimated Initial Density')
ax1.set_xlabel('z (m)')
ax1.set_ylabel('r (m)')
ax1.set_zlabel('n (part/m^3)')
ax1.ticklabel_format(style='sci', scilimits=(0,0))

#Now plot the obtained from macro particles
particleDensity = numpy.zeros((int(Nr), int(Nz + 1)))
'''
int indexR{ aRing.getR() };
int indexZ{ (int)floor(aRing.getZ() / hz) };
int indexRHS{ (refTrap.Nz + 1) * indexR + indexZ };
double z{ aRing.getZ() - indexZ * hz };//Distance from left grid point to particle
double weightFactor{ z / hz };//Weight going to grid point on the right (number between 0 and 1)
RHS.coeffRef(indexRHS) += -macroChargeDensity * (1 - weightFactor) / epsilon;//Point to the left
RHS.coeffRef(indexRHS + 1) += -macroChargeDensity * weightFactor / epsilon;//Point to the right
'''
for row in positionsMatrix:
    indexR = int(row[0])
    indexZ = int(numpy.floor(row[1] / hz))
    z = row[1] - indexZ * hz
    weightFactor = z / hz
    particleDensity[indexR,indexZ] += (1 - weightFactor) * macroChargeDensity / chargeSingle
    particleDensity[indexR,indexZ + 1] += weightFactor * macroChargeDensity / chargeSingle
x = numpy.linspace(0, trapLength, num = int(Nz + 1))
y = numpy.linspace(0, trapRadius, num = int(Nr), endpoint = False)
X, Y = numpy.meshgrid(x, y)
ax2 = figInitialDensity.add_subplot(2, 1, 2, projection='3d')
ax2.plot_surface(X, Y, particleDensity, rstride=1, cstride=1,
                cmap='plasma', edgecolor='none')
ax2.set_title('Initial Density Distribution')
ax2.set_xlabel('z (m)')
ax2.set_ylabel('r (m)')
ax2.set_zlabel('n (part/m^3)')
ax2.ticklabel_format(style='sci', scilimits=(0,0))
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(saveFileNameInitialDensity)
        
        
