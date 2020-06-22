"""
Written by: Daniel Duque

This is part of the third step in a Simple Ekick study
Visualise the evolution of electrons or antiprotons at a fixed r position
"""

# Written by Daniel Duque
# Last modified on 04/02/2020
# Diagnostics for a plasma generated from MCP image
#But we will only look at a specific radius
#Delete all other radius

import csv
import matplotlib
import matplotlib.animation
import math
import copy
import numpy
import scipy
from scipy import constants

#----------------------------------------------
#Input
#----------------------------------------------

#How many bins to use in the speed histogram?
velocityBinsNumber = 140 #Try sqrt of number of macro-particles and see how it looks

#Where is the data stored?

pathRead = "./Data Files/"
pathWrite = "./Plots/"

trapParametersFile = pathRead + "Z Trap Parameters.csv"
positionsFile = pathRead + "Z PositionsAntiprotons.csv"
speedsFile = pathRead + "Z SpeedsAntiprotons.csv"
ParametersFile = pathRead + "Z Antiproton Parameters.csv"
densityFile = pathRead + "Z Expected Antiproton Density.csv"
timesFile = pathRead + "Z Times.csv"
rFile = pathRead + "Z rIndex.txt"

#Where do you want to store the animation?
saveFileNameGif = pathWrite + 'Gif Antiproton Evolution.gif'

#----------------------------------------------
#End Input
#----------------------------------------------

matplotlib.pyplot.rcParams['animation.convert_path'] = 'C:\Program Files\ImageMagick-7.0.10-Q16/magick.exe'

#Data unpacking
rIndex = int(list(csv.reader(open(rFile), delimiter=','))[0][0])

trapParametersMatrix = list(csv.reader(open(trapParametersFile), delimiter=','))
for row in trapParametersMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
trapRadius = trapParametersMatrix[0][0]
Nz, Nr = trapParametersMatrix[4]
trapLength = trapParametersMatrix[6][0]
hz = trapLength / Nz
hr = trapRadius / Nr

ParametersMatrix = list(csv.reader(open(ParametersFile), delimiter=','))
for row in ParametersMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
mass = ParametersMatrix[0][0]
charge = ParametersMatrix[1][0]
chargeMacro = ParametersMatrix[2][0]
densityMacro = 4 * chargeMacro / (math.pi * hz * hr * hr)
if(rIndex != 0):
    chargeMacro = 8 * rIndex * chargeMacro
massMacro = mass * chargeMacro / charge;
weightMacro = chargeMacro / charge

positionsMatrix = list(csv.reader(open(positionsFile), delimiter=','))
for row in positionsMatrix:
    del row[0] #This is the radius index
    del row[-1] #The last position doesnt have a speed (they are shifted by feltaT / 2)
    for i in range(len(row)):
        row[i] = float(row[i])
        
speedsMatrix = list(csv.reader(open(speedsFile), delimiter=','))
#Leap frog. Need linear interpolation
for row in speedsMatrix:
    for i in range(len(row) - 1):
        row[i] = (float(row[i]) + float(row[i + 1])) / 2
    del row[len(row) - 1]
    
kineticMatrix = copy.deepcopy(speedsMatrix)
for row in kineticMatrix:
    for i in range(len(row)):
        row[i] = massMacro * row[i] * row[i] / 2
        
timesMatrix = list(csv.reader(open(timesFile), delimiter=','))[0]
del timesMatrix[-1] #Same as positions, one extra time
for i in range(len(timesMatrix)):
    timesMatrix[i] = float(timesMatrix[i])
        
densitiesMatrix = list(csv.reader(open(densityFile), delimiter=','))

#Make the animation:

#Calculate the average temperature
temperatures = []
#Add total kinetic energy of all particles every time step
for i in range(len(kineticMatrix[0])):
    KE = 0
    for row in kineticMatrix:
        KE += row[i]
    temperatures.append(2 * KE / (len(kineticMatrix) * weightMacro * scipy.constants.Boltzmann))
averageTemperature = numpy.mean(temperatures)

#These are the ranges for the plots. The max and min value of z between all the particles at all times
zMax = 0 
zMin = trapLength
for row in positionsMatrix:
    for element in row:
        if element > zMax:
            zMax = element
        if element < zMin:
            zMin = element
deltaZ = zMax - zMin #Length of axis

#This is the equilibrium expected at this radius
expectedDensities = []
theZExpected = []
anIndex = 0
while anIndex * hz <= zMax:
    if anIndex * hz >= zMin:
        theZExpected.append(anIndex * hz)
        expectedDensities.append(float(densitiesMatrix[rIndex * (int(Nz) + 1) + anIndex][0]) / charge)
    anIndex += 1


printDensities = []
printZExpected = []
printIndex = 0
while printIndex * hz <= 0.037:
    if printIndex * hz >= 0.031:
        printZExpected.append(printIndex * hz)
        printDensities.append(float(densitiesMatrix[rIndex * (int(Nz) + 1) + printIndex][0]) / charge)
    printIndex += 1
'''    
f = open("equilibrium150.txt", "w")
for x, y in zip(printZExpected, printDensities):
    f.write(str(x) + " " + str(y) + "\n")
f.close()
'''
#These are the ranges for the y in Phase Space. The max and min value of v between all the particles at all times
vMax = 0
vMin = 0
for row in speedsMatrix:
    for element in row:
        if element > vMax:
            vMax = element
        if element < vMin:
            vMin = element
deltaV = vMax - vMin #Length of axis

#Determine how many gridpoints are within the z ranges of the simulation and what z are them in
gridPointsZ = []
for i in range(int(Nz) + 1):
    dummy = i * hz
    if dummy >= zMin:
        if len(gridPointsZ) == 0:
            gridPointsZ.append(dummy - hz)
        gridPointsZ.append(dummy)
        if dummy > zMax:
            break
        
#Calculate the max charge density at anytime through the simulation.
#This will just determine the length to be used in the axis later
maxDensity = -1
for i in range(len(timesMatrix)):
    densities = [0 for l in gridPointsZ]
    for row in positionsMatrix:
        for index in range(len(gridPointsZ)):
            if gridPointsZ[index] + hz >= row[i]:
                weight = (row[i] - gridPointsZ[index])/hz
                densities[index + 1] += weight * densityMacro / charge
                densities[index] += (1 - weight) * densityMacro / charge
                break
    candidate = max(densities)
    if candidate > maxDensity:
        maxDensity = candidate
        
#Determine the speed bin positions for the histogram
speedBins = []
binLength = deltaV / velocityBinsNumber
currentCentre = vMin + binLength / 2
while True:
    speedBins.append(currentCentre)
    currentCentre += binLength
    if currentCentre > vMax:
        break

#Now calculate the max bin height ever. Use it as axis y axis length
halfL = binLength / 2
maxFrequency = 0;
for i in range(len(timesMatrix)):
    frequencies = [0 for a in speedBins]
    for row in speedsMatrix:
        for index in range(len(speedBins)):
            if row[i] >= speedBins[index] - halfL and row[i] <= speedBins[index] + halfL:
                frequencies[index] += 1
                break
    candidate = max(frequencies)
    if candidate > maxFrequency:
        maxFrequency = candidate               
maxFrequency /= (binLength * len(positionsMatrix)) #Normalize to probability density

#Three subplots, one for Density, another one for Phase Space and other one for speed distribution
#----------
#Give Format to axis and graphs
#----------
fig = matplotlib.pyplot.figure(figsize=(6, 6))
grid = matplotlib.pyplot.GridSpec(4, 4)
ax2 = fig.add_subplot(grid[2:, :2])#Phase Space
ax3 = fig.add_subplot(grid[2:, 2:], sharey=ax2)#speed
ax1 = fig.add_subplot(grid[:2, :2], sharex=ax2)#density
ax4 = fig.add_subplot(grid[:2, 2:])#Text
ax4.axis("off")

deltaBar = gridPointsZ[-1] - gridPointsZ[0] + hz
densities = [0 for i in gridPointsZ]

ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(True)
ax1.set_ylim(0, maxDensity)
rects = ax1.bar(gridPointsZ, densities, width = hz)
ax1.set_ylabel("Particles / m$^3$")
ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
print('z=' + str(zMin - 0.15 * deltaZ) + ',' + str(zMax + 0.15 * deltaZ))
print('v=' + str(vMin - 0.15 * deltaV) + ',' + str(vMax + 0.15 * deltaV))
ax2.set_xlim(zMin - 0.15 * deltaZ, zMax + 0.15 * deltaZ)
ax2.set_ylim(vMin - 0.15 * deltaV, vMax + 0.15 * deltaV)
line, = ax2.plot([], [], 'ro', markersize=1)
ax2.set_ylabel("v (m/s)")
ax2.set_xlabel("z (m)")
ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)  

frequencies = [0 for i in speedBins]
ax3.yaxis.set_visible(False)
ax3.set_xlim(0, maxFrequency)
histBars = ax3.barh(speedBins, frequencies, height = binLength)
ax3.set_xlabel("Prob. Density")
ax3.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)

timeText = ax4.text(0.1, 0.8, 'Time')
ax4.text(0.1, 0.65, f'Fixed r = {rIndex * hr} m')
ax4.text(0.1, 0.5, f'No. macro-particles at r = {len(positionsMatrix)}')
ax4.text(0.1, 0.35, 'Avg. T = %.2f K' % averageTemperature)
patches = list(rects) + [line,] + [timeText] + list(histBars)

#Plot expected equilibrium
ax1.plot(theZExpected, expectedDensities, 'r')

#Plot Maxwell-Boltzmann fit into the speeds histogram
a = math.sqrt(scipy.constants.Boltzmann * numpy.mean(temperatures) / mass)
v = []
step = deltaV / 50
for i in range(50):
    v.append(vMin + i * step)
prob = []
for aV in v:
    prob.append(math.sqrt(1 / (2 * math.pi)) * pow(a, -1) * math.exp(- 0.5 * pow(aV, 2) * pow(a, -2)))
ax3.plot(prob, v, 'r')

matplotlib.pyplot.tight_layout()
'''
f = open("EquilibriumSpeed150.txt", "w")
for x, y in zip(prob, v):
    f.write(str(x) + " " + str(y) + "\n")
f.close()
'''
#----------
#Animation functions
#----------
def animate(i):
    #Density Evolution
    densities = [0 for l in gridPointsZ]
    for row in positionsMatrix:
        for index in range(len(gridPointsZ)):
            if gridPointsZ[index] + hz >= row[i]:
                weight = (row[i] - gridPointsZ[index])/hz
                densities[index + 1] += weight * densityMacro / charge
                densities[index] += (1 - weight) * densityMacro / charge
                break
    for bar,h in zip(rects, densities):
        bar.set_height(h)
    '''
    f = open("Charge_Density-" + str(i) + ".txt", "w")
    for x, y in zip(gridPointsZ, densities):
        f.write(str(x) + " " + str(y) + "\n")
    f.close()
    '''
    
    #Speed Distribution
    frequencies = [0 for i in speedBins]
    for row in speedsMatrix:
        for index in range(len(speedBins)):
            if row[i] >= speedBins[index] - halfL and row[i] <= speedBins[index] + halfL:
                frequencies[index] += 1
                break
    for bar,h in zip(histBars, frequencies):
        bar.set_width(h / (binLength * len(positionsMatrix)))
    '''
    f = open("Speed_Density-" + str(i) + ".txt", "w")
    for x, y in zip(frequencies, speedBins):
        f.write(str(x / (binLength * len(positionsMatrix))) + " " + str(y) + "\n")
    f.close()
    '''
    #Phase Space
    x = []
    for row in positionsMatrix:
        x.append(row[i])
    y = []
    for row in speedsMatrix:
        y.append(row[i])
    line.set_data(x, y)
    '''
    f = open("Phase_Space-" + str(i) + ".txt", "w")
    for xi, yi in zip(x, y):
        f.write(str(xi) + " " + str(yi) + "\n")
    f.close()
    '''
    #Title
    timeText.set_text('t = %.5E s' % timesMatrix[i])
    '''
    f = open("time-" + str(i) + ".txt", "w")
    f.write(str("{:.1f}".format(timesMatrix[i] * 1e9)))
    f.close()
    '''
    #matplotlib.pyplot.savefig(saveImageName + str(i) + '.png')
    return patches

def init():
    for bar in rects:
        bar.set_height(0)
    for bar in histBars:
        bar.set_width(0)            
    line.set_data([], [])
    return line,

anim = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, frames=len(timesMatrix), interval=30, blit=True)
anim.save(saveFileNameGif, writer='imagemagick', extra_args="convert")