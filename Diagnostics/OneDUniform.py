# Written by Daniel Duque
# Last modified on 04/02/0202
# Plot a Phase space animation of an arbitrary number of particles in a 1d plasma at fixed radius

import csv
import matplotlib
import matplotlib.animation
import math
import copy
import numpy
import scipy.fftpack
from scipy import signal
from scipy.constants import Boltzmann
#----------------------------------------------
#INPUT:
#Name of the files where the data from the simulation was stored:
positionsFile = "PositionsElectrons.csv"
speedsFile = "SpeedsElectrons.csv"
timesFile = "Times.csv"
potentialsFile = "PotentialEnergies.csv" #total potential energy of the system at each time step
trapParametersFile = "TrapParameters.csv"
plasmaParametersFile = "ElectronsParameters.csv"
trapPotentialFile = 'TrapPotential.csv'
selfPotentialFile = 'SelfPotential.csv'

#Output parameters

#Do you want to plot the electrostatic potential of the trap
plotTrapPotential = True
saveFileNameTrap = '3 Electrodes Trap.png'

#Do you want to plot the energy of the system through the simulation?
computeEnergy = True
saveFileNameEnergy = 'MP Energy Early.png'

#Do you want to plot the central density as a function of time?
#This will create an imaginary grid point at the centre of the trap, and compute the density there
#The size of this grid would be exactly the same as the real grid points
computeRhoCentre = True
saveFileNameCentralDensity = 'MP Density Early.png'

#Do you want to analyse the temperature?
#This will make a plot of temperature vs time obtained from the average kinetic energy
#And plot in the gif how the Maxwell-Boltzmann should look like for this temperature
computeTemperature = True
saveFileNameTemperature = 'MP Temperature Early.png'
#Do you want to compute gifs?
#Phase Space, charge density, and speed distribution evolution
#This can be a bit slow, more than 10 times the simulation time
computeGifs = True
saveFileNameGif = 'MP Evolution Early.gif'
velocityBinsNumber = 32
#----------------------------------------------
#Data unpacking
#----------------------------------------------
#Unpack parameters of the trap. Look at the C++ output for the extract method to see how the parameters are formatted in the file
trapParametersMatrix = list(csv.reader(open(trapParametersFile), delimiter=','))
for row in trapParametersMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
#See in C++ extract what each number represents
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
massMacro = plasmaParametersMatrix[2][0]
chargeMacro = plasmaParametersMatrix[3][0]
weightMacro = massMacro / massSingle
#Read the speed of all macroparticles
#These are at different times than positions (Leap frog method)
#Linear interpolate to estimate between given values, delete the last one
speedsMatrix = list(csv.reader(open(speedsFile), delimiter=','))
for row in speedsMatrix:
    for i in range(len(row) - 1):
        row[i] = (float(row[i]) + float(row[i + 1])) / 2
    del row[len(row) - 1]
#Compute the kinetic energy
kineticMatrix = copy.deepcopy(speedsMatrix)
for row in kineticMatrix:
    for i in range(len(row)):
        row[i] = massMacro * row[i] * row[i] / 2
#Read the positions, the first element corresponds to the radius index in the trap grid
#Remove the last element. We have one less speed from before, then make them all the same length as the speeds matrix
positionsMatrix = list(csv.reader(open(positionsFile), delimiter=','))
radiusIndex = int(positionsMatrix[0][0])
for row in positionsMatrix:
    del row[0]
    del row[len(row) - 1]
    for i in range(len(row)):
        row[i] = float(row[i])
#Read times, delete last element for the same reason
timesMatrix = list(csv.reader(open(timesFile), delimiter=','))
for row in timesMatrix:
    del row[len(row) - 1]
    for i in range(len(row)):
        row[i] = float(row[i])
#Read potentials and delete last element
potentialsMatrix = list(csv.reader(open(potentialsFile), delimiter=','))
for row in potentialsMatrix:
    del row[len(row) - 1]
    for i in range(len(row)):
        row[i] = float(row[i])
#Read potential of the trap, remember it is embeded as a 1d array
transitionMatrix = list(csv.reader(open(trapPotentialFile), delimiter=','))
for row in transitionMatrix:
    for i in range(len(row)):
        row[i] = float(row[i])
#Now arrange the 1D array as an appropriate 2D array as the trap grid
trapPotentialMatrix = []
for i in range(len(transitionMatrix)):
    if i % (Nz + 1) == 0:
        aRow = []
        aRow.append(transitionMatrix[i][0])
        trapPotentialMatrix.append(aRow)
    else:
        trapPotentialMatrix[-1].append(transitionMatrix[i][0])
#------------------------------------------------
#Different things you can get:
#------------------------------------------------   
#Trap Potential
#------------------------------------------------
#Here we will only plot the electrostatic well in 1D in which the particles are (depending on r index)
if plotTrapPotential:
    zAxis = [0]
    for i in range(int(Nz)):
        zAxis.append(zAxis[-1] + hz)
    #Now lets plot the field in the range in which particles move, and fit a parabola to it
    zMax = 0 
    zMin = trapParametersMatrix[5][0]
    for row in positionsMatrix:
        for element in row:
            if element > zMax:
                zMax = element
            if element < zMin:
                zMin = element
    shortZAxis = []
    shortTrap = []
    for i in range(len(zAxis)):
        if zAxis[i] >= zMin and zAxis[i] <= zMax:
            shortZAxis.append(zAxis[i])
            shortTrap.append(trapPotentialMatrix[radiusIndex][i])
    #Make a quadratic fit of the are of interest
    coefs = numpy.polynomial.polynomial.polyfit(shortZAxis, shortTrap, 2)
    ffit = numpy.polynomial.polynomial.polyval(shortZAxis, coefs)
    #Plot them
    figTrap, (anAxis1, anAxis2) = matplotlib.pyplot.subplots(2)
    anAxis1.set_title('Trap Electrostatic Potential')
    anAxis1.set_xlabel('z (m)')
    anAxis1.set_ylabel('Potential (V)')
    anAxis1.plot(zAxis, trapPotentialMatrix[radiusIndex]) 
    anAxis1.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    anAxis2.set_title('Zone of Interest')
    anAxis2.set_xlabel('z (m)')
    anAxis2.set_ylabel('Potential (V)')
    anAxis2.plot(shortZAxis, shortTrap)
    anAxis2.plot(shortZAxis, ffit, color = 'red', linewidth = 0.6)
    anAxis2.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    expected = (math.sqrt(2 * coefs[2] * chargeMacro / massMacro)) / (math.pi * 2)
    matplotlib.pyplot.figtext(0.5, 0.01, 'Single particle natural frequency of %.2E Hz in the well' % expected, wrap=True, horizontalalignment='center')
    figTrap.align_labels()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(saveFileNameTrap)
#------------------------------------------------
#ENERGY
#------------------------------------------------
if computeEnergy:
    totalKinetic = []
    #Add total kinetic energy of all particles every time step
    for i in range(len(kineticMatrix[0])):
        KE = 0
        for row in kineticMatrix:
            KE += row[i]
        totalKinetic.append(KE)
    totalEnergy = [x + y for x, y in zip(totalKinetic, potentialsMatrix[0])]
    figEnergy, (axOne, axTwo, axThree) = matplotlib.pyplot.subplots(3, sharex = True)
    axOne.set_title('Potential Energy')
    axOne.set_ylabel('U (V)')
    axOne.plot(timesMatrix[0], potentialsMatrix[0]) 
    axOne.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    axTwo.set_title('Kinetic Energy Energy')
    axTwo.set_ylabel('KE (V)')
    axTwo.plot(timesMatrix[0], totalKinetic)
    axTwo.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    axThree.set_title('Total Energy')
    axThree.set_xlabel('time (s)')
    axThree.set_ylabel('U + KE (V)')
    axThree.plot(timesMatrix[0], totalEnergy)
    axThree.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    figEnergy.align_labels()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(saveFileNameEnergy)

if computeRhoCentre:
    zInterest = trapLength / 2
    #Volume of macro particles at fixed radial index (they are all at the same index in r)
    volume = 0 
    if radiusIndex == 0:
        volume = math.pi * hz * hr * hr / 4
    else:
        volume = math.pi * hz * hr * 2 * radiusIndex * hr
    #Central density at every time
    centralRho = []
    for t in range(len(timesMatrix[0])):
        aDensity = 0
        for row in positionsMatrix:
            if row[t] <= zInterest + hz and row[t] >= zInterest - hz:
                dz = abs(row[t] - zInterest)
                weight = dz / hz
                aDensity += (1 - weight) * chargeMacro / volume
        centralRho.append(aDensity)
    figDensity, (axis1, axis2) = matplotlib.pyplot.subplots(2)
    axis1.set_title('Central Charge Density')
    axis1.set_ylabel('ρ (C/m$^3$)')
    axis1.set_xlabel('time (s)')
    axis1.plot(timesMatrix[0], centralRho)
    axis1.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    #Fourier transform of charge density evolution
    # Number of samplepoints
    N = len(centralRho)
    # sample spacing
    T = timesMatrix[0][1] - timesMatrix[0][0]
    clean = scipy.signal.detrend(centralRho) #Remove the 0th frequency obtained from the average
    yf = scipy.fftpack.rfft(clean)
    xf = numpy.linspace(0.0, 1.0/(2.0*T), N/2)    
    axis2.set_title('Density Fourier Transform')
    axis2.set_ylabel('Amplitude')
    axis2.set_xlabel('Frequency (Hz)')
    axis2.plot(xf, 2.0/N * numpy.abs(yf[:N//2]))
    axis2.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(saveFileNameCentralDensity)
#------------------------------------------------
#Temperature analysis
#------------------------------------------------
if computeTemperature:
    temperatures = []
    #Add total kinetic energy of all particles every time step
    for i in range(len(kineticMatrix[0])):
        KE = 0
        for row in kineticMatrix:
            KE += row[i]
        temperatures.append(2 * KE / (len(kineticMatrix) * weightMacro * scipy.constants.Boltzmann))
    
    figTemp, (oneAxis) = matplotlib.pyplot.subplots(1)
    oneAxis.set_title('Temperature Evolution')
    oneAxis.set_xlabel('time (s)')
    oneAxis.set_ylabel('Temperature (K)')
    oneAxis.plot(timesMatrix[0], temperatures) 
    oneAxis.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    matplotlib.pyplot.figtext(0.5, 0.01, 'Average temperature of %.2E K, obtained from average KE' % numpy.mean(temperatures), wrap=True, horizontalalignment='center')
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(saveFileNameTemperature)
#------------------------------------------------
#Density, Phase Space, and Speed distribution Gifs
#------------------------------------------------
if computeGifs:
    #These are the ranges for the plots. The max and min value of z between all the particles at all times
    zMax = 0 
    zMin = trapParametersMatrix[5][0]
    for row in positionsMatrix:
        for element in row:
            if element > zMax:
                zMax = element
            if element < zMin:
                zMin = element
    deltaZ = zMax - zMin #Length of axis
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
    #Volume of macro particles at fixed radial zone (they are all at the same zone)
    volume = 0 
    if radiusIndex == 0:
        volume = math.pi * hz * hr * hr / 4
    else:
        volume = math.pi * hz * hr * 2 * radiusIndex * hr
    
    #Determine how many gridpoints are within the z ranges of the simulation and what z are them in
    gridPointsZ = []
    for i in range(int(trapParametersMatrix[3][0]) + 1):
        dummy = i * hz
        if dummy >= zMin:
            if len(gridPointsZ) == 0:
                gridPointsZ.append(dummy - hz)
            gridPointsZ.append(dummy)
            if dummy > zMax:
                break
    #What sign is the charge of macro particle
    if chargeMacro < 0:
        sign = -1
        maxDensity = 1
    else:
        sign = 1
        maxDensity = -1
    #Calculate the max charge density at anytime through the simulation.
    #This will just determine the length to be used in the axis later
    for i in range(len(positionsMatrix[0])):
        densities = [0 for l in gridPointsZ]
        for row in positionsMatrix:
            for index in range(len(gridPointsZ)):
                if gridPointsZ[index] + hz >= row[i]:
                    weight = (row[i] - gridPointsZ[index])/hz
                    densities[index + 1] += weight * chargeMacro / volume
                    densities[index] += (1 - weight) * chargeMacro / volume
                    break
        if sign == 1:
            candidate = max(densities)
            if candidate > maxDensity:
                maxDensity = candidate
        else:
            candidate = min(densities)
            if candidate < maxDensity:
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
    for i in range(len(speedsMatrix[0])):
        frequencies = [0 for a in speedBins]
        for row in speedsMatrix:
            for index in range(len(speedBins)):
                if row[i] >= speedBins[index] - halfL and row[i] <= speedBins[index] + halfL:
                    frequencies[index] += 1
                    break
        candidate = max(frequencies)
        if candidate > maxFrequency:
            maxFrequency = candidate               
    maxFrequency = maxFrequency / (binLength * len(positionsMatrix))
    #Three subplots, one for Density, another one for Phase Space and other one for speed distribution
    #----------
    #Give Format to axis and graphs
    #----------
    fig, (ax1, ax2, ax3) = matplotlib.pyplot.subplots(1,3, figsize = (12, 5))
    deltaBar = gridPointsZ[-1] - gridPointsZ[0] + hz
    densities = [0 for i in gridPointsZ]
    
    ax1.set_xlim(gridPointsZ[0] - hz/2 - 0.15 * deltaBar, gridPointsZ[-1] + hz/2 + 0.15 * deltaZ)
    ax1.set_ylim(0, maxDensity)
    rects = ax1.bar(gridPointsZ, densities, width = hz)
    ax1.set_title("Density Evolution", y = 1.08)
    ax1.set_ylabel("Charge Density (C/m$^3$)")
    ax1.set_xlabel("z (m)")
    ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    x0,x1 = ax1.get_xlim()
    y0,y1 = ax1.get_ylim()
    ax1.set_aspect(abs(x1-x0)/abs(y1-y0))
    
    ax2.set_xlim(zMin - 0.15 * deltaZ, zMax + 0.15 * deltaZ)
    ax2.set_ylim(vMin - 0.15 * deltaV, vMax + 0.15 * deltaV)
    line, = ax2.plot([], [], 'ro')
    ax2.set_title("Phase Space", y = 1.08)
    ax2.set_ylabel("v (m/s)")
    ax2.set_xlabel("z (m)")
    ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)  
    x0,x1 = ax2.get_xlim()
    y0,y1 = ax2.get_ylim()
    ax2.set_aspect(abs(x1-x0)/abs(y1-y0))
    
    frequencies = [0 for i in speedBins]
    ax3.set_xlim(vMin - 0.15 * deltaV, vMax + 0.15 * deltaV)
    ax3.set_ylim(0, maxFrequency)
    
    #Parenthesis to plot Maxwell-Boltzmann
    if computeTemperature:
        a = math.sqrt(scipy.constants.Boltzmann * numpy.mean(temperatures) / massSingle)
        v = []
        step = deltaV / 100
        for i in range(100):
            v.append(vMin + i * step)
        prob = []
        for aV in v:
            prob.append(math.sqrt(1 / (2 * math.pi)) * pow(a, -1) * math.exp(- 0.5 * pow(aV, 2) * pow(a, -2)))
        ax3.plot(v, prob, 'r')
    
    histBars = ax3.bar(speedBins, frequencies, width = binLength)
    ax3.set_title("Speed Distribution", y = 1.08)
    ax3.set_ylabel("Density")
    ax3.set_xlabel("v (m/s)")
    ax3.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
    x0,x1 = ax3.get_xlim()
    y0,y1 = ax3.get_ylim()
    ax3.set_aspect(abs(x1-x0)/abs(y1-y0))
    
    matplotlib.pyplot.tight_layout()
    
    timeText = ax1.text(0.52, 0.97, 'Time', horizontalalignment='center', verticalalignment='center', transform=fig.transFigure)
    matplotlib.pyplot.figtext(0.5, 0.01, f'Simulation running {len(positionsMatrix)} macro-particles at a radius index of {radiusIndex}. Macro-particles with mass and charge of {plasmaParametersMatrix[2][0]} kg and {plasmaParametersMatrix[3][0]} C respectively.', wrap=True, horizontalalignment='center')
    patches = list(rects) + [line,] + [timeText] + list(histBars)
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
                    densities[index + 1] += weight * chargeMacro / volume
                    densities[index] += (1 - weight) * chargeMacro / volume
                    break
        for bar,h in zip(rects, densities):
            bar.set_height(h)
        #Speed Distribution
        frequencies = [0 for i in speedBins]
        for row in speedsMatrix:
            for index in range(len(speedBins)):
                if row[i] >= speedBins[index] - halfL and row[i] <= speedBins[index] + halfL:
                    frequencies[index] += 1
                    break
        for bar,h in zip(histBars, frequencies):
            bar.set_height(h / (binLength * len(positionsMatrix)))
        
        #Phase Space
        x = []
        for row in positionsMatrix:
            x.append(row[i])
        y = []
        for row in speedsMatrix:
            y.append(row[i])
        line.set_data(x, y)
        #Title
        timeText.set_text('t = %.5E s' % timesMatrix[0][i])
        return patches
    
    def init():
        for bar in rects:
            bar.set_height(0)
        for bar in histBars:
            bar.set_height(0)            
        line.set_data([], [])
        return line,
    
    anim = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, frames=len(positionsMatrix[0]), interval=30, blit=True)
    anim.save(saveFileNameGif, writer='imagemagick')

