"""
Written by: Daniel Duque

Open a csv file with 2 elements per row
Make an X vs Y plot
"""
import csv
import matplotlib

'''-------------------------------
 Input 
-------------------------------'''
 
pathRead = "./Data Files/"
pathWrite = "./Plots/"

inputFile = pathRead + "Z Electron Temperature Evolution.csv"
imageOutput = pathWrite + "Plot Electron Equilibrium Temperature Evolution.png"
title = "Electron Temperature Evolution"
xlabel = "time (s)"
ylabel = "Temperature (K)"

'''-------------------------------
 End of Input 
-------------------------------'''

x = []
y = []

data = list(csv.reader(open(inputFile), delimiter=','))
for pair in data:
    x.append(float(pair[0]))
    y.append(float(pair[1]))

matplotlib.pyplot.plot(x,y)
matplotlib.pyplot.title(title)
matplotlib.pyplot.xlabel(xlabel)
matplotlib.pyplot.ylabel(ylabel)
matplotlib.pyplot.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(imageOutput)