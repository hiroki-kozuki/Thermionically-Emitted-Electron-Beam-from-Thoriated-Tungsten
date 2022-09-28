# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:15:23 2022

@author: Hiroki
"""

import numpy as np
from numpy import exp
from numpy import polyfit
from numpy import log
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties as unc
from uncertainties import umath
from uncertainties import ufloat
from uncertainties import unumpy # NEW !! Includes inverse trig functions and more!

titleFont = {'fontname':'Bodoni 72','size':10}
axesFont = {'fontname':'CMU Sans Serif','size':8}
ticksFont = {'fontname':'DM Mono','size':7}
errorStyle = {'mew':1,'ms':3,'capsize':3,'color':'blue','ls':''}
pointStyle = {'mew':1,'ms':3,'color':'blue'}
lineStyle = {'linewidth':0.5}
histStyle = {'facecolor':'blue','alpha':0.5,'edgecolor':'black'}

#%%
#Import csv data file.
vCathode, iCathode, iFaraday = np.loadtxt(r"C:\Users\ginta\OneDrive\HIROKI\Imperial College London\Year 1\Y1 Term 3\Y1 Group Project\Cathode and Faraday cup csv data.csv", delimiter=",", skiprows=1, unpack=True)

vCathode1, vCathode2, iCathode1, iCathode2, iFaraday1, iFaraday2 = vCathode[0:23], vCathode[24:42], iCathode[0:23], iCathode[24:42], iFaraday[0:23], iFaraday[24:42]
#52

# Random Uncertainties in input voltage, input current, and output current due to fluctuations in multimetre.
vCathodeUnc = 0.01 # V # 0.01
iCathodeUnc = 0.05 # A # 0.05
iFaradayUnc = 0.019 # nA 0.01
RCathode1, RCathode2 = vCathode1/iCathode1, vCathode2/iCathode2 # Cathode Power
RCathode1Unc, RCathode2Unc = np.sqrt((vCathodeUnc/iCathode1)**2+(vCathode1*iCathodeUnc/(iCathode1)**2)**2), np.sqrt((vCathodeUnc/iCathode2)**2+(vCathode2*iCathodeUnc/(iCathode2)**2)**2) # Uncertainties for input power propagated in quadrature.

SA = np.pi*(0.00058/2)**2 
# cross sectional surface area of Faraday cup wire (m^2).
SAUnc = 2*np.pi*(0.00058/2)*0.00001 # Abs uncertainty in Faraday cup wire area.
jFaraday1, jFaraday2 = iFaraday1/SA, iFaraday2/SA # output current densities
jFaraday1Unc, jFaraday2Unc = np.sqrt((iFaradayUnc/SA)**2+(iFaraday1*SAUnc/(SA**2))), np.sqrt((iFaradayUnc/SA)**2+(iFaraday2*SAUnc/(SA**2))) # Abs uncretianties in output current densities.

#Calculate cathode temperature:

# a (Temperature coefficient of resistance) and reference values at 20 °C (273.15 K) taken from online and empirically. 
a = 0.0045
TRef = 20+273.15
r = 5.65E-8 # resistivity
l = 0.115 # lenght of wire (m)
lUnc = 0.001 # uncertainty in length #0.03
A = np.pi*(0.125E-3/2)**2
# add rho, l and And their unc, propagate error for RRef
RRef = r*l/A
RRefUnc = r*lUnc/A
 
TCathode1 = TRef+((RCathode1/RRef)-1)/a # approx 900 - 1650
TCathode2 = TRef+((RCathode2/RRef)-1)/a # spprox 850 - 1350 
TCathode1Unc = np.sqrt((RCathode1Unc/(a*RRef))**2 + (RCathode1*RRefUnc/(a*RRef**2))**2)
TCathode2Unc = np.sqrt((RCathode2Unc/(a*RRef))**2 + (RCathode2*RRefUnc/(a*RRef**2))**2)

print(RRefUnc)
#%%

# Curve fit 

k = 1.380649E-23 # Boltzmann constant

def exponentialFit1(TCathode1, AG1, W):   # Functions found in lab manuals (equation 1.1, 1.10)
    return AG1*TCathode1**2*np.e**(-W/(k*TCathode1))

def exponentialFit2(TCathode2, AG2, W):   
    return AG2*TCathode2**2*np.e**(-W/(k*TCathode2))

# Initial guesses
guessAG1 = 60 # arbitrary guess of scale factor for trial 1 # 0.05
guessAG2 = 400 # arbitrary guess of scale factor for trial 2 # 0.075
guessW = 2.6*1.602E-19 # reference work function of thoriated tungsten in joules.  
x1 = np.linspace(900, 1650, num=10000) # x values for plotting exponential fit
x2 = np.linspace(820, 1380, num=10000) # x values for plotting exponential fit

#Trial 1 Graph:
curveFit1, covCurveFit1 = curve_fit(exponentialFit1, TCathode1, jFaraday1, p0=[guessAG1, guessW], sigma = jFaraday1Unc, absolute_sigma=True)   # curve_fit finds line of best fit
plt.xlabel("Cathode Temperature (K)", **axesFont)
plt.ylabel("Output Current Density from Faraday Cup (nA $m^{-2}$)", **axesFont)
plt.xticks(**ticksFont)
plt.yticks(**ticksFont)
plt.title('Trial 1: Output Current Density vs. Cathode Temperature', **titleFont)
plt.errorbar(TCathode1, jFaraday1, ls='', xerr = TCathode1Unc, yerr = jFaraday1Unc,  mew=1, ms=3, capsize=3) # plot of data points
plt.plot(x1, exponentialFit1(x1, *curveFit1)) # plot of exponential fit.

plt.savefig(r'C:\Users\ginta\OneDrive\HIROKI\Imperial College London\Year 1\Y1 Term 3\Y1 Group Project\Group Project Graphs\Trial 1.jpeg', dpi = 1000)
plt.show()

#%%

#Trial 2 Graph
curveFit2, covCurveFit2 = curve_fit(exponentialFit2, TCathode2, jFaraday2, p0=[guessAG2, guessW], sigma = jFaraday2Unc, absolute_sigma = True)   # curve_fit finds line of best fit
plt.xlabel("Cathode Temperature (K)", **axesFont)
plt.ylabel("Output Current Density from Faraday Cup (nA $m^{-2}$)", **axesFont)
plt.xticks(**ticksFont)
plt.yticks(**ticksFont)
plt.title('Trial 2: Output Current Density vs. Cathode Temperature', **titleFont)
plt.errorbar(TCathode2, jFaraday2, ls='', xerr = TCathode2Unc, yerr = jFaraday2Unc,  mew=1, ms=3, capsize=2)
plt.plot(x2, exponentialFit2(x2, *curveFit2)) # plot of exponential fit.

plt.savefig(r'C:\Users\ginta\OneDrive\HIROKI\Imperial College London\Year 1\Y1 Term 3\Y1 Group Project\Group Project Graphs\Trial 2.jpeg', dpi = 1000)
plt.show()
#%%
print("Trial 1: AG = %.3e +/- %.3e, W = %.3e +/- %.3e" %(curveFit1[0], np.sqrt(covCurveFit1[0,0]), curveFit1[1], np.sqrt(covCurveFit1[1,1])))
print("Trial 2: AG = %.3e +/- %.3e, W = %.3e +/- %.3e" %(curveFit2[0], np.sqrt(covCurveFit2[0,0]), curveFit2[1], np.sqrt(covCurveFit2[1,1])))


