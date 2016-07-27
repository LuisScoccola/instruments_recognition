"""
plotting.py
"""

import scipy
import numpy as np
import matplotlib.pyplot as plt


# plot (1/4 of the) spectrum 
def plotspec(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        plt.plot(freqs, ex, label=name,c=c)
    
# logarithmic scale for x axis [log(frequency) = note]
def plotspec_log(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        plt.plot(np.log2(freqs), ex, label=name,c=c)
    

# logarithmic scale for both axes [log(frequency) = note, log(amplitude) = volume]
def plotspec_loglog(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        ex = np.log2(ex)
        c = next(color)
        plt.plot(np.log2(freqs), ex, label=name,c=c)
 
    
def plotspec_normalized(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        #normalize
        ex /= np.max(np.abs(ex),axis=0)
        plt.plot(freqs, ex, label=name,c=c)

def plotspec_and_harmonics(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        # normalize
        ex /= np.max(np.abs(ex),axis=0)
        #ex /= rms(ex)
        plt.plot(freqs, ex, label=name,c=c)
        plt.scatter(x,y , label=name,c=c)
#        plt.ylim([0,1])

    
def plotspec_and_harmonics_log(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        # in decibels
        ex = np.log2(ex)
        # normalize
        ex /= np.max(np.abs(ex),axis=0)
        plt.plot(freqs, ex, label=name,c=c)
        plt.scatter(x, y, label=name,c=c)
#        plt.ylim([0,1])
    
# no color version
def plotharmonics_nc(names_and_graphs_and_harm_number) :
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        xx = range(len(y))
        plt.scatter(xx, y,label=name)
    
def plotharmonics(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        xx = range(len(y))
        plt.scatter(xx, y,label=name,c=c)
    

def plotharmonics2(names_and_graphs_and_harm_number) :
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(names_and_graphs_and_harm_number))))
    for (name, ex, freqs), (x,y) in names_and_graphs_and_harm_number :
        c = next(color)
        plt.plot(y,label=name,c=c)
