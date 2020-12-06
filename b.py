from numpy import zeros, arange, zeros
from math import pi, e
import random
import matplotlib. pyplot as plt

class MBdist():
    
    def maxwell_boltzmann(self, temperature, mass, speed):
        mass = mass * 1.6726219*10**(-27)
        k = 1.380648*10**(-23)
        p = 4*pi*speed**2 * (mass / (2*pi*k*temperature))**(3.0/2.0) * \
        e**(-mass*speed**2 /(2*k*temperature))
        return p
    
    #prepared the cumulative weights
    def prepare(self):
        self.z[0] = self.maxwell_boltzmann(self.T, self.m, 0)
        for i in range(1, self.L):
            self.z[i] = self.z[i-1] + self.maxwell_boltzmann(self.T, self.m, i)
        
    def __init__(self, T, m, L):
        self.T = T
        self.m = m
        self.L = L
        self.z = zeros(L)
        self.prepare()

    #takes bins, returns em full with n nummbers 
    def sample_n(self, bins, n):
        #choose rand num from dist
        ran = random.choices(arange(self.L), cum_weights=self.z, k = n)
        #bin it
        for s in ran:
            ind = int(s * (len(bins) - 1)/self.L)
            bins[ind] = bins[ind] + 1
            
    def set_temperature(self, new_t):
        self.T = new_t
        self.prepare()
        
    def set_mass(self, new_m):
        self.m = new_m
        self.prepare()


r=1 #how many subplots
l = 35 #x number of bars in histogram
xdim=1 #number of subplots in a row
T = 300; m = 4; L = 3000
o = MBdist(T, m, L)
for j in range(1, r+1): #produce r subplots
    bins = zeros (l ,int) #bins
    #for some reason it won't print subplots
    o.set_temperature(T+10*j)
    o.sample_n(bins, 50000)
    plt.subplot(r/xdim, xdim, j)
    plt.xticks(arange(0,L,2*L/10))
    plt.bar(arange(0,L,L/l), bins, width =L/l, align ='center', color='blue', edgecolor='cyan')
    plt.xlabel('speed (m/s)')
    plt.ylabel('particle count')
    plt.subplots_adjust(wspace=0.4, hspace=0.55)
title = 'Distribution of the speeds of Helium particles (m=4, T=300)'
plt.title(title)
plt.savefig('a_bonus_plot.png',dpi=200)
del o





