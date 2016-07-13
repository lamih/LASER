# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import math


nc = 0;# probes number defined later
# probes locations array : (x1 y1 z1 x2 y2 ...xnc ync znc) size (probes number)x3
x = [];
# velocity array : (u1(t1) v1(t1) w1(t1) u2(t1) ...vnc(t1) wnc(t1) u1(t2) v1(t2)...vnc(tn) wnc(tn)) size (probes number)x3x(iterations number)
u = [];
taille = 0;

with open("U") as f:
# probes locations reading
    while taille is 0:
        line = f.readline();
        if line[2] is 'P':
            i = 3;
            while line[i] is not "(" :
                i = i + 1;
            j = i+1;
            while line[j] is not " " :
                j = j + 1;
            x.append(float(line[i+1:j]));
            i = j+2;
            while line[i] is not " " :
                i = i + 1;
            x.append(float(line[j+2:i]));
            j = i + 2;
            while line[j] is not ")" :
                j = j + 1;
            x.append(float(line[i+2:j]));
            nc = nc + 1;
        else:
            taille = 1;
    next(f);
# velocities reading
    for line in f:    
        i = 0;
        for l in range (0,nc):
            while line[i] is not "(" :
                i = i + 1;
            j = i+1;
            while line[j] is not " " :
                j = j + 1;
            u.append(float(line[i+1:j]))
            i = j+1;
            while line[i] is not " " :
                i = i + 1;
            u.append(float(line[j+1:i]))
            j = i + 1;
            while line[j] is not ")" :
                j = j + 1;
            u.append(float(line[i+1:j]))

ul = np.array(u);
n = ul.shape[0];

time = int(n / 3 / nc);

# temporal averaging of 3 components of velocity   
# suppose an uniform time step 
averageV = np.zeros((3*nc));
for t in range (0,time): 
    for l in range (0,nc):
        for i in range (0,3):
            averageV[3*l+i] = averageV[3*l+i] + ul[3*t*nc+3*l+i];
            
averageV = averageV / time;

# fluctuating velocity
for t in range (0,time): 
    for l in range (0,nc): 
        ul[3*t*nc+3*l:3*t*nc+3*l+2] = ul[3*t*nc+3*l:3*t*nc+3*l+2] - averageV[3*l:3*l+2];

# component of velocity correlated    
v = np.zeros((nc*time));

for t in range (0,time): 
    for l in range (0,nc):
        v[t*nc+l] = v[t*nc+l] + ul[3*nc*t+3*l];   # u
#        v[t*nc+l] = v[t*nc+l] + ul[3*nc*t+3*l+1]; # v
#        v[t*nc+l] = v[t*nc+l] + ul[3*nc*t+3*l+2]; # w

# calulation of \overline{v'(x0)v'(x0+dx) and rms for each probe}        
average = np.zeros((nc));
rms = np.zeros((nc));
for t in range (0,time): 
    for l in range (0,nc):
        average[l] = average[l] + v[t*nc]*v[t*nc+l];
        rms[l] = rms[l] + v[t*nc+l]*v[t*nc+l];

# normalized two points correlations       
average = average / time;
rms = rms / time;
for l in range (0,nc):
    average[l] = average[l] / np.sqrt(rms[0]) / np.sqrt(rms[l]);

# distance        
d = np.zeros((nc));
for l in range (0,nc):
    d[l] = np.sqrt((x[3*l]-x[0])*(x[3*l]-x[0])+(x[3*l+1]-x[1])*(x[3*l+1]-x[1])+(x[3*l+2]-x[2])*(x[3*l+2]-x[2]));
    
plt.subplot(2,1,1)    
plt.plot(d,average,'+') 


      
