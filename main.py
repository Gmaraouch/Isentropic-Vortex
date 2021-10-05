from math import pi
import numpy as np
import functions
import matplotlib.pyplot as plt
from drawnow import *
from copy import deepcopy


def draw_fig():
    #figure, axis=plt.subplots(2,2)
    #plt.gca().set_title('t=' + str(t_i))
    #axis[0,0].contourf(x,y,w[0])
    #axis[0,1].contourf(x,y,w[1])
    #axis[1,0].contourf(x,y,w[2])
    #axis[1,1].contourf(x,y,w[3])
    plt.contourf(x, y, w[0])
    plt.gca().set_title('t=' + str(t_i))
    plt.colorbar()
    #plt.show()


#initial conditions
S = 13.5 #the strength of the vortex
Ma = 0.4 #the mach number
gamma = 1.4 #the heat capacity ratio
R = 1.5 #the radius of the vortex
x_c = 10 #the center of the vortex in x
y_c = 10 #the center of the vortex in y
Lx = 20 #the length of the domain in x
Ly = 20 #the length of the domain in y
Nx = 60 #number of grid points in x
Ny = 60 #number of grid points in y
tf = 20
dt = 0.01
dx = Lx/Nx
dy = Ly/Ny

x = np.arange(0, Lx, dx)
y = np.arange(0, Ly, dy)
t = np.arange(0, tf, dt)

x, y = np.meshgrid(x, y)

f = (1-(x-x_c)**2-(y-y_c)**2)/(2*R**2)

rho = (1-(S**2*Ma**2*(gamma-1)*np.exp(2*f))/(8*pi**2))**(1/(gamma-1))
u = (S*(y-y_c)*np.exp(f))/(2*pi*R)
v = 1-(S*(x-x_c)*np.exp(f))/(2*pi*R)
p = rho**gamma/(gamma*Ma**2)
E = p/(rho*(gamma-1))+1/2*(u**2+v**2)

w0 = np.array([rho, rho*u, rho*v, rho*E])
w=np.copy(w0)
#F = functions.xflux(w, gamma)
#G = functions.yflux(w, gamma)

for t_i in t:
    wt = functions.fRK44_2(w, dx, dy, dt, gamma)
    w = wt
    drawnow(draw_fig)

# figure, axis=plt.subplots(2,2)
# axis[0,0].contourf(x,y,w[0])
# axis[0,1].contourf(x,y,w[1])
# axis[1,0].contourf(x,y,w[2])
# axis[1,1].contourf(x,y,w[3])
# plt.show()



