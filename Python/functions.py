from copy import deepcopy
import numpy as np


def Second_Order_Central_1(u1, u_1, dx):

    const = 1/(2*dx)
    der = const*(u1-u_1)
    return der


def Fourth_Order_Central_1(u2, u1, u_1, u_2, dx):

    der = (-1/12*u2+2/3*u1-2/3*u_1+1/12*u_2)/dx
    return der


def RK44_Der(F, G, dx, dy):

    n, m = F[0].shape
    Ru = np.zeros_like(F)
    Rux = np.zeros_like(F)
    Ruy = np.zeros_like(F)

    for i in range(len(F)):

        for x in range(m):
            if x == 0:
                Rux[i][:, x] = Second_Order_Central_1(F[i][:, x+1], F[i][:, m-1], dx)
            elif x == (m-1):
                Rux[i][:, x] = Second_Order_Central_1(F[i][:, 0], F[i][:, x-1], dx)
            else:
                Rux[i][:, x] = Second_Order_Central_1(F[i][:, x+1], F[i][:, x-1], dx)

        for y in range(n):
            if y == 0:
                Ruy[i][y, :] = Second_Order_Central_1(G[i][y+1, :], G[i][n-1, :], dy)
            elif y == (n-1):
                Ruy[i][y, :] = Second_Order_Central_1(G[i][0, :], G[i][y-1, :], dy)
            else:
                Ruy[i][y, :] = Second_Order_Central_1(G[i][y+1, :], G[i][y-1, :], dy)

        #dw/dt=-dF/dx-dG/dy
        Ru[i] = -Rux[i]-Ruy[i]

    return Ru


def xflux(w, gamma):

    rho = w[0]
    rhou = w[1]
    rhov = w[2]
    rhoE = w[3]

    u = rhou/rho
    v = rhov/rho

    p = (gamma-1)*(rhoE-(1/2)*rho*(u**2+v**2))

    F = []
    F.append(rhou)
    F.append(rhou*u+p)
    F.append(rhou*v)
    F.append(u*(rhoE+p))

    return F


def yflux(w, gamma):

    rho = w[0]
    rhou = w[1]
    rhov = w[2]
    rhoE = w[3]

    u = rhou/rho
    v = rhov/rho

    p = (gamma-1)*(rhoE-(1/2)*rho*(u**2+v**2))

    G = []
    G.append(rhov)
    G.append(rhov*u)
    G.append(rhov*v+p)
    G.append(v*(rhoE+p))

    return G


def fRK44_2(w, dx, dy, dt, gamma):

    w2 = np.zeros_like(w)
    w3 = np.zeros_like(w)
    w4 = np.zeros_like(w)
    wt = np.zeros_like(w)

    w1 = np.copy(w)
    F1 = xflux(w1, gamma)
    G1 = yflux(w1, gamma)
    Rw1 = RK44_Der(F1, G1, dx, dy)

    for i in range(len(w)):
        w2[i] = w[i]+(dt/2)*Rw1[i]

    F2 = xflux(w2, gamma)
    G2 = yflux(w2, gamma)
    Rw2 = RK44_Der(F2, G2, dx, dy)

    for i in range(len(w)):
        w3[i] = w[i]+(dt/2)*Rw2[i]

    F3 = xflux(w3, gamma)
    G3 = yflux(w3, gamma)
    Rw3 = RK44_Der(F3, G3, dx, dy)

    for i in range(len(w)):
        w4[i] = w[i] + dt * Rw3[i]

    F4 = xflux(w4, gamma)
    G4 = yflux(w4, gamma)
    Rw4 = RK44_Der(F4, G4, dx, dy)

    for i in range(len(w)):
        wt[i] = w[i]+(dt/6)*(Rw1[i]+2*Rw2[i]+2*Rw3[i]+Rw4[i])

    return wt

