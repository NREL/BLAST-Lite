"""Functions for updating time-varying states"""

import numpy as np

def update_power_state(y0, dx, k, p):
    if y0 == 0:
        if dx == 0:
            dydx = 0
        else:
            y0 = k*(dx**p)
            dydx = y0/dx
    else:
        if dx == 0:
            dydx = 0
        else:
            dydx = k*p*((y0/k)**((p-1)/p))
    return dydx * dx

def update_power_B_state(y0, dx, k, p):
    if y0 == 0:
        if dx == 0:
            dydx = 0
        else:
            y0 = (k*dx)**p
            dydx = y0/dx
    else:
        if dx == 0:
            dydx = 0
        else:
            z = (y0 ** (1/p)) / k
            dydx = (p * (k*z)**p)/z
    return dydx * dx

def update_sigmoid_state(y0, dx, y_inf, k, p):
    if y0 == 0:
        if dx == 0:
            dydx = 0
        else:
            dy = 2 * y_inf * (1/2 - 1 / (1 + np.exp((k * dx) ** p)))
            dydx = dy / dx
    else:
        if dx == 0:
            dydx = 0
        else:
            x_inv = (1 / k) * ((np.log(-(2 * y_inf/(y0-y_inf)) - 1)) ** (1 / p) )
            z = (k * x_inv) ** p
            dydx = (2 * y_inf * p * np.exp(z) * z) / (x_inv * (np.exp(z) + 1) ** 2)
    return dydx * dx