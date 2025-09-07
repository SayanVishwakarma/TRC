"""
z is vertically upwards against gravity





"""
#Importing required packages
import sys
import os
import  csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize, differential_evolution, root, root_scalar
import math
import json
import pandas as pd
from math import pi,ceil,isnan
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter
from tabulate import tabulate
import pandas as pd
import tqdm
import time
import copy
from PIL import Image, ImageDraw, ImageFont
from scipy.interpolate import LinearNDInterpolator

#Defining some required constants
g=9.8 #m/s^2
R_earth=6371e3 #meter
dt=10 #seconds
total_time=60*60
latitude= 19.5 / 180 * np.pi #radians
longitude= 41.5 / 180 * np.pi #radians
w_earth=2*np.pi/24/60/60 #rad/s
G_force_limit=3


#importing our classes
import delta_v
import rocket
import stage_optimizer
import trajectory
import tank
import utils


#Creating a function to interpolate combustion data

def interpolate_combustion_data(query_pressure, query_of):
    # Load the data
    data = np.loadtxt("data files/LOXMETHANE.txt", skiprows=1)  # Skip header
    
    pressures = data[:, 0]
    ofs = data[:, 1]
    temps = data[:, 2]
    isps = data[:, 3]
    mws = data[:, 4]
    cstars = data[:, 5]
    gammas = data[:, 6]

    # Stack (pressure, OF) pairs as interpolation input
    points = np.column_stack((pressures, ofs))

    # Create interpolators for each output
    temp_interp = LinearNDInterpolator(points, temps)
    isp_interp = LinearNDInterpolator(points, isps)
    mw_interp = LinearNDInterpolator(points, mws)
    cstar_interp = LinearNDInterpolator(points, cstars)
    gamma_interp = LinearNDInterpolator(points, gammas)

    # Interpolate at desired (pressure, O/F)
    T = temp_interp(query_pressure, query_of)
    Isp = isp_interp(query_pressure, query_of)
    MW = mw_interp(query_pressure, query_of)
    Cstar = cstar_interp(query_pressure, query_of)
    Gamma = gamma_interp(query_pressure, query_of)

    if None in (T, Isp, MW, Cstar, Gamma) or any(np.isnan(x) for x in (T, Isp, MW, Cstar, Gamma)):
        raise ValueError("Interpolation point outside data range or invalid. Pressure= {}, O/F={}".format(query_pressure, query_of))

    return {
        'Pressure': query_pressure,
        'O/F': query_of,
        'Temp (K)': T[()],
        'Isp (m/s)': Isp[()],
        'Mol. Weight': MW[()],
        'C* (m/s)': Cstar[()],
        'Gamma': Gamma[()]
    }

