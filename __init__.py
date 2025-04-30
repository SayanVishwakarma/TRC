"""
z is vertically upwards against gravity





"""
#Importing required packages
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize, differential_evolution
import math
import json
import pandas as pd
from math import pi,ceil
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter
from tabulate import tabulate

#Defining some required constants
g=9.8 #m/s^2
R_earth=6371e3 #meter
dt=0.1 #seconds
latitude=29 / 180 * np.pi #radians
w_earth=2*np.pi/24/60/60 #rad/s


#importing our classes
import rocket
import delta_v
import stage_optimizer
import trajectory