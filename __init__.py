"""
x is vertically upwards against gravity





"""
#Importing required packages
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math
import json
import pandas as pd
from math import pi,ceil

#Defining some required constants
g=9.8 #m/s^2
R_earth=6.4e6 #meter
dt=0.01 #seconds
latitude=29 #degrees
w_earth=2*np.pi/24/60/60 #rad/s


#importing our classes
import rocket
import delta_v
import stage_optimizer