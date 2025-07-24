"""
z is vertically upwards against gravity





"""
#Importing required packages
import sys
import os
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
import utils
import time
import copy
from PIL import Image, ImageDraw, ImageFont
#Defining some required constants
g=9.8 #m/s^2
R_earth=6371e3 #meter
dt=10 #seconds
total_time=60*60
latitude= 29 / 180 * np.pi #radians
w_earth=2*np.pi/24/60/60 #rad/s
G_force_limit=4


#importing our classes
import delta_v
import rocket
import stage_optimizer
import trajectory
import tank
#import path_planning
