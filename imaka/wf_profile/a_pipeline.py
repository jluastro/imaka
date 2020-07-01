# Eden McEwen
# June 29th 2020
# A class for the analysis pipeline

import math
import imageio
import time
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Self-witten code
from corr_code import *
from graph_code import *
from data_pipeline import *


class APipe(DataPipe):
    def __init__(self, fits_file, t_file=""):
        self.fits_pull(fits_file)
        self.set_target(t_file)
        
    def track_max(self):
        #TODO: this function will return the list of max points in a list as it evolves over time
        return 0
    
    
    def calc_pix_speed(self):
        #TODO: uses given telescope diam / subap to figure out speed per pixel
        return 0
    
    def get_windspeed(self):
        #TODO: looks up date in windspeed dir, saves windspeed from start of circular buffer
        return 0