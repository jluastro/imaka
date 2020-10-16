# Eden McEwen
# July 23rd 2020
# A class for pulling velocity

import math
import time
import os
import sys
import numpy as np

##################################### Point class helps us out
                
class TPoint(object):
    def __init__(self, t, x, y, val, srate=153):
        '''Defines x and y variables'''
        self.y = y
        self.x = x
        self.t = t
        self.val = val
        self.srate = srate
    ##### REPRESENTATION #####
    def __repr__(self):
        return "(%s, (%s,%s))"%(self.t, self.x, self.y)    
    def __str__(self):
        return "TPoint(%s,%s,%s)=%s"%(self.t, self.x, self.y, self.val)    
    def get_pt(self):
        return (self.x, self.y)
    ##### DISTANCE #####
    def dist(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        return math.sqrt(dx**2 + dy**2)    
    def tdist(self, other):
        return abs(self.t - other.t)
    ##### VELOCITY #####
    def vel(self):
        if self.t == 0: return 0
        return self.dist(TPoint(0,7,7,0))/self.t    
    def vel_ms(self):
        plt_v = self.vel()
        return plt_v * (2.2/8) * self.srate
    def vel_knots(self):
        plt_v = self.vel()
        return plt_v * (2.2/8) * self.srate * 1.943844
    ##### ANGLE FROM CENTER #####
    def dir_r(self):
        return math.atan2(self.x-7, self.y-7)   
    def dir_d(self):
        dirr = (math.degrees(self.dir_r()) + 180) % 360
        return dirr if dirr >=0 else dirr + 360
    ##### COMPARE #####
    def close(self, other):
        dist = self.dist(other)
        tdist = self.tdist(other)
        return dist < 2 and tdist <=2
    def closer(self, pt1, pt2):
        dist1 = self.dist(pt1)
        dist2 = self.dist(pt2)
        if dist1 == dist2:
            tdist1 = self.tdist(pt1)
            tdist2 = self.tdist(pt2)
            if tdist1 == tdist2:
                return pt1 if pt1.val > pt2.val else pt2
            else:
                return pt1 if tdist1 < tdist2 else pt2
        else:
            return pt1 if dist1 < dist2 else pt2
    def run_closest(self, runs, idxs):
        c_index = idxs[0]
        closest = runs[c_index][-1]
        for i in idxs:
            tmp = self.closer(runs[i][-1], closest)
            if closest != tmp:
                closest = tmp
                c_index = i
        return c_index
    def run_comp(self, run):
        if self.t == run[-1].t:
            return False
        else:
            return self.close(run[-1])