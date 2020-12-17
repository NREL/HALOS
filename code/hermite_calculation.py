# -*- coding: utf-8 -*-
"""

HALOS miscellaneous scripts

This module calculates a Hermite expansion 

"""
import math
import numpy as np
import scipy
import scipy.stats


class HermiteCalculation(object):
    def __init__(self):
        """
        This class calculates a Hermite expansion for a collection of 
        measurement-based locations.  
        """
        self.hermite_terms = np.array([
            [1,0,0,0,0,0,0],
            [0,1,0,0,0,0,0],
            [-1,0,1,0,0,0,0],
            [0,-3,0,1,0,0,0],
            [3,0,-6,0,1,0,0],
            [0,15,0,-10,0,1,0],
            [-15,0,45,0,-15,0,1]
            ],dtype=float)
        self.gaussian_moments = np.array([1,0,1,0,3,0,15],dtype=float)
        
    def calculate_hermite_const(self,moments,i,j):
        return (moments*self.hermite_terms[i]).sum() * (moments*self.hermite_terms[j]).sum()

    def orthonormal_const(self,n):
        return (math.factorial(n))**-0.5

    def evaluate_hermite(self,x,herm_idx):
    #    print(x,herm_idx)
        Hx = 0
        for i in range(7):
            Hx += self.hermite_terms[herm_idx,i]*(x**i)
    #        print(i, hermite_terms[herm_idx,i], Hx)
        return Hx

    def full_hermite_eval(self,x,y,mx,my,sx,sy,moments):
        dx = (x-mx)/sx
        dy = (y-my)/sy
        pdf = (1./np.sqrt(2*np.pi*sx*sy))*np.exp(-0.5*(dx**2+dy**2))
        
        herm_sum = 0.
        for i in range(7):
            for j in range(7):
                herm_sum +=  self.calculate_hermite_const(moments,i,j) * self.evaluate_hermite(dx,i) * self.evaluate_hermite(dy,j) / (math.factorial(i) * math.factorial(j))
        return pdf * herm_sum
    


def calculate_moment_xy(power_x,power_y,func,pts_per_dim,lb1,ub1,lb2,ub2,fname=None):
    """
    calcuate \mu_{xy} across a known support.
    """
    xs = lb1 + (0.5 + np.arange(pts_per_dim,dtype=float))*(ub1-lb1)/pts_per_dim
    ys = (lb2 + (0.5 + np.arange(pts_per_dim,dtype=float))*(ub2-lb2)).reshape([pts_per_dim,1])/pts_per_dim
    evals = np.zeros([pts_per_dim,pts_per_dim],dtype=float)
    for i in range(pts_per_dim):
        evals[i] = func(xs,ys[i])
    if fname != None: 
#        normal_tests.plot_obj_heatmap(evals,fname)
        pass
    if power_x == 0: 
        xs = 1.
    if power_y == 0:
        ys = 1.
    return ((xs**power_x) * evals * (ys**power_y)).sum() * ((ub1-lb1)*(ub2-lb2)) / ((pts_per_dim)*(pts_per_dim))


def calculate_1d_moment(xs,mx,sx,lp,power):
    p = scipy.stats.norm.pdf((xs-mx)/sx)
    px = (xs**power) * p * lp
    return px.sum()
    
    
def gaussian_2d_func(x,y,mx=0,my=0,sx=1,sy=1):
    return scipy.stats.norm.pdf((x-mx)/sx)*scipy.stats.norm.pdf((y-my)/sy)/(sx*sy)


if __name__ == "__main__":
    pts_per_dim = 50
    dist_per_pt = 10 / 50.
    xs = -5 + (0.5 + np.arange(pts_per_dim,dtype=float))*(10)/pts_per_dim
    print(xs)
    #lp = 
    
    map1 = scipy.zeros([pts_per_dim,pts_per_dim],dtype=float)
    map2 = scipy.zeros([pts_per_dim,pts_per_dim],dtype=float)
    
    for i in range(pts_per_dim):
        print (i,"working")
        for j in range(pts_per_dim):
            map1[i,j] = gaussian_2d_func(xs[i],xs[j],mx=0,my=0,sx=1,sy=1)
            map2[i,j] = full_hermite_eval(xs[i],xs[j],0,0,1,1,gaussian_moments)
    print(map1, map2)

#normal_tests.plot_obj_heatmap(map1,"gaussian_map.pdf")
#normal_tests.plot_obj_heatmap(map2,"hermite_map.pdf")


