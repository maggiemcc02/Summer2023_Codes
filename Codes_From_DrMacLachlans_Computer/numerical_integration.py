

# In this file I have numerical integration methods
# In particular, I have trapezoidal and Riemann sum methods


# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# composite, nonuniform, trapezoidal rule

def trap_nonuni(f_vals, x_vals):


  num = len(x_vals)


  # define the mesh spacing 

  h_vals = []

  for i in range(1, num):

    h = x_vals[i][0] - x_vals[i-1][0]

    h_vals.append(h)


  # trapezoidal sum

  sum_trap = 0

  for i in range(1, num):

    val = ((f_vals[i][0] + f_vals[i-1][0]) / 2) * h_vals[i-1]

    sum_trap += val 


  
  return sum_trap




# A 'Riemann Integral'
# SUM(i) f(xi) * hi

def int_constant(f_vals, x_vals):



  num = len(f_vals)



  # define the mesh spacing 

  h_vals = []

  for i in range(1, num):

    h = x_vals[i][0] - x_vals[i-1][0]

    h_vals.append(h)


  int_final = 0

  for i in range(1, num):

    val = f_vals[i][0] * h_vals[i-1]

    int_final += val

  return int_final
