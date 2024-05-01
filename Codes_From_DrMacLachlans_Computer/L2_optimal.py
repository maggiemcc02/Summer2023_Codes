


# In this file I have the codes needed to compute
# the mesh density function that minimizes the H1 norm


# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


from numerical_integration import *


def interior_u_xx(Uc, Ub, Ua, xc, xb, xa):

  """
  Uc : U_i+1 
  Ua : U_i-1
  xc, xb, xa : corresponding grid points
  """

  U0 = Ua # U_i-1
  U1 = Ub # Ui
  U2 = Uc # U_i+1

  h1 = xb - xa 
  h2 = xc - xb


  # compute the derivative u_xx

  u_xx_top = h1 * Uc + h2 * Ua - (h1 + h2) * Ub

  u_xx_bottom = (1/2)* (h1*h2**2 + h2 * h1**2)

  u_xx = u_xx_top / u_xx_bottom


  return u_xx



def endpoint_u_xx(Ua, Ub, Uc, xa, xb, xc):
  """
  Uc = U_i+2 or U_i-2
  Ub : U_i 
  Ua : Ui+1 or U_i-1
  xc, xb, xa : corresponding grid points
  """

  hb = xb - xa
  hc = xc - xb



  term_a = 2 / (hb**2 + hb*hc)

  term_b = -2 / (hb*hc)

  term_c = 2 / (hc**2 + hb*hc)


  u_xx_a = term_a * Ua

  u_xx_b = term_b * Ub

  u_xx_c = term_c * Uc

  u_xx = u_xx_a + u_xx_b + u_xx_c


  return u_xx





# Mesh density function for minimizing L2 norm







# code to generate M values on a given grid

def M_calc_optimal_L2(U, grid):

  """
  U : Vector of Approximations

  grid : mesh

  """




  u_dpr = [] # empty list to save M values


  u_dpr.append(endpoint_u_xx(U[0][0], U[1][0], U[2][0], grid[0][0],  grid[1][0], grid[2][0])) # add M value at first grid point using forward approx

  for i in range(len(grid)-2): # add M values for interior nodes

    # print("finding the derivative of term", i+2)
    # print()


    # set values of U and y to pass to mesh density function

    U1 = U[i+1][0]
    U0 = U[i][0]
    U2 = U[i+2][0]
    x0 = grid[i][0]
    x1 = grid[i+1][0]
    x2 = grid[i+2][0]

    val = interior_u_xx(U2, U1, U0, x2, x1, x0) # compute mesh density
    
    u_dpr.append(val)



  u_dpr.append(endpoint_u_xx(U[-3][0], U[-2][0], U[-1][0],  grid[-3][0], grid[-2][0], grid[-1][0])) # add M value for last grid point

  u_dpr = np.array([u_dpr]).T


  # compute the desired integral 

  u_dpr_power = abs(u_dpr)**(2/5)

  # int_u_dpr_power = trap_nonuni(u_dpr_power, grid)

  int_u_dpr_power = int_constant(u_dpr_power, grid)

  alpha_int = ((1 / (grid[-1][0] - grid[0][0])) * int_u_dpr_power )**(5)

  # compute the mesh density function 

  rho = ( 1 + (1/alpha_int) * (abs(u_dpr))**2)**(1/5)

  return rho




