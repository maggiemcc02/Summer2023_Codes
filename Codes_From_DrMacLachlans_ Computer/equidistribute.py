
# In this file I define my equidistribution code

# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



def equidistribute(oldmesh, rho, uni_grid, x0, xn, n):
  """
  oldmesh : current mesh we are updating with this code

  P : function defined in Russell and Huang which represents the integral of the linearized mesh density function

  rho : value of the mesh density function at each point on the oldmesh

  uni_grid : a uniform grid with the same number of points and on the same interval the oldmesh

  x0, xn : first and last grid points

  n : number of grid points
  """

  
  # compute P(y) for each y in oldmesh

  P_list = [0 for i in range(len(oldmesh))] # P(a) = 0


  for q in range(1, n): # q = second grid point, ..., last grid point

    y1 = oldmesh[q-1][0]
    y2 = oldmesh[q][0]
    M1 = rho[q-1][0]
    M2 = rho[q][0]

    val = (y2 - y1) * ( (M2 + M1) / 2)

    result = P_list[q-1] + val

    P_list[q] = result

  # main code to compute the new grid points

  newmesh = [0 for i in range(n)] # add in boundary condition x(0) = 0
  newmesh[0] = x0

  last_min = P_list[0]
  last_max = P_list[1]
  i = 0
  k = 1

  for j in range(1, n): # j = second grid point, ...,last grid point 

    # compute xi * P(b)

    xi = uni_grid[j][0] # THIS COULD BE REPLACED WITH j/N

    P_b = P_list[-1] # P(b) is last P value in P_list

    xi_pb = xi * P_b # compute xi * P(b)


    # find what interval on the old mesh xi*P(b) is in

    found = False
    i += 1 
    q = k
    l = 0 # iter counter

    min = last_min 
    max = last_max 

    while not found:


      if ((xi_pb > min) and (xi_pb <= max)):

        k = q
        found = True
        last_min = min
        last_max = max



      else:

        l += 1
        q += 1

        min = P_list[q-1]
        max = P_list[q]


    
    # compute the new mesh point


    a = oldmesh[k-1][0] # a = y_k-1

    b = 2 * ( xi_pb - P_list[k-1] ) # b = 2 * (xi*P(b) - P(y_k-1))

    c = rho[k-1][0] + rho[k][0] # c = M(y_k-1) + M(y_k)

    new_x = a + b / c # new_x = y_k-1 + ( 2 * (xi*P(b) - P(y_k-1))) / (M(y_k-1) + M(y_k))

    newmesh[j] = new_x # save new grid point to mesh


  

  # set the last endpoint to be 1 

  newmesh[-1] = xn

  newmesh = np.array([newmesh]).T # make mesh a column vector 

  
  return newmesh



