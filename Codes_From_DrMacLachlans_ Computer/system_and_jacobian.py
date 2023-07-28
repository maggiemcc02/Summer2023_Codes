# In this notebook I will define the system and jacobian functions
# I will also define the test problem


# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d




# Test problem




def k(u, t):
  return u*(u-1)*(u-t-(3/2))

def k_prime(u, t):
  return 3*u**2 - 2*u*t + (-10*u + 3) / 2 + t



# tridiagonal matrix code



def tridiag(a, b, c, k1 = -1, k2 = 0, k3 = 1):
  """
  a = lower diag
  b = main diag
  c = upper diag

  """

  return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)





# Our System






def F(x, mesh, mesh_spacing, k, u_0, u_n, E):

  """
  x : array of guesses
  mesh : vector of grid points
  mesh_spacing : vector of mesh widths
  k : our function in the ODE
  u_0, u_n : boundary conditions
  E : epsilon
  solution : found solution
  sol_mesh : mesh we found r on

  """
  g_list = []

  for i in range(len(x)):

    # set the grid point
    point = mesh[i][0]
    # set mesh width
    h1 = mesh_spacing[i][0]
    h2 = mesh_spacing[i+1][0]

    # create entries
    if (i == 0): # u_i-1 is u_0

      u_1 = x[i][0]
      u_2 = x[i+1][0]
      u_minus = u_0

      g_list.append( [ ( ( h1*u_2 + h2*u_minus - (h1 + h2)*u_1 ) \
                        / ( (1/2) * (h1*(h2**2) + h2*(h1**2)) ) \
                        - ( (k(u_1, point)) / (E**2) ) ) ] )
   
      
    elif i == (len(x) - 1): # U_i+1 is u_n

      u_1 = x[i][0]
      u_2 = u_n
      u_minus = x[i-1][0]  

      g_list.append( [ ( ( h1*u_2 + h2*u_minus - (h1 + h2)*u_1 ) \
                        / ( (1/2) * (h1*(h2**2) + h2*(h1**2)) ) \
                        - ( (k(u_1, point)) / (E**2) ) ) ] )

    else:

      u_1 = x[i][0]
      u_2 = x[i+1][0]
      u_minus = x[i-1][0]  

      g_list.append( [ ( ( h1*u_2 + h2*u_minus - (h1 + h2)*u_1 ) \
                        / ( (1/2) * (h1*(h2**2) + h2*(h1**2)) ) \
                        - ( (k(u_1, point)) / (E**2) ) ) ] )

      
    

  # create array 
  G = np.array(g_list)


  return G









# Our Jacobian













def JF(x, mesh, mesh_spacing, k, k_prime, u_0, u_n, E):

  """
  x : array of guesses
  mesh : vector of grid points
  mesh_spacing : vector of mesh widths
  k : function in ODE
  k_prime : derivaive of function in ODE
  E : epsilon
  solution : previously found solution
  sol_mesh : mesh we found r on

  """


  main  = []
  upper = []
  lower = []
 

  for i in range(len(x)):

    # set grid point
    point = mesh[i][0]
    # set mesh width
    h1 = mesh_spacing[i][0]
    h2 = mesh_spacing[i+1][0]



    # create diagonals 

    if i == 0: # first row

      main.append( ( (-1*(h1+h2)) / ( (1/2)*(h1*(h2**2) + h2*(h1**2)) ) ) - ( (1 / (E**2)) * ( k_prime(x[i][0], point)) ) )
      #upper.append(0) # for banded matrix
      upper.append( h1 / ( (1/2)* (h1*(h2**2) + h2*(h1**2)) ) )

    elif i == (len(x)-1): # last row

      main.append( ( (-1*(h1+h2)) / ( (1/2)*(h1*(h2**2) + h2*(h1**2)) ) ) - ( (1 / (E**2)) * ( k_prime(x[i][0], point)) ) )
      #lower.append(0) # for banded matrix 
      lower.append( h2 / ( (1/2)* (h1*(h2**2) + h2*(h1**2)) ) )


    else: 

      upper.append( h1 / ( (1/2)* (h1*(h2**2) + h2*(h1**2)) ) )
      main.append( ( (-1*(h1+h2)) / ( (1/2)*(h1*(h2**2) + h2*(h1**2)) ) ) - ( (1 / (E**2)) * ( k_prime(x[i][0], point)) ) )
      lower.append( h2 / ( (1/2)* (h1*(h2**2) + h2*(h1**2)) ) )

    

  # create banded array
  #jacob = np.array([upper, main, lower])

  # create tridiagonal matrix

  jacob = tridiag(lower, main, upper)



  return jacob
