


# In this file I have the optimal arclength mesh density function



from numerical_integration import *



# mesh density function for interior nodes 

def interior_u_x(Uc, Ua, xc, xb, xa):

  """
  Uc : U_i+1 
  Ua : U_i-1
  xc, xb, xa : corresponding grid points
  """

  U0 = Ua # U_i-1
  U2 = Uc # U_i+1

  h1 = xb - xa
  h2 = xc - xb

  ux = (U2 - U0) / (h2 + h1) # centered approx of u'

  # approx = np.sqrt(1 + (ux)**2)

  return ux







# mesh density function for endpoints 

def endpoint_u_x(Ub, Ua, xb, xa):
  """
  Ub : U_i+1 or Ui
  Ua : Ui or U_i-1
  xb, xa : corresponding grid points
  """
  ux = (Ub - Ua) / (xb - xa) # one sided approx of u' for endpoints


  # approx = np.sqrt( 1 + (ux)**2)

  return ux

# code to generate M values on a given grid









def M_calc_optimal_arc(U, grid):

  """
  U : Vector of Approximations

  grid : mesh

  """


  u_pr = [] # empty list to save M values

  u_pr.append(endpoint_u_x(U[1][0], U[0][0], grid[1][0], grid[0][0])) # add M value at first grid point using forward approx

  for i in range(len(grid)-2): # add M values for interior nodes


    # set values of U and y to pass to mesh density function

    U0 = U[i][0]
    U2 = U[i+2][0]
    x0 = grid[i][0]
    x1 = grid[i+1][0]
    x2 = grid[i+2][0]

    val = interior_u_x(U2, U0, x2, x1, x0) # compute mesh density
    
    u_pr.append(val)


  u_pr.append(endpoint_u_x(U[-1][0], U[-2][0], grid[-1][0], grid[-2][0])) # add M value for last grid point

  u_pr = np.array([u_pr]).T

  # compute the integral 

  u_pr_power = abs(u_pr)**(2/3)

  # int_u_pr_power = trap_nonuni(u_pr_power, grid)
  
  int_u_pr_power = int_constant(u_pr_power, grid)

  alpha_int = ((1 / (grid[-1][0] - grid[0][0])) * int_u_pr_power )**3

  # result 

  rho = ( 1 + (1/alpha_int) * ((abs(u_pr))**2) )**(1/3)

  return rho
