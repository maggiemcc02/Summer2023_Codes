
# In this file I have the codes I use to deflate


# In particular I have

# interpolation matrix code
# eta abd eta' code
# smiley face code
# grumpy face code



# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d






# Interpolation matrix




def interp_matrix(A, B):

  """
  Forms interpolation matrix from mesh A and mesh B to mesh C = AUB

  U_A : Known function values on A

  A : Mesh With Known Values

  B : Mesh with Unknown values 

  """


  matrix_list = []



  for i in range(len(A)-1): # for each A

    # add an identity row
    lst = [0 for i in range(len(A))]
    lst[i] = 1
    matrix_list.append(lst)



    # set values if we need to interpolate 

    x0 = A[i][0]
    x1 = A[i+1][0]



    counter = 0 # count how many B values are in the current interval
  
    for j in range(len(B)): # for each B

      xj = B[j][0]

      if (xj >= x1) :
        # move on to next interval
        break

      if ((xj > x0) and (xj < x1)):

        # B = xj  is in the interval so Interpolate !


        lst = [ 0 for i in range(len(A))]

        lst[i] = (x1 - xj) / (x1 - x0)

        lst[i+1] = (xj - x0) / (x1 - x0)

        matrix_list.append(lst)


  # add last row 

  lst = [0 for i in range(len(A))]
  lst[-1] = 1
  matrix_list.append(lst)


  # make the list a matrix

  result  = np.array(matrix_list) 

  return result








# eta







def eta(U_A, r, A, B, power): # GRID 2 NORM 

  """
  U_A : Approximation of solution on A

  r : Approximation of solution on B

  A : Mesh related to U_A

  B : Mesh related to r

  """


  # compute P from A to C

  P_AC = interp_matrix(A, B)


  # Multiply P_AC by U_A to interpolate U_A onto C = AUB

  
  U_C = np.matmul(P_AC, U_A)

  # compute P from B to C

  P_BC = interp_matrix( B, A)

  # Multiply P_BC by r to interpolate r onto C = AUB

  r_C = np.matmul(P_BC, r)

  # Create Combined Mesh C

  C_list = []
  for i in range(len(A)):
     C_list.append(A[i][0])
     C_list.append(B[i][0])
  C_set = set(C_list)
  C_list2 = list(C_set)
  C = sorted(C_list2)

  C = np.array([ C ]).T

  # Compute mesh spacing of C

  space = [0 for i in range(len(C))]

  for i in range(1, len(C)):

    space[i] = ( C[i][0] - C[i-1][0] )


  # compute P_AC * x_A - P_BC * r_B

  val = U_C - r_C

  # construct the vector whose ith entry is |sqrt(h_i) * xi |

  norm_vect = []

  for i in range(len(C)):

    h = space[i]
    xi = val[i][0]

    norm_vect.append(abs((np.sqrt(h))*xi))

  norm_vect = np.array( [ norm_vect ] ).T


  # take numpys norm

  result = norm(norm_vect)

  result = result ** power

  return(result)




# eta'







def eta_prime(U_A, r, A, B, power): # BASED OF GRID 2 NORM 

  """
  U_A : Approximation of solution on A

  r : Approximation of solution on B

  A : Mesh related to U_A

  B : Mesh related to r

  P_AC : interpolation matrix from A to C

  P_BC : interpolation matrix from B to C

  """
  # compute P from A to C

  P_AC = interp_matrix( A, B) # could eliminate this step and pass matrix in

  # Multiply P_AC by U_A to interpolate U_A onto C = AUB

  U_C = np.matmul(P_AC, U_A)

  # compute P from B to C

  P_BC = interp_matrix( B, A) # could eliminate step and pass matrix in

  # Multiply P_BC by r to interpolate r onto C = AUB

  r_C = np.matmul(P_BC, r)

  # Create Combined Mesh C

  C_list = []
  for i in range(len(A)):
    C_list.append(A[i][0])
    C_list.append(B[i][0])
  C_set = set(C_list)
  C_list2 = list(C_set)
  C = sorted(C_list2)

  C = np.array([ C ]).T

  # Compute mesh spacing of C

  space = [0 for i in range(len(C))] # letting first h be 0

  for i in range(1, len(C)):

    space[i] = ( C[i][0] - C[i-1][0] )


  # compute the Z = P_AC*x_a - P_Bc * r_B

  Z = U_C - r_C

  # construct the vector whose ith entry is |sqrt(h_i) * xi |

  norm_vect = []

  for i in range(len(C)):

    h = space[i]
    xi = Z[i][0]

    norm_vect.append(abs((np.sqrt(h))*xi))

  norm_vect = np.array( [ norm_vect ] ).T

  # compute norm : || Uc - rc || 

  euclid = norm(norm_vect)

  # compute transpose of P(A to C)

  P_AC_T = np.transpose(P_AC)

  # multiply space by Z
  
  hZ = []

  for i in range(len(C)):

    h = space[i]
    xi = Z[i][0]

    hZ.append(h*xi)  

  hZ = np.array([hZ]).T

  # compute outside derivative for chain rule : power/2 ( || ||^2 )^(power/2 - 1 )

  a = power * ( euclid ** 2 ) ** ( (power/2) - 1 )

  # compute inside derivative for chain rule : derivative of || ||^2 = 2(P_A^C)^T(sub) (2 cancels with 2 in 'a')

  b = np.matmul(P_AC_T, hZ)

  # result is a * b

  result = a * b

  return result








# SMILEY FACE code







def smiley_face(U, mesh, mesh_spacing, u_0, u_n, solutions, sol_mesh, alpha, power):


  """
   U : current approximation in physical solver

   mesh : current mesh from de Boor

   mesh_spacing : mesh width of mesh

   u_0, u_n : boundary conditions

   solutions : vector of previously found solutions

   sol_mesh : corresponding meshes for the found solutions

   alpha, power : deflation parameters

  """

  # if we have no solutions we just solve the original system
  
  if len(solutions) == 0 :

    M = 1
    eta_i = []
    full_U = []
    full_mesh = []

    return [M, [eta_i, full_U, full_mesh]]
   
  # consider the boundaries as well as the approximations for U and the mesh
    
  full_U = [u_0]
  full_mesh = [0]

  for i in range(len(U)):
    full_U.append(U[i][0])
    full_mesh.append(mesh[i][0])
  full_U.append(u_n)
  full_mesh.append(1)

  full_U = np.array([full_U]).T
  full_mesh = np.array([full_mesh]).T


  # compute the eta's for each of the solutions found

  eta_i = []

  for i in range(len(solutions)):

    n_i = eta(full_U, solutions[i], full_mesh, sol_mesh[i], power)

    eta_i.append(n_i)


  # construct product of M's

  M = 1

  for i in range(len(eta_i)):

    n = eta_i[i]

    M = M * ( ( 1 / (n) ) + alpha )


  # we return this M value 

  return [M, [eta_i, full_U, full_mesh]]











# GRUMPY FACE code








def grumpy_face( U, mesh, mesh_spacing, smiley_face_results, solutions,\
                sol_mesh, alpha, power):
  

  """
   U : current approximation in physical solver

   mesh : current mesh from de Boor

   mesh_spacing : mesh width of mesh

   u_0, u_n : boundary conditions

   solutions : vector of previously found solutions

   sol_mesh : corresponding meshes for the found solutions

   alpha, power : deflation parameters

  """

  # if we have no solutions we just use the original jacobian 

  if len(solutions) == 0 : 

    return np.zeros((1, (len(U))))
  

  # from smiley face results pull out eta_i, full_U, full_mesh

  eta_i = smiley_face_results[0]

  full_U = smiley_face_results[1]

  full_mesh = smiley_face_results[2]


  # construct the eta primes for each solution given

  eta_pr = []

  for i in range(len(eta_i)):


    d = eta_prime(full_U, solutions[i], full_mesh, sol_mesh[i], power)

    # take endpoints out

    d2 = []
    for j in range(len(d)):
      if j != 0:
        if j != (len(d)-1):
          d2.append(d[j][0])
    d = np.array([d2]).T

    # save to list 

    eta_pr.append(d)

  # construct the weird sum

  my_sum = 0


  for i in range(len(eta_i)): # each term of sum

    prod = 1
    

    for j in range(len(eta_i)): # product for ith term of sum


      if i != j : 


        eta_j = eta_i[j]

        prod *= (( 1 / (eta_j) ) + alpha) 

    # transpose of ith eta prime

    d_T = np.transpose(eta_pr[i])

    # create that derivative 

    deriv = (-1 / (eta_i[i]**(2))) * d_T

    # multiply it by prod to get the term for the sum

    term = prod * deriv


    # add it to the sum 

    my_sum  = my_sum +  term
    

  # return that sum as grumpy face 

  return my_sum








