import numpy as np


X = np.array([[9,1],
             [5,3],
             [1,2]])
X
col1 = X[:,0] 
col2 = X[:,1] 


# https://matplotlib.org/tutorials/introductory/pyplot.html

# n points in p-dimensional space
# - each column of X is a dimension
# - variability occurs in more than one direction
# - sample mean is center of gravity
# - determinant of sample covmatrix provides measure of variability


import matplotlib.pyplot as plt
plt.scatter(col1,col2,label="observations")
plt.scatter(np.mean(col1), np.mean(col2), color="r", label="sample mean")
plt.xlabel('column 1')
plt.ylabel('column 2')
plt.suptitle("n points in p-dimensional space")
plt.title("Each row is an observation, each column a variable. n=3, p=2", size="small")
plt.legend()



# p VECTORS in n-dimensional space

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
#ax.autoscale(enable=False)
ax.set_xlim(0,9)
ax.set_ylim(0,6)
ax.set_zlim(0,6)
origin = [0], [0], [0]
# In Python3.5+, the asterisk unpacks a list
ax.quiver(*origin, X[0,:],X[1,:], X[2,:], arrow_length_ratio=0.1 )
plt.title("Quiver plot")

# as scatterplot: 
ax = fig.add_subplot(122, projection='3d')
ax.set_xlim(0,9)
ax.set_ylim(0,6)
ax.set_zlim(0,6)
ax.scatter(X[0,:],X[1,:], X[2,:])
plt.suptitle("p vectors in n-dimensional space. p=2, n=3")
plt.title("Scatter plot")


#####################

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
#ax.autoscale(enable=False)
ax.set_xlim(0,9)
ax.set_ylim(0,6)
ax.set_zlim(0,6)
origin = [0], [0], [0]
# In Python3.5+, the asterisk unpacks a list
ax.quiver(*origin, X[0,:],X[1,:], X[2,:], arrow_length_ratio=0.1 )
plt.title("Quiver plot")


xbar = np.mean(X, axis=0) 
nrow = X.shape[0] # in R: nrow(X2)
Xbar = np.outer(np.ones(nrow),xbar)
Z1 = X - Xbar # deviation matrix

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.autoscale(enable=False)
ax.set_xlim(0,9)
ax.set_ylim(0,6)
ax.set_zlim(0,6)
origin = [0], [0], [0]
# In Python3.5+, the asterisk unpacks a list
ax.quiver(*origin, X[0,:],X[1,:], X[2,:], arrow_length_ratio=0.1, label="x_i", color="b" )
ax.quiver(*origin, Z1[0,:],Z1[1,:], Z1[2,:], arrow_length_ratio=0.1 , label="deviation", color = "r")
ax.quiver(*origin, Xbar[0,:],Xbar[1,:], Xbar[2,:], arrow_length_ratio=0.1 , label="xbar", color = "g")
plt.title("Quiver plot")






#####################


X2 = np.array([[12,17,29],
               [18,20,38],
               [14,16,30],
               [20,18,38],
               [16,19,35]])

# X2 is a multivariate matrix with 5 observations with 3 dimensions(columns)

#axis = 0 means along the column and axis = 1 means working along the row.
xbar = np.mean(X2, axis=0) # same as colMeans(X2) in R 
type(xbar)

# visualize in 3d:

fig = plt.figure()
# Alternative form for add_subplot(111) is add_subplot(1, 1, 1)
# These are subplot grid parameters encoded as a single integer. For example, "111" means "1x1 grid
#(numrows, numcols, position)
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(121, projection='3d')
ax.scatter(xs= X2[:,0], ys = X2[:,1], zs = X2[:,2])
ax.scatter(xs=xbar[0], ys=xbar[1], zs = xbar[2] , color="r", label = "mean vector")
ax.legend()


nrow = X2.shape[0] # in R: nrow(X2)
# convert column means as n by p matrix:
Xbar = np.outer(np.ones(nrow),xbar) # dot product with 1d arrays in np
#see  https://stackoverflow.com/questions/23566515/multiplication-of-1d-arrays-in-numpy 
Xbar

# the mean corrected matrix / the deviation matrix:
Z = X2 - Xbar
Z

ax = fig.add_subplot(122, projection='3d')
ax.scatter(xs= Z[:,0], ys = Z[:,1], zs = Z[:,2])

# notice that col1+col2=col3
# therefore col1+col2 - col3 = 0
# ie columns are linearly dependent
a = np.array([1,1,-1])
Z@a # array([0., 0., 0., 0., 0.])



def my_sample_covmatrix(X):
    """
    Parameters
    ----------
    X : n by p matrix

    Returns
    -------
    Covariance matrix of sample (p by p)
    Unbiased (n-1).
    """
    
    n = X.shape[0]
    I = np.identity(n) 
    J = np.outer(np.ones(n),np.ones(n))
    S = (1/(n-1))*X.T@(I-(1/n)*J)@X
    return(S)

S = my_sample_covmatrix(X2)
S

np.linalg.det(S) # the generalized sample variance is basically zero 

S@a # = 0 , ie a can be rescaled to be an eigenvector corresponding to eigenvalue zero

np.linalg.eigvals(S)



