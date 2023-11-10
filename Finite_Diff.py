import numpy as np

m = 5.68  # special boy units
h_bar = 0.6582  # eV * fs
c = 300 # nm/fs
width = 1  # nm
delta_x = 0.01  # nm
delta_t = 0.01  # fs
n = int((width/delta_x))  # number of steps


# The potential function V is defined for all time. (Taken to be constant)
# for now, so we will not need to worry about passing A any time information.
# 

# B requires the solutions from our previous matrix. The best bet seems to
# be saving the previous solution vector and simply tell B which "i" is
# required, and then just index out from the previous solutions.

def A(diagCount):
    i = diagCount

    Ai = -2 + ((4*i*m*delta_x**2)/(h_bar*delta_t)) - ((2*m*delta_x**2)/(h_bar)) # times potential function at n + 1

    return (Ai)
    
def Bi(solutions,row):
    i = row
    psi_iplus = solutions[i+1]
    psi_iminus = solutions[i-1]
    psi_i = solutions[i]

    B = -psi_iplus - psi_iminus + psi_i * (2+((4*i*m*delta_x**2)/(h_bar*delta_t))+((2*m*delta_x**2)/(h_bar**2))) # Times potential at previous time)
    return(B)



def array_populator(size,solutions):
    M = np.zeros((size, size))
    B = np.zeroes((size,))
    diagCounter = 0


    # i and j should start from 1

    for i in range(size):
        for j in range(size):

            if i == j:  # on the diagonal
                diagCounter += 1
                M[i,j] = A(diagCounter)

            elif i == (diagCounter + 1) and (j == diagCounter):
                M[i,j] = 1
        # assign value based on A_i
            elif (i == diagCounter) and j == (diagCounter + 1):
                M[i,j] = 1
        # assign value based on A_i




        # if i == 0:
         # = function of psi

            # ask about zero index
            # assign boundary value

        if (i > 0) and (i < n):
            B[i] = Bi(solutions,i)
            # assign values based on B_n

        elif i == n:
            pass
            # assign other boundary value

    # return matrix and vector

for t in range(n):  # time iteration loop

    pass

