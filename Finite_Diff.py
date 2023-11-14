# TODO
# Set the initial conditions properly, (0 on some infinite barrier or something)
# Set the i = 0 conditions properly as well, right now I believe it's also unconstrained.




import numpy as np
import matplotlib.pyplot as plt


m = 5.68  # special boy units
h_bar = 0.6582  # eV * fs
c = 300 # nm/fs
width = 1  # nm
delta_x = 0.01  # nm
delta_t = 0.01  # fs
n = int((width/delta_x))  # number of time steps


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

    if i != len(solutions):
        psi_iplus = solutions[i+1]
        psi_iminus = solutions[i-1]
        psi_i = solutions[i]
    else:
        psi_iplus = 0
        psi_iminus = solutions[i-1]
        psi_i = solutions[i]

    B = -psi_iplus - psi_iminus + psi_i * (2+((4*i*m*delta_x**2)/(h_bar*delta_t))+((2*m*delta_x**2)/(h_bar**2))) # Times potential at previous time)
    return(B)



def array_populator(size,solutions):
    M = np.zeros((size, size))
    B = np.zeros(size)
    diagCounter = 0


    # i and j should start from 1

    for i in range(size):
        for j in range(size):

            if i == j:  # on the diagonal
                diagCounter += 1
                M[i,j] = A(diagCounter)

            elif (i == diagCounter - 1) and (j == diagCounter): # fixed iterator in these steps, seemed built for 1 index
                M[i,j] = 1
        # assign value based on A_i
            elif (i == diagCounter) and j == (diagCounter - 1):
                M[i,j] = 1
        # assign value based on A_i




        # if i == 0:
         # = function of psi

            # ask about zero index
            # assign boundary value

        if (i > 0) and (i < size - 1):
            B[i] = Bi(solutions,i)
            # assign values based on B_n
        
        # elif i == n:
        #     pass
        #     # assign other boundary value

    return (M,B)




num_space = 50
size = np.linspace(0,1,num_space)

psi0 = np.exp((-size-1.5)**2) # 100 space steps given, this value should be the size passed to the populator function
initial_condition = psi0

plt.subplot(211)
plt.plot(size,initial_condition)

full_solution_vector = np.empty((n,num_space))
full_solution_vector[0] = initial_condition

M, B = array_populator(num_space,initial_condition)
solutions = np.linalg.solve(M,B)
print(solutions)


plt.subplot(212)
plt.plot(size,solutions)

plt.tight_layout()
plt.show()

for t in range(n):  # time iteration loop
    
    pass

