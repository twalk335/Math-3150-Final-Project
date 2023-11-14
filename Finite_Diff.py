# TODO
# Set the initial conditions properly, (0 on some infinite barrier or something)
# Set the i = 0 conditions properly as well, right now I believe it's also unconstrained.


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

m = 5.68  # special boy units
h_bar = 0.6582  # eV * fs
c = 300  # nm/fs
width = 1  # nm
delta_x = 0.01  # nm
delta_t = 0.0001  # fs
time = 1  # fs
n = int((width/delta_x))  # number of space steps
t = int((time/delta_t))  # number of time steps


# The potential function V is defined for all time. (Taken to be constant)
# for now, so we will not need to worry about passing A any time information.
#

# B requires the solutions from our previous matrix. The best bet seems to
# be saving the previous solution vector and simply tell B which "i" is
# required, and then just index out from the previous solutions.


def A(diagCount):

    i = np.imag(1j)  # fixed by adding complex numbers

    Ai = -2 + ((4*i*m*delta_x**2)/(h_bar*delta_t)) - ((2*m*delta_x**2)/(h_bar**2))  # times potential function at n + 1

    return Ai


def Bi(solutions, row):

    i = row  # 'i' is an index here, not an imaginary number

    if i != len(solutions):
        psi_iplus = solutions[i+1]
        psi_iminus = solutions[i-1]
        psi_i = solutions[i]
    else:
        psi_iplus = 0
        psi_iminus = solutions[i-1]
        psi_i = solutions[i]

    j = np.imag(1j)
    B = -psi_iplus - psi_iminus + psi_i * (2 + ((4*j*m*delta_x**2)/(h_bar*delta_t)) + ((2*m*delta_x**2)/(h_bar**2)))
    # Times potential at previous time

    return B


def array_populator(size, solutions):

    M = np.zeros((size, size))
    B = np.zeros(size)
    diagCounter = 0

    for i in range(size):
        for j in range(size):

            if i == j:  # on the diagonal
                diagCounter += 1
                M[i, j] = A(diagCounter)

            elif (i == diagCounter - 1) and (j == diagCounter):
                M[i, j] = 1

            elif (i == diagCounter) and j == (diagCounter - 1):
                M[i, j] = 1

        if (i > 0) and (i < size - 1):
            B[i] = Bi(solutions, i)

    return M, B


domain = np.linspace(0, width, n)

psi0 = np.exp(-(domain-0.5)**2)  # n space steps given, this value should be the size passed to the populator function
normalization_constant = (1/1.2533)  # value that normalizes wave function squared of initial conditions
psi0_pdf = psi0**2 * normalization_constant  # normalized probability density function for initial conditions

solution_matrix = np.empty((t, n))  # each new row corresponds to the wave function at a different time
solution_matrix[0, :] = psi0  # fixed to assign whole row to initial conditions

"""
plt.subplot(211)
plt.plot(size, normalization_constant * psi0_pdf)


M, B = array_populator(n, psi0)
solutions = np.linalg.solve(M, B)
# print(solutions)


plt.subplot(212)
plt.plot(size, np.absolute(solutions))

plt.tight_layout()
plt.show()

"""

for i in range(t-1):  # time iteration loop

    M, B = array_populator(n, solution_matrix[i, :])
    psi_n = np.linalg.solve(M, B)
    solution_matrix[i + 1, :] = psi_n

    print(f"Progress: {np.round(i/t * 100, 2)}%")


fig, ax = plt.subplots()


def update(frame):
    ax.clear()
    ax.plot(domain, np.absolute(solution_matrix[frame, :])**2)  # plotting magnitude of wave function squared
    ax.set_title(f'Time Step: {frame}')
    ax.set_xlabel('Space')
    ax.set_ylabel('Function Value')


# Create the animation
num_frames = solution_matrix.shape[0]
animation = FuncAnimation(fig, update, frames=num_frames, interval=100)

# Display the animation
plt.show()


