# TODO
# Set the initial conditions properly, (0 on some infinite barrier or something)
# Set the i = 0 conditions properly as well, right now I believe it's also unconstrained.


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation

m = 5.68  # special boy units
h_bar = 0.6582  # eV * fs
c = 300  # nm/fs
width = 1  # nm
delta_x = 0.005  # nm
delta_t = 0.1  # fs
time = 1  # fs
n = int((width/delta_x))  # number of space steps
t = int((time/delta_t))  # number of time steps

# The potential function V is defined for all time. (Taken to be constant)
# for now, so we will not need to worry about passing A any time information.
#

# B requires the solutions from our previous matrix. The best bet seems to
# be saving the previous solution vector and simply tell B which "i" is
# required, and then just index out from the previous solutions.


def A(diagCount, V):

    i = np.imag(1j)  # fixed by adding complex numbers

    Ai = -2 + ((4*i*m*delta_x**2)/(h_bar*delta_t)) - ((2*m*delta_x**2)/(h_bar**2)) * V[diagCount - 1]

    return Ai


def Bi(solutions, row, V):

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
    B = -psi_iplus - psi_iminus + psi_i * (2 + ((4*j*m*delta_x**2)/(h_bar*delta_t)) + ((2*m*delta_x**2)/(h_bar**2)) * V[i])
    # Times potential at previous time

    return B


def array_populator(size, solutions, V):

    M = np.zeros((size, size))
    B = np.zeros(size)
    diagCounter = 0

    for i in range(size):
        for j in range(size):

            if i == j:  # on the diagonal
                diagCounter += 1
                M[i, j] = A(diagCounter, V)

            elif (i == diagCounter - 1) and (j == diagCounter):
                M[i, j] = 1

            elif (i == diagCounter) and j == (diagCounter - 1):
                M[i, j] = 1

        if (i > 0) and (i < size - 1):
            B[i] = Bi(solutions, i, V)

    return M, B


"""
def initial_conditions(x, x0, sigma):
    i = np.imag(1j)
    p0 = 0.00001
    A = 1/(np.pi**(1/4) * np.sqrt(sigma))
    B = -(((x - x0)**2)/(2 * sigma**2))
    C = (i * p0 * x)/h_bar

    output = A * np.exp(B) * np.exp(C)
    return output
"""

"""
def initial_conditions(x, x0, sigma):
    output = np.exp(-((x - x0)**2)/(2 * sigma**2))
    normalization_constant = 1/(sigma * np.sqrt(2 * np.pi))
    return output * normalization_constant
"""


def initial_conditions(x, mode):
    return np.sqrt(2) * np.sin(mode * np.pi * x)


def potential(x, mu, sigma):
    V = 700 * np.exp((-(x - mu)**2)/(2 * sigma**2))
    return V


"""
def analytical_solution(x, t, x0, sigma):
    i = np.imag(1j)
    p0 = 20
    E = m * c**2
    A = 1/(np.pi**(1/4) * np.sqrt(sigma * (1 + (i * h_bar * t)/(m * sigma**2))))
    B = (-(x - (x0 + (p0 * t)/m)**2))
    C = 2 * sigma**2 * (1 + (i * h_bar * t)/(m * sigma**2))
    D = i * (p0 * x - E * t)/h_bar

    output = A * np.exp(B/C) * np.exp(D)
    return output
"""

domain_x = np.linspace(0, width, n)
domain_t = np.linspace(0, time, t)

psi0 = initial_conditions(domain_x, 3)

solution_matrix = np.empty((t, n))  # each new row corresponds to the wave function at a different time
solution_matrix[0, :] = psi0  # fixed to assign whole row to initial conditions
V_x = potential(domain_x, 0.5, 0.1)

# plt.plot(domain_x, psi0)
# plt.plot(domain_x, V_x)
# plt.show()
for i in range(t-1):  # time iteration loop

    M, B = array_populator(n, solution_matrix[i, :], V_x)
    psi_n = np.linalg.solve(M, B)
    psi_n = psi_n/np.sqrt(np.sum(psi_n**2) * delta_x)
    solution_matrix[i + 1, :] = psi_n

    print(f"Progress: {np.round(i/t * 100, 2)}%")


true_solution_matrix = np.zeros((t, n))
# true_solution_matrix[0, :] = initial_conditions(domain_x, -0.5, 0.01)  # x, x0, sigma

"""
for i in range(t - 1):
    true_solution_matrix[i + 1, :] = analytical_solution(domain_x, domain_t[i], -0.5, 0.1)
    print(f"2nd Progress: {np.round(i/t * 100, 2)}%")
"""

fig, ax = plt.subplots()


def update(frame):
    ax.clear()
    ax.plot(domain_x, np.absolute(solution_matrix[frame, :])**2)  # plotting magnitude of wave function squared
    # ax.plot(domain_x, np.absolute(true_solution_matrix[frame, :])**2)
    area = np.trapz(np.absolute(solution_matrix[frame, :])**2, domain_x)
    ax.set_title(f'Area: {area}')
    ax.set_xlabel('Space')
    ax.set_ylabel(r"$|\psi(x,t)|^2$")


# Create the animation
num_frames = solution_matrix.shape[0]
animation = FuncAnimation(fig, update, frames=num_frames, interval=100)

# Display the animation
plt.show()
