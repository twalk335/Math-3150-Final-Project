import numpy as np

m = 5.68  # special boy units
h_bar = 0.06582  # eV * fs

width = 1  # nm
delta_x = 0.01  # nm
delta_t = 0.01  # fs
n = int((width/delta_x))  # number of steps


def array_populator(size):
    np.zeros((size, size))
    diagCounter = 0

    for i in range(size):
        for j in range(size):

            if i == j:  # on the diagonal
                diagCounter += 1

            elif i == (diagCounter + 1) and (j == diagCounter):
        # assign value based on A_i
            elif (i == diagCounter) and j == (diagCounter + 1):
        # assign value based on A_i

        if i == 1:
            # assign boundary value
        elif (i > 1) and (i < n):
            # assign values based on B_n
        elif i == n:
            # assign other boundary value

    # return matrix and vector

for t in range(n):  # time iteration loop
