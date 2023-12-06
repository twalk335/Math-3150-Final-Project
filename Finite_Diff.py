import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

i = np.imag(1j)
hbar = 1
m = 1
x = np.linspace(0, 1, 100)
dx = x[1]-x[0]
t = np.linspace(0, 1, 1000)
dt = t[1]-t[0]
alpha = (i*dt)/(2*hbar)
beta = hbar/(2*m)
A = (alpha * beta) / dx**2
mode = 5  # number of anti-nodes in psi0
print(f"Stability Ratio: {dt/dx**2}")

# Terms used in matrix algebra
Kt = lambda x, ti: (1 + 2*A + Vt(ti)[x] * alpha) / A
Gt = lambda x, ti: (1 - 2*A - Vt(ti)[x] * alpha) / A


### POTENTIAL FUNCTIONS
# Vt = lambda ti: 1e2 * np.exp(-((x-1+2*ti*1/(len(t)))**2)*100)  # moving gaussian
Vt = lambda ti: 1e2 * (x - ti/len(t))**2  # moving harmonic oscillator
# Vt = lambda ti: 1e2 * (x - 0.5)**2 + ti*0  # static harmonic oscillator
# Vt = lambda ti: 0 * ti * x  # infinite square well
###

### INITIAL CONDITIONS
L = x[-1] - x[0]  # length of box
psi0 = (np.sqrt(2/L) * np.sin(mode * np.pi/L * x))  # standing wave initial conditions
psi0 = psi0/np.sqrt(np.sum(psi0**2 * dx))  # normalizing
###


def analytical_solution(x_input, t_input):  # analytical solution for infintie square well under given ICs

    sin_part = np.sqrt(2/L) * np.sin(mode * np.pi/L * x_input)
    exp_part = np.exp(((-1j * mode**2 * np.pi**2 * hbar)/(2 * m * L**2)) * t_input)
    return sin_part * exp_part


plt.plot(x, psi0**2, label='psi0 for approximation')
plt.title(f"psi0**2 area under graph: {np.round(np.trapz(psi0**2, x), 3)}")
plt.plot(x, analytical_solution(x, t[0])**2, 'r.', label='exact psi0')
plt.legend()
plt.show()

true_sol_mat = np.empty((len(t), len(x)), dtype=np.complex128)
for k in range(len(t)):
    true_sol_mat[k, :] = analytical_solution(x, t[k])


def MBarrs(solutions, t=None):
    diag_count = 0

    M = np.zeros((len(x), len(x)))
    B = np.zeros((len(x),))
    for i in range(len(x)):
        for j in range(len(x)):

            if i == j:
                # M[i,j] = K (i)
                M[i, j] = Kt(i, t)
                diag_count += 1

            elif (i == diag_count - 1) and (j == diag_count):
                M[i, j] = -1

            elif (i == diag_count) and j == (diag_count - 1):
                M[i, j] = -1
        
        if 1 < i < len(x) - 1:
            # B[i] = solutions[i] * G(i) + solutions[i+1] + solutions[i-1]
            B[i] = solutions[i] * Gt(i, t) + solutions[i+1] + solutions[i-1]
        elif i == 0:
            B[i] = 0
        elif i == len(x):
            B[i] = 0

    return M, B


M1, B1 = MBarrs(psi0, 0)

V_overtime = np.empty((len(t), len(x)))

sols = np.linalg.solve(M1, B1)

solution_matrix = np.empty((len(t), len(x)))
solution_matrix[0, :] = psi0
solution_matrix[1, :] = sols

for i in range(2, len(t)):
    if np.mod(i, 40) == 0:
        print(f"Progress: {np.round(i/len(t) * 100, 3)}%")
    Mx, Bx = MBarrs(sols, i)
    sols = np.linalg.solve(Mx, Bx)
    sols = sols / np.sqrt(np.sum(sols**2 * dx))
    solution_matrix[i] = sols
    V_overtime[i] = Vt(i)

fig, ax = plt.subplots()


def update(frame):
    ax.clear()
    # plotting magnitude of wave function squared
    ax.plot(x, np.absolute(solution_matrix[frame, :])**2)
    # ax.plot(x, np.absolute(true_sol_mat[frame, :])**2, 'r.')  # analytical solution for standing wave
    ax.plot(x, V_overtime[frame, :] * 1/5, 'r-')
    area = np.trapz(np.absolute(solution_matrix[frame, :])**2, x)
    ax.set_title(f'Area Under Graph: {np.round(area, 4)}, time-step: {frame}')
    ax.set_xlabel('Distance')
    ax.set_ylabel(r"$|\Psi(x,t)|^2$")
    ax.set_ylim(0, 6)


# Create the animation
num_frames = solution_matrix.shape[0]
animation = FuncAnimation(fig, update, frames=num_frames, interval=30)
# animation.save("test.gif", writer="Pillow")
animation.save("test.mp4", writer="ffmpeg", fps=60)

