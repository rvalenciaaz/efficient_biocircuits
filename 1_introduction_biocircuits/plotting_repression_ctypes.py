import ctypes
import numpy as np
import matplotlib.pyplot as plt

# Load the shared library
lib = ctypes.CDLL('./repression_intervals_sundials_ctypes.so')

# Define the function signature
lib.solve_repression.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double]
lib.solve_repression.restype = None

# Parameters
n_steps = 500
dt = 0.1
results = np.zeros(n_steps, dtype=np.float64)

# Call the C function
lib.solve_repression(results.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_steps, dt)

# Plot the results
t = np.linspace(0, dt * (n_steps - 1), n_steps)
plt.plot(t, results)
plt.xlabel('Time')
plt.ylabel('Protein Concentration')
plt.title('Repression Model with Intervals')
plt.show()
