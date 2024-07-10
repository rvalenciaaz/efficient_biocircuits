import ctypes
import numpy as np
import matplotlib.pyplot as plt

# Load the shared library
lib = ctypes.CDLL('./dichotomous_feedback_sundials.so')

# Define the function signature
lib.solve_dichotomous_feedback.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_double]
lib.solve_dichotomous_feedback.restype = None

# Parameters
n_steps = 500
dt = 0.1
I = 1.0
results = np.zeros((n_steps, 8), dtype=np.float64)

# Call the C function
lib.solve_dichotomous_feedback(results.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_steps, dt, I)

# Plot the results
t = np.linspace(0, dt * (n_steps - 1), n_steps)
plt.figure(figsize=(12, 8))

# Plot each variable
labels = ["HK", "HKp", "RR", "RRp", "SR", "SRp", "PH", "Output"]
for i in range(8):
    plt.plot(t, results[:, i], label=labels[i])

plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Dichotomous Feedback Model')
plt.legend()
plt.show()
