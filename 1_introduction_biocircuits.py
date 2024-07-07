import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the ODE models from the document

# Simple gene expression model
def simple_gene_expression(x, t, beta, gamma):
    dxdt = beta - gamma * x
    return dxdt

# Model with separate transcription and translation steps
def transcription_translation(y, t, beta_m, gamma_m, beta_p, gamma_p):
    m, p = y
    dmdt = beta_m - gamma_m * m
    dpdt = beta_p * m - gamma_p * p
    return [dmdt, dpdt]

# Parameters for the simple gene expression model
beta = 1.0
gamma = 0.5
x0 = 0  # initial condition
t = np.linspace(0, 50, 100)  # time points

# Solve the simple gene expression model
x = odeint(simple_gene_expression, x0, t, args=(beta, gamma))

# Parameters for the transcription and translation model
beta_m = 1.0
gamma_m = 0.5
beta_p = 1.0
gamma_p = 0.5
y0 = [0, 0]  # initial conditions for mRNA and protein

# Solve the transcription and translation model
y = odeint(transcription_translation, y0, t, args=(beta_m, gamma_m, beta_p, gamma_p))
m, p = y.T

# Plotting the results
plt.figure(figsize=(12, 6))

# Plot for the simple gene expression model
plt.subplot(1, 2, 1)
plt.plot(t, x, label='Protein concentration (x)')
plt.title('Simple Gene Expression Model')
plt.xlabel('Time')
plt.ylabel('Protein concentration')
plt.legend()

# Plot for the transcription and translation model
plt.subplot(1, 2, 2)
plt.plot(t, m, label='mRNA concentration (m)')
plt.plot(t, p, label='Protein concentration (p)')
plt.title('Transcription and Translation Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.tight_layout()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the ODE models for other topics from the document

# Repression model
def repression(y, t, beta0, alpha0, Kd, gamma):
    p, r = y
    dpdt = beta0 / (1 + (r / Kd)) + alpha0 - gamma * p
    drdt = 0  # assuming r is constant for simplicity
    return [dpdt, drdt]

# Activation model
def activation(y, t, beta0, Kd, gamma):
    p, a = y
    dpdt = beta0 * (a / Kd) / (1 + (a / Kd)) - gamma * p
    dadt = 0  # assuming a is constant for simplicity
    return [dpdt, dadt]

# Parameters for the repression model
beta0_repression = 1.0
alpha0_repression = 0.1
Kd_repression = 0.5
gamma_repression = 0.5
y0_repression = [0, 0]  # initial conditions for protein and repressor

# Parameters for the activation model
beta0_activation = 1.0
Kd_activation = 0.5
gamma_activation = 0.5
y0_activation = [0, 0]  # initial conditions for protein and activator

# Time points
t = np.linspace(0, 50, 100)

# Solve the repression model
y_repression = odeint(repression, y0_repression, t, args=(beta0_repression, alpha0_repression, Kd_repression, gamma_repression))
p_repression, r_repression = y_repression.T

# Solve the activation model
y_activation = odeint(activation, y0_activation, t, args=(beta0_activation, Kd_activation, gamma_activation))
p_activation, a_activation = y_activation.T

# Plotting the results
plt.figure(figsize=(12, 6))

# Plot for the repression model
plt.subplot(1, 2, 1)
plt.plot(t, p_repression, label='Protein concentration (p)')
plt.title('Repression Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

# Plot for the activation model
plt.subplot(1, 2, 2)
plt.plot(t, p_activation, label='Protein concentration (p)')
plt.title('Activation Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.tight_layout()
plt.show()


# Hill function models for ultrasensitivity

# Activating Hill function model
def hill_activation(y, t, beta0, Kd, n, gamma):
    p, a = y
    dpdt = beta0 * (a**n / (Kd**n + a**n)) - gamma * p
    dadt = 0  # assuming a is constant for simplicity
    return [dpdt, dadt]

# Repressing Hill function model
def hill_repression(y, t, beta0, Kd, n, gamma):
    p, r = y
    dpdt = beta0 / (1 + (r / Kd)**n) - gamma * p
    drdt = 0  # assuming r is constant for simplicity
    return [dpdt, drdt]

# Parameters for the Hill function models
beta0_hill = 1.0
Kd_hill = 0.5
n_hill = 2  # Hill coefficient for ultrasensitivity
gamma_hill = 0.5
y0_hill_activation = [0, 0]  # initial conditions for protein and activator
y0_hill_repression = [0, 0]  # initial conditions for protein and repressor

# Time points
t = np.linspace(0, 50, 100)

# Solve the Hill function activation model
y_hill_activation = odeint(hill_activation, y0_hill_activation, t, args=(beta0_hill, Kd_hill, n_hill, gamma_hill))
p_hill_activation, a_hill_activation = y_hill_activation.T

# Solve the Hill function repression model
y_hill_repression = odeint(hill_repression, y0_hill_repression, t, args=(beta0_hill, Kd_hill, n_hill, gamma_hill))
p_hill_repression, r_hill_repression = y_hill_repression.T

# Plotting the results
plt.figure(figsize=(12, 6))

# Plot for the Hill function activation model
plt.subplot(1, 2, 1)
plt.plot(t, p_hill_activation, label='Protein concentration (p)')
plt.title('Hill Function Activation Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

# Plot for the Hill function repression model
plt.subplot(1, 2, 2)
plt.plot(t, p_hill_repression, label='Protein concentration (p)')
plt.title('Hill Function Repression Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

plt.tight_layout()
plt.show()
