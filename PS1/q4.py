import numpy as np
import random
import matplotlib.pyplot as plt


# Discrete SIS calculation using Euler's method
def run_mathematical_SIS_euler(Si, Ii, beta, gamma, steps, h): 

    # Initialize populations 
    S = Si
    I = Ii

    # Generate result for each time step of size h
    history = []
    for step in range(steps):

        # Calculate the change in the number of infected and susceptible individuals using Euler's method
        dI_dt = ((beta * S * I) - (gamma * I))
        dS_dt = -((beta * S * I) - (gamma * I))

        next_I = I + dI_dt * h    
        next_S = S + dS_dt * h     

        # Update the number of infected and susceptible individuals
        I = next_I 
        S = next_S

        # Append the current state (S, I) to the history
        history.append((S, I))

    # Convert history to NumPy arrays for easy analysis and plotting
    St, It = map(np.array, zip(*history))
    return St, It

# Discrete SIS calculation using Heun's method
def run_mathematical_SIS_heun(Si, Ii, beta, gamma, steps, h):

    # Initialize populations 
    S = Si
    I = Ii

    # Generate result for each time step of size h
    history = []
    for step in range(steps):

  
        dI_dt = ((beta * S * I) - (gamma * I))
        dS_dt = -((beta * S * I) - (gamma * I))
        
        I_e = I + dI_dt * h    
        S_e = S + dS_dt * h    

        dI_dt_e = ((beta * S_e * I_e) - (gamma * I_e))
        dS_dt_e = -((beta * S_e * I_e) - (gamma * I_e))

        next_I = I + h * (dI_dt + dI_dt_e) / 2
        next_S = S + h * (dS_dt + dS_dt_e) / 2

        # Update the number of infected and susceptible individuals
        I = next_I 
        S = next_S

        # Append the current state (S, I) to the history
        history.append((S, I))

    # Convert history to NumPy arrays for easy analysis and plotting
    St, It = map(np.array, zip(*history))
    return St, It


# Define hyperparameters
Si = 90
Ii = 10
beta_values = [0.03, 0.06, 0.1]
gamma = 0.25
steps = 50
h_values = [0.01, 0.5, 2.0] # step sizes



fig, ax = plt.subplots(3, 3, figsize=(12, 12))  
for i, beta in enumerate(beta_values):
    for j, h in enumerate(h_values):
        ax[i, j].plot(run_mathematical_SIS_euler(Si, Ii, beta, gamma, steps, h)[0], label="Susceptible")
        ax[i, j].plot(run_mathematical_SIS_euler(Si, Ii, beta, gamma, steps, h)[1], label="Infected")
        ax[i, j].set_title(f"h: {h}, Beta: {beta}", fontsize=10)
        ax[i, j].set_xlabel("Time Steps")
        ax[i, j].set_ylabel("Population")
        ax[i, j].legend()
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
plt.savefig('Figures/euler_plots.png')
plt.show()


fig, ax = plt.subplots(3, 3, figsize=(12, 12))  
for i, beta in enumerate(beta_values):
    for j, h in enumerate(h_values):
        ax[i, j].plot(run_mathematical_SIS_heun(Si, Ii, beta, gamma, steps, h)[0], label="Susceptible")
        ax[i, j].plot(run_mathematical_SIS_heun(Si, Ii, beta, gamma, steps, h)[1], label="Infected")
        ax[i, j].set_title(f"h: {h}, Beta: {beta}", fontsize=10)
        ax[i, j].set_xlabel("Time Steps")
        ax[i, j].set_ylabel("Population")
        ax[i, j].legend()
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
plt.savefig('Figures/heun_plots.png')
plt.show()


