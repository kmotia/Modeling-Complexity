import numpy as np
import matplotlib.pyplot as plt
import random

def prey_predator_hunter(alpha, beta, gamma, delta, steps, h):
    history = []
    
    Pi = 50
    Di = 100
    Hi = 40  

    P = Pi
    D = Di
    H = Hi

    for step in range(steps):
        # Every [kill_step] time steps
        # if step % 10 == 0:

        predators_killed = H * 0.001*D                         
        D -= predators_killed
        hunters_killed = H * 0.001
        H -= hunters_killed

        dP_dt = (alpha * P - beta * P * D)       # change in prey
        dD_dt = (delta * P * D - gamma * D)      # change in predators
        P_next = P + dP_dt * h
        D_next = D + dD_dt * h
        H_next = H  # Hunters don't reproduce or die naturally
        
        P = P_next  
        D = D_next
        H = H_next

        # Append the current state (P, D, H) to the history
        history.append((P, D, H))

    # Convert history to NumPy arrays
    P, D, H = map(np.array, zip(*history))
    return P, D, H

# Set parameters
alpha = 0.1                                     # Prey birth rate
beta_values = [0.002, 0.005, 0.008]             # Predation rate
gamma = 0.3                                     # Predator reproduction rate
delta = 0.01                                    # Predator death rate                                                
steps = 1000                                    # Number of Steps
h_values = [0.1, 0.2, 0.3]                      # Step Sizes

fig, ax = plt.subplots(3, 3, figsize=(12, 12))  
for i, beta in enumerate(beta_values):
    for j, h in enumerate(h_values):
        prey, predators, hunters = prey_predator_hunter(alpha, beta, gamma, delta, steps, h)
        ax[i, j].plot(prey, label="Prey")
        ax[i, j].plot(predators, label="Predators")
        ax[i, j].plot(hunters, label="Hunters")
        ax[i, j].set_title(f"h: {h} Beta: {beta}")
        ax[i, j].set_xlabel("Time Steps")
        ax[i, j].set_ylabel("Population")
        ax[i, j].legend()
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.4)
plt.savefig('Figures/hunter_plots.png')
plt.show()

