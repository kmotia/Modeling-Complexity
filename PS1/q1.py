import numpy as np
import matplotlib.pyplot as plt

def lotka_volterra_phase_space(K1, K2, alpha1, alpha2):

    # define largest possible N1 and N2 values
    largest_N1_value = max(K1, K2/alpha1) 
    largest_N2_value = max(K2, K1/alpha2) 

    # define values of N1 and N2
    all_n1, all_n2 = np.meshgrid(np.arange(0, largest_N1_value, 1), np.arange(0, largest_N2_value, 1))

    # define time-derivatives of N1 and N2
    r1 = 1
    r2 = 1
    dn1_dt = r1 * all_n1 * (K1 - all_n1 - alpha2 * all_n2) / K1
    dn2_dt = r2 * all_n2 * (K2 - all_n2 - alpha1 * all_n1) / K2
 
    # define fixed points
    n1_fixed = (K1 - alpha2 * K2) / (1 - alpha1 * alpha2)
    n2_fixed = (K2 - alpha1 * K1) / (1 - alpha1 * alpha2)

    # define equilibrium points dependent on case
    if K1/alpha2 >= K2 and K1 >= K2/alpha1:
        equilibrium_points = [[0, K1]]
        title = '$K_1/\\alpha_2 \geq K_2 \quad K_1 \geq K_2/\\alpha_1$'

    if K2 >= K1/alpha2 and K2/alpha1 >= K1:
        equilibrium_points = [[0, K2]]
        title = '$K_2 \geq K_1/\\alpha_2 \quad K_2/\\alpha_1 \geq K_1$'

    if K1/alpha2 >= K2 and K2/alpha1 >= K1:
        equilibrium_points = [[n1_fixed, n2_fixed]]
        title = '$K_1/\\alpha2 \geq K_2 \quad K_2/\\alpha_1 \geq K_1$'

    if K2 >= K1/alpha2 and K1 >= K2/alpha1:
        equilibrium_points = [[K1, 0], [0, K2], [n1_fixed, n2_fixed]] 
        title = '$K_2 \geq K_1 /\\alpha_2 \quad K_1 \geq K_2 /\\alpha_1$'

    # plot figure
    plt.figure(figsize=(6, 6))

    # plot nullclines
    plt.plot([0, K2/alpha1], [K2, 0], 'r', label='$N_2$ Isocline') 
    plt.plot([0, K1], [K1/alpha2, 0], 'b', label='$N_1$ Isocline')  

    # plot the perfect diagonal
    x_values = np.linspace(0, largest_N1_value, 2) 
    y_values = ((K1**2)*alpha1)/((K2**2)*alpha2) * x_values  
    plt.plot(x_values, y_values, 'purple', label='$N_2= \\frac{K_1^2\\alpha_1}{K_2^2\\alpha_2} N_1$', linewidth=2.0)  

    # plot equilibrium points
    equilibrium_points = np.array(equilibrium_points)
    plt.scatter(equilibrium_points[:, 0], equilibrium_points[:, 1], color='black', s=100, label='Equilibrium Points')

    # plot vector field 
    plt.streamplot(all_n1, all_n2, dn1_dt, dn2_dt, color='gray', linewidth=1, arrowsize=1.5, density=1.5, cmap='viridis')

    # annotate teh x and y intercepts
    plt.annotate(f'$K_2/\\alpha_1$', xy=(K2/alpha1, 0), xytext=(0, -30), textcoords='offset points', fontsize=10) 
    plt.annotate(f'$K_2$', xy=(0, K2), xytext=(-33, 0), textcoords='offset points', fontsize=10)
    plt.annotate(f'$K_1/\\alpha_2$', xy=(0, K1/alpha2), xytext=(-43, 0), textcoords='offset points', fontsize=10)
    plt.annotate(f'$K_1$', xy=(K1, 0), xytext=(0, -30), textcoords='offset points', fontsize=10) 

    # figure details
    plt.xlabel("$N_1$")
    plt.ylabel("$N_2$")
    plt.title(f'Lotka-Volterra Model \n {title}')
    plt.legend(loc='upper right')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    # show figure
    plt.grid()
    plt.savefig('Figures/phase_space.png')
    plt.show()



# K1, K2, alpha1, alpha2
# Question 1
lotka_volterra_phase_space(100, 100, 3, 3)  # case4


# lotka_volterra_phase_space(100, 200, 1, 2)  # case2