import matplotlib.pyplot as plt
import numpy as np

def plot_potential(potential_info, image_file):
    potential = potential_info[0]
    plt.figure(1, figsize=(8, 8))
    plt.imshow(potential, cmap="gray")
    plt.axis("off")
    # plt.colorbar()
    plt.title("Plot of potential", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"./plots/{image_file}_plots/potential_plot.png")


def plot_eigenfunctions(eigenvectors, potential_info, max_level, image_file):
    N = potential_info[1]
    edges = potential_info[2]
    for state in range(max_level):
        plt.figure(state, figsize=(8, 8))
        plt.imshow(edges, cmap="gray", alpha = 0.3)
        eig = plt.imshow(eigenvectors.T[state].reshape((N,N)), cmap="viridis", alpha=0.7)
        plt.title(f"State {state} eigenfunction", fontsize=16)
        plt.axis("off")
        # plt.colorbar(eig)
        plt.tight_layout()
        plt.savefig(f"./plots/{image_file}_plots/{state}/{state}_eigenfunction.png")

def plot_prob_densities(eigenvectors, potential_info, max_level, image_file):
    N = potential_info[1]
    edges = potential_info[2]
    for state in range(max_level):
        plt.figure(state, figsize=(8, 8))
        plt.imshow(edges, cmap="gray", alpha = 0.3)
        prob_density = (np.abs(eigenvectors.T[state].reshape((N,N))))**2
        eig = plt.imshow(prob_density, cmap="viridis", alpha=0.7)
        plt.title(f"State {state} probability density", fontsize=16)
        plt.axis("off")
        # plt.colorbar(eig)
        plt.tight_layout()
        plt.savefig(f"./plots/{image_file}_plots/{state}/{state}_prob_density.png")

def plot_eigenfunction_zero_crossings(eigenvectors, potential_info, max_level, image_file):
    N = potential_info[1]
    edges = potential_info[2]
    for state in range(max_level):
        plt.figure(state, figsize=(8, 8))
        print(f"state {state} zero's calculating")
        eigenfunction = eigenvectors.T[state].reshape((N, N))
        
        zero_crossings = np.isclose(eigenfunction, 0, atol=1e-5)
        plt.imshow(edges, cmap="gray", alpha=0.3)
        
        zero_points = np.argwhere(zero_crossings)
        for y, x in zero_points:
            plt.scatter(x, y, color='red', s=10, zorder=5)  # Plot each point

        plt.title(f"State {state} zero crossings", fontsize=16)
        plt.axis("off")
        plt.tight_layout()
        
        # Save the figure
        plt.savefig(f"./plots/{image_file}_plots/{state}/{state}_zero_crossings.png")
