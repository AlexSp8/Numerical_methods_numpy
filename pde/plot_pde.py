
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def plot_contour(mesh, u: np.ndarray[tuple[int]],
    u_exact: np.ndarray[tuple[int]]):
    """
    Plots the numerical solution, analytical solution, and absolute error 
    for a 2D structured or curvilinear mesh.
    """
    # 1. Extract coordinates from your mesh class
    # Assumes mesh.x_mesh is an array of shape (nodes_total, 2)
    x_coords = mesh.x_mesh[:, 0]
    y_coords = mesh.x_mesh[:, 1]
    
    # 3. Compute the local absolute error
    abs_error = np.abs(u - u_exact)
    
    # 4. Create a uniform evaluation grid for smooth contour plotting
    xi = np.linspace(0.0, 1.0, 200)
    yi = np.linspace(0.0, 1.0, 200)
    X, Y = np.meshgrid(xi, yi)
    
    # Interpolate the flat unstructured data onto the uniform plotting grid
    U_num = griddata((x_coords, y_coords), u, (X, Y), method='cubic')
    U_ana = griddata((x_coords, y_coords), u_exact, (X, Y), method='cubic')
    U_err = griddata((x_coords, y_coords), abs_error, (X, Y), method='cubic')
    
    # 5. Initialize the matplotlib figure
    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    
    # --- Plot 1: Numerical Solution ---
    contour1 = axes[0].contourf(X, Y, U_num, levels=10, cmap='plasma')
    fig.colorbar(contour1, ax=axes[0], label='Temperature')
    axes[0].set_title('Numerical Solution ($u_{num}$)')
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    
    # --- Plot 2: Analytical Solution ---
    contour2 = axes[1].contourf(X, Y, U_ana, levels=10, cmap='plasma')
    fig.colorbar(contour2, ax=axes[1], label='Temperature')
    axes[1].set_title('Analytical Solution ($u_{exact}$)')
    axes[1].set_xlabel('X')
    
    # --- Plot 3: Absolute Error Distribution ---
    # Uses a log scale colorbar if errors span orders of magnitude
    contour3 = axes[2].contourf(X, Y, U_err, levels=10, cmap='viridis')
    fig.colorbar(contour3, ax=axes[2], label='Absolute Error')
    axes[2].set_title('Spatial Error Distribution ($|u_{num} - u_{exact}|$)')
    axes[2].set_xlabel('X')
    
    # Optional: Overlay your actual node points as tiny dots to see grid density
    # axes[0].scatter(x_coords, y_coords, color='white', s=1, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def plot_3d_volume(mesh, u_numerical, u_analytical):

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), subplot_kw={'projection': '3d'})
    
    # Extract structural spatial trajectories
    xm = mesh.x_mesh[:, 0]
    ym = mesh.x_mesh[:, 1]
    zm = mesh.x_mesh[:, 2]
    
    # Draw every node as a colored sphere
    sc = axes[0].scatter(xm, ym, zm, c=u_numerical, cmap='jet', s=60, edgecolors='black', alpha=0.8)
    
    fig.colorbar(sc, ax=axes[0], label='u_numerical')
    axes[0].set_title("3D Field Distribution Vector")
    axes[0].set_xlabel('X Axis')
    axes[0].set_ylabel('Y Axis')
    axes[0].set_zlabel('Z Axis')
    
    sc2 = axes[1].scatter(xm, ym, zm, c=u_analytical, cmap='jet', s=60, edgecolors='black', alpha=0.8)
    
    fig.colorbar(sc2, ax=axes[1], label='u_analytical')
    axes[1].set_title("3D Field Distribution Vector")
    axes[1].set_xlabel('X Axis')
    axes[1].set_ylabel('Y Axis')
    axes[1].set_zlabel('Z Axis')
    
    du = np.abs(u_numerical-u_analytical)
    sc3 = axes[2].scatter(xm, ym, zm, c=du, cmap='jet', s=60, edgecolors='black', alpha=0.8)
    
    fig.colorbar(sc3, ax=axes[2], label='Error')
    axes[2].set_title("3D Field Distribution Vector")
    axes[2].set_xlabel('X Axis')
    axes[2].set_ylabel('Y Axis')
    axes[2].set_zlabel('Z Axis')
    
    plt.show()