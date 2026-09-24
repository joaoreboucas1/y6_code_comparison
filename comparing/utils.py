import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ARCMIN_TO_RAD = 2.90888208665721580e-4
RAD_TO_ARCMIN = 1/ARCMIN_TO_RAD
NUM_SOURCE_BINS = 4
NUM_LENS_BINS = 6
NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2
NUM_ANG_BINS = 26
BLOCK_SIZE = NUM_ANG_BINS

COSMOSIS_DATA_PATH = "../cosmosis/lcdm_datavector_run/"

def get_theta(theta_min, theta_max, num_theta):
    """
    Angle values used in Cosmolike
    """

    THETA_MIN_RAD = theta_min * ARCMIN_TO_RAD
    THETA_MAX_RAD = theta_max * ARCMIN_TO_RAD
    DLOG_THETA = (np.log(THETA_MAX_RAD) - np.log(THETA_MIN_RAD))/num_theta
    theta = np.zeros(num_theta)
    for i in range(num_theta):
        THETA_MIN_BIN = np.exp(np.log(THETA_MIN_RAD) + i * DLOG_THETA)
        THETA_MAX_BIN = np.exp(np.log(THETA_MIN_RAD) + (i + 1) * DLOG_THETA)
        theta[i] = (2/3) * (THETA_MAX_BIN**3 - THETA_MIN_BIN**3) / (THETA_MAX_BIN**2 - THETA_MIN_BIN**2)

    theta *= RAD_TO_ARCMIN
    return theta

def load_cosmosis_shear_data():
    # Load Cosmosis data
    shear_cosmosis = []
    theta_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/shear_xi_plus/theta.txt")*RAD_TO_ARCMIN
    theta_edges_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/shear_xi_plus/theta_edges.txt")*RAD_TO_ARCMIN
    for corr_func in ["plus", "minus"]:
        for i in range(1, NUM_SOURCE_BINS + 1):
            for j in range(i, NUM_SOURCE_BINS + 1):
                xi_ij = np.loadtxt(f"{COSMOSIS_DATA_PATH}/shear_xi_{corr_func}/bin_{j}_{i}.txt")

                # theta_cosmosis is different than our theta, we interpolate
                # xi_ij_interp = interp1d(theta_cosmosis, xi_ij, bounds_error=False, fill_value=0.0)
                # xi_ij = xi_ij_interp(theta)
                shear_cosmosis.append(xi_ij)

    return theta_cosmosis, np.array(shear_cosmosis).flatten()

def load_cosmosis_ggl_data():
    # Load Cosmosis data
    ggl_cosmosis = []
    theta_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_shear_xi/theta.txt")*RAD_TO_ARCMIN
    theta_edges_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_shear_xi/theta_edges.txt")*RAD_TO_ARCMIN
    for j in range(1, NUM_LENS_BINS + 1):
        for i in range(1, NUM_SOURCE_BINS + 1):
            gamma_t_ij = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_shear_xi/bin_{j}_{i}.txt")
            ggl_cosmosis.append(gamma_t_ij)

    return theta_cosmosis, np.array(ggl_cosmosis).flatten()

def load_cosmosis_wg_data():
    # Load Cosmosis data
    wg_cosmosis = []
    theta_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/theta.txt")*RAD_TO_ARCMIN
    theta_edges_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/theta_edges.txt")*RAD_TO_ARCMIN
    for i in range(1, NUM_LENS_BINS + 1):
        wg_i = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/bin_{i}_{i}.txt")
        wg_cosmosis.append(wg_i)

    return theta_cosmosis, np.array(wg_cosmosis).flatten()

def plot_shear_datavectors(theta, shear_cosmosis, shear_cosmolike):
    num_cols = NUM_SOURCE_BINS+1

    xi_plus_cosmolike  = shear_cosmolike[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cosmolike = shear_cosmolike[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_plus_cosmosis  = shear_cosmosis[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cosmosis = shear_cosmosis[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]

    fig, axes = plt.subplots(NUM_SOURCE_BINS, num_cols, figsize=(13, 13), sharex=True, sharey=True, constrained_layout=True)
    plus_idx = 0
    minus_idx = 0

    for i in range(NUM_SOURCE_BINS):
        for col in range(num_cols):
            ax = axes[i, col]

            if col <= i:
                j = col
                start = plus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                ax.loglog(theta, xi_plus_cosmolike[start:end], color='C0', marker='o', linestyle='-')
                ax.loglog(theta, xi_plus_cosmosis[start:end], color='C1', marker='P', linestyle='-')
                line_idx = plus_idx
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                ax.loglog(theta, xi_minus_cosmolike[start:end], color='C0', marker='o', linestyle='-')
                ax.loglog(theta, xi_minus_cosmosis[start:end], color='C1', marker='P', linestyle='-')
                line_idx = NUM_SHEAR_BLOCKS + minus_idx
                minus_idx += 1
                label = r'$\xi_-$'

            if i == NUM_SOURCE_BINS - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\xi$', fontsize=20)

            ax.set_title(f'{label} bins {i+1},{j+1}', fontsize=12)
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
        Line2D([], [], marker="o", color="C0", markersize=20, label="Cosmolike"),
        Line2D([], [], marker="P", color="C1", markersize=20, label="Cosmosis"),
    ]
    fig.legend(handles=handles, bbox_to_anchor=(1, 1.05), fontsize=20)

    fig.suptitle('DES-Y6 Simulated Data Vectors', fontsize=24)
    plt.savefig("comparison_shear.pdf", bbox_inches="tight")
    plt.show()

def plot_shear_relative_errors(theta, shear_cosmosis, shear_cosmolike):
    num_cols = NUM_SOURCE_BINS+1
    NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2

    xi_plus_cosmolike  = shear_cosmolike[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cosmolike = shear_cosmolike[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_plus_cosmosis  = shear_cosmosis[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cosmosis = shear_cosmosis[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]

    fig, axes = plt.subplots(NUM_SOURCE_BINS, num_cols, figsize=(13, 13), sharex=True, sharey=True, constrained_layout=True)
    plus_idx = 0
    minus_idx = 0

    for i in range(NUM_SOURCE_BINS):
        for col in range(num_cols):
            ax = axes[i, col]

            if col <= i:
                j = col
                start = plus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                rel_error = (xi_plus_cosmolike[start:end] - xi_plus_cosmosis[start:end])/xi_plus_cosmosis[start:end]
                ax.semilogx(theta, rel_error, color='C3', marker='o', linestyle='-')
                line_idx = plus_idx
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                rel_error = (xi_minus_cosmolike[start:end] - xi_minus_cosmosis[start:end])/xi_minus_cosmosis[start:end]
                ax.semilogx(theta, rel_error, color='C3', marker='o', linestyle='-')
                line_idx = NUM_SHEAR_BLOCKS + minus_idx
                minus_idx += 1
                label = r'$\xi_-$'

            if i == NUM_SOURCE_BINS - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\frac{\\xi_\\mathrm{Cosmolike} - \\xi_\\mathrm{Cosmosis}}{\\xi_\\mathrm{Cosmosis}}$', fontsize=20)

            ax.set_title(f'{label} bins {i+1},{j+1}', fontsize=12)
            ax.set_ylim([-0.1, 0.1])
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors between Cosmolike and Cosmosis', fontsize=24)
    plt.savefig("errors_shear.pdf", bbox_inches="tight")
    plt.show()

def plot_ggl_datavectors(theta, ggl_cosmosis, ggl_cosmolike):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            ax.loglog(theta, ggl_cosmolike[start:end], color='C0', marker='o', linestyle='-')
            ax.loglog(theta, ggl_cosmosis[start:end], color='C1', marker='P', linestyle='-')
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\gamma_t$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
            Line2D([], [], marker="o", color="C0", markersize=20, label="Cosmolike"),
            Line2D([], [], marker="P", color="C1", markersize=20, label="Cosmosis"),
        ]
    fig.legend(handles=handles, bbox_to_anchor=(1, 1.05), fontsize=20)
    

    fig.suptitle('DES-Y6 Simulated Data Vectors: GGL comparison between Cosmolike and Cosmosis', fontsize=20)
    plt.savefig("comparison_ggl.pdf", bbox_inches="tight")
    plt.show()

def plot_ggl_relative_errors(theta, ggl_cosmosis, ggl_cosmolike):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            rel_error = (ggl_cosmolike[start:end] - ggl_cosmosis[start:end])/ggl_cosmosis[start:end]
            ax.semilogx(theta, rel_error, color='C3', marker='o', linestyle='-')
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\frac{\\gamma_t^\\mathrm{Cosmolike} - \\gamma_t^\\mathrm{Cosmosis}}{\\gamma_t^\\mathrm{Cosmosis}}$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            ax.set_ylim([-0.1, 0.1])
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors in GGL between Cosmolike and Cosmosis', fontsize=20)
    plt.savefig("errors_ggl.pdf", bbox_inches="tight")
    plt.show()

def plot_wg_datavectors(theta, wg_cosmosis, wg_cosmolike):
    num_cols = NUM_LENS_BINS

    fig, axes = plt.subplots(1, num_cols, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    
    for col in range(num_cols):
        ax = axes[col]
        start = col * BLOCK_SIZE
        end    = start + BLOCK_SIZE
        ax.loglog(theta, wg_cosmolike[start:end], color='C0', marker='o', linestyle='-')
        ax.loglog(theta, wg_cosmosis[start:end], color='C1', marker='P', linestyle='-')
        label = '$w_g$'

        ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
        if col == 0:
            ax.set_ylabel("$w_g$", fontsize=20)

        ax.set_title(f'$w^{col+1}$', fontsize=12)
        ax.set_xlim([2, 255])
        ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
            Line2D([], [], marker="o", color="C0", markersize=20, label="Cosmolike"),
            Line2D([], [], marker="P", color="C1", markersize=20, label="Cosmosis"),
        ]
    fig.legend(handles=handles, bbox_to_anchor=(1, 1.05), fontsize=20)
    

    fig.suptitle('DES-Y6 Simulated Data Vectors: Clustering comparison between Cosmolike and Cosmosis', fontsize=20)
    plt.savefig("comparison_wg.pdf", bbox_inches="tight")
    plt.show()
