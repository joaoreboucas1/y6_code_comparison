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
MIN_ANGLE = 2.5
MAX_ANGLE = 995.2679263837432

COSMOSIS_DATA_PATH = "../cosmosis/lcdm_datavector_run/"

handles_comparison_plots = [
    Line2D([], [], marker="o", color="C0", markersize=20, label="Cocoa"),
    Line2D([], [], marker="P", color="C1", markersize=20, label="Cosmosis"),
    Line2D([], [], marker="^", color="C2", markersize=20, label="Cosmolike Lighthouse"),
]
handles_error_plots = [
    Line2D([], [], marker="o", color="C3", markersize=20, label="A: Cocoa, B: Cosmosis"),
    Line2D([], [], marker="P", color="C4", markersize=20, label="A: Cocoa, B: Cosmolike"),
    Line2D([], [], marker="^", color="C5", markersize=20, label="A: Cosmolike, B: Cosmosis"),
]

def get_theta(theta_min, theta_max, num_theta):
    """
    Angle values used in cocoa
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

theta = get_theta(MIN_ANGLE, MAX_ANGLE, NUM_ANG_BINS)

def load_cosmolike_formatted_data(modelvector_path):
    _, dv   = np.loadtxt(modelvector_path, unpack=True)
    shear = dv[:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE]
    ggl = dv[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE]
    gc = dv[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE:]
    return shear, ggl, gc

def load_cocoa_data():
    return load_cosmolike_formatted_data("../cocoa/LCDM.modelvector")

def load_cosmolike_data():
    return load_cosmolike_formatted_data("../cosmolike_lighthouse/LCDM_lighthouse.modelvector")

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

def load_cosmosis_gc_data():
    # Load Cosmosis data
    gc_cosmosis = []
    theta_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/theta.txt")*RAD_TO_ARCMIN
    theta_edges_cosmosis = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/theta_edges.txt")*RAD_TO_ARCMIN
    for i in range(1, NUM_LENS_BINS + 1):
        gc_i = np.loadtxt(f"{COSMOSIS_DATA_PATH}/galaxy_xi/bin_{i}_{i}.txt")
        gc_cosmosis.append(gc_i)

    return theta_cosmosis, np.array(gc_cosmosis).flatten()

def load_cosmosis_data():
    theta_cosmosis, shear_cosmosis = load_cosmosis_shear_data()
    _, ggl_cosmosis = load_cosmosis_ggl_data()
    _, gc_cosmosis = load_cosmosis_gc_data()
    return theta_cosmosis, shear_cosmosis, ggl_cosmosis, gc_cosmosis

def plot_shear_datavectors(theta, shear_cosmosis, shear_cocoa, shear_cosmolike):
    num_cols = NUM_SOURCE_BINS+1

    xi_plus_cocoa  = shear_cocoa[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cocoa = shear_cocoa[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_plus_cosmosis  = shear_cosmosis[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
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
                ax.loglog(theta, xi_plus_cocoa[start:end], color='C0', marker='o', linestyle='-')
                ax.loglog(theta, xi_plus_cosmosis[start:end], color='C1', marker='P', linestyle='-')
                ax.loglog(theta, xi_plus_cosmolike[start:end], color='C2', marker='^', linestyle='-')
                line_idx = plus_idx
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                ax.loglog(theta, xi_minus_cocoa[start:end], color='C0', marker='o', linestyle='-')
                ax.loglog(theta, xi_minus_cosmosis[start:end], color='C1', marker='P', linestyle='-')
                ax.loglog(theta, xi_minus_cosmolike[start:end], color='C2', marker='^', linestyle='-')
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

    fig.legend(handles=handles_comparison_plots, bbox_to_anchor=(0.9, 1.07), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated Shear', fontsize=24)
    plt.savefig("comparison_shear.pdf", bbox_inches="tight")
    plt.show()

def plot_shear_relative_errors(theta, shear_cosmosis, shear_cocoa, shear_cosmolike):
    num_cols = NUM_SOURCE_BINS+1
    NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2

    xi_plus_cocoa  = shear_cocoa[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    xi_minus_cocoa = shear_cocoa[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
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
                rel_error_cocoa_cosmosis     = (xi_plus_cocoa[start:end] - xi_plus_cosmosis[start:end])/xi_plus_cosmosis[start:end]
                rel_error_cocoa_cosmolike    = (xi_plus_cocoa[start:end] - xi_plus_cosmolike[start:end])/xi_plus_cosmolike[start:end]
                rel_error_cosmolike_cosmosis = (xi_plus_cosmolike[start:end] - xi_plus_cosmosis[start:end])/xi_plus_cosmosis[start:end]
                ax.semilogx(theta, rel_error_cocoa_cosmosis, color='C3', marker='o', linestyle='-')
                ax.semilogx(theta, rel_error_cocoa_cosmolike, color='C4', marker='P', linestyle='-')
                ax.semilogx(theta, rel_error_cosmolike_cosmosis, color='C5', marker='^', linestyle='-')
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                rel_error_cocoa_cosmosis     = (xi_minus_cocoa[start:end] - xi_minus_cosmosis[start:end])/xi_minus_cosmosis[start:end]
                rel_error_cocoa_cosmolike    = (xi_minus_cocoa[start:end] - xi_minus_cosmolike[start:end])/xi_minus_cosmolike[start:end]
                rel_error_cosmolike_cosmosis = (xi_minus_cosmolike[start:end] - xi_minus_cosmosis[start:end])/xi_minus_cosmosis[start:end]
                ax.semilogx(theta, rel_error_cocoa_cosmosis, color='C3', marker='o', linestyle='-')
                ax.semilogx(theta, rel_error_cocoa_cosmolike, color='C4', marker='P', linestyle='-')
                ax.semilogx(theta, rel_error_cosmolike_cosmosis, color='C5', marker='^', linestyle='-')
                minus_idx += 1
                label = r'$\xi_-$'

            if i == NUM_SOURCE_BINS - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\frac{\\xi_A - \\xi_B}{\\xi_B}$', fontsize=20)

            ax.set_title(f'{label} bins {i+1},{j+1}', fontsize=12)
            ax.set_ylim([-0.1, 0.1])
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.legend(handles=handles_error_plots, bbox_to_anchor=(1, 1.07), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors between Cocoa and Cosmosis', fontsize=24)
    plt.savefig("errors_shear.pdf", bbox_inches="tight")
    plt.show()

def plot_ggl_datavectors(theta, ggl_cosmosis, ggl_cocoa, ggl_cosmolike):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            ax.loglog(theta, ggl_cocoa[start:end], color='C0', marker='o', linestyle='-')
            ax.loglog(theta, ggl_cosmosis[start:end], color='C1', marker='P', linestyle='-')
            ax.loglog(theta, ggl_cosmolike[start:end], color='C2', marker='^', linestyle='-')
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\gamma_t$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.legend(handles=handles_comparison_plots, bbox_to_anchor=(0.9, 1.1), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated GGL', fontsize=20, )
    plt.savefig("comparison_ggl.pdf", bbox_inches="tight")
    plt.show()

def plot_ggl_relative_errors(theta, ggl_cosmosis, ggl_cocoa, ggl_cosmolike):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            rel_error_cocoa_cosmosis     = (ggl_cocoa[start:end] - ggl_cosmosis[start:end])/ggl_cosmosis[start:end]
            rel_error_cocoa_cosmolike    = (ggl_cocoa[start:end] - ggl_cosmolike[start:end])/ggl_cosmolike[start:end]
            rel_error_cosmolike_cosmosis = (ggl_cosmolike[start:end] - ggl_cosmosis[start:end])/ggl_cosmosis[start:end]
            ax.semilogx(theta, rel_error_cocoa_cosmosis, color='C3', marker='o', linestyle='-')
            ax.semilogx(theta, rel_error_cocoa_cosmolike, color='C4', marker='P', linestyle='-')
            ax.semilogx(theta, rel_error_cosmolike_cosmosis, color='C5', marker='^', linestyle='-')
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\frac{\\gamma_t^A - \\gamma_t^B}{\\gamma_t^B}$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            # ax.set_ylim([-0.1, 0.1])
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.legend(handles=handles_error_plots, bbox_to_anchor=(1.0, 1.1), fontsize=20, ncol=3)
    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors in GGL between cocoa and Cosmosis', fontsize=20)
    plt.savefig("errors_ggl.pdf", bbox_inches="tight")
    plt.show()

def plot_gc_datavectors(theta, gc_cosmosis, gc_cocoa, gc_cosmolike):
    num_cols = NUM_LENS_BINS

    fig, axes = plt.subplots(1, num_cols, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    
    for col in range(num_cols):
        ax = axes[col]
        start = col * BLOCK_SIZE
        end    = start + BLOCK_SIZE
        ax.loglog(theta, gc_cocoa[start:end], color='C0', marker='o', linestyle='-')
        ax.loglog(theta, gc_cosmosis[start:end], color='C1', marker='P', linestyle='-')
        ax.loglog(theta, gc_cosmolike[start:end], color='C2', marker='^', linestyle='-')

        ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
        if col == 0:
            ax.set_ylabel("$w$", fontsize=20)

        ax.set_title(f'$w^{col+1}$', fontsize=12)
        ax.set_xlim([2, 255])
        ax.grid(True, linewidth=0.5, alpha=0.5)

    fig.legend(handles=handles_comparison_plots, bbox_to_anchor=(0.85, 1.2), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated Galaxy Clustering', fontsize=20)
    plt.savefig("comparison_gc.pdf", bbox_inches="tight")
    plt.show()

def plot_gc_relative_errors(theta, gc_cosmosis, gc_cocoa, gc_cosmolike):
    num_cols = NUM_LENS_BINS

    fig, axes = plt.subplots(1, num_cols, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    
    for col in range(num_cols):
        ax = axes[col]
        start = col * BLOCK_SIZE
        end    = start + BLOCK_SIZE
        rel_error_cocoa_cosmosis     = (gc_cocoa[start:end] - gc_cosmosis[start:end])/gc_cosmosis[start:end]
        rel_error_cocoa_cosmolike    = (gc_cocoa[start:end] - gc_cosmolike[start:end])/gc_cosmolike[start:end]
        rel_error_cosmolike_cosmosis = (gc_cosmolike[start:end] - gc_cosmosis[start:end])/gc_cosmosis[start:end]
        ax.semilogx(theta, rel_error_cocoa_cosmosis, color='C3', marker='o', linestyle='-')
        ax.semilogx(theta, rel_error_cocoa_cosmolike, color='C4', marker='P', linestyle='-')
        ax.semilogx(theta, rel_error_cosmolike_cosmosis, color='C5', marker='^', linestyle='-')

        ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
        if col == 0:
            ax.set_ylabel(r"$\frac{ w_A - w_B }{w_B}$", fontsize=20)

        ax.set_title(f'$w^{col+1}$', fontsize=12)
        ax.set_xlim([2, 255])
        ax.set_ylim([-0.07, 0.07])
        ax.grid(True, linewidth=0.5, alpha=0.5)
    
    fig.legend(handles=handles_error_plots, fontsize=20, bbox_to_anchor=(1, 1.2), ncol=3)
    fig.suptitle('DES-Y6 Simulated Data Vectors: Clustering relative errors between cocoa and Cosmosis', fontsize=20)
    plt.savefig("errors_gc.pdf", bbox_inches="tight")
    plt.show()
