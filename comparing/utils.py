from itertools import combinations
import numpy as np
import matplotlib as mpl
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

colors_comparison_plots = mpl.colors.TABLEAU_COLORS
colors_error_plots = mpl.colors.XKCD_COLORS

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

def plot_shear_datavectors(theta, shear_dvs, labels):
    num_cols = NUM_SOURCE_BINS+1

    xi_pluses = [shear_dv[:NUM_SHEAR_BLOCKS*BLOCK_SIZE] for shear_dv in shear_dvs]
    xi_minuses = [shear_dv[NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE] for shear_dv in shear_dvs]
    
    # xi_plus_cocoa  = shear_cocoa[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_minus_cocoa = shear_cocoa[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_plus_cosmosis  = shear_cosmosis[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_plus_cosmolike  = shear_cosmolike[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_minus_cosmolike = shear_cosmolike[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_plus_cosmosis  = shear_cosmosis[:NUM_SHEAR_BLOCKS * BLOCK_SIZE]
    # xi_minus_cosmosis = shear_cosmosis[NUM_SHEAR_BLOCKS * BLOCK_SIZE:2 * NUM_SHEAR_BLOCKS * BLOCK_SIZE]

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
                for xi_plus, color in zip(xi_pluses, colors_comparison_plots):
                    ax.loglog(theta, xi_plus[start:end], color=color, marker="o")
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                for xi_minus, color in zip(xi_minuses, colors_comparison_plots):
                    ax.loglog(theta, xi_minus[start:end], color=color, marker="o")
                minus_idx += 1
                label = r'$\xi_-$'

            if i == NUM_SOURCE_BINS - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\xi$', fontsize=20)

            ax.set_title(f'{label} bins {i+1},{j+1}', fontsize=12)
            ax.set_xlim([2, 255])
            ax.set_ylim([1e-7, 2e-4])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
        Line2D([], [], color=color, marker="o", label=label) for color, label in zip(colors_comparison_plots, labels)
    ]
    fig.legend(handles=handles, bbox_to_anchor=(0.9, 1.07), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated Shear', fontsize=24)
    plt.savefig("comparison_shear.pdf", bbox_inches="tight")
    plt.show()

def plot_shear_relative_errors(theta, shear_dvs, labels, save_fig_filename="errors_shear.pdf"):
    num_cols = NUM_SOURCE_BINS+1
    NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2

    xi_pluses = [shear_dv[:NUM_SHEAR_BLOCKS*BLOCK_SIZE] for shear_dv in shear_dvs]
    xi_minuses = [shear_dv[NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE] for shear_dv in shear_dvs]
    
    fig, axes = plt.subplots(NUM_SOURCE_BINS, num_cols, figsize=(13, 13), sharex=True, sharey=True, constrained_layout=True)
    plus_idx = 0
    minus_idx = 0

    for i in range(NUM_SOURCE_BINS):
        for col in range(num_cols):
            ax = axes[i, col]
            # JVR NOTE: this weird next line calculates all possible pairs between the provided shear_dvs
            # JVR NOTE 2: this needs to be calculated at every iteration since the iterator gets exhausted
            comparison_pairs = combinations(range(len(shear_dvs)), 2)
            if col <= i:
                j = col
                start = plus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                ax.axhline(color="black", ls="--")
                for (a, b), color in zip(comparison_pairs, colors_error_plots):
                    rel_error = (xi_pluses[a][start:end] - xi_pluses[b][start:end])/xi_pluses[b][start:end]
                    ax.semilogx(theta, rel_error, color=color, marker='o')
                plus_idx += 1
                label = r'$\xi_+$'
            else:
                j = col - 1
                start = minus_idx * BLOCK_SIZE
                end = start + BLOCK_SIZE
                ax.axhline(color="black", ls="--")
                for (a, b), color in zip(comparison_pairs, colors_error_plots):
                    rel_error = (xi_minuses[a][start:end] - xi_minuses[b][start:end])/xi_minuses[b][start:end]
                    ax.semilogx(theta, rel_error, color=color, marker='o')
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

    handles = [Line2D([], [], label=f"A: {labels[a]}, B: {labels[b]}", color=color, marker="o") for (a, b), color in zip(combinations(range(len(shear_dvs)), 2), colors_error_plots)]
    fig.legend(handles=handles, bbox_to_anchor=(1, 0), fontsize=20, ncol=1)

    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors in Shear', fontsize=24)
    plt.savefig(save_fig_filename, bbox_inches="tight")
    plt.show()

def plot_ggl_datavectors(theta, ggl_dvs, labels, save_fig_filename="comparison_ggl.pdf"):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            for ggl_dv, color in zip(ggl_dvs, colors_comparison_plots):
                ax.loglog(theta, ggl_dv[start:end], color=color, marker='o')
            
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\gamma_t$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
        Line2D([], [], color=color, marker="o", label=label) for color, label in zip(colors_comparison_plots, labels)
    ]
    fig.legend(handles=handles, bbox_to_anchor=(0.9, 1.1), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated GGL', fontsize=20, )
    plt.savefig(save_fig_filename, bbox_inches="tight")
    plt.show()

def plot_ggl_relative_errors(theta, ggl_dvs, labels, save_fig_filename="errors_ggl.pdf"):
    num_cols = NUM_LENS_BINS
    num_rows = NUM_SOURCE_BINS

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(13, 10), sharex=True, sharey=True, constrained_layout=True)
    
    for row in range(num_rows):
        for col in range(num_cols):
            ax = axes[row, col]
            start = (row*NUM_SOURCE_BINS + col) * BLOCK_SIZE
            end = start + BLOCK_SIZE
            comparison_pairs = combinations(range(len(ggl_dvs)), 2)
            for (a, b), color in zip(comparison_pairs, colors_error_plots):
                rel_error = (ggl_dvs[a][start:end] - ggl_dvs[b][start:end])/ggl_dvs[b][start:end]
                ax.semilogx(theta, rel_error, color=color, marker='o')
            
            ax.axhline(color="black", ls="--")
            
            label = r'$\gamma_t$'

            if row == num_rows - 1:
                ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
            if col == 0:
                ax.set_ylabel('$\\frac{\\gamma_t^A - \\gamma_t^B}{\\gamma_t^B}$', fontsize=20)

            ax.set_title(f'{label} bins {row+1},{col+1}', fontsize=12)
            ax.set_ylim([-0.1, 0.1])
            ax.set_xlim([2, 255])
            ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [Line2D([], [], label=f"A: {labels[a]}, B: {labels[b]}", color=color, marker="o") for (a, b), color in zip(combinations(range(len(ggl_dvs)), 2), colors_error_plots)]
    fig.legend(handles=handles, bbox_to_anchor=(1, 0), fontsize=20, ncol=1)
    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative Errors in GGL', fontsize=20)
    plt.savefig(save_fig_filename, bbox_inches="tight")
    plt.show()

def plot_gc_datavectors(theta, gc_dvs, labels, save_fig_filename="comparison_gc.pdf"):
    num_cols = NUM_LENS_BINS

    fig, axes = plt.subplots(1, num_cols, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    
    for col in range(num_cols):
        ax = axes[col]
        start = col * BLOCK_SIZE
        end    = start + BLOCK_SIZE
        for gc_dv, color in zip(gc_dvs, colors_comparison_plots):
            ax.loglog(theta, gc_dv[start:end], color=color, marker='o')

        ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
        if col == 0:
            ax.set_ylabel("$w$", fontsize=20)

        ax.set_title(f'$w^{col+1}$', fontsize=12)
        ax.set_xlim([2, 255])
        ax.grid(True, linewidth=0.5, alpha=0.5)

    handles = [
        Line2D([], [], color=color, marker="o", label=label) for color, label in zip(colors_comparison_plots, labels)
    ]
    fig.legend(handles=handles, bbox_to_anchor=(0.9, 1.1), fontsize=20, ncol=3)

    fig.suptitle('DES-Y6 Simulated Galaxy Clustering', fontsize=20)
    plt.savefig(save_fig_filename, bbox_inches="tight")
    plt.show()

def plot_gc_relative_errors(theta, gc_dvs, labels, save_fig_filename="errors_gc.pdf"):
    num_cols = NUM_LENS_BINS

    fig, axes = plt.subplots(1, num_cols, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    
    for col in range(num_cols):
        ax = axes[col]
        start = col * BLOCK_SIZE
        end    = start + BLOCK_SIZE
        ax.axhline(color="black", ls="--")
        comparison_pairs = combinations(range(len(gc_dvs)), 2)
        for (a, b), color in zip(comparison_pairs, colors_error_plots):
            rel_error = (gc_dvs[a][start:end] - gc_dvs[b][start:end])/gc_dvs[b][start:end]
            ax.semilogx(theta, rel_error, color=color, marker='o')

        ax.set_xlabel('$\\theta$ (arcmin)', fontsize=20)
        if col == 0:
            ax.set_ylabel(r"$\frac{ w_A - w_B }{w_B}$", fontsize=20)

        ax.set_title(f'$w^{col+1}$', fontsize=12)
        ax.set_xlim([2, 255])
        ax.set_ylim([-0.07, 0.07])
        ax.grid(True, linewidth=0.5, alpha=0.5)
    
    handles = [Line2D([], [], label=f"A: {labels[a]}, B: {labels[b]}", color=color, marker="o") for (a, b), color in zip(combinations(range(len(gc_dvs)), 2), colors_error_plots)]
    fig.legend(handles=handles, bbox_to_anchor=(1, 0), fontsize=20, ncol=1)
    fig.suptitle('DES-Y6 Simulated Data Vectors: Relative errors in Clustering', fontsize=20)
    plt.savefig(save_fig_filename, bbox_inches="tight")
    plt.show()
