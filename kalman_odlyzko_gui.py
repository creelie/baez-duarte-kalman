
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
import numpy as np
from scipy.linalg import lstsq
from scipy.ndimage import uniform_filter1d
import time

# -----------------------------
# Load real Odlyzko zeros from single-column file
# -----------------------------
def load_odlyzko_zeros_from_file(path, max_lines=None):
    gammas = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            try:
                gamma = float(line.strip())
                gammas.append(gamma)
            except:
                continue
    zeros = 0.5 + 1j * np.array(gammas)
    print(f"Loaded {len(zeros)} Odlyzko zeros.")
    print("First 5 zeros:", zeros[:5])
    return zeros

odlyzko_zeros = load_odlyzko_zeros_from_file("C://zeros6.txt", max_lines=50000)

def rk_sequence(k, N):
    n = np.arange(1, N + 1)
    return k * (n / k - np.floor(n / k))

def get_weights(N):
    n = np.arange(1, N + 1)
    w = 1.0 / np.log(n + 2)
    return w / np.max(w)

def kalman_filter(signal, R_kalman=1e-6):
    N = len(signal)
    x = 1.0
    P = 0.0
    Q = 0.0
    output = np.zeros(N)
    for i in range(N):
        x_pred = x
        P_pred = P + Q
        K = P_pred / (P_pred + R_kalman) if (P_pred + R_kalman) != 0 else 0.0
        x = x_pred + K * (signal[i] - x_pred)
        P = (1 - K) * P_pred
        output[i] = x
    return output

# Global values to update across zoom panel
approx = None
weighted = None
kalman = None

# -----------------------------
# Update function for main GUI
# -----------------------------
def update(val):
    global approx, weighted, kalman

    start = time.time()
    epsilon = epsilon_slider.val
    M = int(M_slider.val)
    N = int(N_slider.val)
    window = 9

    one_seq = np.ones(N)
    R = np.vstack([rk_sequence(k, N) for k in range(1, M + 1)]).T
    w = get_weights(N)
    W = np.diag(np.sqrt(w))
    Rw = W @ R
    bw = W @ one_seq
    coeffs, *_ = lstsq(Rw, bw)
    approx = R @ coeffs
    weighted = W @ approx
    smoothed = uniform_filter1d(weighted, size=window, mode='nearest')
    kalman = kalman_filter(smoothed)

    ax1.clear()
    ax2.clear()

    ax1.plot(one_seq, '--', color='gray', label='Target (1)')
    ax1.plot(approx, label='Raw LS')
    ax1.plot(weighted, label='Log-Weighted', color='orange')
    ax1.plot(smoothed, label='Smoothed', color='blue')
    ax1.plot(kalman, label='Kalman Output', color='green')
    ax1.set_title("Kalman-Filtered Symbolic Approximants")
    ax1.set_ylabel("Function Value")
    ax1.legend()
    ax1.grid(True)

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 50)
    ax2.axvline(0.5 - epsilon, color='purple', linestyle='--')
    ax2.axvline(0.5 + epsilon, color='purple', linestyle='--')
    ax2.fill_betweenx(np.linspace(0, 50, 500), 0.5 - epsilon, 0.5 + epsilon,
                      color='purple', alpha=0.1, label=r'$\delta_\epsilon$')

    if odlyzko_toggle.get_status()[0]:
        mask = odlyzko_zeros.imag <= 50
        jitter = 1e-3 * np.random.randn(np.sum(mask))
        jittered_real = odlyzko_zeros.real[mask] + jitter
        ax2.scatter(jittered_real, odlyzko_zeros.imag[mask],
                    color='blue', s=15, alpha=0.8, marker='x', label='Odlyzko zeros')

    if density_toggle.get_status()[0]:
        ax2_twin = ax2.twinx()
        all_real_parts = odlyzko_zeros.real
        hist, bins = np.histogram(all_real_parts, bins=50, range=(0, 1))
        ax2_twin.plot((bins[:-1] + bins[1:]) / 2, hist, 'k--', alpha=0.5, label='Density')
        ax2_twin.set_ylabel("Zero Density", color='black')
        ax2_twin.tick_params(axis='y', labelcolor='black')

    ax2.set_title("Zero Distribution in Critical Strip")
    ax2.set_xlabel("Re(s)")
    ax2.set_ylabel("Im(s)")
    ax2.legend()
    ax2.grid(True)

    fig.canvas.draw_idle()
    print(f"Execution time: {time.time() - start:.4f} seconds")

# -----------------------------
# GUI setup
# -----------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
plt.subplots_adjust(bottom=0.35)

ax_eps = plt.axes([0.2, 0.2, 0.65, 0.03])
ax_M = plt.axes([0.2, 0.15, 0.65, 0.03])
ax_N = plt.axes([0.2, 0.10, 0.65, 0.03])

epsilon_slider = Slider(ax_eps, 'ε (band half-width)', 0.01, 0.24, valinit=0.1, valstep=0.005)
M_slider = Slider(ax_M, 'M (basis terms)', 10, 500, valinit=100, valstep=10)
N_slider = Slider(ax_N, 'N (sample size)', 200, 2000, valinit=1000, valstep=100)

check_odlyzko = plt.axes([0.01, 0.3, 0.15, 0.05])
check_density = plt.axes([0.01, 0.25, 0.15, 0.05])
odlyzko_toggle = CheckButtons(check_odlyzko, ['Show Odlyzko'], [True])
density_toggle = CheckButtons(check_density, ['Show Density'], [False])

epsilon_slider.on_changed(update)
M_slider.on_changed(update)
N_slider.on_changed(update)
odlyzko_toggle.on_clicked(update)
density_toggle.on_clicked(update)

update(None)
plt.show()

# -----------------------------
# Zoomed view with stage selector
# -----------------------------
import matplotlib.pyplot as plt2

zoom_fig, zoom_ax = plt2.subplots(figsize=(10, 4))
plt2.subplots_adjust(bottom=0.25)

zoom_slider_ax = plt2.axes([0.2, 0.1, 0.65, 0.03])
zoom_slider = Slider(zoom_slider_ax, 'Zoom Position', 0, 950, valinit=0, valstep=10)

stage_check_ax = plt2.axes([0.01, 0.25, 0.15, 0.15])
stage_checkboxes = CheckButtons(stage_check_ax, ['Raw LS', 'Log-Weighted', 'Kalman'], [True, True, True])

def update_zoom(val):
    index_start = int(zoom_slider.val)
    index_end = index_start + 50
    if index_end > len(kalman):
        index_end = len(kalman)

    zoom_ax.clear()
    zoom_ax.set_title("Zoomed-in View of Approximation Convergence")
    zoom_ax.set_xlim(index_start, index_end)
    zoom_ax.grid(True)

    x = np.arange(index_start, index_end)
    show = stage_checkboxes.get_status()

    if show[0]:
        zoom_ax.plot(x, approx[index_start:index_end], label='Raw LS')
    if show[1]:
        zoom_ax.plot(x, weighted[index_start:index_end], label='Log-Weighted', color='orange')
    if show[2]:
        zoom_ax.plot(x, kalman[index_start:index_end], label='Kalman', color='green')

    zoom_ax.legend()
    zoom_fig.canvas.draw_idle()

zoom_slider.on_changed(update_zoom)
stage_checkboxes.on_clicked(update_zoom)

update_zoom(None)
plt2.show()
