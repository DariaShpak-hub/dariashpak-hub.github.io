import numpy as np
from scipy.ndimage import laplace
import matplotlib.pyplot as plt

# ==========================================
# 1. THE HINGE TEST (Coherence Resonance)
# ==========================================
def sweep_coherence_resonance(D_values, dt=0.01, T=5000, rng=None):
    """
    Sweeps noise intensity D to find non-monotonic 1/CV.
    Evaluates a single stochastic FHN node in the excitable regime.
    """
    print("Running Hinge Test: 1/CV Sweep...")
    # Standard excitable FHN parameters
    epsilon, a, b = 0.08, 0.7, 0.8
    I_ext = 0.3 # Safely below deterministic oscillation threshold

    steps = int(T / dt)
    sqrtdt = np.sqrt(dt)
    rng = np.random.default_rng() if rng is None else rng

    regularity_metrics = []

    for D in D_values:
        v, w = -1.2, -0.62 # Approximate resting state
        spike_times = []
        is_spiking = False

        # Draw the whole noise trace up front instead of calling
        # np.random.normal() once per step: same distribution, far fewer
        # Python-level calls over hundreds of thousands of steps.
        noise = rng.normal(size=steps) * np.sqrt(2 * D) * sqrtdt

        for t_idx in range(steps):
            # Euler-Maruyama integration
            dv = (v - (v**3)/3 - w + I_ext) * dt + noise[t_idx]
            dw = epsilon * (v + a - b * w) * dt
            v += dv
            w += dw

            # Simple threshold crossing for spike detection
            if v > 0.0 and not is_spiking:
                spike_times.append(t_idx * dt)
                is_spiking = True
            elif v < -0.5:
                is_spiking = False

        # Calculate 1/CV
        if len(spike_times) > 2:
            isis = np.diff(spike_times)
            cv = np.std(isis) / np.mean(isis)
            regularity = 1.0 / cv if cv > 0 else 0
        else:
            regularity = 0

        regularity_metrics.append(regularity)

    return D_values, regularity_metrics

# ==========================================
# 2. THE R_c SCALING TEST
# ==========================================
def check_droplet_survival(R, g, delta_I, N=50, dt=0.05, T=150):
    """
    Seeds a droplet of radius R in a deterministic 2D network.
    Returns True if the droplet consumes the network, False if it collapses.
    """
    epsilon, a, b = 0.08, 0.7, 0.8
    I_th = 0.33 # Approximate theoretical threshold
    I_ext = I_th - delta_I

    # Explicit-Euler + 5-point Laplacian on a unit grid is only stable for
    # g*dt < 0.25; past that the diffusion term blows up and silently
    # corrupts the survival verdict instead of raising.
    if g * dt >= 0.25:
        raise ValueError(
            f"Unstable step size for explicit Euler diffusion: "
            f"g*dt={g * dt:.3f} >= 0.25 (g={g}, dt={dt}). Reduce dt or g."
        )

    steps = int(T / dt)

    # Initialize resting state grid
    v = np.full((N, N), -1.2)
    w = np.full((N, N), -0.62)

    # Seed the central droplet
    center = N // 2
    y, x = np.ogrid[-center:N-center, -center:N-center]
    mask = x**2 + y**2 <= R**2
    v[mask] = 2.0 # Force into active state
    w[mask] = 0.0

    for t in range(steps):
        # Calculate spatial coupling (Laplacian for nearest neighbor)
        coupling = g * laplace(v)

        dv = (v - (v**3)/3 - w + I_ext + coupling) * dt
        dw = epsilon * (v + a - b * w) * dt

        v += dv
        w += dw

    # Check final state: did the activation spread or die?
    active_fraction = np.sum(v > 0) / (N * N)
    return active_fraction > 0.5

def find_critical_radius(g_values, delta_I_values, R_max=15):
    """
    Uses binary search to find Rc across different coupling and barrier
    parameters. Cells where no radius up to R_max survives are reported
    as NaN, distinguishing "no critical radius found in range" from
    "Rc happens to equal R_max".
    """
    print("Running Rc Scaling Test...")
    results = np.full((len(g_values), len(delta_I_values)), np.nan)

    for i, g in enumerate(g_values):
        for j, delta_I in enumerate(delta_I_values):
            # Binary search for Rc
            low_R, high_R = 1, R_max
            Rc = None

            while low_R <= high_R:
                mid_R = (low_R + high_R) // 2
                survives = check_droplet_survival(mid_R, g, delta_I)

                if survives:
                    Rc = mid_R
                    high_R = mid_R - 1 # Try a smaller droplet
                else:
                    low_R = mid_R + 1  # Need a bigger droplet

            results[i, j] = Rc if Rc is not None else np.nan
    return results

# ==========================================
# 3. NULL BATTERY GENERATOR
# ==========================================
def generate_null_battery(mode, N=20, dt=0.05, T=1000, seed=42):
    """
    Generates NxT matrix of time series for blind pipeline testing.
    Modes: 'A' (True Handover), 'B' (Decoupled Drift), 'C' (No Collapse)
    """
    print(f"Generating Null Battery Mode {mode}...")
    epsilon, a, b = 0.08, 0.7, 0.8
    D = 0.05
    sqrtdt = np.sqrt(dt)
    steps = int(T / dt)

    # Configure modes
    g = 0.5 if mode in ['A', 'C'] else 0.0 # B is decoupled
    I_start = 0.1
    I_end = 0.32 if mode in ['A', 'B'] else 0.2 # C stops short of collapse

    # Drive ramp
    I_ramp = np.linspace(I_start, I_end, steps)

    v = np.full((N, N), -1.2)
    w = np.full((N, N), -0.62)

    # We will record the 1D flattened array over time
    history = np.zeros((N*N, steps))

    # Local RNG seeded for direct causal comparison across modes - using
    # a local Generator instead of np.random.seed() keeps this from
    # clobbering the global RNG state for any other code in the process.
    rng = np.random.default_rng(seed)

    for t in range(steps):
        noise = rng.normal(size=(N, N)) * np.sqrt(2 * D) * sqrtdt
        coupling = g * laplace(v)

        dv = (v - (v**3)/3 - w + I_ramp[t] + coupling) * dt + noise
        dw = epsilon * (v + a - b * w) * dt

        v += dv
        w += dw
        history[:, t] = v.flatten()

    return history

# ==========================================
# EXECUTION SCRIPT
# ==========================================
if __name__ == "__main__":
    # 1. Run the Hinge Test
    D_vals = np.logspace(-4, -0.5, 20)
    Ds, regularities = sweep_coherence_resonance(D_vals)

    plt.figure(figsize=(8, 4))
    plt.plot(Ds, regularities, marker='o')
    plt.xscale('log')
    plt.xlabel('Noise Intensity (D)')
    plt.ylabel('Regularity (1/CV)')
    plt.title('Hinge Test: FHN Coherence Resonance')
    plt.grid(True)
    plt.show()

    # 2. Generate the Blind Validation Matrices
    data_A = generate_null_battery('A') # Feed this to your pipeline (should pass)
    data_B = generate_null_battery('B') # Feed this to your pipeline (should fail)
    data_C = generate_null_battery('C') # Feed this to your pipeline (should fail)
