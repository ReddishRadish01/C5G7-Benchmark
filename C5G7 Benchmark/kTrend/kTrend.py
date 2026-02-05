import os
import re
import numpy as np
import matplotlib.pyplot as plt

PAIR_RE = re.compile(r"^\s*(\d+)\s+([0-9eE\.\+\-]+)\s*$")
KMEAN_RE = re.compile(r"^\s*k_mean:\s*([0-9eE\.\+\-]+)\s*$")

def read_k_history(path: str):
    cycles = []
    keffs = []
    k_mean = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            m = KMEAN_RE.match(line)
            if m:
                try:
                    k_mean = float(m.group(1))
                except ValueError:
                    pass
                continue

            m = PAIR_RE.match(line)
            if m:
                c = int(m.group(1))
                k = float(m.group(2))
                cycles.append(c)
                keffs.append(k)

    if not cycles:
        raise RuntimeError(f"No (cycle, keff) pairs parsed from {path}")

    cycles = np.asarray(cycles, dtype=int)
    keffs = np.asarray(keffs, dtype=float)

    # ensure sorted by cycle
    order = np.argsort(cycles)
    cycles = cycles[order]
    keffs = keffs[order]

    return cycles, keffs, k_mean

def rolling_mean(x: np.ndarray, window: int):
    if window <= 1:
        return x.copy()
    w = np.ones(window, dtype=float) / float(window)
    return np.convolve(x, w, mode="valid")

def plot_k_trend(infile: str, out_png: str = "k_trend.png", window: int = 50):
    cycles, keffs, k_mean = read_k_history(infile)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(cycles, keffs, linewidth=1.0, label="keff per cycle")

    if len(keffs) >= window:
        rm = rolling_mean(keffs, window)
        ax.plot(cycles[window - 1 :], rm, linewidth=2.0, label=f"rolling mean (window={window})")

    if k_mean is not None:
        ax.axhline(k_mean, linewidth=1.5, label=f"k_mean = {k_mean:.6f}")

    ax.set_xlabel("Cycle")
    ax.set_ylabel("keff")
    ax.set_title("keff history")
    ax.grid(True)
    ax.legend()

    out_dir = os.path.dirname(out_png)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    plot_k_trend("k_history_20260205_161314.txt", out_png="plots/k_trend.png", window=50)
