import os
import re
import numpy as np
import matplotlib.pyplot as plt

PAIR_RE = re.compile(r"^\s*(\d+)\s+([0-9eE\.\+\-]+)\s*$")

def read_k_history(path: str):
    cycles = []
    keffs = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = PAIR_RE.match(line)
            if m:
                cycles.append(int(m.group(1)))
                keffs.append(float(m.group(2)))

    if not cycles:
        raise RuntimeError(f"No (cycle, keff) pairs parsed from {path}")

    cycles = np.asarray(cycles, dtype=int)
    keffs = np.asarray(keffs, dtype=float)

    order = np.argsort(cycles)
    return cycles[order], keffs[order]

def rolling_mean(x: np.ndarray, window: int):
    if window <= 1 or len(x) < window:
        return None
    w = np.ones(window, dtype=float) / float(window)
    return np.convolve(x, w, mode="valid")

def plot_k_trend(
    infile: str,
    out_png: str = "plots/k_trend.png",
    window: int = 50,
    inactive_cut: int = 1500,
):
    cycles, keffs = read_k_history(infile)

    active_mask = cycles >= inactive_cut
    if not np.any(active_mask):
        raise RuntimeError(f"No active cycles found with cut={inactive_cut}")

    k_active_mean = float(np.mean(keffs[active_mask]))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(cycles, keffs, linewidth=1.0, label="keff per cycle")

    rm = rolling_mean(keffs, window)
    if rm is not None:
        ax.plot(cycles[window - 1:], rm, linewidth=2.0, label=f"rolling mean (window={window})")

    # Vertical cut line
    ax.axvline(inactive_cut, linewidth=2.0, linestyle="--", label=f"cut = {inactive_cut}")

    # Active mean line
    ax.axhline(k_active_mean, linewidth=2.0, linestyle="-.", label=f"active mean = {k_active_mean:.6f}")

    # Region labels + explicit mean text near the cut line
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    y_text = ymax - 0.05 * (ymax - ymin)

    ax.text(inactive_cut - 0.02 * (cycles.max() - cycles.min()), y_text, "Inactive", ha="right", va="top")
    ax.text(inactive_cut + 0.02 * (cycles.max() - cycles.min()), y_text, "Active", ha="left", va="top")

    ax.text(
        xmin + 0.01 * (xmax - xmin),
        k_active_mean,
        f"  mean={k_active_mean:.6f}",
        ha="left",
        va="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7),
    )

    ax.set_xlabel("Cycle")
    ax.set_ylabel("keff")
    ax.set_title("keff history")
    ax.grid(True)
    ax.legend()

    out_dir = os.path.dirname(out_png)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    plot_k_trend(
        "k_history_20260205_165504.txt",
        out_png="plots/k_trend.png",
        window=50,
        inactive_cut=1500,
    )
