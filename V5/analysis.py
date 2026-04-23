"""
analysis.py

Generate comparison plots for the three-case bacterial transport study:
    A1: flow + slip pillars
    A2: flow + noSlip pillars
    B:  no flow, run-and-tumble with pillar collisions

Run this AFTER `gui_control.py`'s "RUN COMPARISON STUDY" finishes. It reads
the three HDF5 trajectory files and produces:
    - msd_comparison.png      : mean squared displacement curves
    - density_heatmaps.png    : top-down particle density for each case
    - pillar_proximity.png    : distribution of distances to nearest pillar
    - escape_fraction.png     : fraction of particles still in domain vs time
    - summary_stats.txt       : numerical summary
"""

import os
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# Consistent styling
mpl.rcParams['figure.dpi'] = 120
mpl.rcParams['savefig.dpi'] = 150
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.alpha'] = 0.3

CASES = {
    "A1": {
        "dir": "case_A1/particle_output",
        "label": "A1: Flow + Slip Pillars",
        "color": "tab:blue",
    },
    "A2": {
        "dir": "case_A2/particle_output",
        "label": "A2: Flow + No-Slip Pillars",
        "color": "tab:orange",
    },
    "B": {
        "dir": "case_B/particle_output",
        "label": "B: No Flow + Collisions",
        "color": "tab:green",
    },
}


def load_case(case_key):
    """Load trajectories.h5 for one case. Returns None if missing."""
    h5_path = Path(CASES[case_key]["dir"]) / "trajectories.h5"
    if not h5_path.exists():
        print(f"[WARN] {h5_path} not found - skipping case {case_key}")
        return None
    with h5py.File(h5_path, "r") as f:
        data = {
            "time": f["time"][:],
            "positions": f["positions"][:],  # (T, N, 3)
            "headings": f["headings"][:],
            "active": f["active"][:],
            "ids": f["ids"][:],
        }
    print(f"[LOAD] Case {case_key}: {data['positions'].shape[1]} particles, "
          f"{len(data['time'])} snapshots, t_end={data['time'][-1]:.2f}s")
    return data


def plot_msd(cases_data, out_path):
    """Mean squared displacement for each case on a log-log axis.
    For active particles only (deactivated ones mess up the statistics)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for key, data in cases_data.items():
        if data is None:
            continue
        t = data["time"]
        pos = data["positions"]
        active = data["active"]

        # MSD at each timestep, averaged over ACTIVE particles at that time
        # Using initial position as reference for each particle
        msd = np.zeros(len(t))
        for i in range(len(t)):
            mask = active[i] & active[0]  # particles active at both t and t=0
            if mask.sum() > 0:
                displacements = pos[i, mask] - pos[0, mask]
                msd[i] = np.mean(np.sum(displacements ** 2, axis=1))
            else:
                msd[i] = np.nan

        # Avoid log(0) for the first point
        valid = (t > 0) & np.isfinite(msd) & (msd > 0)
        ax.loglog(t[valid], msd[valid],
                  label=CASES[key]["label"],
                  color=CASES[key]["color"],
                  linewidth=2)

    # Reference slopes
    t_ref = np.logspace(-3, 1, 100)
    # Ballistic: MSD ~ t^2
    ax.loglog(t_ref, 1e-11 * t_ref**2, "k--", alpha=0.3, label="Ballistic (slope 2)")
    # Diffusive: MSD ~ t
    ax.loglog(t_ref, 1e-10 * t_ref, "k:", alpha=0.3, label="Diffusive (slope 1)")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("MSD (m²)")
    ax.set_title("Mean Squared Displacement: All Three Cases")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[PLOT] Saved {out_path}")


def plot_density_heatmaps(cases_data, out_path):
    """Top-down (XY-plane) 2D histogram of particle positions, aggregated
    over all time. Shows WHERE bacteria spend their time."""
    n_cases = sum(1 for d in cases_data.values() if d is not None)
    if n_cases == 0:
        return

    fig, axes = plt.subplots(1, n_cases, figsize=(5 * n_cases, 4.5),
                             squeeze=False)
    axes = axes[0]

    # Use a consistent domain extent across panels
    all_pos = np.concatenate([d["positions"][d["active"]].reshape(-1, 3)
                              for d in cases_data.values() if d is not None])
    x_bounds = (all_pos[:, 0].min(), all_pos[:, 0].max())
    y_bounds = (all_pos[:, 1].min(), all_pos[:, 1].max())

    idx = 0
    for key, data in cases_data.items():
        if data is None:
            continue
        ax = axes[idx]
        pos = data["positions"]
        active = data["active"]

        # Stack all active positions across all times
        xs = pos[active][:, 0]
        ys = pos[active][:, 1]

        # 2D histogram, log-normalized for visibility
        h, xe, ye = np.histogram2d(
            xs, ys,
            bins=80,
            range=[x_bounds, y_bounds],
        )
        # Log scale (add 1 to avoid log(0))
        im = ax.imshow(np.log10(h.T + 1),
                       origin="lower",
                       extent=[xe[0] * 1e6, xe[-1] * 1e6, ye[0] * 1e6, ye[-1] * 1e6],
                       cmap="viridis",
                       aspect="auto")
        ax.set_title(CASES[key]["label"])
        ax.set_xlabel("x (μm)")
        ax.set_ylabel("y (μm)")
        fig.colorbar(im, ax=ax, label="log10(count+1)")
        idx += 1

    fig.suptitle("Particle Density (integrated over full trajectory)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[PLOT] Saved {out_path}")


def plot_escape_fraction(cases_data, out_path):
    """What fraction of particles are still in the domain at each time?
    Case A: advected out by flow, escapes fast.
    Case B: no flow, escapes slowly if at all."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, data in cases_data.items():
        if data is None:
            continue
        t = data["time"]
        active = data["active"]
        frac_active = active.sum(axis=1) / active.shape[1]
        ax.plot(t, frac_active,
                label=CASES[key]["label"],
                color=CASES[key]["color"],
                linewidth=2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Fraction of particles in domain")
    ax.set_ylim(0, 1.05)
    ax.set_title("Retention: How long do particles stay in the domain?")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[PLOT] Saved {out_path}")


def plot_pillar_proximity(cases_data, out_path, stl_path=None):
    """Distribution of distances from bacteria to nearest pillar surface.
    If pillars attract bacteria (no-slip, boundary layer trapping), we expect
    the A2 distribution to be shifted toward smaller distances vs A1."""
    if stl_path is None:
        stl_path = "constant/triSurface/pillars.stl"
    if not os.path.exists(stl_path):
        print(f"[WARN] STL not found at {stl_path}, skipping pillar proximity plot")
        return

    try:
        import pyvista as pv
    except ImportError:
        print("[WARN] pyvista not available, skipping pillar proximity plot")
        return

    stl = pv.read(stl_path)
    # STL is in mm by default, OpenFOAM is in m - need to match units
    # Check scale: if bounding box > 1, probably mm; scale down
    bbox = stl.bounds
    max_extent = max(abs(bbox[1] - bbox[0]),
                     abs(bbox[3] - bbox[2]),
                     abs(bbox[5] - bbox[4]))
    if max_extent > 0.01:  # More than 1 cm -> probably mm
        print(f"[INFO] Scaling STL by 0.001 (mm -> m)")
        stl = stl.scale(0.001, inplace=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    for key, data in cases_data.items():
        if data is None:
            continue
        # Sample up to 50000 positions to keep it fast
        pos = data["positions"][data["active"]]
        if len(pos) > 50000:
            idx = np.random.choice(len(pos), 50000, replace=False)
            pos = pos[idx]

        # Compute distance from each position to nearest point on STL
        points = pv.PolyData(pos)
        try:
            dists = points.compute_implicit_distance(stl)["implicit_distance"]
        except Exception:
            # Older pyvista API
            from scipy.spatial import cKDTree
            tree = cKDTree(stl.points)
            dists, _ = tree.query(pos)

        # Absolute distance (we don't care if inside or outside)
        dists = np.abs(dists)
        # Convert to micrometers for readability
        dists_um = dists * 1e6

        ax.hist(dists_um, bins=50, alpha=0.5,
                label=CASES[key]["label"],
                color=CASES[key]["color"],
                density=True)
    ax.set_xlabel("Distance to nearest pillar surface (μm)")
    ax.set_ylabel("Density (normalized)")
    ax.set_title("Pillar Proximity Distribution")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[PLOT] Saved {out_path}")


def write_summary(cases_data, out_path):
    """Text file with numerical statistics for the lab notebook."""
    lines = ["BACTERIAL TRANSPORT STUDY - SUMMARY STATISTICS", "=" * 55, ""]
    for key, data in cases_data.items():
        if data is None:
            lines.append(f"Case {key}: NOT RUN (missing trajectories.h5)")
            continue
        t = data["time"]
        pos = data["positions"]
        active = data["active"]

        n_initial = active[0].sum()
        n_final = active[-1].sum()
        escape_frac = 1.0 - n_final / n_initial

        # Final MSD
        mask = active[-1] & active[0]
        if mask.sum() > 0:
            disp = pos[-1, mask] - pos[0, mask]
            final_msd = np.mean(np.sum(disp ** 2, axis=1))
            final_rms = np.sqrt(final_msd)
        else:
            final_msd = float("nan")
            final_rms = float("nan")

        lines.append(f"Case {key}: {CASES[key]['label']}")
        lines.append(f"  t_end           = {t[-1]:.3f} s")
        lines.append(f"  particles       = {n_initial} initial, {n_final} final")
        lines.append(f"  escape fraction = {escape_frac:.1%}")
        lines.append(f"  final MSD       = {final_msd:.3e} m²")
        lines.append(f"  final RMS disp. = {final_rms * 1e6:.2f} μm")
        lines.append("")

    # Slip vs no-slip comparison (the scientific question)
    if cases_data.get("A1") is not None and cases_data.get("A2") is not None:
        lines.append("SCIENTIFIC COMPARISON: A1 vs A2 (slip vs no-slip pillars)")
        lines.append("-" * 55)
        a1 = cases_data["A1"]
        a2 = cases_data["A2"]
        esc_a1 = 1.0 - a1["active"][-1].sum() / a1["active"][0].sum()
        esc_a2 = 1.0 - a2["active"][-1].sum() / a2["active"][0].sum()
        lines.append(f"  Escape fraction (A1 slip):    {esc_a1:.1%}")
        lines.append(f"  Escape fraction (A2 noSlip):  {esc_a2:.1%}")
        if esc_a1 > esc_a2:
            lines.append(f"  -> Slip pillars allow more particles to escape downstream.")
            lines.append(f"     No-slip boundary layer is retaining particles.")
        elif esc_a2 > esc_a1:
            lines.append(f"  -> No-slip pillars allow more particles to escape.")
            lines.append(f"     Unexpected: investigate further.")
        else:
            lines.append(f"  -> No measurable difference in retention.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[SUMMARY] Saved {out_path}")
    print()
    print("\n".join(lines))


def main():
    # Load all available cases
    cases_data = {key: load_case(key) for key in CASES}

    if all(v is None for v in cases_data.values()):
        print("[ERROR] No cases found. Run the comparison study from the GUI first.")
        sys.exit(1)

    out_dir = Path("analysis_output")
    out_dir.mkdir(exist_ok=True)

    plot_msd(cases_data, out_dir / "msd_comparison.png")
    plot_density_heatmaps(cases_data, out_dir / "density_heatmaps.png")
    plot_escape_fraction(cases_data, out_dir / "escape_fraction.png")
    plot_pillar_proximity(cases_data, out_dir / "pillar_proximity.png")
    write_summary(cases_data, out_dir / "summary_stats.txt")

    print(f"\n[DONE] All plots saved to {out_dir}/")


if __name__ == "__main__":
    main()
