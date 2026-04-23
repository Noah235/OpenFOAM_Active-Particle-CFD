"""
particle_tracker.py

Run-and-tumble particle tracker for E. coli in a pre-computed OpenFOAM fluid field.

Physics:
    u_particle(t) = u_fluid(x_particle(t)) + v_swim * heading(t)

The particle advects with the local fluid velocity AND swims at v_swim in its
current heading direction. Heading changes at random 'tumble' events.

Run-and-tumble rules (Berg model for E. coli):
    - Between tumbles, heading is constant (straight 'run').
    - Time between tumbles is exponentially distributed with mean tau_run.
    - At each tumble, heading rotates by a random angle drawn from a
      gamma distribution peaked at ~68 deg (Berg & Brown 1972).
    - Tumble duration is treated as instantaneous (dt_tumble << dt_run).

Inputs:
    - A VTK/OpenFOAM case containing a 3D fluid velocity field 'U'.
    - We take the LAST (converged) timestep - assumes steady-state fluid.

Outputs:
    - VTK polydata file per timestep (for ParaView visualization)
    - HDF5 file with full trajectory data (for analysis)
"""

import numpy as np
import pyvista as pv
import argparse
import json
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Default E. coli parameters (Berg, "E. coli in Motion", 2004)
# ---------------------------------------------------------------------------
DEFAULTS = {
    "n_particles": 2000,
    "v_swim": 20e-6,          # m/s, typical E. coli swim speed
    "tau_run": 1.0,           # s, mean run duration (exponential)
    "tumble_angle_mean": 68.0,  # degrees, mean tumble angle
    "tumble_angle_std": 35.0,   # degrees, spread
    "dt": 1e-3,               # s, tracker timestep (1 ms is fine for this problem)
    "t_end": 2.0,             # s, total simulation time
    "write_interval": 0.02,   # s, how often to write a VTK snapshot
    "seed": 42,

    # Mode selector: "openfoam" or "zero_fluid"
    # - "openfoam": load fluid from case, interpolate velocities (advection + swim + tumble)
    # - "zero_fluid": no fluid, pure run-and-tumble in free space
    "mode": "openfoam",

    # Injection region (used for both modes; defaults match near-inlet zone)
    "inject_xmin": -0.00019,
    "inject_xmax": -0.00017,
    "inject_ymin": -0.0001,
    "inject_ymax": 0.0001,
    "inject_zmin": -4e-5,
    "inject_zmax": 4e-5,

    # Domain bounds for zero_fluid mode (particles that leave this box deactivate)
    # Only used when mode == "zero_fluid". For "openfoam" mode we use mesh bounds.
    "domain_xmin": -0.0002,
    "domain_xmax": 0.0002,
    "domain_ymin": -0.000125,
    "domain_ymax": 0.000125,
    "domain_zmin": -5e-5,
    "domain_zmax": 5e-5,
}


# ---------------------------------------------------------------------------
# Fluid field loader
# ---------------------------------------------------------------------------
def load_steady_fluid(case_dir, field_name="U"):
    """Load the last timestep from an OpenFOAM case as a pyvista dataset.

    We use pyvista's OpenFOAMReader which wraps VTK's reader.

    Returns:
        mesh: pyvista.UnstructuredGrid with 'U' as a point array.
    """
    foam_file = Path(case_dir) / "case.foam"
    if not foam_file.exists():
        foam_file.touch()

    reader = pv.POpenFOAMReader(str(foam_file))
    # Use last timestep (steady state)
    times = reader.time_values
    if not times:
        raise RuntimeError(f"No timesteps found in {case_dir}")
    reader.set_active_time_value(times[-1])
    print(f"[fluid] Loading t={times[-1]} from {case_dir}")

    # The reader returns a MultiBlock; for an OpenFOAM case the
    # 'internalMesh' block holds the volume field.
    multi = reader.read()
    internal = multi["internalMesh"]

    # pyvista loads cell-centered data; we need point data for interpolation.
    if field_name in internal.cell_data:
        internal = internal.cell_data_to_point_data()
    if field_name not in internal.point_data:
        raise RuntimeError(
            f"Field '{field_name}' not found. "
            f"Available: {list(internal.point_data.keys())}"
        )

    print(f"[fluid] Loaded {internal.n_points} points, {internal.n_cells} cells")
    U = internal.point_data[field_name]
    print(f"[fluid] U magnitude: min={np.linalg.norm(U, axis=1).min():.3e}, "
          f"max={np.linalg.norm(U, axis=1).max():.3e} m/s")
    return internal


# ---------------------------------------------------------------------------
# Velocity interpolator
# ---------------------------------------------------------------------------
class FluidSampler:
    """Wraps VTK's probe filter for interpolating U at arbitrary points.

    VTK handles the unstructured mesh, cell-finding, and trilinear
    interpolation. We just hand it a point array and get velocities back.

    Points outside the fluid domain (in pillars, or outside the box) return
    zero velocity - which is what we want. A particle that finds itself
    inside a pillar won't be advected further into it.
    """

    def __init__(self, mesh, field_name="U"):
        self.mesh = mesh
        self.field_name = field_name

    def sample(self, points):
        """points: (N, 3) array. Returns (N, 3) velocity array."""
        probe = pv.PolyData(points)
        probed = probe.sample(self.mesh)
        U = np.asarray(probed.point_data[self.field_name])
        # pyvista marks missed points; zero them to be safe
        valid = np.asarray(probed.point_data.get("vtkValidPointMask", np.ones(len(points), dtype=bool)))
        U[~valid.astype(bool)] = 0.0
        return U

    def in_fluid(self, points):
        """Return boolean array: True if point is inside the fluid mesh."""
        probe = pv.PolyData(points)
        probed = probe.sample(self.mesh)
        valid = np.asarray(probed.point_data.get("vtkValidPointMask", np.ones(len(points), dtype=bool)))
        return valid.astype(bool)


class CollisionChecker:
    """Mesh-only variant: checks whether points are inside the fluid mesh.
    Used in no_flow_with_pillars mode - particles outside fluid mesh have
    entered a pillar."""

    def __init__(self, mesh):
        self.mesh = mesh

    def sample(self, points):
        """No advection in this mode - return zeros."""
        return np.zeros_like(points)

    def in_fluid(self, points):
        """Return boolean array: True if point is inside the fluid mesh."""
        probe = pv.PolyData(points)
        probed = probe.sample(self.mesh)
        valid = np.asarray(probed.point_data.get("vtkValidPointMask", np.ones(len(points), dtype=bool)))
        return valid.astype(bool)


def load_mesh_only(case_dir):
    """Load the OpenFOAM mesh without requiring a fluid solve.

    We use the 0/ timestep (which exists from case setup) so we have a
    geometry to query even if icoFoam was never run.
    """
    foam_file = Path(case_dir) / "case.foam"
    if not foam_file.exists():
        foam_file.touch()

    reader = pv.POpenFOAMReader(str(foam_file))
    times = reader.time_values
    if not times:
        raise RuntimeError(
            f"No timesteps found in {case_dir}. Run blockMesh + snappyHexMesh first."
        )
    # Use earliest timestep - we only need the geometry
    reader.set_active_time_value(times[0])
    print(f"[mesh] Loading geometry from t={times[0]}")

    multi = reader.read()
    internal = multi["internalMesh"]
    print(f"[mesh] Loaded {internal.n_points} points, {internal.n_cells} cells")
    print(f"[mesh] Bounds: {internal.bounds}")
    return internal


# ---------------------------------------------------------------------------
# Tumble angle sampler
# ---------------------------------------------------------------------------
def sample_tumble_angles(n, mean_deg, std_deg, rng):
    """Sample tumble angles from a gamma distribution peaked at mean_deg.

    Gamma is a reasonable fit to Berg & Brown's empirical distribution.
    Shape k and scale theta from method of moments: k = (mean/std)^2,
    theta = std^2/mean. Returns angles in radians.
    """
    mean = np.deg2rad(mean_deg)
    std = np.deg2rad(std_deg)
    k = (mean / std) ** 2
    theta = std ** 2 / mean
    angles = rng.gamma(k, theta, size=n)
    # Clip to [0, pi] since angles beyond pi are nonsensical on a sphere
    return np.clip(angles, 0.0, np.pi)


def rotate_headings(headings, angles, rng):
    """Rotate each heading by its corresponding angle about a random axis
    perpendicular to itself.

    For each heading h, we pick a random unit vector u perpendicular to h,
    then rotate h about u by angle.
    Uses Rodrigues' rotation formula.
    """
    n = len(headings)
    # Build a random axis perpendicular to each heading
    # Trick: take a random 3-vector, subtract its projection onto heading
    random_vec = rng.normal(size=(n, 3))
    proj = np.sum(random_vec * headings, axis=1, keepdims=True) * headings
    axis = random_vec - proj
    axis_norm = np.linalg.norm(axis, axis=1, keepdims=True)
    # Safety: if we got unlucky and axis is near-zero, fall back to any perp vector
    bad = (axis_norm < 1e-10).flatten()
    if bad.any():
        # A vector not parallel to heading: swap the smallest component
        fallback = np.zeros((bad.sum(), 3))
        fallback[:, 0] = 1.0
        parallel = np.abs(headings[bad, 0]) > 0.9
        fallback[parallel, 0] = 0.0
        fallback[parallel, 1] = 1.0
        axis[bad] = np.cross(headings[bad], fallback)
        axis_norm[bad] = np.linalg.norm(axis[bad], axis=1, keepdims=True)
    axis = axis / axis_norm

    # Rodrigues: h_new = h cos(a) + (axis x h) sin(a) + axis (axis . h) (1-cos(a))
    # The last term vanishes because axis is perpendicular to h by construction.
    cos_a = np.cos(angles)[:, None]
    sin_a = np.sin(angles)[:, None]
    cross = np.cross(axis, headings)
    new_headings = headings * cos_a + cross * sin_a
    # Renormalize to kill floating-point drift
    new_headings /= np.linalg.norm(new_headings, axis=1, keepdims=True)
    return new_headings


def random_unit_sphere(n, rng):
    """Draw n unit vectors uniformly on the sphere."""
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# Main tracker
# ---------------------------------------------------------------------------
def run_tracker(case_dir, params, output_dir):
    rng = np.random.default_rng(params["seed"])

    mode = params.get("mode", "openfoam")
    print(f"[tracker] mode = {mode}")

    # --- Load fluid (or not) ---
    if mode == "openfoam":
        mesh = load_steady_fluid(case_dir, "U")
        sampler = FluidSampler(mesh, "U")
        bounds = mesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
        collision_check = True  # Particles outside fluid mesh = in pillar
    elif mode == "no_flow_with_pillars":
        # Load mesh for geometry (pillar detection) but NO fluid advection.
        # This isolates pure run-and-tumble behavior within your pillar array.
        print("[tracker] Mesh-only mode: pillar collisions on, no advection")
        mesh = load_mesh_only(case_dir)
        sampler = CollisionChecker(mesh)  # Used only for 'is this in a pillar?'
        bounds = mesh.bounds
        collision_check = True
    elif mode == "zero_fluid":
        print("[tracker] Zero-fluid mode: no advection, no geometry, free space")
        sampler = None
        bounds = (
            params["domain_xmin"], params["domain_xmax"],
            params["domain_ymin"], params["domain_ymax"],
            params["domain_zmin"], params["domain_zmax"],
        )
        collision_check = False
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'openfoam', 'no_flow_with_pillars', or 'zero_fluid'.")

    # --- Initialize particles ---
    n = params["n_particles"]
    x = np.column_stack([
        rng.uniform(params["inject_xmin"], params["inject_xmax"], n),
        rng.uniform(params["inject_ymin"], params["inject_ymax"], n),
        rng.uniform(params["inject_zmin"], params["inject_zmax"], n),
    ])
    heading = random_unit_sphere(n, rng)

    # Next-tumble time for each particle, exponentially distributed
    time_to_tumble = rng.exponential(params["tau_run"], size=n)

    # Particle IDs for tracking
    ids = np.arange(n)
    active = np.ones(n, dtype=bool)  # flag: still in the domain?

    # --- Storage: full trajectory as a list of arrays ---
    # We write per-timestep VTK for ParaView and keep trajectories in memory
    # for a final HDF5 dump.
    dt = params["dt"]
    t_end = params["t_end"]
    write_dt = params["write_interval"]
    n_steps = int(t_end / dt)
    steps_per_write = max(1, int(write_dt / dt))

    # Domain bounds set above (either from mesh or config)
    print(f"[domain] bounds: {bounds}")

    trajectory_snapshots = []

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "vtk"), exist_ok=True)

    print(f"[tracker] n_particles={n}, dt={dt}, t_end={t_end}")
    print(f"[tracker] v_swim={params['v_swim']} m/s, tau_run={params['tau_run']} s")

    write_count = 0
    t = 0.0
    for step in range(n_steps):
        # Only advance active particles
        if not active.any():
            print(f"[tracker] All particles escaped at t={t:.3f}s; stopping.")
            break

        idx = np.where(active)[0]
        x_active = x[idx]

        # Sample fluid velocity at active particle positions (or zero it)
        if sampler is not None:
            u_fluid = sampler.sample(x_active)
        else:
            u_fluid = np.zeros_like(x_active)

        # Total velocity = advection + self-propulsion
        u_total = u_fluid + params["v_swim"] * heading[idx]

        # Tentative Euler step
        x_new = x_active + u_total * dt

        # Pillar collision check: if the proposed position is outside the
        # fluid mesh, the particle just swam into a pillar.
        # Action: revert to the old position AND force a tumble.
        # This models bacteria reorienting upon contact with a surface.
        if collision_check:
            in_fluid = sampler.in_fluid(x_new)
            collided = ~in_fluid
            if collided.any():
                # Revert positions for collided particles
                x_new[collided] = x_active[collided]
                # Force immediate tumble for collided particles
                # (Berg's "contact tumble" - reorient on surface interaction)
                collided_idx = idx[collided]
                if len(collided_idx) > 0:
                    # Use uniform-on-sphere for contact tumbles (full reorientation,
                    # not the narrower gamma distribution used for normal tumbles)
                    heading[collided_idx] = random_unit_sphere(len(collided_idx), rng)
                    time_to_tumble[collided_idx] = rng.exponential(
                        params["tau_run"], size=len(collided_idx)
                    )

        # Commit the (possibly reverted) positions
        x[idx] = x_new

        # Update tumble clock (normal tumbles)
        time_to_tumble[idx] -= dt
        tumbling = idx[time_to_tumble[idx] <= 0]
        if len(tumbling) > 0:
            angles = sample_tumble_angles(
                len(tumbling),
                params["tumble_angle_mean"],
                params["tumble_angle_std"],
                rng,
            )
            heading[tumbling] = rotate_headings(heading[tumbling], angles, rng)
            # Reset tumble clock for these particles
            time_to_tumble[tumbling] = rng.exponential(
                params["tau_run"], size=len(tumbling)
            )

        # Deactivate particles that left the domain
        out_of_box = (
            (x[:, 0] < bounds[0]) | (x[:, 0] > bounds[1]) |
            (x[:, 1] < bounds[2]) | (x[:, 1] > bounds[3]) |
            (x[:, 2] < bounds[4]) | (x[:, 2] > bounds[5])
        )
        active &= ~out_of_box

        t += dt

        # Write snapshot
        if step % steps_per_write == 0:
            write_vtk_snapshot(
                os.path.join(output_dir, "vtk", f"particles_{write_count:05d}.vtp"),
                x, heading, ids, active, t
            )
            trajectory_snapshots.append({
                "t": t,
                "x": x.copy(),
                "heading": heading.copy(),
                "active": active.copy(),
            })
            write_count += 1
            n_active = active.sum()
            print(f"[tracker] t={t:.3f}s, active={n_active}/{n}")

    # Final write
    write_vtk_snapshot(
        os.path.join(output_dir, "vtk", f"particles_{write_count:05d}.vtp"),
        x, heading, ids, active, t
    )
    trajectory_snapshots.append({
        "t": t, "x": x.copy(), "heading": heading.copy(), "active": active.copy(),
    })

    # Write HDF5 for analysis
    write_trajectory_hdf5(
        os.path.join(output_dir, "trajectories.h5"),
        trajectory_snapshots,
        ids,
    )

    # Write ParaView time-series file
    write_pvd(os.path.join(output_dir, "particles.pvd"), write_count + 1, write_dt)

    print(f"[tracker] Done. Wrote {write_count + 1} snapshots to {output_dir}")


def write_vtk_snapshot(filepath, x, heading, ids, active, t):
    """Write one timestep of particles as VTK polydata.

    ParaView can read a series of these as an animated particle field.
    We only write ACTIVE particles (the ones still in the domain).
    """
    mask = active
    if not mask.any():
        # Write an empty polydata so ParaView doesn't choke
        poly = pv.PolyData(np.zeros((0, 3)))
    else:
        poly = pv.PolyData(x[mask])
        poly.point_data["heading"] = heading[mask]
        poly.point_data["id"] = ids[mask]
        poly.point_data["time"] = np.full(mask.sum(), t)
    poly.save(filepath)


def write_pvd(filepath, n_snapshots, write_dt):
    """Write a ParaView .pvd index file pointing to the VTP series."""
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
             '  <Collection>']
    for i in range(n_snapshots):
        t = i * write_dt
        lines.append(f'    <DataSet timestep="{t}" group="" part="0" file="vtk/particles_{i:05d}.vtp"/>')
    lines.append('  </Collection>')
    lines.append('</VTKFile>')
    with open(filepath, "w") as f:
        f.write("\n".join(lines))


def write_trajectory_hdf5(filepath, snapshots, ids):
    """Dump all trajectory snapshots to HDF5 for offline analysis."""
    try:
        import h5py
    except ImportError:
        print("[warn] h5py not installed; skipping HDF5 output.")
        return
    with h5py.File(filepath, "w") as f:
        n_snap = len(snapshots)
        n_part = len(ids)
        times = np.array([s["t"] for s in snapshots])
        positions = np.stack([s["x"] for s in snapshots], axis=0)  # (T, N, 3)
        headings = np.stack([s["heading"] for s in snapshots], axis=0)
        active = np.stack([s["active"] for s in snapshots], axis=0)
        f.create_dataset("time", data=times)
        f.create_dataset("positions", data=positions, compression="gzip")
        f.create_dataset("headings", data=headings, compression="gzip")
        f.create_dataset("active", data=active, compression="gzip")
        f.create_dataset("ids", data=ids)
    print(f"[tracker] Wrote HDF5: {filepath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, help="OpenFOAM case directory")
    parser.add_argument("--output", default="particle_output",
                        help="Output directory for tracks")
    parser.add_argument("--config", default=None,
                        help="JSON config overriding defaults")
    args = parser.parse_args()

    params = DEFAULTS.copy()
    if args.config and os.path.exists(args.config):
        with open(args.config) as f:
            params.update(json.load(f))

    run_tracker(args.case, params, args.output)


if __name__ == "__main__":
    main()
