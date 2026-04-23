"""
E. coli Lab: OpenFOAM Fluid + Python Run-and-Tumble Pipeline
=============================================================

Architecture:
    [GUI] -> writes OpenFOAM dicts + particle_config.json
    [Allrun] -> icoFoam (fluid) -> particle_tracker.py (E. coli)
    [ParaView] -> view case.foam + particle_output/particles.pvd together
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import subprocess
import os
import re
import shutil
import json


class OpenFOAMGui:
    def __init__(self, root):
        self.root = root
        self.root.title("E. coli Lab: Fluid + Run-and-Tumble (v2512 + Python)")
        self.root.geometry("620x1000")

        # All parameters
        self.params = {
            # Geometry
            "stl_path": tk.StringVar(value="Select STL..."),
            "L": tk.StringVar(value="0.0004"),
            "W": tk.StringVar(value="0.00025"),
            "H": tk.StringVar(value="0.0001"),
            "dx": tk.StringVar(value="5e-6"),

            # Fluid
            "uFluid": tk.StringVar(value="0.001"),
            "rhoFluid": tk.StringVar(value="1000"),
            "nuFluid": tk.StringVar(value="1e-06"),
            "endTimeFluid": tk.StringVar(value="0.5"),  # Long enough for steady state

            # E. coli (run-and-tumble)
            "n_particles": tk.StringVar(value="2000"),
            "v_swim": tk.StringVar(value="20e-6"),
            "tau_run": tk.StringVar(value="1.0"),
            "tumble_angle_mean": tk.StringVar(value="68.0"),
            "tumble_angle_std": tk.StringVar(value="35.0"),
            "tracker_dt": tk.StringVar(value="1e-3"),
            "tracker_t_end": tk.StringVar(value="2.0"),
            "write_interval": tk.StringVar(value="0.02"),
            "seed": tk.StringVar(value="42"),

            "force_mesh": tk.BooleanVar(value=False),
            "skip_particles": tk.BooleanVar(value=False),
            "zero_fluid": tk.BooleanVar(value=False),
            "no_flow_with_pillars": tk.BooleanVar(value=False),
        }

        self.last_mesh_state = ""
        self.load_settings()
        self.create_widgets()

    def create_widgets(self):
        ttk.Label(self.root,
                  text="Fluid (OpenFOAM) + E. coli (Python)",
                  font=("Arial", 14, "bold")).pack(pady=10)

        # Geometry & Mesh
        geo = ttk.LabelFrame(self.root, text=" Geometry & Mesh ")
        geo.pack(padx=20, pady=5, fill="x")
        self.add_entry(geo, "Length (L):", self.params["L"])
        self.add_entry(geo, "Width (W):", self.params["W"])
        self.add_entry(geo, "Height (H):", self.params["H"])
        self.add_entry(geo, "Cell Size (dx):", self.params["dx"])

        # Fluid
        fluid = ttk.LabelFrame(self.root, text=" Fluid (icoFoam) ")
        fluid.pack(padx=20, pady=5, fill="x")
        self.add_entry(fluid, "Inlet Velocity (m/s):", self.params["uFluid"])
        self.add_entry(fluid, "Fluid Density (kg/m3):", self.params["rhoFluid"])
        self.add_entry(fluid, "Viscosity nu (m2/s):", self.params["nuFluid"])
        self.add_entry(fluid, "Fluid End Time (s):", self.params["endTimeFluid"])

        # E. coli
        coli = ttk.LabelFrame(self.root, text=" E. coli (Run-and-Tumble) ")
        coli.pack(padx=20, pady=5, fill="x")
        self.add_entry(coli, "Num Particles:", self.params["n_particles"])
        self.add_entry(coli, "Swim Speed (m/s):", self.params["v_swim"])
        self.add_entry(coli, "Mean Run Duration (s):", self.params["tau_run"])
        self.add_entry(coli, "Tumble Angle Mean (deg):", self.params["tumble_angle_mean"])
        self.add_entry(coli, "Tumble Angle Std (deg):", self.params["tumble_angle_std"])
        self.add_entry(coli, "Tracker dt (s):", self.params["tracker_dt"])
        self.add_entry(coli, "Tracker End Time (s):", self.params["tracker_t_end"])
        self.add_entry(coli, "Write Interval (s):", self.params["write_interval"])
        self.add_entry(coli, "Random Seed:", self.params["seed"])

        # STL
        stl = ttk.LabelFrame(self.root, text=" STL Geometry ")
        stl.pack(padx=20, pady=5, fill="x")
        ttk.Entry(stl, textvariable=self.params["stl_path"]).pack(
            side="left", expand=True, fill="x", padx=5)
        ttk.Button(stl, text="Browse", command=self.browse_stl).pack(side="right")

        # Buttons
        btn = ttk.Frame(self.root)
        btn.pack(pady=10)
        ttk.Button(btn, text="SAVE SETTINGS", command=self.save_settings).grid(row=0, column=0, padx=5)
        ttk.Button(btn, text="RUN PIPELINE", command=self.run_simulation).grid(row=0, column=1, padx=5)
        ttk.Button(btn, text="ParaView", command=self.open_paraview).grid(row=0, column=2, padx=5)
        ttk.Button(btn, text="RUN COMPARISON STUDY (A1+A2+B)",
                   command=self.run_comparison_study).grid(row=1, column=0, columnspan=3, pady=5)
        ttk.Checkbutton(btn, text="Force Remesh",
                        variable=self.params["force_mesh"]).grid(row=2, column=0, pady=5)
        ttk.Checkbutton(btn, text="Skip Particles (fluid only)",
                        variable=self.params["skip_particles"]).grid(row=2, column=1, pady=5)
        ttk.Checkbutton(btn, text="Zero Fluid (pure run-and-tumble, no pillars)",
                        variable=self.params["zero_fluid"]).grid(row=3, column=0, columnspan=3, pady=2)
        ttk.Checkbutton(btn, text="No Flow + Pillars (tumble with collisions, no advection)",
                        variable=self.params["no_flow_with_pillars"]).grid(row=4, column=0, columnspan=3, pady=2)

        self.log_text = tk.Text(self.root, height=14, width=75,
                                font=("Courier", 9), bg="#f0f0f0")
        self.log_text.pack(padx=20, pady=10)

    def add_entry(self, parent, label, var):
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=10, pady=2)
        ttk.Label(frame, text=label, width=25).pack(side="left")
        ttk.Entry(frame, textvariable=var).pack(side="right", expand=True, fill="x")

    def browse_stl(self):
        fp = filedialog.askopenfilename(initialdir=os.getcwd(),
                                        filetypes=[("STL files", "*.stl")])
        if fp:
            self.params["stl_path"].set(fp)

    def save_settings(self):
        data = {k: str(v.get()) for k, v in self.params.items()}
        data["last_mesh_state"] = self.last_mesh_state
        with open("gui_settings.json", "w") as f:
            json.dump(data, f)
        self.log(">>> [SYSTEM] Settings saved.")

    def load_settings(self):
        self.last_mesh_state = ""
        if os.path.exists("gui_settings.json"):
            with open("gui_settings.json") as f:
                data = json.load(f)
            for k, v in data.items():
                if k in self.params:
                    if isinstance(self.params[k], tk.BooleanVar):
                        self.params[k].set(v == 'True')
                    else:
                        self.params[k].set(v)
            self.last_mesh_state = data.get("last_mesh_state", "")

    def log(self, message):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    # -----------------------------------------------------------------------
    # Config writers
    # -----------------------------------------------------------------------
    def update_configs(self):
        try:
            L = float(self.params["L"].get())
            W = float(self.params["W"].get())
            H = float(self.params["H"].get())
            dx = float(self.params["dx"].get())
            uFluid = self.params["uFluid"].get()

            nx, ny, nz = int(L/dx), int(W/dx), int(max(1, H/(dx*2)))
            xmin, xmax = -L/2, L/2
            ymin, ymax = -W/2, W/2
            zmin, zmax = -H/2, H/2

            safe_x = xmin + (L * 0.05)

            # 1. blockMesh
            self.update_file("system/blockMeshDict", {
                r"nx\s+\d+;": f"nx  {nx};",
                r"ny\s+\d+;": f"ny  {ny};",
                r"nz\s+\d+;": f"nz  {nz};",
                r"xmin\s+[^\n]+;": f"xmin  {xmin};",
                r"xmax\s+[^\n]+;": f"xmax  {xmax};",
                r"ymin\s+[^\n]+;": f"ymin  {ymin};",
                r"ymax\s+[^\n]+;": f"ymax  {ymax};",
                r"zmin\s+[^\n]+;": f"zmin  {zmin};",
                r"zmax\s+[^\n]+;": f"zmax  {zmax};",
            })

            # 2. snappyHexMesh center
            self.update_file("system/snappyHexMeshDict",
                             {r"locationInMesh\s+\(.*\);":
                              f"locationInMesh ({safe_x} 0.0 0.0);"})

            # 3. Write 0/U and 0/p from scratch (correct 3D patch types)
            self.write_fluid_U("0/U", uFluid)
            self.write_fluid_p("0/p")

            # 4. Transport properties
            self.update_file("constant/transportProperties", {
                r"rhoInf\s+\[.*\]\s+[\d.e-]+;":
                    f"rhoInf          [1 -3 0 0 0 0 0] {self.params['rhoFluid'].get()};",
                r"nu\s+\[.*\]\s+[\d.e-]+;":
                    f"nu              [0 2 -1 0 0 0 0] {self.params['nuFluid'].get()};",
            })

            # 5. controlDict - switch to icoFoam, set endTime
            self.write_control_dict(
                application="icoFoam",
                endTime=self.params["endTimeFluid"].get(),
            )

            # 6. Particle tracker config
            self.write_particle_config(xmin, xmax, ymin, ymax, zmin, zmax)

            # 7. STL
            stl_src = self.params["stl_path"].get()
            if os.path.exists(stl_src):
                os.makedirs("constant/triSurface", exist_ok=True)
                shutil.copy(stl_src, "constant/triSurface/pillars.stl")

            # 8. Clean up Lagrangian cruft (no longer needed)
            for f in ["constant/kinematicCloudProperties",
                      "constant/kinematicCloudPositions"]:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                        self.log(f">>> [CLEAN] Removed obsolete {f}")
                    except OSError:
                        pass

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_file(self, filepath, replacements):
        if not os.path.exists(filepath):
            return
        with open(filepath) as f:
            content = f.read()
        for pat, repl in replacements.items():
            content = re.sub(pat, repl, content)
        with open(filepath, "w") as f:
            f.write(content)

    def write_fluid_U(self, filepath, uFluid, pillar_bc="noSlip"):
        """Write 0/U with configurable pillar boundary condition.
        pillar_bc: either 'slip' or 'noSlip'. This is the core variable in
        our slip-vs-noslip comparison study."""
        content = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  2512                                  |
|   \  /    A nd           | Website:  www.openfoam.com                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({uFluid} 0 0);
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    walls
    {{
        type            noSlip;
    }}
    "pillars.*"
    {{
        type            {pillar_bc};
    }}
}}

// ************************************************************************* //
""".format(uFluid=uFluid, pillar_bc=pillar_bc)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            f.write(content)

    def write_fluid_p(self, filepath):
        content = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  2512                                  |
|   \  /    A nd           | Website:  www.openfoam.com                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    walls
    {
        type            zeroGradient;
    }
    "pillars.*"
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            f.write(content)

    def write_control_dict(self, application, endTime):
        content = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}

application     {app};

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {endTime};

deltaT          1e-6;

// Adaptive timestepping: solver auto-shrinks dt to keep Courant stable.
// This prevents numerical explosion when flow accelerates around pillars.
adjustTimeStep  yes;
maxCo           0.5;
maxDeltaT       1e-5;

writeControl    adjustableRunTime;
writeInterval   {writeInt};

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
compression     off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""".format(app=application,
           endTime=endTime,
           writeInt=max(float(endTime) / 20.0, 1e-3))
        with open("system/controlDict", "w") as f:
            f.write(content)

    def write_particle_config(self, xmin, xmax, ymin, ymax, zmin, zmax,
                              mode=None, t_end=None):
        """Write JSON config for particle_tracker.py.

        Injection zone: a thin slab just downstream of the inlet.

        If mode or t_end are None, use GUI-set values. The comparison runner
        passes explicit values to override per-case.
        """
        inject_xmin = xmin + (xmax - xmin) * 0.025
        inject_xmax = xmin + (xmax - xmin) * 0.075
        inject_ymargin = (ymax - ymin) * 0.1
        inject_zmargin = (zmax - zmin) * 0.1

        if mode is None:
            mode = "openfoam"
            if self.params["zero_fluid"].get():
                mode = "zero_fluid"
            elif self.params["no_flow_with_pillars"].get():
                mode = "no_flow_with_pillars"

        if t_end is None:
            t_end = float(self.params["tracker_t_end"].get())

        config = {
            "mode": mode,
            "n_particles": int(self.params["n_particles"].get()),
            "v_swim": float(self.params["v_swim"].get()),
            "tau_run": float(self.params["tau_run"].get()),
            "tumble_angle_mean": float(self.params["tumble_angle_mean"].get()),
            "tumble_angle_std": float(self.params["tumble_angle_std"].get()),
            "dt": float(self.params["tracker_dt"].get()),
            "t_end": float(t_end),
            "write_interval": float(self.params["write_interval"].get()),
            "seed": int(self.params["seed"].get()),
            "inject_xmin": inject_xmin,
            "inject_xmax": inject_xmax,
            "inject_ymin": ymin + inject_ymargin,
            "inject_ymax": ymax - inject_ymargin,
            "inject_zmin": zmin + inject_zmargin,
            "inject_zmax": zmax - inject_zmargin,
            "domain_xmin": xmin,
            "domain_xmax": xmax,
            "domain_ymin": ymin,
            "domain_ymax": ymax,
            "domain_zmin": zmin,
            "domain_zmax": zmax,
        }
        with open("particle_config.json", "w") as f:
            json.dump(config, f, indent=2)

    # -----------------------------------------------------------------------
    # Run
    # -----------------------------------------------------------------------
    def run_simulation(self):
        self.update_configs()

        state = (f"{self.params['L'].get()}_{self.params['W'].get()}_"
                 f"{self.params['H'].get()}_{self.params['dx'].get()}_"
                 f"{self.params['stl_path'].get()}")
        mesh_exists = os.path.exists("constant/polyMesh")
        force = self.params["force_mesh"].get()

        do_mesh = force or not mesh_exists or state != self.last_mesh_state
        if do_mesh:
            self.last_mesh_state = state
            self.params["force_mesh"].set(False)
            self.save_settings()
            self.log(">>> [INIT] Meshing ENABLED.")
        else:
            self.log(">>> [INIT] Meshing SKIPPED.")

        cmd = ["bash", "./Allrun"]
        if do_mesh:
            cmd.append("--mesh")
        if self.params["skip_particles"].get():
            cmd.append("--no-particles")
        if self.params["zero_fluid"].get():
            cmd.append("--zero-fluid")
        if self.params["no_flow_with_pillars"].get():
            cmd.append("--no-flow")

        self.log(">>> [RUN] Starting pipeline...")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            self.log(line.rstrip())
        proc.wait()
        self.log(">>> [FINISH] Pipeline complete.")

    def open_paraview(self):
        try:
            subprocess.Popen(["paraFoam"])
        except Exception:
            subprocess.Popen(["paraview", "case.foam"])

    # -----------------------------------------------------------------------
    # Comparison study runner (the main scientific pipeline)
    # -----------------------------------------------------------------------
    def setup_case_dir(self, case_name, pillar_bc, mode, t_end,
                       xmin, xmax, ymin, ymax, zmin, zmax):
        """Create case_NAME/ as a full OpenFOAM case with symlinked mesh/system.

        This gives each case a self-contained directory that ParaView can
        open independently, while avoiding mesh duplication.
        """
        case_dir = case_name  # Relative to working dir (cwd)

        # Blow away existing case dir (user asked for clean slate on rerun)
        if os.path.exists(case_dir):
            shutil.rmtree(case_dir)
        os.makedirs(case_dir)
        os.makedirs(os.path.join(case_dir, "0"))

        # Symlink constant/ and system/ so the case can find the mesh and configs
        # Use absolute paths for symlinks so they work regardless of cwd
        for shared in ("constant", "system"):
            src = os.path.abspath(shared)
            dst = os.path.join(case_dir, shared)
            if os.path.exists(src):
                os.symlink(src, dst)
            else:
                self.log(f">>> [STUDY] WARNING: {src} does not exist, symlink skipped")

        # Write per-case 0/U with the right pillar BC
        self.write_fluid_U(
            os.path.join(case_dir, "0", "U"),
            self.params["uFluid"].get(),
            pillar_bc=pillar_bc,
        )
        # Write 0/p (same for all cases)
        self.write_fluid_p(os.path.join(case_dir, "0", "p"))

        # Write particle config inside the case directory
        # (Allrun's tracker invocation uses --config particle_config.json
        # relative to the case dir)
        config_path_saved = "particle_config.json"
        self.write_particle_config(
            xmin, xmax, ymin, ymax, zmin, zmax,
            mode=mode, t_end=t_end,
        )
        # write_particle_config writes to cwd/particle_config.json;
        # move it into the case dir
        shutil.move(config_path_saved,
                    os.path.join(case_dir, "particle_config.json"))

        return case_dir

    def verify_pillar_bc(self, case_dir, expected_bc):
        """Read the 0/U we just wrote and confirm the pillar BC matches expectation."""
        u_path = os.path.join(case_dir, "0", "U")
        try:
            with open(u_path) as f:
                content = f.read()
            import re as re_local
            match = re_local.search(
                r'"pillars\.\*".*?type\s+(\w+);', content, re_local.DOTALL
            )
            if match:
                actual = match.group(1)
                self.log(f">>> [STUDY] Verified {case_dir}/0/U: pillar BC = '{actual}' "
                         f"(expected '{expected_bc}')")
                return actual == expected_bc
            else:
                self.log(f">>> [STUDY] WARNING: could not find pillar BC in {u_path}")
                return False
        except Exception as e:
            self.log(f">>> [STUDY] ERROR reading {u_path}: {e}")
            return False

    def run_comparison_study(self):
        """Run the full three-case comparison: A1 (slip pillars + flow),
        A2 (noSlip pillars + flow), B (no flow + pillar collisions).

        Each case runs in its own subdirectory (case_A1/, case_A2/, case_B/)
        with symlinked mesh and system configs. This lets you open any case
        in ParaView independently and preserves each case's full history.

        All three share mesh, particle count/seed/parameters, injection zone.
        Only the fluid pillar BC (A1 vs A2) and mode (A vs B) differ.
        """
        # Compute shared geometry values
        try:
            L = float(self.params["L"].get())
            W = float(self.params["W"].get())
            H = float(self.params["H"].get())
            xmin, xmax = -L/2, L/2
            ymin, ymax = -W/2, W/2
            zmin, zmax = -H/2, H/2
        except Exception as e:
            messagebox.showerror("Error", f"Bad geometry params: {e}")
            return

        # Set up shared configs in the ROOT working directory.
        # This writes: system/blockMeshDict, system/snappyHexMeshDict,
        # constant/transportProperties, system/controlDict, etc.
        # Each case_X/ will symlink to these.
        self.update_configs()

        # Clean up any stray 0/U and 0/p in the root (only case dirs should
        # have initial conditions going forward)
        for stray in ("0/U", "0/p"):
            if os.path.exists(stray):
                try:
                    os.remove(stray)
                except OSError:
                    pass

        # Decide if we need to mesh
        state = (f"{self.params['L'].get()}_{self.params['W'].get()}_"
                 f"{self.params['H'].get()}_{self.params['dx'].get()}_"
                 f"{self.params['stl_path'].get()}")
        mesh_exists = os.path.exists("constant/polyMesh")
        force = self.params["force_mesh"].get()
        do_mesh = force or not mesh_exists or state != self.last_mesh_state
        if do_mesh:
            self.last_mesh_state = state
            self.params["force_mesh"].set(False)
            self.save_settings()
            self.log(">>> [STUDY] Mesh will be generated (once, in root constant/).")
        else:
            self.log(">>> [STUDY] Reusing existing mesh.")

        # Define the three cases
        # NOTE: With fluid velocity ~0.1 um/s and swim speed ~20 um/s, swimming
        # dominates advection. Case A needs long t_end for meaningful exploration.
        cases = [
            {
                "name": "case_A1",
                "description": "Flow ON, slip pillars",
                "pillar_bc": "slip",
                "mode": "openfoam",
                "t_end": 30.0,
                "flags": [],
            },
            {
                "name": "case_A2",
                "description": "Flow ON, noSlip pillars",
                "pillar_bc": "noSlip",
                "mode": "openfoam",
                "t_end": 30.0,
                "flags": [],
            },
            {
                "name": "case_B",
                "description": "Flow OFF, particles with pillar collisions",
                "pillar_bc": "noSlip",  # irrelevant (no fluid solve)
                "mode": "no_flow_with_pillars",
                "t_end": 30.0,
                "flags": ["--no-flow"],
            },
        ]

        for i, case in enumerate(cases):
            self.log(f"\n>>> [STUDY] =============================================")
            self.log(f">>> [STUDY] {case['name']}: {case['description']}")
            self.log(f">>> [STUDY] =============================================")

            # 1. Set up the case subdirectory
            case_dir = self.setup_case_dir(
                case["name"],
                case["pillar_bc"],
                case["mode"],
                case["t_end"],
                xmin, xmax, ymin, ymax, zmin, zmax,
            )

            # 2. Verify the BC actually got written correctly
            if not self.verify_pillar_bc(case_dir, case["pillar_bc"]):
                if case["mode"] == "openfoam":  # Only care for fluid cases
                    self.log(f">>> [STUDY] ERROR: BC verification failed for {case_dir}. Aborting.")
                    return

            # 3. Run Allrun pointing at this case directory
            cmd = ["bash", "./Allrun", "--case", case["name"]]
            # Mesh only on first case (if needed)
            if do_mesh and i == 0:
                cmd.append("--mesh")
            cmd.extend(case["flags"])

            self.log(f">>> [STUDY] Command: {' '.join(cmd)}")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
            for line in proc.stdout:
                self.log(line.rstrip())
            proc.wait()

            if proc.returncode != 0:
                self.log(f">>> [STUDY] {case['name']} FAILED (exit {proc.returncode}). Aborting study.")
                return

            self.log(f">>> [STUDY] {case['name']} complete. Results in {case_dir}/")

            # After first case, no need to remesh even if do_mesh was set
            do_mesh = False

        self.log(f"\n>>> [STUDY] =============================================")
        self.log(f">>> [STUDY] Comparison study complete.")
        self.log(f">>> [STUDY] Per-case results:")
        self.log(f">>> [STUDY]   case_A1/  - Flow + slip pillars")
        self.log(f">>> [STUDY]   case_A2/  - Flow + noSlip pillars")
        self.log(f">>> [STUDY]   case_B/   - No flow + particle-pillar collisions")
        self.log(f">>> [STUDY] Run `python3 analysis.py` for comparison plots.")
        self.log(f">>> [STUDY] =============================================")


if __name__ == "__main__":
    root = tk.Tk()
    app = OpenFOAMGui(root)
    root.mainloop()
