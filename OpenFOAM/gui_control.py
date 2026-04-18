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
        self.root.title("E. coli Lab: Centered Origin Pipeline (v2512)")
        self.root.geometry("600x950")

        # MASTER PARAMETERS LIST
        self.params = {
            "stl_path": tk.StringVar(value="Select STL..."),
            "L": tk.StringVar(value="0.0004"),
            "W": tk.StringVar(value="0.00025"),
            "H": tk.StringVar(value="0.0001"),
            "dx": tk.StringVar(value="5e-6"),
            
            "uFluid": tk.StringVar(value="0.001"),
            "rhoFluid": tk.StringVar(value="1000"),
            "nuFluid": tk.StringVar(value="1e-06"),
            
            "dPart": tk.StringVar(value="1e-06"),
            "rhoPart": tk.StringVar(value="1100"),
            "sigma": tk.StringVar(value="1e-08"),
            "slip": tk.StringVar(value="0.8"),
            
            "endTime": tk.StringVar(value="1.0")
        }
        
        self.load_settings()
        self.create_widgets()

    def create_widgets(self):
        ttk.Label(self.root, text="OpenFOAM 2512: Master Control", font=("Arial", 14, "bold")).pack(pady=10)

        # 1. Geometry & Mesh
        geo_frame = ttk.LabelFrame(self.root, text=" Geometry & Mesh Calibration (Centered) ")
        geo_frame.pack(padx=20, pady=5, fill="x")
        self.add_entry(geo_frame, "Length (L):", self.params["L"])
        self.add_entry(geo_frame, "Width (W):", self.params["W"])
        self.add_entry(geo_frame, "Height (H):", self.params["H"])
        self.add_entry(geo_frame, "Cell Size (dx):", self.params["dx"])

        # 2. Fluid Properties
        fluid_frame = ttk.LabelFrame(self.root, text=" Fluid Properties ")
        fluid_frame.pack(padx=20, pady=5, fill="x")
        self.add_entry(fluid_frame, "Inlet Velocity (U):", self.params["uFluid"])
        self.add_entry(fluid_frame, "Fluid Density (rhoInf):", self.params["rhoFluid"])
        self.add_entry(fluid_frame, "Viscosity (nu):", self.params["nuFluid"])

        # 3. Particle Physics
        phys_frame = ttk.LabelFrame(self.root, text=" E. coli Physics ")
        phys_frame.pack(padx=20, pady=5, fill="x")
        self.add_entry(phys_frame, "Particle Dia (dPart):", self.params["dPart"])
        self.add_entry(phys_frame, "Particle Dens (rho0):", self.params["rhoPart"])
        self.add_entry(phys_frame, "Active Diffusion (Sigma):", self.params["sigma"])
        self.add_entry(phys_frame, "Pillar Slip (0-1):", self.params["slip"])

        # 4. STL Selection
        stl_frame = ttk.LabelFrame(self.root, text=" STL Geometry ")
        stl_frame.pack(padx=20, pady=5, fill="x")
        ttk.Entry(stl_frame, textvariable=self.params["stl_path"]).pack(side="left", expand=True, fill="x", padx=5)
        ttk.Button(stl_frame, text="Browse", command=self.browse_stl).pack(side="right")

        # 5. Buttons (Save Settings Retained)
        btn_frame = ttk.Frame(self.root)
        btn_frame.pack(pady=10)
        ttk.Button(btn_frame, text="SAVE SETTINGS", command=self.save_settings).grid(row=0, column=0, padx=5)
        ttk.Button(btn_frame, text="RUN PIPELINE", command=self.run_simulation).grid(row=0, column=1, padx=5)
        ttk.Button(btn_frame, text="ParaView", command=self.open_paraview).grid(row=0, column=2, padx=5)

        self.log_text = tk.Text(self.root, height=12, width=70, font=("Courier", 9), bg="#f0f0f0")
        self.log_text.pack(padx=20, pady=10)

    def add_entry(self, parent, label, var):
        frame = ttk.Frame(parent)
        frame.pack(fill="x", padx=10, pady=2)
        ttk.Label(frame, text=label, width=25).pack(side="left")
        ttk.Entry(frame, textvariable=var).pack(side="right", expand=True, fill="x")

    def browse_stl(self):
        file_path = filedialog.askopenfilename(initialdir=os.getcwd(), filetypes=[("STL files", "*.stl")])
        if file_path: self.params["stl_path"].set(file_path)
        
    def save_settings(self):
        data = {k: v.get() for k, v in self.params.items()}
        with open("gui_settings.json", "w") as f:
            json.dump(data, f)
        self.log(">>> [SYSTEM] Settings saved locally.")

    def load_settings(self):
        if os.path.exists("gui_settings.json"):
            with open("gui_settings.json", "r") as f:
                data = json.load(f)
                for k, v in data.items():
                    if k in self.params: self.params[k].set(v)

    def log(self, message):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def update_configs(self):
        try:
            L, W, H, dx = float(self.params["L"].get()), float(self.params["W"].get()), float(self.params["H"].get()), float(self.params["dx"].get())
            nx, ny, nz = int(L/dx), int(W/dx), int(max(1, H/(dx*2)))

            # Origin Centered Coordinates
            xmin, xmax = -L/2, L/2
            ymin, ymax = -W/2, W/2
            zmin, zmax = -H/2, H/2
            
            # Safe Point: 5% into the channel from the left inlet, dead center in Y and Z
            safe_x = xmin + (L * 0.05)
            safe_y = 0.0
            safe_z = 0.0

            # 1. Update blockMesh (Centered Bounds)
            self.update_file("system/blockMeshDict", {
                r"nx\s+\d+;": f"nx  {nx};", r"ny\s+\d+;": f"ny  {ny};", r"nz\s+\d+;": f"nz  {nz};",
                r"xmin\s+[^\n]+;": f"xmin  {xmin};", r"xmax\s+[^\n]+;": f"xmax  {xmax};",
                r"ymin\s+[^\n]+;": f"ymin  {ymin};", r"ymax\s+[^\n]+;": f"ymax  {ymax};",
                r"zmin\s+[^\n]+;": f"zmin  {zmin};", r"zmax\s+[^\n]+;": f"zmax  {zmax};"
            })

            # 2. Update snappyCenter
            self.update_file("system/snappyHexMeshDict", {r"locationInMesh\s+\(.*\);": f"locationInMesh ({safe_x} {safe_y} {safe_z});"})
            
            # 3. Write Strict Header kinematicCloudPositions
            positions_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2512                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant";
    object      kinematicCloudPositions;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

(
({safe_x} {safe_y} {safe_z})
)
"""
            with open("constant/kinematicCloudPositions", "w") as f: f.write(positions_content)

            # 4. Update Fluid Properties
            self.update_file("constant/transportProperties", {
                r"rhoInf\s+[\d.e-]+;": f"rhoInf          {self.params['rhoFluid'].get()};",
                r"nu\s+nu\s+\[.*\]\s+[\d.e-]+;": f"nu              nu [0 2 -1 0 0 0 0] {self.params['nuFluid'].get()};"
            })

            # 5. Update Particle Cloud Parameters
            self.update_file("constant/kinematicCloudProperties", {
                r"rho0Val\s+[\d.e-]+;": f"rho0Val     {self.params['rhoPart'].get()};",
                r"dPart\s+[\d.e-]+;": f"dPart       {self.params['dPart'].get()};",
                r"sigmaVal\s+[\d.e-]+;": f"sigmaVal    {self.params['sigma'].get()};",
                r"eRest\s+[\d.e-]+;": f"eRest       {self.params['slip'].get()};"
            })

            # 6. Local STL Copy
            stl_src = self.params["stl_path"].get()
            if os.path.exists(stl_src): shutil.copy(stl_src, "constant/triSurface/pillars.stl")

            self.log(">>> [CONFIG] Centered geometry math applied.")
        except Exception as e: messagebox.showerror("Error", str(e))

    def update_file(self, filepath, replacements):
        if not os.path.exists(filepath): return
        with open(filepath, 'r') as f: content = f.read()
        for pattern, replacement in replacements.items(): content = re.sub(pattern, replacement, content)
        with open(filepath, 'w') as f: f.write(content)

    def run_simulation(self):
        self.update_configs()
        self.log(">>> [RUN] Starting Pipeline...")
        process = subprocess.Popen(['bash', './Allrun'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout: self.log(line.strip())
        process.wait()
        self.log(">>> [FINISH] Done.")

    def open_paraview(self):
        try: subprocess.Popen(['paraFoam'])
        except: subprocess.Popen(['paraview', 'case.foam'])

if __name__ == "__main__":
    root = tk.Tk()
    app = OpenFOAMGui(root)
    root.mainloop()