import os

# --- INOTIFY CRASH FIX ---
# Create the Streamlit config file to permanently disable the file watcher.
# This prevents the OSError: [Errno 28] inotify watch limit reached.
os.makedirs(".streamlit", exist_ok=True)
config_path = os.path.join(".streamlit", "config.toml")
if not os.path.exists(config_path):
    with open(config_path, "w") as f:
        f.write("[server]\nfileWatcherType = \"none\"\n")
    print("CREATED STREAMLIT CONFIG. PLEASE RESTART APP IF INOTIFY ERROR PERSISTS.")

import streamlit as st
import pandas as pd
import time
import shutil
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import uuid
import re
import warnings
import glob
import math
import hashlib
from datetime import datetime

# --- SAFE 3D RENDERER IMPORT ---
try:
    import pyvista as pv
    from stpyvista import stpyvista
    import platform
    os.environ["PYVISTA_OFF_SCREEN"] = "true"
    pv.OFF_SCREEN = True
    if platform.system() == "Linux":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
    PYVISTA_AVAILABLE = True
    PV_ERROR_MSG = ""
except Exception as e:
    PYVISTA_AVAILABLE = False
    PV_ERROR_MSG = str(e)

# Import directly from the root folder
from case_builder import CaseBuilder
from runner import SimulationRunner
from stl_helper import STLHelper
import setup_template_v13_multicomponent

# --- CONFIGURATION ---
TEMPLATE_PATH = os.path.join("templates", "base_case_multicomp")
SIM_ROOT = "simulation_runs"
ARCHIVE_DIR = "Scenario_Results_Zip_Files"
STATE_FILE = "sim_state.json"
SCENARIO_BACKUP = "current_scenarios.csv"
DEFAULTS_FILE = "user_defaults.json"

# --- HAZARD LIMITS (PPM) FOR PLOTTING (Graded Alarms) ---
HAZARD_LIMITS = {
    "CH4": {"type": "flammable", "visible": 100, "low": 5000, "high": 25000, "emergency": 50000},
    "H2": {"type": "flammable", "visible": 100, "low": 4000, "high": 20000, "emergency": 40000},
    "C3H8": {"type": "flammable", "visible": 100, "low": 2100, "high": 10500, "emergency": 21000},
    "C2H6": {"type": "flammable", "visible": 100, "low": 3000, "high": 15000, "emergency": 30000},
    "C4H10": {"type": "flammable", "visible": 100, "low": 1800, "high": 9000, "emergency": 18000},
    "C2H4": {"type": "flammable", "visible": 100, "low": 2700, "high": 13500, "emergency": 27000},
    "C2H2": {"type": "flammable", "visible": 100, "low": 2500, "high": 12500, "emergency": 25000},
    "C5H12": {"type": "flammable", "visible": 100, "low": 1400, "high": 7000, "emergency": 14000},
    "C6H14": {"type": "flammable", "visible": 100, "low": 1100, "high": 5500, "emergency": 11000},
    
    "H2S": {"type": "toxic", "visible": 5, "low": 20, "high": 40, "emergency": 50},
    "CO": {"type": "toxic", "visible": 10, "low": 50, "high": 100, "emergency": 600},
    "NH3": {"type": "toxic", "visible": 10, "low": 50, "high": 100, "emergency": 150},
    "Cl2": {"type": "toxic", "visible": 0.5, "low": 1, "high": 2, "emergency": 5},
    "SO2": {"type": "toxic", "visible": 1, "low": 5, "high": 10, "emergency": 50},
    "NO2": {"type": "toxic", "visible": 1, "low": 5, "high": 10, "emergency": 10},
    "NO": {"type": "toxic", "visible": 5, "low": 25, "high": 50, "emergency": 50},
    "HCN": {"type": "toxic", "visible": 2, "low": 10, "high": 20, "emergency": 25},
    "C6H6": {"type": "toxic", "visible": 0.5, "low": 1, "high": 5, "emergency": 250},
    "C7H8": {"type": "toxic", "visible": 50, "low": 200, "high": 400, "emergency": 250},
    "CH3OH": {"type": "toxic", "visible": 50, "low": 200, "high": 400, "emergency": 3000},
}

# --- EXTENDED COLOR MAP ---
RESIDUAL_COLOR_MAP = {
    "p_rgh": "#1f77b4", "p": "#1f77b4",
    "Ux": "#d62728", "Uy": "#ff7f0e", "Uz": "#9467bd",
    "k": "#2ca02c", "epsilon": "#8c564b", "omega": "#8c564b",
    "h": "#e377c2", "T": "#e377c2",
    "Air": "#000000", "CH4": "#bcbd22", "C3H8": "#17becf",
    "H2": "#ff9896", "NH3": "#98df8a", "CO2": "#c5b0d5",
    "C2H6": "#a6cee3", "C4H10": "#b2df8a", "C2H4": "#fb9a99", "C2H2": "#fdbf6f",
    "C5H12": "#cab2d6", "C6H14": "#ffff99", "H2S": "#b15928", "CO": "#8dd3c7",
    "Cl2": "#bebada", "SO2": "#fb8072", "NO2": "#80b1d3", "NO": "#fccde5",
    "HCN": "#d9d9d9", "C6H6": "#bc80bd", "C7H8": "#ccebc5", "CH3OH": "#ffed6f"
}

for d in [SIM_ROOT, ARCHIVE_DIR]:
    if not os.path.exists(d): os.makedirs(d)

builder = CaseBuilder(TEMPLATE_PATH)
runner = SimulationRunner(solver="foamRun")
stl_helper = STLHelper()

st.set_page_config(page_title="OpenFOAM Auto-CFD", layout="wide", page_icon="🌊")
st.title("🌊 OpenFOAM v13: Gas Dispersion Workflow")

# --- UI DEFAULTS MANAGEMENT ---
DEFAULT_SETTINGS = {
    "sim_mode": "Closed Room CFD",
    "hvac_cfm": 500.0,
    "hvac_temp": 72.0,
    "use_manual_area": False,
    "hvac_area": 1.0,
    "hvac_press": 0.01,
    "wind_speed": 5.0,
    "wind_dir": 90.0,
    "stability_class": "D (Neutral)",
    "z0": 0.1,
    "sensors_input": "",
    "auto_stop": True,
    "sim_duration": 200,
    "mesh_res": "Fast",
    "refine_lvl": 2,
    "max_courant": 1.0,
    "starting_delta_t": 0.0001,
    "num_cores": 1
}

user_defaults = {}
if os.path.exists(DEFAULTS_FILE):
    try:
        with open(DEFAULTS_FILE, "r") as f:
            user_defaults = json.load(f)
    except: pass

for key, default_val in DEFAULT_SETTINGS.items():
    if key not in st.session_state:
        st.session_state[key] = user_defaults.get(key, default_val)

if "wind_speed" not in st.session_state: st.session_state.wind_speed = 5.0
if "wind_dir" not in st.session_state: st.session_state.wind_dir = 90.0
if "stability_class" not in st.session_state: st.session_state.stability_class = "D (Neutral)"
if "z0" not in st.session_state: st.session_state.z0 = 0.1

# --- HELPER FUNCS ---
def load_and_clean_csv(file_or_path):
    if isinstance(file_or_path, str) and not os.path.exists(file_or_path):
        return None
    df = pd.read_csv(file_or_path)
    df.columns = df.columns.str.strip().str.upper()
    
    col_map = {}
    for col in df.columns:
        if 'SCEN' in col: col_map[col] = 'scenario'
        elif 'COMP' in col or 'GAS' in col or 'SPEC' in col or 'MIX' in col: col_map[col] = 'COMPOSITION'
        elif 'RATE' in col or 'FLOW' in col or 'MASS' in col: col_map[col] = 'RATE'
        elif 'INV' in col and 'LB' in col: col_map[col] = 'INVENTORY_LB'
        elif 'INV' in col and 'KMOL' in col: col_map[col] = 'INVENTORY_KMOL'
        elif 'X' == col: col_map[col] = 'X'
        elif 'Y' == col: col_map[col] = 'Y'
        elif 'Z' == col: col_map[col] = 'Z'
        
    df.rename(columns=col_map, inplace=True)
    return df

def save_state(status_data, current_idx, pid, case_name, log_path, is_running, inlet, outlet, area):
    with open(STATE_FILE, "w") as f:
        json.dump({
            "status_data": status_data,
            "current_idx": current_idx,
            "pid": pid,
            "case_name": case_name,
            "log_path": log_path,
            "is_active": True,
            "is_running": is_running,
            "vent_config": {
                "inlet": inlet,
                "outlet": outlet,
                "area": area
            }
        }, f)

def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE, "r") as f: return json.load(f)
    return None

def clear_state():
    if os.path.exists(STATE_FILE): os.remove(STATE_FILE)
    if os.path.exists(SCENARIO_BACKUP): os.remove(SCENARIO_BACKUP)

def check_logs_for_errors(case_path):
    logs = ["log.blockMesh", "log.topoSet", "log.createPatch", "log.decompose", "log.run"]
    for log_file in logs:
        full_path = os.path.join(case_path, log_file)
        if os.path.exists(full_path):
            with open(full_path, "r") as f:
                lines = f.readlines()
                tail_content = "".join(lines[-2000:])
                if "FOAM FATAL" in tail_content or "FOAM FATAL IO ERROR" in tail_content or "command not found" in tail_content or "No such file or directory" in tail_content:
                    return log_file, tail_content[-1500:]
    return None, None

def check_steady_state(df_out_timeseries, inventory_remaining, is_blowdown, window_size=20, threshold=0.05):
    if df_out_timeseries is None or len(df_out_timeseries) < window_size:
        return False
        
    latest_time = df_out_timeseries['Time'].max()
    
    latest_data = df_out_timeseries[df_out_timeseries['Time'] == latest_time]
    eval_data = latest_data[~latest_data['Species'].isin(['Air'])]
    if eval_data.empty: return False
    
    primary_spec = eval_data.loc[eval_data['Concentration (ppm)'].idxmax()]['Species']
    limit_data = HAZARD_LIMITS.get(primary_spec, {"low": 5000})
    
    spec_max = df_out_timeseries[df_out_timeseries['Species'] == primary_spec]['Concentration (ppm)'].max()
    if spec_max < limit_data["low"] or latest_time < 20.0: return False
    
    spec_data = df_out_timeseries[df_out_timeseries['Species'] == primary_spec]['Concentration (ppm)'].values[-window_size:]
    
    if len(spec_data) == window_size:
        val_mean = np.mean(spec_data)
        
        if is_blowdown:
            if inventory_remaining < 0.01 and val_mean < 10.0: 
                return True
        else:
            if val_mean > limit_data["low"]:
                deviation = (np.max(spec_data) - np.min(spec_data)) / val_mean
                if deviation <= threshold:
                    return True
                
    return False

# --- SESSION STATE INITIALIZATION ---
if "sim_status" not in st.session_state: st.session_state.sim_status = []
if "saved_files" not in st.session_state: st.session_state.saved_files = []
if "is_running" not in st.session_state: st.session_state.is_running = False
if "current_idx" not in st.session_state: st.session_state.current_idx = 0
if "vents" not in st.session_state: st.session_state.vents = []
if "inlet_center" not in st.session_state: st.session_state.inlet_center = None
if "outlet_center" not in st.session_state: st.session_state.outlet_center = None
if "inlet_area_sqft" not in st.session_state: st.session_state.inlet_area_sqft = 1.0
if "error_logs" not in st.session_state: st.session_state.error_logs = []
if "scenario_start_time" not in st.session_state: st.session_state.scenario_start_time = None
if "error_found" not in st.session_state: st.session_state.error_found = False

pid = None
active_run = load_state()
is_resuming = False

if active_run and active_run.get("is_running", False):
    st.session_state.is_running = True
    st.session_state.sim_status = active_run.get("status_data", [])
    st.session_state.current_idx = active_run.get("current_idx", 0)
    pid = active_run.get("pid")
    
    if "vent_config" in active_run:
        vc = active_run["vent_config"]
        if st.session_state.inlet_center is None: st.session_state.inlet_center = vc.get("inlet")
        if st.session_state.outlet_center is None: st.session_state.outlet_center = vc.get("outlet")
        if st.session_state.inlet_area_sqft == 1.0: st.session_state.inlet_area_sqft = vc.get("area", 1.0)

tab1, tab2, tab3 = st.tabs(["🚀 Setup & Run", "📈 Live Monitor", "🔍 3D Analysis & Download"])

# ================= TAB 3: ANALYSIS =================
with tab3:
    st.header("🔍 Post-Processing & Visualization")
    
    col_refresh, _ = st.columns([1, 4])
    with col_refresh:
        if st.button("🔄 Refresh Scenarios List"):
            st.rerun()

    if os.path.exists(SIM_ROOT):
        cases_in_dir = [d for d in os.listdir(SIM_ROOT) if d.startswith("case_")]
        cases = []
        df_backup = load_and_clean_csv(SCENARIO_BACKUP)
        if df_backup is not None and not df_backup.empty:
            ordered_scen_ids = df_backup['scenario'].unique()
            for sid in ordered_scen_ids:
                c_name = f"case_{sid}"
                if c_name in cases_in_dir:
                    cases.append(c_name)
        for c in cases_in_dir:
            if c not in cases:
                cases.append(c)

        if cases:
            sel_case = st.selectbox("Select Scenario:", cases)
            case_path = os.path.join(SIM_ROOT, sel_case)
            
            col_a, col_b = st.columns(2)
            with col_a:
                if st.button(f"📦 Archive & Download {sel_case}"):
                    archive_base = os.path.join(ARCHIVE_DIR, f"{sel_case}_results")
                    runner.generate_outlet_report(case_path)
                    runner.generate_probe_report(case_path)
                    zip_path = shutil.make_archive(archive_base, 'zip', case_path)
                    with open(zip_path, "rb") as f: st.download_button("⬇️ Download Zip", f, file_name=f"{sel_case}.zip")
            
            with col_b:
                if st.button(f"🖥️ Launch ParaView for {sel_case}"):
                    import subprocess
                    foam_file = os.path.join(case_path, f"{sel_case}.foam")
                    if not os.path.exists(foam_file):
                        open(foam_file, 'a').close()
                    try:
                        subprocess.Popen(
                            ["paraview", f"{sel_case}.foam"], 
                            cwd=case_path,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            start_new_session=True
                        )
                        st.success("ParaView launched in a separate window! You can continue using the app.")
                    except FileNotFoundError:
                        st.error("ParaView executable not found. Ensure it is installed and in your system PATH.")

            st.divider()
            st.markdown("#### 🧊 In-App 3D Viewer (Animated Timeline)")
            
            if PYVISTA_AVAILABLE:
                vtk_dir = os.path.join(case_path, "postProcessing", "surfaces1")
                streamlines_hvac_dir = os.path.join(case_path, "postProcessing", "streamlines_hvac")
                streamlines_leak_dir = os.path.join(case_path, "postProcessing", "streamlines_leak")
                
                if os.path.exists(vtk_dir):
                    times = sorted([d for d in os.listdir(vtk_dir) if d.replace('.','').isdigit()], key=float)
                    if times:
                        
                        col_vis1, col_vis2, col_vis3 = st.columns(3)
                        with col_vis1:
                            show_room = st.checkbox("Show Room/Obstacles", value=True)
                            show_coord_grid = st.checkbox("Show Origin Planes (X=0, Y=0, Z=0) & Coordinate Outline (ft)", value=True)
                            show_markers = st.checkbox("Show Component Markers", value=True)
                            show_vis_cloud = st.checkbox("Show Visible Plume", value=True)
                        with col_vis2:
                            show_low_cloud = st.checkbox("Show Low Alarm Cloud", value=True)
                            show_high_cloud = st.checkbox("Show High Alarm Cloud", value=False)
                            show_emerg_cloud = st.checkbox("Show Emergency Cloud", value=False)
                        with col_vis3:
                            show_leak_lines = st.checkbox("Show Leak Streamlines", value=True)
                            show_hvac_lines = st.checkbox("Show HVAC Flowlines", value=False)
                            cloud_opacity = st.slider("Cloud Opacity", 0.1, 1.0, 0.5)

                        selected_time_float = st.select_slider("Scrub Timeline (Seconds):", options=[float(t) for t in times], value=float(times[-1]))
                        latest_time = next(t for t in times if float(t) == selected_time_float)

                        if st.button(f"Render 3D Model for {sel_case}"):
                            with st.spinner("Rendering 3D Model..."):
                                try:
                                    plotter = pv.Plotter(window_size=[900, 700])
                                    scale_ft = 3.28084 
                                    
                                    all_meshes_for_bounds = []
                                    stl_file_path = os.path.join(case_path, "constant", "triSurface", "room.stl")
                                    if show_room and os.path.exists(stl_file_path):
                                        room = pv.read(stl_file_path)
                                        room.points *= scale_ft
                                        all_meshes_for_bounds.append(room)
                                        
                                    cfg_path = os.path.join(case_path, "case_config.json")
                                    if stl_file_path and show_room and os.path.exists(cfg_path):
                                        with open(cfg_path, 'r') as f:
                                            c_cfg = json.load(f)
                                        inl = c_cfg.get("inlet")
                                        outl = c_cfg.get("outlet")
                                        leak_orig = c_cfg.get("leak")
                                        
                                        if inl and show_markers:
                                            plotter.add_point_labels([np.array(inl)*scale_ft], ["HVAC Inlet\nEntry Point"], point_color='blue', text_color='blue', shape_color='white', shape_opacity=0.8, point_size=30, bold=True)
                                        if outl and show_markers:
                                            plotter.add_point_labels([np.array(outl)*scale_ft], ["HVAC Outlet\nExit Point"], point_color='green', text_color='green', shape_color='white', shape_opacity=0.8, point_size=30, bold=True)
                                        if leak_orig and show_markers:
                                            plotter.add_point_labels([np.array(leak_orig)*scale_ft], ["Leak Origin"], point_color='yellow', text_color='yellow', shape_color='black', shape_opacity=0.8, point_size=30, bold=True)

                                    if show_vis_cloud:
                                        if not show_room:
                                            show_room = True 
                                            
                                    if all_meshes_for_bounds:
                                        combined_context = all_meshes_for_bounds[0]
                                        for m in all_meshes_for_bounds[1:]: combined_context += m
                                        bounds = combined_context.bounds
                                        bounds_min = [bounds[0], bounds[2], bounds[4]]
                                        bounds_max = [bounds[1], bounds[3], bounds[5]]
                                        plotter.add_mesh(combined_context, style='wireframe', color='black', opacity=0.3, label='Room Geometry context')
                                    else:
                                        bounds_min, bounds_max = [-10,-10,-10], [10,10,10]
                                        
                                    if show_coord_grid:
                                        p_center = [(bounds_max[i]+bounds_min[i])/2.0 for i in range(3)]
                                        p_size = [max(bounds_max[i]-bounds_min[i], 20.0)*1.2 for i in range(3)] 
                                        
                                        plane_yz = pv.Plane(center=(0.0, p_center[1], p_center[2]), direction=(1, 0, 0), i_size=p_size[1], j_size=p_size[2])
                                        plane_xz = pv.Plane(center=(p_center[0], 0.0, p_center[2]), direction=(0, 1, 0), i_size=p_size[0], j_size=p_size[2])
                                        plane_xy = pv.Plane(center=(p_center[0], p_center[1], 0.0), direction=(0, 0, 1), i_size=p_size[0], j_size=p_size[1])
                                        
                                        plotter.add_mesh(plane_yz, color='red', opacity=0.2, show_edges=True, edge_color='darkred', label='X=0 (YZ Plane)')
                                        plotter.add_mesh(plane_xz, color='green', opacity=0.2, show_edges=True, edge_color='darkgreen', label='Y=0 (XZ Plane)')
                                        plotter.add_mesh(plane_xy, color='blue', opacity=0.2, show_edges=True, edge_color='darkblue', label='Z=0 (XY Plane)')

                                        with warnings.catch_warnings():
                                            warnings.simplefilter("ignore")
                                            try:
                                                plotter.show_bounds(grid=False, location='outer', ticks='outside', xtitle='X (ft)', ytitle='Y (ft)', ztitle='Z (ft)', color="black")
                                            except TypeError:
                                                plotter.show_bounds(grid=False, location='outer', ticks='outside', xlabel='X (ft)', ylabel='Y (ft)', zlabel='Z (ft)', color="black")
                                            plotter.add_axes() 
                                        
                                    df_backup = load_and_clean_csv(SCENARIO_BACKUP)
                                    primary_gas = "CH4" 
                                    scen_data = pd.DataFrame()
                                    if df_backup is not None:
                                        scen_id = sel_case.replace("case_", "")
                                        scen_data = df_backup[df_backup['scenario'] == scen_id]
                                        if not scen_data.empty:
                                            comp_raw = str(scen_data.iloc[0].get('COMPOSITION', "CH4:1.0")).split(';')[0]
                                            primary_gas = comp_raw.split(':')[0].strip().upper()
                                            primary_gas = builder.NAME_MAP.get(primary_gas, primary_gas)
                                    
                                    limits = HAZARD_LIMITS.get(primary_gas, {"visible": 100, "low": 5000, "high": 25000, "emergency": 50000})
                                    limit_type = "LEL" if limits.get("type", "flammable") == "flammable" else "PEL"
                                    
                                    cloud_path = os.path.join(vtk_dir, latest_time)
                                    
                                    def render_cloud(cloud_suffix, color, label_text):
                                        if os.path.exists(cloud_path):
                                            vtks = glob.glob(os.path.join(cloud_path, f"*Hazard_Cloud_{cloud_suffix}*.vtk"))
                                            if vtks:
                                                mesh = pv.read(vtks[0])
                                                if mesh.n_points > 0:
                                                    mesh.points *= scale_ft
                                                    mesh = mesh.compute_normals()
                                                    plotter.add_mesh(mesh, color=color, opacity=cloud_opacity, smooth_shading=True, label=label_text)

                                    if show_vis_cloud:
                                        render_cloud("Visible", 'lightblue', f"{primary_gas} Visible Plume (> {limits['visible']} ppm)")
                                    if show_low_cloud:
                                        render_cloud("Low", 'orange', f"{primary_gas} Low Alarm (> {limits['low']} ppm)")
                                    if show_high_cloud:
                                        render_cloud("High", 'red', f"{primary_gas} High Alarm (> {limits['high']} ppm)")
                                    if show_emerg_cloud:
                                        render_cloud("Emergency", 'darkred', f"{primary_gas} Emergency Alarm (> {limits['emergency']} ppm)")
                                    
                                    def render_streamlines(s_dir, label_text, cmap, solid_color=None, tube_radius=0.015, glyph_factor=0.3):
                                        s_path = os.path.join(s_dir, latest_time)
                                        if os.path.exists(s_path):
                                            s_vtks = glob.glob(os.path.join(s_path, "*.vtk"))
                                            if s_vtks:
                                                full_mesh = pv.read(s_vtks[0])
                                                combined = full_mesh
                                                if isinstance(full_mesh, pv.MultiBlock):
                                                    combined = full_mesh.combine()

                                                if combined is not None and combined.n_points > 0 and combined.n_cells > 0:
                                                    combined.points *= scale_ft
                                                    try:
                                                        tube = combined.tube(radius=tube_radius * scale_ft) 
                                                        if tube.n_points > 0 and "U" in tube.array_names:
                                                            tube["U"] *= scale_ft
                                                            if solid_color:
                                                                plotter.add_mesh(tube, color=solid_color, opacity=0.8, label=label_text)
                                                                try:
                                                                    arrows = combined.glyph(orient="U", scale=False, factor=glyph_factor * scale_ft, tolerance=0.05, geom=pv.Arrow())
                                                                    plotter.add_mesh(arrows, color=solid_color)
                                                                except: pass
                                                            else:
                                                                plotter.add_mesh(tube, scalars="U", cmap=cmap, opacity=0.8, label=label_text, scalar_bar_args={'title': 'Velocity U (ft/s)'})
                                                                try:
                                                                    arrows = combined.glyph(orient="U", scale=False, factor=glyph_factor * scale_ft, tolerance=0.05, geom=pv.Arrow())
                                                                    plotter.add_mesh(arrows, color="black")
                                                                except: pass
                                                    except Exception as e:
                                                        pass

                                    if show_leak_lines:
                                        render_streamlines(streamlines_leak_dir, 'Leak Streamlines', "jet", tube_radius=0.015, glyph_factor=0.3)

                                    if show_hvac_lines:
                                        render_streamlines(streamlines_hvac_dir, 'HVAC Flowlines', None, solid_color="cyan", tube_radius=0.005, glyph_factor=0.2)
                                    
                                    plotter.add_legend()
                                    stpyvista(plotter, key=f"pv_{sel_case}_{latest_time}_{uuid.uuid4().hex}")
                                except Exception as render_err:
                                    st.error(f"Failed to render 3D view. Missing dependencies? Error: {render_err}")
                    else:
                        st.info("No surface data found in postProcessing.")
                else:
                    st.info("Simulation hasn't produced 3D surfaces yet. Run the solver first.")
            else:
                st.warning("PyVista is not installed or Xvfb is missing. Cannot render 3D inline.")
        else:
            st.info("No simulation data found yet.")

# ================= TAB 1: SETUP =================
with tab1:
    st.subheader("⚙️ Simulation Mode")
    sim_mode = st.radio("Choose the operational mode:", ["Closed Room CFD", "Open Air Gas Mapping"], horizontal=True, key="sim_mode")

    col1, col2 = st.columns([1, 2])
    with col1:
        st.subheader("1. Geometry & Topology")
        uploaded_stl = st.file_uploader("Geometry (STL)", type=["stl"])
        
        if uploaded_stl is None:
            st.session_state.vents = []
            st.session_state.inlet_center = None
            st.session_state.outlet_center = None
            if os.path.exists("temp.stl"):
                os.remove("temp.stl")
        else:
            with open("temp.stl", "wb") as f: f.write(uploaded_stl.getvalue())

        if sim_mode == "Closed Room CFD":
            if uploaded_stl is not None:
                if not st.session_state.vents:
                    with st.spinner("Analyzing Geometry for Openings..."):
                        st.session_state.vents = stl_helper.find_openings("temp.stl")
                if st.session_state.vents:
                    st.success(f"✅ Detected {len(st.session_state.vents)} Openings/Vents")
                    vent_opts = {f"Vent {i} (Center: {v['center']})": i for i, v in enumerate(st.session_state.vents)}
                    sel_inlet = st.selectbox("Define Inlet Vent (Entry)", options=vent_opts.keys(), key="in_sel")
                    sel_outlet = st.selectbox("Define Outlet Vent (Exit)", options=vent_opts.keys(), index=1 if len(vent_opts)>1 else 0, key="out_sel")
                    inlet_idx = vent_opts[sel_inlet]
                    outlet_idx = vent_opts[sel_outlet]
                    st.session_state.inlet_center = st.session_state.vents[inlet_idx]['center']
                    st.session_state.outlet_center = st.session_state.vents[outlet_idx]['center']
                    detected_area = st.session_state.vents[inlet_idx].get('area_m2', 0.25) * 10.7639
                    if st.session_state.inlet_area_sqft == 1.0: st.session_state.inlet_area_sqft = detected_area
                else:
                    st.warning("⚠️ No holes detected. Using default centers.")
                    st.session_state.inlet_center = [0,0,0]
                    st.session_state.outlet_center = [5,0,0]

        uploaded_csv = st.file_uploader("Scenarios (CSV)", type=["csv"])
        st.info("ℹ️ Support added for Inventory-based blowdown releases. Ensure your CSV contains an 'INVENTORY_LB' or 'INVENTORY_KMOL' column. If specified, the leak rate will decay over time and the simulation will automatically run until the inventory is exhausted. If omitted, the leak is assumed to be a constant source.")

        with st.popover("🧪 View Supported Gas Species & Graded Alarm Limits"):
            st.markdown("### Gas Species Alarm Setpoints (PEL & LEL based)")
            st.markdown("Detector networks are typically configured with graded alarms rather than single IDLH limits to ensure early intervention. *Note: Concentrations are internally calculated in precise Mass Fractions based on your mix.*")
            
            table_data = []
            for sp, limits in HAZARD_LIMITS.items():
                table_data.append({"Species": sp, "Type": limits["type"].title(), "Low Alarm (ppm)": limits["low"], "High Alarm (ppm)": limits["high"], "Emergency (ppm)": limits["emergency"]})
            st.dataframe(pd.DataFrame(table_data), hide_index=True)

        st.divider()

        if sim_mode == "Closed Room CFD":
            st.subheader("2. HVAC Properties (Imperial)")
            with st.expander("💨 Vent Inlet/Outlet Settings", expanded=True):
                col_hvac1, col_hvac2 = st.columns(2)
                with col_hvac1:
                    hvac_cfm = st.number_input("Inlet Air Flow (CFM)", step=10.0, key="hvac_cfm")
                    hvac_temp = st.number_input("Inlet Temperature (°F)", step=1.0, key="hvac_temp")
                with col_hvac2:
                    use_manual_area = st.checkbox("Override Vent Area?", key="use_manual_area")
                    if use_manual_area:
                        hvac_area = st.number_input("Vent Area (sq. ft)", step=0.1, key="hvac_area")
                    else:
                        detected_area = st.session_state.inlet_area_sqft if st.session_state.inlet_area_sqft > 0 else 1.0
                        if uploaded_stl and len(st.session_state.vents) > 0:
                            st.info(f"Detected Area: {detected_area:.2f} sq. ft (Calculated from the uploaded 3D Model)")
                        elif uploaded_stl and len(st.session_state.vents) == 0:
                            st.warning(f"No vents detected in the uploaded 3D Model (enclosed space). Using default: {detected_area:.2f} sq. ft")
                        else:
                            st.info(f"Default Area: {detected_area:.2f} sq. ft")
                        hvac_area = detected_area
                    hvac_press = st.number_input("Outlet Pressure (in. WC)", step=0.001, format="%.4f", key="hvac_press")
            
            wind_speed = st.session_state.wind_speed
            wind_dir = st.session_state.wind_dir
            stability_class = st.session_state.stability_class
            z0 = st.session_state.z0

        elif sim_mode == "Open Air Gas Mapping":
            st.subheader("2. Atmospheric Boundary Layer (ABL)")
            with st.expander("🌍 Wind Rose & Stability Settings", expanded=True):
                col_abl1, col_abl2 = st.columns(2)
                with col_abl1:
                    wind_speed = st.number_input("Wind Speed (m/s)", step=0.5, key="wind_speed")
                    wind_dir = st.number_input("Wind Direction (deg, 0=N)", step=10.0, key="wind_dir")
                with col_abl2:
                    stability_class = st.selectbox("Pasquill-Gifford Stability", ["A (Very Unstable)", "B", "C", "D (Neutral)", "E", "F (Very Stable)"], key="stability_class")
                    z0 = st.number_input("Surface Roughness z0 (m)", step=0.05, key="z0")
            
            hvac_cfm = st.session_state.hvac_cfm
            hvac_temp = st.session_state.hvac_temp
            hvac_press = st.session_state.hvac_press
            hvac_area = st.session_state.hvac_area

        st.divider()

        st.subheader("3. Sensor/Probe Locations (Optional)")
        st.info("📝 How to enter sensor coordinates:\nEnter the X, Y, Z coordinates enclosed in parentheses separated by spaces. For multiple sensors, separate each group with a comma.\n\nExample: (1.5 2.0 1.0), (3.0 1.5 2.5)")
        sensors_input = st.text_input("Probe Coordinates", key="sensors_input")
        probe_locations = []
        if sensors_input.strip():
            matches = re.findall(r'\((.*?)\)', sensors_input)
            for m in matches:
                try:
                    coords = [float(c) for c in m.split()]
                    if len(coords) == 3: probe_locations.append(coords)
                except ValueError: pass

        st.divider()

        st.subheader("4. Simulation Control")
        st.info("ℹ️ If a simulation status says 'Completed (Reached Max Time)', it means the scenario ran for the full 'Simulation Duration' specified below without the gas completely clearing out. Increase the duration if you want to watch the gas fully dissipate.")
        auto_stop = st.checkbox("Enable Auto-Stop (Inventory Released + decaying concentrations)", key="auto_stop")
        col_sim1, col_sim2 = st.columns(2)
        with col_sim1:
            sim_duration = st.number_input("Simulation Duration (s)", step=50, key="sim_duration")
            mesh_res = st.radio("Mesh Resolution", ["Fast", "Accurate"], key="mesh_res")
            refine_lvl = st.slider("Refinement Level", 0, 3, key="refine_lvl")
        with col_sim2:
            max_courant = st.number_input("Max Courant Number", step=0.1, key="max_courant")
            starting_delta_t = st.number_input("Starting Delta T (s)", format="%.5f", key="starting_delta_t")
            max_cores = os.cpu_count() or 1
            safe_cores = min(st.session_state.num_cores, max_cores)
            st.session_state.num_cores = safe_cores
            num_cores = st.slider("CPU Cores", 1, max_cores, key="num_cores")

    with col2:
        st.subheader("Job Preview")
        if uploaded_csv:
            uploaded_csv.seek(0)
            raw_df = load_and_clean_csv(uploaded_csv)
            if raw_df is None or 'scenario' not in raw_df.columns or 'COMPOSITION' not in raw_df.columns:
                st.error("⚠️ Invalid CSV Format. Ensure you have 'Scenario' and 'Composition' (or 'Gas'/'Species') columns.")
                df = None
            else:
                all_scenarios = raw_df['scenario'].unique().tolist()
                st.markdown("🎯 Active Scenarios (Remove scenarios to exclude them from the batch)")
                selected_scenarios = st.multiselect("Active Scenarios:", all_scenarios, default=all_scenarios, label_visibility="collapsed")
                if selected_scenarios:
                    df = raw_df[raw_df['scenario'].isin(selected_scenarios)].copy()
                    st.info("ℹ️ Only the first 50 entries will be shown in this preview for performance.")
                    df_display = df.head(50).reset_index(drop=True)
                    df_display.index = df_display.index + 1
                    
                    st.dataframe(df_display, width="stretch")
                else:
                    st.warning("⚠️ No scenarios selected.")
                    df = None

            if st.button("💾 Load Configuration") and df is not None:
                if not uploaded_stl:
                    st.error("Upload STL first.")
                elif not selected_scenarios:
                    st.error("Select at least one scenario.")
                else:
                    df.to_csv(SCENARIO_BACKUP, index=False)
                    st.session_state.sim_status = []
                    st.session_state.error_logs = []
                    for scen_id in selected_scenarios:
                        st.session_state.sim_status.append({"Scenario": scen_id, "Status": "Ready", "Leaks": 0, "Log Path": ""})
                    st.session_state.current_idx = 0
                    save_state(st.session_state.sim_status, 0, None, "", "", False, st.session_state.inlet_center, st.session_state.outlet_center, st.session_state.inlet_area_sqft)
                    st.success("Loaded! Go to Live Monitor.")
        
        st.divider()
        st.subheader("⚙️ App Preferences")
        if st.button("💾 Save Current Inputs as Default"):
            defaults_to_save = {}
            for k in DEFAULT_SETTINGS.keys():
                defaults_to_save[k] = st.session_state.get(k, DEFAULT_SETTINGS[k])
            with open(DEFAULTS_FILE, "w") as f:
                json.dump(defaults_to_save, f)
            st.success("Default settings saved successfully! They will be applied next time you launch the app.")
            
        if st.button("🔄 Reset to Factory Defaults"):
            if os.path.exists(DEFAULTS_FILE):
                os.remove(DEFAULTS_FILE)
            for k in DEFAULT_SETTINGS.keys():
                st.session_state[k] = DEFAULT_SETTINGS[k]
            st.rerun()

# ================= TAB 2: MONITOR =================
with tab2:
    st.subheader("Batch Execution")

    st.markdown("#### 📊 Overnight Run Statistics")
    m_col1, m_col2, m_col3 = st.columns(3)
    queue_metric = m_col1.empty()
    success_metric = m_col2.empty()
    failed_metric = m_col3.empty() 

    def update_status_table(data):
        if not data: return
        df_stat = pd.DataFrame(data)
        df_stat.index = df_stat.index + 1
        status_table.dataframe(df_stat, width="stretch")
        
        c_count = sum(1 for s in data if "✅" in s["Status"])
        f_count = sum(1 for s in data if "❌" in s["Status"])
        queue_metric.metric("Total Scenarios Queue", len(data))
        success_metric.metric("Successfully Completed", c_count)
        failed_metric.metric("Failed / Aborted", f_count)

    st.divider()

    if st.session_state.error_logs:
        st.error("⚠️ The following critical errors occurred during the batch run:")
        for err in st.session_state.error_logs:
            st.code(err)
        st.divider()

    mon_col1, mon_col2 = st.columns([1, 1])

    with mon_col1:
        col_prog1, col_prog2, col_prog3 = st.columns([2, 1, 1])
        with col_prog1: progress_bar = st.empty()
        with col_prog2: timer_placeholder = st.empty()
        with col_prog3: inv_placeholder = st.empty() 
        
        status_table = st.empty()
        if st.session_state.sim_status:
            update_status_table(st.session_state.sim_status)

        st.divider()
        st.markdown("#### 🚨 Alarm Trigger Log (Outlet/Vents)")
        alarm_table = st.empty()

        st.divider()
        st.markdown("#### 🧪 Outlet Monitor")
        flow_unit = st.selectbox("Flow Rate Unit", ["Mass (kg/s)", "Molar (kmol/s)", "Volumetric (m3/s)"], key="flow_unit_select")
        outlet_table = st.empty()
        outlet_chart_placeholder = st.empty()

    with mon_col2:
        st.markdown("### 📡 Live Residuals")
        st.caption("🔍 Variables: Ux/Uy/Uz (Velocity), p/p_rgh (Pressure), k/epsilon/omega (Turbulence), h (Enthalpy), Species (Mass Fractions)")
        chart_placeholder = st.empty()
        st.divider()
        st.caption("Terminal Log Tail:")
        live_log_tail = st.empty()
        error_placeholder = st.empty()

    if st.session_state.is_running:
        st.button("⏳ Processing Batch...", disabled=True)
    else:
        if st.button("▶️ Start Simulation Loop"):
            st.session_state.is_running = True
            st.session_state.error_found = False
            active_run = load_state()
            if active_run:
                active_run["is_running"] = True
                with open(STATE_FILE, "w") as f: json.dump(active_run, f)
            st.rerun()

    if st.session_state.is_running:
        df = load_and_clean_csv(SCENARIO_BACKUP) if os.path.exists(SCENARIO_BACKUP) else None
        if df is None or df.empty:
            st.error("❌ Data missing. Please configure scenarios in Setup.")
            st.session_state.is_running = False
            st.stop()

        status_data = st.session_state.sim_status
        total_scenarios = len(status_data)

        # Main Loop
        while st.session_state.current_idx < total_scenarios:
            current_idx = st.session_state.current_idx
            scen_id = status_data[current_idx]["Scenario"]
            case_name = f"case_{scen_id}"
            case_path = os.path.join(SIM_ROOT, case_name)
            log_path = os.path.join(case_path, "log.run")

            progress_bar.progress(current_idx / total_scenarios, text=f"Processing Scenario {current_idx + 1} of {total_scenarios} (ID: {scen_id})...")

            # Dynamic Configuration Hash generation
            group_df = df[df["scenario"] == scen_id]
            hvac_area_local = st.session_state.inlet_area_sqft
            case_run_config = {
                "csv_data": group_df.to_dict('records'),
                "sim_mode": st.session_state.sim_mode,
                "hvac_cfm": st.session_state.hvac_cfm,
                "hvac_temp": st.session_state.hvac_temp,
                "hvac_press": st.session_state.hvac_press,
                "hvac_area": hvac_area_local, # evaluated locally
                "wind_speed": st.session_state.wind_speed,
                "wind_dir": st.session_state.wind_dir,
                "stability_class": st.session_state.stability_class,
                "z0": st.session_state.z0,
                "sim_duration": st.session_state.sim_duration,
                "max_courant": st.session_state.max_courant,
                "mesh_res": st.session_state.mesh_res,
                "refine_lvl": st.session_state.refine_lvl
            }
            config_str = json.dumps(case_run_config, sort_keys=True)
            current_hash = hashlib.md5(config_str.encode('utf-8')).hexdigest()

            # Check if this exact configuration has been run successfully before
            hash_file = os.path.join(case_path, "scenario_hash.txt")
            hash_matches = False
            if os.path.exists(hash_file):
                with open(hash_file, "r") as f:
                    saved_hash = f.read().strip()
                if saved_hash == current_hash:
                    hash_matches = True

            is_completed = runner.check_if_completed(case_path)

            if is_completed and hash_matches:
                status_data[current_idx]["Status"] = "✅ Already Done"
                update_status_table(status_data)
                st.session_state.current_idx += 1
                save_state(status_data, st.session_state.current_idx, None, "", "", True, st.session_state.inlet_center, st.session_state.outlet_center, st.session_state.inlet_area_sqft)
                continue

            if pid is None or not runner.is_process_running(pid):
                case_exists = os.path.exists(case_path)
                
                # Initialize if case doesnt exist OR settings have changed (hash mismatch)
                if not case_exists or not hash_matches:
                        
                    status_data[current_idx]["Status"] = "🔄 Initializing..."
                    update_status_table(status_data)
                    stl_obj = open("temp.stl", "rb") if os.path.exists("temp.stl") else None
                    
                    builder.create_case_from_group(
                        case_name, stl_obj, group_df, st.session_state.mesh_res, st.session_state.refine_lvl, st.session_state.num_cores,
                        sim_mode=st.session_state.sim_mode, hvac_cfm=st.session_state.hvac_cfm, hvac_temp_f=st.session_state.hvac_temp, hvac_press_inwc=st.session_state.hvac_press, hvac_area_sqft=hvac_area_local,
                        vent_inlet_center=st.session_state.inlet_center, vent_outlet_center=st.session_state.outlet_center,
                        wind_speed=st.session_state.wind_speed, wind_dir=st.session_state.wind_dir, stability_class=st.session_state.stability_class, z0=st.session_state.z0,
                        probe_locations=probe_locations, sim_duration=st.session_state.sim_duration, max_courant=st.session_state.max_courant, delta_t=st.session_state.starting_delta_t
                    )
                    
                    # Save hash immediately so future runs know this configuration was built
                    with open(os.path.join(case_path, "scenario_hash.txt"), "w") as f:
                        f.write(current_hash)
                else:
                    status_data[current_idx]["Status"] = "🔄 Resuming..."
                    update_status_table(status_data)

                runner._optimize_control_dict(case_path)
                pid, _ = runner.run_case_detached(case_path, st.session_state.num_cores) 
                
                status_data[current_idx]["Status"] = "🔄 Running..."
                update_status_table(status_data)
                
                save_state(status_data, current_idx, pid, case_name, log_path, True, st.session_state.inlet_center, st.session_state.outlet_center, st.session_state.inlet_area_sqft)
                st.session_state.scenario_start_time = time.time()

            latest_outlet_fig = None
            alarm_status_dict = {}

            chart_placeholder.empty()
            outlet_table.empty()
            alarm_table.empty()
            outlet_chart_placeholder.empty()
            live_log_tail.empty()

            # Active Monitor Loop
            while runner.is_process_running(pid):
                if st.session_state.scenario_start_time:
                    elapsed = int(time.time() - st.session_state.scenario_start_time)
                    mins, secs = divmod(elapsed, 60)
                    timer_placeholder.markdown(f"⏱️ Elapsed Time: {mins:02d}:{secs:02d}")

                failed_log, error_msg = check_logs_for_errors(case_path)
                if failed_log:
                    st.session_state.error_found = True
                    err_txt = f"[{datetime.now().strftime('%H:%M:%S')}] Fatal Error in {case_name} ({failed_log}):\n{error_msg}"
                    st.session_state.error_logs.append(err_txt)
                    error_placeholder.error(err_txt)
                    status_data[current_idx]["Status"] = "❌ Aborted"
                    update_status_table(status_data)
                    try: os.kill(pid, 15)
                    except: pass
                    break

                if os.path.exists(log_path):
                    with open(log_path, "r") as f:
                        lines = f.readlines()
                        if lines:
                            live_log_tail.code(f"--- Last Updated: {datetime.now().strftime('%H:%M:%S')} ---\n" + "".join(lines[-10:]))
                            
                            # ROBUST END DETECTION TO BREAK LOOP
                            tail_content = "".join(lines[-50:])
                            if "\nEnd\n" in tail_content:
                                break

                current_time = 0.0
                rem_inv_kg = 0.0
                is_blowdown = False
                df_resid = runner.get_residuals(log_path)
                if df_resid is not None and not df_resid.empty and len(df_resid['Time'].unique()) > 1:
                        
                    df_resid.sort_values(by="Variable", inplace=True)
                    current_time = df_resid['Time'].max()
                    
                    unique_chart_key = f"resid_{scen_id}_{uuid.uuid4().hex}"
                    
                    # Disabled markers for massive performance gain when rendering full log histories
                    fig = px.line(df_resid, x="Time", y="Residual", color="Variable", log_y=True, color_discrete_map=RESIDUAL_COLOR_MAP, markers=False)
                    fig.update_xaxes(range=[0, max(10, current_time * 1.05)]) # Rigidly lock X-axis to 0
                    chart_placeholder.plotly_chart(fig, key=unique_chart_key)
                else:
                    chart_placeholder.info(f"⏳ Accumulating residual data for {case_name}...")
                    
                # Update Live Inventory Countdown if meta exists
                meta_path = os.path.join(case_path, "leak_meta.json")
                if os.path.exists(meta_path):
                    try:
                        with open(meta_path, "r") as mf:
                            leak_meta = json.load(mf)
                        if leak_meta:
                            is_blowdown = True # Inform auto-stop that this has an inventory limit
                            for lk in leak_meta:
                                m0 = lk.get("M0", 0)
                                tau = lk.get("tau", 1)
                                if tau > 0 and current_time < 5.0 * tau:
                                    rem_inv_kg += m0 * math.exp(-current_time / tau)
                            
                            # Convert internal mass back to pounds for the UI display
                            rem_inv_lb = rem_inv_kg / 0.453592
                            inv_placeholder.metric("⏳ Remaining Mass Inventory (lb)", f"{rem_inv_lb:.3f}")
                    except: pass

                outlet_df = runner.get_outlet_data(case_path)
                if outlet_df is not None:
                    disp_df = outlet_df.copy().reset_index(drop=True)
                    disp_df.index = disp_df.index + 1
                    disp_df["Concentration (ppm)"] = disp_df["Concentration (ppm)"].map(lambda x: "{:.3e}".format(x))
                    
                    # Dynamic Flow Rate Column
                    if flow_unit == "Mass (kg/s)":
                        disp_df["Flow Rate"] = disp_df["Mass Flow (kg/s)"].map(lambda x: "{:.3e}".format(x))
                    elif flow_unit == "Molar (kmol/s)":
                        disp_df["Flow Rate"] = disp_df["Molar Flow (kmol/s)"].map(lambda x: "{:.3e}".format(x))
                    else:
                        # Approx Volumetric (m3/s) = Mass (kg/s) / Density of Air (~1.2 kg/m3)
                        disp_df["Flow Rate"] = (disp_df["Mass Flow (kg/s)"] / 1.2).map(lambda x: "{:.3e}".format(x))
                        
                    outlet_table.dataframe(disp_df[["Species", "Concentration (ppm)", "Flow Rate"]], width="stretch")
                else:
                    outlet_table.info(f"⏳ Awaiting outlet flow data for {case_name}...")

                raw_timeseries = runner.get_outlet_timeseries(case_path)
                if raw_timeseries is not None and not raw_timeseries.empty:
                    df_out_timeseries = raw_timeseries
                    
                    if not df_out_timeseries.empty and len(df_out_timeseries['Time'].unique()) > 1:
                        
                        # Use make_subplots to put Air on a secondary Y-axis so trace gases dont get squashed
                        fig_out = make_subplots(specs=[[{"secondary_y": True}]])
                        fig_out.update_layout(title_text="Live Outlet Gas Concentrations", legend_title_text="Species", margin=dict(t=50, l=25, r=25, b=25))
                        
                        for spec in df_out_timeseries["Species"].unique():
                            spec_df = df_out_timeseries[df_out_timeseries['Species'] == spec]
                            is_air = (spec == "Air")
                            
                            fig_out.add_trace(
                                go.Scatter(
                                    x=spec_df["Time"], 
                                    y=spec_df["Concentration (ppm)"], 
                                    name=spec, 
                                    mode='lines', 
                                    # Fallback explicitly specified to avoid implicit grey fallback 
                                    line=dict(color=RESIDUAL_COLOR_MAP.get(spec, "#ff0000"))
                                ),
                                secondary_y=is_air
                            )

                        fig_out.update_yaxes(title_text="Trace Gases (ppm)", secondary_y=False, autorange=True)
                        fig_out.update_yaxes(title_text="Air (ppm)", secondary_y=True, range=[0, 1050000], showgrid=False)
                        fig_out.update_xaxes(range=[0, max(10, current_time * 1.05)], title_text="Time")
                        
                        alarm_log_data = []
                        
                        for spec in df_out_timeseries["Species"].unique():
                            if spec in HAZARD_LIMITS:
                                limits = HAZARD_LIMITS[spec]
                                spec_df = df_out_timeseries[df_out_timeseries['Species'] == spec]
                                
                                for i, (lvl_name, lvl_key, color) in enumerate([("Low (Warning)", "low", "orange"), ("High (Danger)", "high", "red"), ("Emergency", "emergency", "darkred")]):
                                    threshold = limits[lvl_key]
                                    # Add horizontal line to the primary y-axis
                                    fig_out.add_hline(y=threshold, line_dash="dash", annotation_text=f"{spec} {lvl_name}", annotation_position="top left", line_color=color, secondary_y=False)
                                    
                                    exceed_df = spec_df[spec_df['Concentration (ppm)'] >= threshold]
                                    trigger_time = f"{exceed_df['Time'].min():.2f}s" if not exceed_df.empty else "Not Reached"
                                    max_ppm = f"{spec_df['Concentration (ppm)'].max():.1f}"
                                    
                                    alarm_log_data.append({
                                        "Species": spec,
                                        "Alarm Level": lvl_name,
                                        "Limit (ppm)": threshold,
                                        "Trigger Time": trigger_time,
                                        "Max Detected (ppm)": max_ppm
                                    })
                        
                        if alarm_log_data:
                            alarm_table.dataframe(pd.DataFrame(alarm_log_data), width="stretch", hide_index=True)
                            
                        unique_outlet_key = f"outlet_{scen_id}_{uuid.uuid4().hex}"
                        outlet_chart_placeholder.plotly_chart(fig_out, key=unique_outlet_key)
                        latest_outlet_fig = fig_out

                        # Check Auto-Stop using the smart logic (Inventory vs Continuous)
                        if st.session_state.auto_stop and check_steady_state(df_out_timeseries, rem_inv_kg, is_blowdown):
                            status_data[current_idx]["Status"] = "✅ Full released curve captured" if is_blowdown else "✅ Steady-State Reached"
                            update_status_table(status_data)
                            
                            # GRACEFUL SHUTDOWN FIX: Modify controlDict to tell OpenFOAM to write the final 3D frame and safely exit
                            cd_path = os.path.join(case_path, "system", "controlDict")
                            if os.path.exists(cd_path):
                                with open(cd_path, "r") as f:
                                    cd_data = f.read()
                                cd_data = cd_data.replace("stopAt endTime;", "stopAt writeNow;")
                                with open(cd_path, "w") as f:
                                    f.write(cd_data)
                            
                            # Wait for OpenFOAM to process the writeNow command and exit naturally
                            time.sleep(5)
                            try: os.kill(pid, 15) # Failsafe cleanup
                            except: pass
                            break
                    else:
                        outlet_chart_placeholder.info(f"⏳ Accumulating outlet concentration data for {case_name}...")
                else:
                    outlet_chart_placeholder.info(f"⏳ Accumulating outlet concentration data for {case_name}...")

                time.sleep(2)

            # Process successfully finished or forcefully aborted
            sim_success = False
            if os.path.exists(log_path):
                with open(log_path, "r") as f:
                    lines = f.readlines()
                    tail_content = "".join(lines[-200:])
                    if "\nEnd\n" in tail_content or "released" in status_data[current_idx]["Status"] or "Steady-State" in status_data[current_idx]["Status"]:
                        sim_success = True

            if sim_success and not st.session_state.error_found:
                # If it stopped due to time expiration, use a simpler success message
                if "released" not in status_data[current_idx]["Status"] and "Steady-State" not in status_data[current_idx]["Status"]:
                    status_data[current_idx]["Status"] = "✅ Completed (Reached Max Time)"
                update_status_table(status_data)
                
                try:
                    with open(os.path.join(case_path, "DONE"), "w") as f:
                        f.write(f"Completed on {datetime.now()}")
                    # Save exact configuration hash upon success to prevent redundant re-runs
                    with open(os.path.join(case_path, "scenario_hash.txt"), "w") as f:
                        f.write(current_hash)
                except: pass
                
                try:
                    runner.generate_outlet_report(case_path)
                    runner.generate_probe_report(case_path)
                    if latest_outlet_fig is not None:
                        latest_outlet_fig.write_html(os.path.join(case_path, "postProcessing", "Outlet_Concentration_Chart.html"))
                    archive_base = os.path.join(ARCHIVE_DIR, f"case_{scen_id}_results")
                    shutil.make_archive(archive_base, 'zip', case_path)
                    st.session_state.saved_files.append(f"case_{scen_id}_results.zip")
                except: pass
            else:
                status_data[current_idx]["Status"] = "❌ Failed/Aborted"
                update_status_table(status_data)

            # Advance accurately to the next scenario
            st.session_state.error_found = False
            st.session_state.current_idx += 1
            pid = None
            st.session_state.scenario_start_time = None
            save_state(status_data, st.session_state.current_idx, None, "", "", True, st.session_state.inlet_center, st.session_state.outlet_center, st.session_state.inlet_area_sqft)

        progress_bar.progress(1.0, text="✅ All Scenarios Completed!")
        clear_state()
        st.session_state.is_running = False
        st.balloons()
        st.success("✅ Batch Run Complete!")
        st.rerun()