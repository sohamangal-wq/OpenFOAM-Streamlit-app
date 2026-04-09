import subprocess
import os
import re
import pandas as pd
import shutil
import numpy as np

try:
    import pyvista as pv
    import plotly.graph_objects as go
except ImportError:
    pv = None
    go = None

class SimulationRunner:
    def __init__(self, solver="foamRun"):
        self.solver = solver

        self.MOL_WEIGHTS = {
            "Air": 28.96, "N2": 28.01, "O2": 31.99, "CO2": 44.01, "H2O": 18.02,
            "H2": 2.016, "CH4": 16.04, "C2H6": 30.07, "C3H8": 44.1, "C4H10": 58.12,
            "C2H4": 28.05, "C2H2": 26.04,
            "NH3": 17.03, "CO": 28.01, "H2S": 34.08, "SO2": 64.06, "Cl2": 70.9,
            "He": 4.00, "Ar": 39.95, "C6H6": 78.11, "C7H8": 92.14, "C8H10": 106.17,
            "C6H14": 86.18, "C5H12": 72.15, "CH3OH": 32.04, "NO2": 46.005,
            "NO": 30.01, "HCN": 27.025
        }

    def check_if_completed(self, case_path):
        if not os.path.exists(case_path): return False
        if os.path.exists(os.path.join(case_path, "DONE")): return True
        return False

    def _read_pp_dat(self, pp_dir, obj_name):
        obj_dir = os.path.join(pp_dir, obj_name)
        if not os.path.exists(obj_dir): return None
        
        dfs = []
        for d in os.listdir(obj_dir):
            try:
                float(d) 
                file_path = os.path.join(obj_dir, d, "surfaceFieldValue.dat")
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)
                    if not df.empty and df.shape[1] >= 2:
                        df = df.iloc[:, [0, 1]] 
                        df.columns = ["Time", "Value"]
                        dfs.append(df)
            except ValueError:
                pass
                
        if not dfs: return None
        
        combined = pd.concat(dfs, ignore_index=True)
        combined["Time"] = combined["Time"].astype(float).round(4)
        combined = combined.sort_values(by="Time").drop_duplicates(subset=["Time"], keep="last")
        combined.set_index("Time", inplace=True)
        return combined

    def generate_outlet_report(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        if not os.path.exists(pp_dir): return None

        df_flow = self._read_pp_dat(pp_dir, "outletFlow")
        if df_flow is None: return None
        df_flow.rename(columns={"Value": "Total_Mass_Flow_kg_s"}, inplace=True)
        df_flow["Total_Mass_Flow_kg_s"] = df_flow["Total_Mass_Flow_kg_s"].abs()

        species_dfs = []
        found_species = []

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet_") and folder not in ["outletFlow", "outlet_Air"]:
                spec = folder.replace("outlet_", "")
                df_s = self._read_pp_dat(pp_dir, folder)
                if df_s is not None:
                    df_s.rename(columns={"Value": f"MassFrac_{spec}"}, inplace=True)
                    species_dfs.append(df_s)
                    found_species.append(spec)

        if not species_dfs: return None

        full_df = df_flow.join(species_dfs, how='outer').fillna(0.0)
        
        sum_mass_frac = sum(full_df[f"MassFrac_{s}"] for s in found_species) if found_species else 0.0
        full_df["MassFrac_Air"] = (1.0 - sum_mass_frac).clip(lower=0.0)

        full_df["Molar_Air"] = full_df["MassFrac_Air"] / 28.96
        total_molar = full_df["Molar_Air"].copy()

        for spec in found_species:
            mw = self.MOL_WEIGHTS.get(spec, 28.96)
            full_df[f"Molar_{spec}"] = full_df[f"MassFrac_{spec}"] / mw
            total_molar += full_df[f"Molar_{spec}"]

        safe_total = total_molar.replace(0, 1.0)

        for spec in found_species:
            full_df[f"Concentration_{spec}ppm"] = (full_df[f"Molar_{spec}"] / safe_total) * 1e6
            full_df[f"MassFlow_{spec}kg_s"] = full_df["Total_Mass_Flow_kg_s"] * full_df[f"MassFrac_{spec}"]

        output_file = os.path.join(pp_dir, "EXIT_DATA_SUMMARY.csv")
        full_df.to_csv(output_file)

        return output_file

    def get_outlet_timeseries(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        if not os.path.exists(pp_dir): return None

        df_flow = self._read_pp_dat(pp_dir, "outletFlow")
        if df_flow is None: return None
        df_flow.rename(columns={"Value": "Total_Mass_Flow_kg_s"}, inplace=True)
        df_flow["Total_Mass_Flow_kg_s"] = df_flow["Total_Mass_Flow_kg_s"].abs()

        species_dfs = []
        found_species = []

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet_") and folder not in ["outletFlow", "outlet_Air"]:
                spec = folder.replace("outlet_", "")
                df_s = self._read_pp_dat(pp_dir, folder)
                if df_s is not None:
                    df_s.rename(columns={"Value": f"MassFrac_{spec}"}, inplace=True)
                    species_dfs.append(df_s)
                    found_species.append(spec)

        if not species_dfs: return None

        full_df = df_flow.join(species_dfs, how='outer').fillna(0.0)

        sum_mass_frac = sum(full_df[f"MassFrac_{s}"] for s in found_species) if found_species else 0.0
        full_df["MassFrac_Air"] = (1.0 - sum_mass_frac).clip(lower=0.0)

        full_df["Molar_Air"] = full_df["MassFrac_Air"] / 28.96
        total_molar = full_df["Molar_Air"].copy()
        
        for spec in found_species:
            mw = self.MOL_WEIGHTS.get(spec, 28.96)
            full_df[f"Molar_{spec}"] = full_df[f"MassFrac_{spec}"] / mw
            total_molar += full_df[f"Molar_{spec}"]

        safe_total = total_molar.replace(0, 1.0)

        plot_data = []
        for spec in found_species:
            conc_ppm = (full_df[f"Molar_{spec}"] / safe_total) * 1e6
            temp_df = pd.DataFrame({
                "Time": full_df.index,
                "Species": spec,
                "Concentration (ppm)": conc_ppm.values
            })
            plot_data.append(temp_df)
            
        conc_ppm_air = (full_df["Molar_Air"] / safe_total) * 1e6
        plot_data.append(pd.DataFrame({
            "Time": full_df.index,
            "Species": "Air",
            "Concentration (ppm)": conc_ppm_air.values
        }))

        if plot_data:
            return pd.concat(plot_data, ignore_index=True)

        return None

    def generate_probe_report(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        sensors_dir = os.path.join(pp_dir, "gasSensors")

        if not os.path.exists(sensors_dir): return None

        time_dirs = []
        for d in os.listdir(sensors_dir):
            try:
                time_dirs.append((float(d), d))
            except ValueError:
                pass

        if not time_dirs: return None
        time_dirs.sort(key=lambda x: x[0])

        file_names = set()
        for _, d_str in time_dirs:
            for f in os.listdir(os.path.join(sensors_dir, d_str)):
                if os.path.isfile(os.path.join(sensors_dir, d_str, f)):
                    file_names.add(f)

        all_probe_data = []

        for file_name in file_names:
            dfs = []
            for _, d_str in time_dirs:
                file_path = os.path.join(sensors_dir, d_str, file_name)
                if os.path.exists(file_path):
                    try:
                        df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)
                        num_probes = df.shape[1] - 1
                        col_names = ["Time"] + [f"Probe_{i}_{file_name}" for i in range(num_probes)]
                        df.columns = col_names
                        df["Time"] = df["Time"].astype(float).round(4)
                        dfs.append(df)
                    except Exception: pass
            
            if dfs:
                combined = pd.concat(dfs, ignore_index=True)
                combined = combined.sort_values(by="Time").drop_duplicates(subset=["Time"], keep="last")
                combined.set_index("Time", inplace=True)
                all_probe_data.append(combined)

        if not all_probe_data: return None

        full_df = pd.concat(all_probe_data, axis=1).ffill()
        output_file = os.path.join(pp_dir, "PROBE_DATA_SUMMARY.csv")
        full_df.to_csv(output_file)

        return output_file

    def get_residuals(self, log_path):
        data = []
        if not os.path.exists(log_path): return pd.DataFrame()

        try:
            with open(log_path, "r") as f:
                current_time = 0.0
                for line in f:
                    if line.startswith("Time = "):
                        try:
                            time_str = line.strip().split()[2].replace('s', '').replace(',', '')
                            current_time = round(float(time_str), 4)
                        except: pass
                    elif "Initial residual" in line and "Solving for" in line:
                        try:
                            parts = line.split(',')
                            var_str = parts[0].split('Solving for ')[1].strip()
                            res_str = parts[1].split('Initial residual = ')[1].strip()
                            data.append({"Time": current_time, "Variable": var_str, "Residual": float(res_str)})
                        except: pass

            if not data: return pd.DataFrame()

            df = pd.DataFrame(data)
            df = df.sort_values(by=["Variable", "Time"]).drop_duplicates(subset=["Time", "Variable"], keep='last')
            
            if len(df) > 50000:
                df = df.iloc[::(len(df)//50000)]
                
            return df
        except: return pd.DataFrame()

    def get_outlet_data(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        if not os.path.exists(pp_dir): return None

        data_frames = {}

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet") or folder.startswith("bounding"):
                sub_dir = os.path.join(pp_dir, folder)
                if os.path.exists(sub_dir):
                    time_dirs = []
                    for d in os.listdir(sub_dir):
                        try:
                            time_dirs.append((float(d), d))
                        except ValueError: pass
                    
                    if time_dirs:
                        time_dirs.sort(key=lambda x: x[0])
                        latest_dir = time_dirs[-1][1]
                        file_path = os.path.join(sub_dir, latest_dir, "surfaceFieldValue.dat")
                        if os.path.exists(file_path):
                            try:
                                df = pd.read_csv(file_path, comment='#', sep=r'\s+', header=None)
                                data_frames[folder] = df.iloc[-1][1]
                            except: pass

        if not data_frames: return None

        total_mass_flow = abs(data_frames.get("outletFlow", abs(data_frames.get("outletTotal", 0.0))))
        safe_total_mass_flow = total_mass_flow if total_mass_flow > 1e-12 else 1e-12

        results = []
        species_molar_flows = {}
        sum_species_fractions = 0.0

        for key, val in data_frames.items():
            if key in ["outletFlow", "outletTotal", "outlet_Air"]: continue
            spec_name = key.replace("outlet_", "")
            mass_fraction = abs(val)
            sum_species_fractions += mass_fraction
            m_dot_spec = mass_fraction 
            species_molar_flows[spec_name] = m_dot_spec / self.MOL_WEIGHTS.get(spec_name, 28.96)

        n_dot_air = (max(0, 1.0 - sum_species_fractions)) / 28.96
        total_molar = sum(species_molar_flows.values()) + n_dot_air
        safe_total_molar = total_molar if total_molar > 1e-12 else 1.0

        for spec, n_dot in species_molar_flows.items():
            results.append({
                "Species": spec,
                "Concentration (ppm)": (n_dot / safe_total_molar) * 1e6,
                "Mass Flow (kg/s)": safe_total_mass_flow * data_frames.get(f"outlet_{spec}", 0),
                "Molar Flow (kmol/s)": n_dot
            })
            
        results.append({
            "Species": "Air",
            "Concentration (ppm)": (n_dot_air / safe_total_molar) * 1e6,
            "Mass Flow (kg/s)": safe_total_mass_flow * max(0, 1.0 - sum_species_fractions),
            "Molar Flow (kmol/s)": n_dot_air
        })

        return pd.DataFrame(results)

    def _cleanup_ghost_time_folders(self, case_path):
        if not os.path.exists(case_path): return

        write_interval = 5.0
        control_dict_path = os.path.join(case_path, "system", "controlDict")
        if os.path.exists(control_dict_path):
            with open(control_dict_path, 'r') as f:
                for line in f:
                    if "writeInterval" in line and not line.strip().startswith(("//", "/*")):
                        try:
                            write_interval = float(line.strip().split()[1].replace(';', ''))
                        except: pass

        time_dirs = []
        for d in os.listdir(case_path):
            full_dir = os.path.join(case_path, d)
            if os.path.isdir(full_dir):
                try:
                    val = float(d)
                    time_dirs.append((val, d, full_dir))
                except ValueError:
                    pass
        
        if not time_dirs: return
        
        time_dirs.sort(key=lambda x: x[0])
        
        gap_threshold = write_interval * 2.5
        last_valid_time = time_dirs[0][0]
        
        for i in range(1, len(time_dirs)):
            prev_time = time_dirs[i-1][0]
            curr_time = time_dirs[i][0]
            
            if (curr_time - prev_time) > gap_threshold:
                for j in range(i, len(time_dirs)):
                    dir_to_remove = time_dirs[j][2]
                    shutil.rmtree(dir_to_remove, ignore_errors=True)
                    if os.path.exists(dir_to_remove):
                        os.system(f"rm -rf {dir_to_remove}")
                break
            else:
                last_valid_time = curr_time

        pp_dir = os.path.join(case_path, "postProcessing")
        if os.path.exists(pp_dir):
            for func_obj in os.listdir(pp_dir):
                func_obj_path = os.path.join(pp_dir, func_obj)
                if os.path.isdir(func_obj_path):
                    for d in os.listdir(func_obj_path):
                        try:
                            time_val = float(d)
                            if time_val > last_valid_time + 1e-5:
                                dir_to_remove = os.path.join(func_obj_path, d)
                                shutil.rmtree(dir_to_remove, ignore_errors=True)
                                if os.path.exists(dir_to_remove):
                                    os.system(f"rm -rf {dir_to_remove}")
                        except ValueError:
                            pass

    def _optimize_control_dict(self, case_path):
        control_dict_path = os.path.join(case_path, "system", "controlDict")
        if not os.path.exists(control_dict_path): return
        with open(control_dict_path, "r") as f:
            content = f.read()

        content = content.replace("writeFormat binary;", "writeFormat ascii;")
        content = content.replace("deltaT 0.1;", "deltaT 0.0001;")
        content = content.replace("stopAt writeNow;", "stopAt endTime;")
        content = content.replace("startFrom startTime;", "startFrom latestTime;")

        with open(control_dict_path, "w") as f:
            f.write(content)

    def _ensure_v13_transport_failsafe(self, case_path):
        constant_dir = os.path.join(case_path, "constant")
        os.makedirs(constant_dir, exist_ok=True)
        
        tt_path = os.path.join(constant_dir, "thermophysicalTransport")
        with open(tt_path, "w") as f:
            f.write("""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object thermophysicalTransport; }

RAS
{
    model   unityLewisEddyDiffusivity;
    Prt     0.85;
    Sct     0.85;
}
""")

        cp_path = os.path.join(constant_dir, "combustionProperties")
        with open(cp_path, "w") as f:
            f.write("""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object combustionProperties; }

combustionModel none;
""")

    def get_available_surface_times(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing", "surfaces1")
        if not os.path.exists(pp_dir): return []
        
        times = []
        for d in os.listdir(pp_dir):
            try:
                t = float(d)
                times.append((t, d))
            except ValueError:
                pass
        times.sort(key=lambda x: x[0])
        return times

    def is_process_running(self, pid):
        if pid is None: return False
        try:
            stat_file = f"/proc/{pid}/stat"
            if os.path.exists(stat_file):
                with open(stat_file, "r") as f:
                    stat_content = f.read().split()
                    if len(stat_content) > 2 and stat_content[2] == 'Z':
                        return False
            os.kill(pid, 0)
            return True
        except OSError: return False
        except Exception: return False

    def run_case_detached(self, case_path, num_cores=1):
        abs_case_path = os.path.abspath(case_path)
        log_path = os.path.join(abs_case_path, "log.run")

        if not os.path.exists(abs_case_path): os.makedirs(abs_case_path)

        self._ensure_v13_transport_failsafe(abs_case_path)
        self._cleanup_ghost_time_folders(abs_case_path)
        self._optimize_control_dict(abs_case_path)

        poly_mesh = os.path.join(abs_case_path, "constant", "polyMesh", "points")
        has_mesh = os.path.exists(poly_mesh)

        foam_source = "source /home/cfd/OpenFOAM/OpenFOAM-13/etc/bashrc"

        if has_mesh:
            setup_cmds = [f"cd {abs_case_path}"]
        else:
            setup_cmds = [
                f"cd {abs_case_path}",
                f"{foam_source} && blockMesh > log.blockMesh 2>&1",
                f"{foam_source} && surfaceFeatures > log.surfaceFeatures 2>&1",
                f"{foam_source} && snappyHexMesh -overwrite > log.snappyHexMesh 2>&1",
                f"{foam_source} && topoSet > log.topoSet 2>&1",
                f"{foam_source} && createPatch -overwrite > log.createPatch 2>&1"
            ]

        # Overwrite the log file (>) so old aborted logs do not cause false "Success" triggers.
        if num_cores > 1:
            setup_cmds.append(f"{foam_source} && decomposePar -force > log.decompose 2>&1")
            solver_cmd = f"{foam_source} && mpirun -np {num_cores} {self.solver} -parallel -fileHandler collated > log.run 2>&1"
        else:
            solver_cmd = f"{foam_source} && {self.solver} > log.run 2>&1"

        full_cmd = " && ".join(setup_cmds + [solver_cmd])

        env_vars = dict(os.environ, MPI_BUFFER_SIZE="20000000", OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")

        process = subprocess.Popen(
            full_cmd, shell=True, executable="/bin/bash",
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True, env=env_vars
        )

        return process.pid, log_path