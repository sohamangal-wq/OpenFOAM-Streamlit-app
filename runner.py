import subprocess
import os
import re
import pandas as pd
import shutil

class SimulationRunner:
    def __init__(self, solver="foamRun"):
        self.solver = solver

        # EXPANDED MOLECULAR WEIGHTS DATABASE
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

        try:
            items = os.listdir(case_path)
            times = [float(x) for x in items if x.replace('.', '', 1).isdigit()]
            if times and max(times) > 190.0: return True
        except: pass

        return False

    def generate_outlet_report(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        if not os.path.exists(pp_dir): return None

        flow_file = os.path.join(pp_dir, "outletFlow", "0", "surfaceFieldValue.dat")
        if not os.path.exists(flow_file): return None

        try:
            df_flow = pd.read_csv(flow_file, sep=r'\s+', comment='#', names=["Time", "Total_Mass_Flow_kg_s"])
            df_flow.set_index("Time", inplace=True)
            df_flow["Total_Mass_Flow_kg_s"] = df_flow["Total_Mass_Flow_kg_s"].abs()
        except: return None

        species_dfs = []
        found_species = []

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet_") and folder != "outletFlow":
                spec = folder.replace("outlet_", "")
                s_file = os.path.join(pp_dir, folder, "0", "surfaceFieldValue.dat")
                if os.path.exists(s_file):
                    try:
                        df_s = pd.read_csv(s_file, sep=r'\s+', comment='#', names=["Time", f"MassFrac_{spec}"])
                        df_s.set_index("Time", inplace=True)
                        species_dfs.append(df_s)
                        found_species.append(spec)
                    except: pass

        if not species_dfs: return None

        full_df = df_flow.join(species_dfs, how='outer').fillna(0.0)
        total_molar_flow = 0.0

        for spec in found_species:
            mw = self.MOL_WEIGHTS.get(spec, 28.96)
            full_df[f"MolarFlow_{spec}kmol_s"] = (full_df["Total_Mass_Flow_kg_s"] * full_df[f"MassFrac_{spec}"]) / mw
            total_molar_flow += full_df[f"MolarFlow_{spec}_kmol_s"]

        sum_mass_frac = sum(full_df[f"MassFrac_{s}"] for s in found_species)
        mass_frac_air = (1.0 - sum_mass_frac).clip(lower=0.0)
        molar_flow_air = (full_df["Total_Mass_Flow_kg_s"] * mass_frac_air) / 28.96
        total_molar_flow += molar_flow_air

        for spec in found_species:
            safe_total = total_molar_flow.replace(0, 1.0)
            full_df[f"Concentration_{spec}ppm"] = (full_df[f"MolarFlow_{spec}_kmol_s"] / safe_total) * 1e6
            full_df[f"MassFlow_{spec}kg_s"] = full_df["Total_Mass_Flow_kg_s"] * full_df[f"MassFrac_{spec}"]

        output_file = os.path.join(pp_dir, "EXIT_DATA_SUMMARY.csv")
        full_df.to_csv(output_file)

        return output_file

    def get_outlet_timeseries(self, case_path):
        """Generates a long-format DataFrame of Concentration over Time for live plotting."""
        pp_dir = os.path.join(case_path, "postProcessing")
        flow_file = os.path.join(pp_dir, "outletFlow", "0", "surfaceFieldValue.dat")
        if not os.path.exists(flow_file): return None

        try:
            df_flow = pd.read_csv(flow_file, sep=r'\s+', comment='#', names=["Time", "Total_Mass_Flow_kg_s"])
            df_flow.set_index("Time", inplace=True)
            df_flow["Total_Mass_Flow_kg_s"] = df_flow["Total_Mass_Flow_kg_s"].abs()
        except: return None

        species_dfs = []
        found_species = []

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet_") and folder != "outletFlow":
                spec = folder.replace("outlet_", "")
                s_file = os.path.join(pp_dir, folder, "0", "surfaceFieldValue.dat")
                if os.path.exists(s_file):
                    try:
                        df_s = pd.read_csv(s_file, sep=r'\s+', comment='#', names=["Time", f"MassFrac_{spec}"])
                        df_s.set_index("Time", inplace=True)
                        species_dfs.append(df_s)
                        found_species.append(spec)
                    except: pass

        if not species_dfs: return None

        full_df = df_flow.join(species_dfs, how='outer').fillna(0.0)
        total_molar_flow = 0.0

        for spec in found_species:
            mw = self.MOL_WEIGHTS.get(spec, 28.96)
            full_df[f"MolarFlow_{spec}"] = (full_df["Total_Mass_Flow_kg_s"] * full_df[f"MassFrac_{spec}"]) / mw
            total_molar_flow += full_df[f"MolarFlow_{spec}"]

        sum_mass_frac = sum(full_df[f"MassFrac_{s}"] for s in found_species)
        mass_frac_air = (1.0 - sum_mass_frac).clip(lower=0.0)
        molar_flow_air = (full_df["Total_Mass_Flow_kg_s"] * mass_frac_air) / 28.96
        total_molar_flow += molar_flow_air

        plot_data = []
        for spec in found_species:
            safe_total = total_molar_flow.replace(0, 1.0)
            conc_ppm = (full_df[f"MolarFlow_{spec}"] / safe_total) * 1e6
            temp_df = pd.DataFrame({
                "Time": full_df.index,
                "Species": spec,
                "Concentration (ppm)": conc_ppm.values
            })
            plot_data.append(temp_df)

        if plot_data:
            return pd.concat(plot_data, ignore_index=True)

        return None

    def generate_probe_report(self, case_path):
        """Parses OpenFOAM probe data files and compiles them into a single clean CSV."""
        pp_dir = os.path.join(case_path, "postProcessing")
        sensors_dir = os.path.join(pp_dir, "gasSensors")

        if not os.path.exists(sensors_dir):
            return None

        time_dirs = [d for d in os.listdir(sensors_dir) if os.path.isdir(os.path.join(sensors_dir, d))]
        if not time_dirs:
            return None

        time_dirs.sort(key=float)
        active_sensor_dir = os.path.join(sensors_dir, time_dirs[0])
        all_probe_data = []

        for file_name in os.listdir(active_sensor_dir):
            file_path = os.path.join(active_sensor_dir, file_name)
            if os.path.isfile(file_path):
                try:
                    df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)
                    num_probes = df.shape[1] - 1
                    col_names = ["Time"] + [f"Probe_{i}_{file_name}" for i in range(num_probes)]
                    df.columns = col_names
                    df.set_index("Time", inplace=True)
                    all_probe_data.append(df)
                except Exception: pass

        if not all_probe_data:
            return None

        full_df = pd.concat(all_probe_data, axis=1).ffill()
        output_file = os.path.join(pp_dir, "PROBE_DATA_SUMMARY.csv")
        full_df.to_csv(output_file)

        return output_file

    def _optimize_control_dict(self, case_path):
        control_dict_path = os.path.join(case_path, "system", "controlDict")
        if not os.path.exists(control_dict_path): return
        with open(control_dict_path, "r") as f:
            content = f.read()

        content = content.replace("writeFormat binary;", "writeFormat ascii;")
        content = content.replace("deltaT 0.1;", "deltaT 0.0001;")

        with open(control_dict_path, "w") as f:
            f.write(content)

    def run_case_detached(self, case_path, num_cores=1):
        abs_case_path = os.path.abspath(case_path)
        log_path = os.path.join(abs_case_path, "log.run")

        if not os.path.exists(abs_case_path): os.makedirs(abs_case_path)

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

        if num_cores > 1:
            solver_cmd = f"{foam_source} && mpirun -np {num_cores} {self.solver} -parallel -fileHandler collated > log.run 2>&1"
        else:
            solver_cmd = f"{foam_source} && {self.solver} > log.run 2>&1"

        full_cmd = " && ".join(setup_cmds + [solver_cmd])

        process = subprocess.Popen(
            full_cmd, shell=True, executable="/bin/bash",
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True, env=dict(os.environ, MPI_BUFFER_SIZE="20000000")
        )

        return process.pid, log_path

    def is_process_running(self, pid):
        if pid is None: return False
        try:
            os.kill(pid, 0)
            return True
        except OSError: return False

    def get_residuals(self, log_path):
        data = []
        if not os.path.exists(log_path): return pd.DataFrame()

        try:
            with open(log_path, "r") as f:
                lines = f.readlines()[-1000:]
            current_time = 0.0

            for line in lines:
                t_match = re.search(r"^Time = ([0-9.eE+-]+)", line)
                if t_match: current_time = float(t_match.group(1))

                if "Solving for" in line and "Initial residual" in line:
                    r_match = re.search(r"Solving for ([a-zA-Z0-9_]+).*, Initial residual = ([0-9.eE+-]+)", line)
                    if r_match:
                        data.append({"Time": current_time, "Variable": r_match.group(1), "Residual": float(r_match.group(2))})

            if not data: return pd.DataFrame()

            return pd.DataFrame(data).drop_duplicates(subset=["Time", "Variable"], keep='last')
        except: return pd.DataFrame()

    def get_outlet_data(self, case_path):
        pp_dir = os.path.join(case_path, "postProcessing")
        if not os.path.exists(pp_dir): return None

        data_frames = {}

        for folder in os.listdir(pp_dir):
            if folder.startswith("outlet") or folder.startswith("bounding"):
                sub_dir = os.path.join(pp_dir, folder)
                if os.path.exists(sub_dir):
                    time_dirs = sorted([d for d in os.listdir(sub_dir) if d.isdigit() or d=='0'])
                    if time_dirs:
                        file_path = os.path.join(sub_dir, time_dirs[-1], "surfaceFieldValue.dat")
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
            if key in ["outletFlow", "outletTotal"]: continue
            spec_name = key.replace("outlet_", "")
            mass_fraction = abs(val)
            sum_species_fractions += mass_fraction
            m_dot_spec = safe_total_mass_flow * mass_fraction
            species_molar_flows[spec_name] = m_dot_spec / self.MOL_WEIGHTS.get(spec_name, 28.96)

        n_dot_air = (safe_total_mass_flow * max(0, 1.0 - sum_species_fractions)) / 28.96
        total_molar = sum(species_molar_flows.values()) + n_dot_air
        safe_total_molar = total_molar if total_molar > 1e-12 else 1.0

        for spec, n_dot in species_molar_flows.items():
            results.append({
                "Species": spec,
                "Concentration (ppm)": (n_dot / safe_total_molar) * 1e6,
                "Mass Flow (kg/s)": n_dot * self.MOL_WEIGHTS.get(spec, 28.96),
                "Molar Flow (kmol/s)": n_dot
            })

        return pd.DataFrame(results)