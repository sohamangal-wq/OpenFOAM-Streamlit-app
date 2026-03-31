import os
import shutil
import pandas as pd
import numpy as np
import math
import json

class CaseBuilder:

    def __init__(self, base_template):
        self.template_path = base_template
        self.output_root = "simulation_runs"
        self.vent_inlet = "inlet"
        self.vent_outlet = "outlet"
        if not os.path.exists(self.output_root): os.makedirs(self.output_root)

        # --- HAZARD LIMITS (PPM) FOR PLOTTING (Graded Alarms) ---
        self.HAZARD_LIMITS = {
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
        
        self.NAME_MAP = {
            "METHANE": "CH4", "PROPANE": "C3H8", "HYDROGEN": "H2",
            "AMMONIA": "NH3", "CHLORINE": "Cl2", "BENZENE": "C6H6",
            "TOLUENE": "C7H8", "HEXANE": "C6H14", "CARBON MONOXIDE": "CO",
            "HYDROGEN SULFIDE": "H2S", "H2S": "H2S", "CO": "CO", "CH4": "CH4", "H2": "H2", "HE": "He", "HELIUM": "He"
        }

        # EXPANDED JANAF COEFFICIENTS DATABASE
        self.COEFFS = {
            "Air": { "molWeight": 28.96, "high": "(3.0879271 0.001245971 -4.237188e-07 6.747207e-11 -3.97077e-15 -995.2627 5.959609)", "low": "(3.568396 -0.0006787294 1.5537e-06 -3.29937e-12 -4.686e-12 -995.2627 3.67482)" },
            "N2": { "molWeight": 28.01, "high": "(2.92664 0.0014879768 -5.68476e-07 1.0097038e-10 -6.753351e-15 -922.7977 5.980528)", "low": "(3.298677 0.0014082404 -3.9632222e-06 5.641515e-09 -2.444854e-12 -1020.8999 3.950372)" },
            "O2": { "molWeight": 31.99, "high": "(3.28253784 0.00148308754 -7.57966669e-07 2.09470555e-10 -2.16717794e-14 -1088.45772 5.45323129)", "low": "(3.78245636 -0.00299673416 9.84730201e-06 -9.68129509e-09 3.24372837e-12 -1063.94356 3.65767573)" },
            "CO2": { "molWeight": 44.01, "high": "(4.4608041 0.0030981719 -1.2393358e-06 2.2741325e-10 -1.552625e-14 -48966.961 -0.9818554)", "low": "(2.3567735 0.0089845967 -7.12356269e-06 2.45919022e-09 -1.43699548e-13 -48371.9698 9.90105222)" },
            "H2O": { "molWeight": 18.02, "high": "(2.6721459 0.003056293 -8.73026e-07 1.2009964e-10 -6.391618e-15 -29899.2 6.862817)", "low": "(3.3868426 0.0034749825 -6.354696e-06 6.9685812e-09 -2.5065885e-12 -30208.11 2.5902328)" },
            "H2": { "molWeight": 2.016, "high": "(2.9328305 0.00082660796 -1.4640233e-07 1.5410035e-11 -6.8880443e-16 -813.06559 -1.0243289)", "low": "(2.3443311 0.0079805208 -1.9478151e-05 2.0157209e-08 -7.3761176e-12 -917.93517 0.68301024)" },
            "CH4": { "molWeight": 16.04, "high": "(0.074851495 0.0133909467 -5.73285809e-06 1.22292535e-09 -1.01815251e-13 -9468.34459 18.437318)", "low": "(5.14987613 -0.0136709788 4.91800599e-05 -4.84743026e-08 1.66693956e-11 -10246.6476 -4.64130376)" },
            "C2H6": { "molWeight": 30.07, "high": "(4.046666 0.01931326 -5.336044e-06 6.559929e-10 -3.0037e-14 -11939.9 2.197968)", "low": "(1.45832 0.0151121 -9.56013e-07 -6.2239e-09 2.85966e-12 -11075.7 16.1266)" },
            "C3H8": { "molWeight": 44.1, "high": "(7.5341368 0.018872239 -6.2718491e-06 9.1475649e-10 -4.7838069e-14 -16467.516 -17.892349)", "low": "(0.93355381 0.033308552 -1.2765203e-05 1.7153859e-09 0.0 -15652.98 13.53663)" },
            "C4H10": { "molWeight": 58.12, "high": "(9.578335 0.02672522 -9.32717e-06 1.419266e-09 -7.67512e-14 -19999.5 -26.9749)", "low": "(1.93608 0.0392336 -1.13982e-05 -4.3496e-09 3.03714e-12 -17537.1 13.8967)" },
            "H2S": { "molWeight": 34.08, "high": "(2.80905 0.0039266 -1.11666e-06 1.39655e-10 -6.1155e-15 -3595.66 7.42655)", "low": "(3.93181 -0.0014299 5.25332e-06 -3.3721e-09 6.27285e-13 -3656.62 2.87539)" },
            "CO": { "molWeight": 28.01, "high": "(3.0484859 0.0013517281 -4.8579405e-07 7.8853644e-11 -4.6980746e-15 -14266.117 6.0170977)", "low": "(3.5795334 -0.00061035368 1.0168143e-06 9.0700588e-10 -9.0442449e-13 -14344.086 3.5084092)" },
            "SO2": { "molWeight": 64.06, "high": "(4.65614 0.00277322 -9.79979e-07 1.54228e-10 -8.82522e-15 -37678.9 1.56561)", "low": "(2.76679 0.00767073 -8.49071e-06 4.3023e-09 -7.86976e-13 -37073.4 11.5872)" },
            "Cl2": { "molWeight": 70.9, "high": "(4.37927 0.00018805 -7.77353e-08 1.44218e-11 -9.57014e-16 -1257.26 1.3323)", "low": "(3.41909 0.0029379 -4.26979e-06 3.03714e-09 -7.9576e-13 -955.517 6.45241)" },
            "NH3": { "molWeight": 17.03, "high": "(2.6344521 0.005666256 -1.7278676e-06 2.3867161e-10 -1.2578786e-14 -6544.6958 -1.005811)", "low": "(4.2862486 -0.0045608234 2.0298253e-05 -1.0372175e-08 1.7511429e-12 -6985.9868 -7.0510864)" },
            "C6H6": { "molWeight": 78.11, "high": "(11.1606 0.0309968 -1.0964e-05 1.7067e-09 -9.739e-14 6269.4 -36.23)", "low": "(-0.5042 0.05777 -3.968e-05 1.341e-08 -1.725e-12 9904.7 26.24)" },
            "C7H8": { "molWeight": 92.14, "high": "(13.73 0.0359 -1.27e-05 1.98e-09 -1.13e-13 1335.0 -47.9)", "low": "(-1.75 0.0765 -5.67e-05 2.05e-08 -2.85e-12 5900.0 32.8)" },
            "C6H14": { "molWeight": 86.18, "high": "(13.68 0.0398 -1.39e-05 2.15e-09 -1.22e-13 -24700.0 -49.0)", "low": "(1.85 0.0583 -3.22e-05 8.35e-09 -3.15e-13 -21200.0 15.5)" },
            "CH3OH": { "molWeight": 32.04, "high": "(3.69 0.0093 -3.07e-06 4.63e-10 -2.60e-14 -25700.0 4.16)", "low": "(1.80 0.0125 -4.05e-06 7.62e-10 -5.00e-14 -25000.0 12.0)" },
            "NO2": { "molWeight": 46.005, "high": "(4.0 0.0025 -9.6e-07 1.6e-10 -1.0e-14 3600.0 3.5)", "low": "(2.8 0.0055 -3.2e-06 1.0e-09 -1.2e-13 3900.0 9.0)" },
            "HCN": { "molWeight": 27.025, "high": "(3.802 0.0029 -1.0e-06 1.6e-10 -9.0e-15 15900.0 1.5)", "low": "(2.25 0.0075 -5.3e-06 2.0e-09 -3.0e-13 16400.0 8.8)" }
        }

    def _get_coeffs(self, species_name):
        return self.COEFFS.get(species_name, self.COEFFS["CH4"])

    def convert_hvac_units(self, cfm, temp_f, press_in_wc, vent_area_sqft):
        flow_m3_s = cfm * 0.0283168 / 60.0
        temp_k = (temp_f - 32.0) * (5.0/9.0) + 273.15
        p_gauge_pa = press_in_wc * 249.089
        p_abs_pa = 101325 + p_gauge_pa
        return flow_m3_s, temp_k, p_abs_pa

    def get_mass_frac_hazard(self, species_name, limit_type="low"):
        limits = self.HAZARD_LIMITS.get(species_name, {"visible": 100, "low": 5000, "high": 25000, "emergency": 50000})
        vol_frac = limits.get(limit_type, 5000) / 1000000.0
        mw_gas = self.MOL_WEIGHTS.get(species_name, 16.04)
        mw_air = 28.96
        mass_frac = (vol_frac * mw_gas) / ((vol_frac * mw_gas) + ((1.0 - vol_frac) * mw_air))
        return mass_frac

    def create_case_from_group(self, case_name, stl_file, leak_group_df, mesh_res="Fast", refine_lvl=1, num_cores=1,
                               sim_mode="Closed Room CFD", hvac_cfm=1, hvac_temp_f=72, hvac_press_inwc=0.01, hvac_area_sqft=1.0,
                               vent_inlet_center=None, vent_outlet_center=None,
                               wind_speed=5.0, wind_dir=0, stability_class="D", z0=0.1,
                               probe_locations=None, sim_duration=200, max_courant=1.0, delta_t=0.0001):

        case_path = os.path.join(self.output_root, case_name)
        if os.path.exists(case_path): 
            shutil.rmtree(case_path, ignore_errors=True)
            
        shutil.copytree(self.template_path, case_path, dirs_exist_ok=True)

        zero_dir = os.path.join(case_path, "0")
        if os.path.exists(zero_dir):
            valid_base_files = ["U", "p", "p_rgh", "T", "k", "epsilon", "omega", "nut", "alphat", "Y_default"]
            for f in os.listdir(zero_dir):
                if os.path.isfile(os.path.join(zero_dir, f)) and f not in valid_base_files:
                    os.remove(os.path.join(zero_dir, f))

        with open(os.path.join(case_path, f"{case_name}.foam"), "w") as f:
            f.write("")

        # Save config meta to position labels accurately in the 3D Viewer later
        cfg_dict = {"inlet": vent_inlet_center, "outlet": vent_outlet_center}
        with open(os.path.join(case_path, "case_config.json"), "w") as f:
            json.dump(cfg_dict, f)

        flow_rate_val, t_val, p_val = self.convert_hvac_units(hvac_cfm, hvac_temp_f, hvac_press_inwc, hvac_area_sqft)
        wind_dir_rad = math.radians(wind_dir)
        dir_x = -math.sin(wind_dir_rad)
        dir_y = -math.cos(wind_dir_rad)

        tri_path = os.path.join(case_path, "constant", "triSurface")
        if not os.path.exists(tri_path): os.makedirs(tri_path)

        if stl_file:
            stl_data = stl_file.getvalue() if hasattr(stl_file, "getvalue") else stl_file.read()
            with open(os.path.join(tri_path, "room.stl"), "wb") as f:
                f.write(stl_data)

        all_scenario_species = set(["Air"])
        
        injections = {}
        is_hydrogen_case = False
        leak_meta = []
        max_exhaust_time = 0.0

        for index, row in leak_group_df.iterrows():
            comp_raw = str(row.get('COMPOSITION', "CH4:1.0")).strip()
            if comp_raw == "nan" or comp_raw == "" or comp_raw.lower() == "none":
                comp_raw = "CH4:1.0"
                
            if ":" not in comp_raw:
                comp_raw = f"{comp_raw}:1.0"
                
            components = []
            try:
                pairs = comp_raw.split(';')
                for p in pairs:
                    spec_raw, frac = p.split(':')
                    spec_clean = spec_raw.strip().upper()
                    spec = self.NAME_MAP.get(spec_clean, spec_raw.strip())
                    components.append((spec, float(frac)))
                    all_scenario_species.add(spec)
                    if spec in ["H2", "He"]: is_hydrogen_case = True
            except:
                components = [("CH4", 1.0)]
                all_scenario_species.add("CH4")

            total_rate_kgs = float(row.get('RATE', 0.0))
            
            inv_lb = float(row.get('INVENTORY_LB', 0.0)) if 'INVENTORY_LB' in leak_group_df.columns else 0.0
            inv_kmol = float(row.get('INVENTORY_KMOL', 0.0)) if 'INVENTORY_KMOL' in leak_group_df.columns else 0.0
            
            x_coord = float(row.get('X', 0.0))
            y_coord = float(row.get('Y', 0.0))
            z_coord = float(row.get('Z', 0.0))

            for spec, frac in components:
                rate = total_rate_kgs * frac
                spec_inv_kg = 0.0
                
                # Convert inventory to kg for OpenFOAM internal calculation
                if inv_lb > 0:
                    spec_inv_kg = inv_lb * 0.453592 * frac
                elif inv_kmol > 0:
                    mw = self.MOL_WEIGHTS.get(spec, 28.96)
                    spec_inv_kg = inv_kmol * frac * mw
                    
                key = f"{spec}_{x_coord}_{y_coord}_{z_coord}"
                
                if key not in injections:
                    injections[key] = {
                        "species": spec, 
                        "x": x_coord, 
                        "y": y_coord, 
                        "z": z_coord, 
                        "rates": [],
                        "is_blowdown": False
                    }
                    
                if spec_inv_kg > 0 and rate > 0:
                    tau = spec_inv_kg / rate
                    table_rates = []
                    for i in range(21):
                        t_i = (i / 20.0) * (5.0 * tau)
                        r_i = rate * math.exp(-t_i / tau)
                        table_rates.append((t_i, r_i))
                    table_rates.append((5.0 * tau + 0.1, 0.0)) 
                    
                    injections[key]["rates"] = table_rates
                    injections[key]["is_blowdown"] = True
                    
                    leak_meta.append({"species": spec, "M0": spec_inv_kg, "tau": tau})
                    if (5.0 * tau) > max_exhaust_time:
                        max_exhaust_time = 5.0 * tau
                else:
                    injections[key]["rates"].append((0.0, rate))
                    
        if max_exhaust_time > 0:
            sim_duration = max_exhaust_time * 1.2 
            
        with open(os.path.join(case_path, "leak_meta.json"), "w") as mf:
            json.dump(leak_meta, mf)

        self._edit_fvOptions_multi(case_path, list(injections.values()))

        primary_gas = list(injections.values())[0]['species'] if injections else "CH4"
        leak_loc = (list(injections.values())[0]['x'], list(injections.values())[0]['y'], list(injections.values())[0]['z']) if injections else (0,0,0)

        self._rewrite_controlDict(
            case_path, self.vent_inlet, self.vent_outlet, list(all_scenario_species),
            is_light_gas=is_hydrogen_case, probe_locations=probe_locations,
            sim_duration=sim_duration, max_courant=max_courant, delta_t=delta_t,
            primary_gas=primary_gas, leak_loc=leak_loc
        )

        self._edit_thermo_multi(case_path, list(all_scenario_species))
        self._edit_fvSolution(case_path, mode=mesh_res, is_light_gas=is_hydrogen_case)

        if sim_mode == "Closed Room CFD":
            if vent_inlet_center and vent_outlet_center:
                self._edit_blockMesh_split(case_path, mesh_res)
                self._write_topoSet_createPatch_subset(case_path, vent_inlet_center, vent_outlet_center)
            else:
                self._edit_blockMesh(case_path, mesh_res)
        else:
            self._edit_blockMesh_outdoor(case_path, mesh_res)

        first_leak_inj = [{'x': list(injections.values())[0]['x'], 'y': list(injections.values())[0]['y'], 'z': list(injections.values())[0]['z']}] if injections else []
        self._edit_snappyHexMesh(case_path, first_leak_inj, refine_lvl, mode=mesh_res, sim_mode=sim_mode)
        self._write_decomposeParDict(case_path, num_cores)

        # WRITE ZERO DIRECTORY ABSOLUTE PHYSICS (ENFORCES IMPERMEABLE SOLID WALLS)
        self._setup_physics_fields(case_path, sim_mode, flow_rate_val, t_val, p_val, wind_speed, dir_x, dir_y, z0)
        self._setup_species_fields_multi(case_path, list(all_scenario_species))

    def _setup_physics_fields(self, case_path, sim_mode, flow_rate_val, t_val, p_val, wind_speed, dir_x, dir_y, z0):
        # Physics Base Header
        header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/\n"""

        # --- 0/U ---
        u_content = header + """FoamFile { version 2.0; format ascii; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
"""
        if sim_mode == "Closed Room CFD":
            u_content += f"""    {self.vent_inlet} {{ type flowRateInletVelocity; volumetricFlowRate constant {flow_rate_val:.6f}; value uniform (0 0 0); }}
    {self.vent_outlet} {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
    ".*" {{ type noSlip; }}
}}"""
        else:
            u_content += f"""    inlet {{ type atmBoundaryLayerInletVelocity; flowDir ({dir_x} {dir_y} 0); zDir (0 0 1); Uref {wind_speed}; Zref 10.0; z0 uniform {z0}; zGround uniform 0.0; value uniform ({dir_x * wind_speed} {dir_y * wind_speed} 0); }}
    outlet {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
    sky {{ type slip; }}
    ground {{ type noSlip; }}
    ".*" {{ type noSlip; }}
}}"""
        with open(os.path.join(case_path, "0", "U"), "w") as f: f.write(u_content)

        # --- 0/p_rgh ---
        p_rgh_content = header + f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object p_rgh; }}
dimensions [1 -1 -2 0 0 0 0];
internalField uniform {p_val:.1f};
boundaryField
{{
"""
        if sim_mode == "Closed Room CFD":
            p_rgh_content += f"""    {self.vent_inlet} {{ type zeroGradient; }}
    {self.vent_outlet} {{ type fixedValue; value uniform {p_val:.1f}; }}
    ".*" {{ type fixedFluxPressure; value uniform {p_val:.1f}; }}
}}"""
        else:
            p_rgh_content += f"""    inlet {{ type zeroGradient; }}
    outlet {{ type fixedValue; value uniform 101325; }}
    sky {{ type fixedValue; value uniform 101325; }}
    ground {{ type fixedFluxPressure; value uniform 101325; }}
    ".*" {{ type fixedFluxPressure; value uniform 101325; }}
}}"""
        with open(os.path.join(case_path, "0", "p_rgh"), "w") as f: f.write(p_rgh_content)

        # --- 0/p ---
        p_content = header + f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object p; }}
dimensions [1 -1 -2 0 0 0 0];
internalField uniform {p_val:.1f};
boundaryField
{{
    ".*" {{ type calculated; value uniform {p_val:.1f}; }}
}}"""
        with open(os.path.join(case_path, "0", "p"), "w") as f: f.write(p_content)

        # --- 0/T ---
        t_content = header + f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object T; }}
dimensions [0 0 0 1 0 0 0];
internalField uniform {t_val:.2f};
boundaryField
{{
"""
        if sim_mode == "Closed Room CFD":
            t_content += f"""    {self.vent_inlet} {{ type fixedValue; value uniform {t_val:.2f}; }}
    {self.vent_outlet} {{ type inletOutlet; inletValue uniform {t_val:.2f}; value uniform {t_val:.2f}; }}
    ".*" {{ type zeroGradient; }}
}}"""
        else:
            t_content += f"""    inlet {{ type fixedValue; value uniform 300; }}
    outlet {{ type inletOutlet; inletValue uniform 300; value uniform 300; }}
    sky {{ type inletOutlet; inletValue uniform 300; value uniform 300; }}
    ground {{ type zeroGradient; }}
    ".*" {{ type zeroGradient; }}
}}"""
        with open(os.path.join(case_path, "0", "T"), "w") as f: f.write(t_content)

        # --- 0/k ---
        k_content = header + """FoamFile { version 2.0; format ascii; class volScalarField; object k; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0.01;
boundaryField
{
"""
        if sim_mode == "Closed Room CFD":
            k_content += f"""    {self.vent_inlet} {{ type fixedValue; value uniform 0.01; }}
    {self.vent_outlet} {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    ".*" {{ type kqRWallFunction; value uniform 0.01; }}
}}"""
        else:
            k_content += f"""    inlet {{ type atmBoundaryLayerInletK; flowDir ({dir_x} {dir_y} 0); zDir (0 0 1); Uref {wind_speed}; Zref 10.0; z0 uniform {z0}; zGround uniform 0.0; value uniform 0.1; }}
    outlet {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    sky {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    ground {{ type kqRWallFunction; value uniform 0.01; }}
    ".*" {{ type kqRWallFunction; value uniform 0.01; }}
}}"""
        with open(os.path.join(case_path, "0", "k"), "w") as f: f.write(k_content)

        # --- 0/epsilon ---
        eps_content = header + """FoamFile { version 2.0; format ascii; class volScalarField; object epsilon; }
dimensions [0 2 -3 0 0 0 0];
internalField uniform 0.01;
boundaryField
{
"""
        if sim_mode == "Closed Room CFD":
            eps_content += f"""    {self.vent_inlet} {{ type fixedValue; value uniform 0.01; }}
    {self.vent_outlet} {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    ".*" {{ type epsilonWallFunction; value uniform 0.01; }}
}}"""
        else:
            eps_content += f"""    inlet {{ type atmBoundaryLayerInletEpsilon; flowDir ({dir_x} {dir_y} 0); zDir (0 0 1); Uref {wind_speed}; Zref 10.0; z0 uniform {z0}; zGround uniform 0.0; value uniform 0.01; }}
    outlet {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    sky {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    ground {{ type epsilonWallFunction; value uniform 0.01; }}
    ".*" {{ type epsilonWallFunction; value uniform 0.01; }}
}}"""
        with open(os.path.join(case_path, "0", "epsilon"), "w") as f: f.write(eps_content)
        
        # --- 0/nut ---
        nut_content = header + """FoamFile { version 2.0; format ascii; class volScalarField; object nut; }
dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
"""
        if sim_mode == "Closed Room CFD":
            nut_content += f"""    {self.vent_inlet} {{ type calculated; value uniform 0; }}
    {self.vent_outlet} {{ type calculated; value uniform 0; }}
    ".*" {{ type nutkWallFunction; value uniform 0; }}
}}"""
        else:
            nut_content += f"""    inlet {{ type calculated; value uniform 0; }}
    outlet {{ type calculated; value uniform 0; }}
    sky {{ type calculated; value uniform 0; }}
    ground {{ type nutkWallFunction; value uniform 0; }}
    ".*" {{ type nutkWallFunction; value uniform 0; }}
}}"""
        with open(os.path.join(case_path, "0", "nut"), "w") as f: f.write(nut_content)

        # --- 0/alphat ---
        alphat_content = header + """FoamFile { version 2.0; format ascii; class volScalarField; object alphat; }
dimensions [1 -1 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
"""
        if sim_mode == "Closed Room CFD":
            alphat_content += f"""    {self.vent_inlet} {{ type calculated; value uniform 0; }}
    {self.vent_outlet} {{ type calculated; value uniform 0; }}
    ".*" {{ type calculated; value uniform 0; }}
}}"""
        else:
            alphat_content += f"""    inlet {{ type calculated; value uniform 0; }}
    outlet {{ type calculated; value uniform 0; }}
    sky {{ type calculated; value uniform 0; }}
    ground {{ type calculated; value uniform 0; }}
    ".*" {{ type calculated; value uniform 0; }}
}}"""
        with open(os.path.join(case_path, "0", "alphat"), "w") as f: f.write(alphat_content)


    def _edit_blockMesh_outdoor(self, case_path, res):
        cells = 30 if res == "Fast" else 50
        block_mesh_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

convertToMeters 1;

vertices
(
    (-100 -100 0) ( 100 -100 0) ( 100 100 0) (-100 100 0)
    (-100 -100 50) ( 100 -100 50) ( 100 100 50) (-100 100 50)
);

blocks ( hex (0 1 2 3 4 5 6 7) ({cells} {cells} {cells}) simpleGrading (1 1 1) );

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces ( (0 4 7 3) (0 1 5 4) );
    }}
    outlet
    {{
        type patch;
        faces ( (1 2 6 5) (2 3 7 6) );
    }}
    ground
    {{
        type wall;
        faces ( (0 3 2 1) );
    }}
    sky
    {{
        type patch;
        faces ( (4 5 6 7) );
    }}
);

mergePatchPairs ();
"""
        with open(os.path.join(case_path, "system", "blockMeshDict"), "w") as f:
            f.write(block_mesh_content)

    def _edit_fvSolution(self, case_path, mode="Fast", is_light_gas=False):
        correctors = "2" if is_light_gas else ("1" if mode == "Fast" else "3")
        path = os.path.join(case_path, "system", "fvSolution")
        if os.path.exists(path):
            with open(path, 'r') as f: content = f.read()
            content = content.replace("nOuterCorrectors 3;", f"nOuterCorrectors {correctors};")
            with open(path, 'w') as f: f.write(content)

    def _edit_blockMesh_split(self, case_path, res):
        cells = 20 if res == "Fast" else 40
        # Replaces generic 35 resolution with requested
        path = os.path.join(case_path, "system", "blockMeshDict")
        if os.path.exists(path):
            with open(path, 'r') as f: content = f.read()
            content = content.replace("(35 35 35)", f"({cells} {cells} {cells})")
            with open(path, 'w') as f: f.write(content)

    def _edit_blockMesh(self, case_path, res):
        cells = 20 if res == "Fast" else 40
        path = os.path.join(case_path, "system", "blockMeshDict")
        if os.path.exists(path):
            with open(path, 'r') as f: content = f.read()
            content = content.replace("(35 35 35)", f"({cells} {cells} {cells})")
            with open(path, 'w') as f: f.write(content)

    def _write_topoSet_createPatch_subset(self, case_path, in_center, out_center):
        box_size = 0.5
        header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/\n"""

        topo_content = header + f"""FoamFile {{ version 2.0; format ascii; class dictionary; object topoSetDict; }}
actions (
{{ name inletFaceSetRaw; type faceSet; action new; source boxToFace; sourceInfo {{ box ({in_center[0]-box_size} {in_center[1]-box_size} {in_center[2]-box_size}) ({in_center[0]+box_size} {in_center[1]+box_size} {in_center[2]+box_size}); }} }}
{{ name inletFaceSet; type faceSet; action new; source faceToFace; set inletFaceSetRaw; }}
{{ name inletFaceSet; type faceSet; action subset; source patchToFace; sourceInfo {{ name room; }} }}
{{ name outletFaceSetRaw; type faceSet; action new; source boxToFace; sourceInfo {{ box ({out_center[0]-box_size} {out_center[1]-box_size} {out_center[2]-box_size}) ({out_center[0]+box_size} {out_center[1]+box_size} {out_center[2]+box_size}); }} }}
{{ name outletFaceSet; type faceSet; action new; source faceToFace; set outletFaceSetRaw; }}
{{ name outletFaceSet; type faceSet; action subset; source patchToFace; sourceInfo {{ name room; }} }}
);
"""
        with open(os.path.join(case_path, "system", "topoSetDict"), "w") as f:
            f.write(topo_content)

        patch_content = header + """FoamFile { version 2.0; format ascii; class dictionary; object createPatchDict; }

pointSync false;

patches (
{ name inlet; patchInfo { type patch; } constructFrom set; set inletFaceSet; }
{ name outlet; patchInfo { type patch; } constructFrom set; set outletFaceSet; }
);
"""
        with open(os.path.join(case_path, "system", "createPatchDict"), "w") as f:
            f.write(patch_content)

    def _edit_snappyHexMesh(self, case_path, injections, refine_lvl, mode="Fast", sim_mode="Closed Room CFD"):
        if not injections: return

        lx, ly, lz = injections[0]['x'], injections[0]['y'], injections[0]['z']
        box_min = f"{lx-0.25} {ly-0.25} {lz-0.25}"
        box_max = f"{lx+0.25} {ly+0.25} {lz+0.25}"

        surface_lvl = "(1 1)" if mode == "Fast" else "(2 2)"
        region_lvl = "1" if mode == "Fast" else str(refine_lvl)

        refinement_block = f"""
geometry
{{
    room.stl {{ type triSurfaceMesh; file "room.stl"; name room; }}
    refineBox {{ type searchableBox; min ({box_min}); max ({box_max}); }}
}}

castellatedMeshControls
{{
    features ( {{ file "room.eMesh"; level 1; }} );
    refinementSurfaces {{ room {{ level {surface_lvl}; patchInfo {{ type wall; }} }} }}
    refinementRegions {{ refineBox {{ mode inside; levels ((1.0 {region_lvl})); }} }}

    resolveFeatureAngle 30; insidePoint (2.0 2.0 2.0);
    maxLocalCells 1000000; maxGlobalCells 2000000; minRefinementCells 10;
    maxLoadUnbalance 0.10; nCellsBetweenLevels 3; allowFreeStandingZoneFaces true;
}}
"""
        path = os.path.join(case_path, "system", "snappyHexMeshDict")
        with open(path, 'r') as f: content = f.read()

        start_idx = content.find("geometry {")
        end_idx = content.find("snapControls")

        if start_idx != -1 and end_idx != -1:
            new_content = content[:start_idx] + refinement_block + "\n" + content[end_idx:]
            with open(path, 'w') as f: f.write(new_content)

    def _write_decomposeParDict(self, case_path, num_cores):
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object decomposeParDict; }}

numberOfSubdomains {num_cores};
method scotch;
"""
        with open(os.path.join(case_path, "system", "decomposeParDict"), "w") as f:
            f.write(content)

    def _edit_fvOptions_multi(self, case_path, injections):
        options_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }
"""
        for i, inj in enumerate(injections):
            rates = sorted(inj['rates'], key=lambda x: x[0])
            if inj.get('is_blowdown', False) or len(rates) > 1:
                table_str = " ".join([f"({t:.3f} {r:.6f})" for t, r in rates])
                rate_str = f"explicit table ({table_str}); implicit 0;"
            else:
                rate_str = f"explicit {rates[0][1]}; implicit 0;"

            options_content += f"""
leakSource_{i}_{inj['species']}
{{
    type semiImplicitSource;
    active true;
    selectionMode cellZone;
    cellZone
    {{
        type containsPoints;
        points ( ({inj['x']} {inj['y']} {inj['z']}) );
    }}
    volumeMode absolute;
    sources
    {{
        {inj['species']} {{ {rate_str} }}
    }}
}}
"""
        with open(os.path.join(case_path, "system", "fvOptions"), "w") as f:
            f.write(options_content)

    def _rewrite_controlDict(self, case_path, inlet_name, outlet_name, species_list,
                             is_light_gas=False, probe_locations=None, sim_duration=200, max_courant=1.0,
                             delta_t=0.0001, primary_gas="CH4", leak_loc=(0,0,0)):

        field_avg_str = f"""
    fieldAverage1
    {{
        type fieldAverage;
        libs ("libfieldFunctionObjects.so");
        writeControl writeTime;
        fields
        (
"""
        for spec in species_list:
            if spec != "Air":
                field_avg_str += f"            {spec} {{ mean on; prime2Mean off; base time; }}\n"

        field_avg_str += """        );
    }"""

        hazard_mass_frac_vis = self.get_mass_frac_hazard(primary_gas, "visible")
        hazard_mass_frac_low = self.get_mass_frac_hazard(primary_gas, "low")
        hazard_mass_frac_high = self.get_mass_frac_hazard(primary_gas, "high")
        hazard_mass_frac_emerg = self.get_mass_frac_hazard(primary_gas, "emergency")

        surfaces_str = f"""
    surfaces1
    {{
        type surfaces;
        libs ("libsampling.so");
        writeControl writeTime;
        surfaceFormat vtk;
        fields ( U p {primary_gas} );
        interpolationScheme cellPoint;
        surfaces
        (
            Hazard_Cloud_Visible
            {{
                type isoSurface;
                isoField {primary_gas};
                isoValue {hazard_mass_frac_vis};
                interpolate true;
            }}
            Hazard_Cloud_Low
            {{
                type isoSurface;
                isoField {primary_gas};
                isoValue {hazard_mass_frac_low};
                interpolate true;
            }}
            Hazard_Cloud_High
            {{
                type isoSurface;
                isoField {primary_gas};
                isoValue {hazard_mass_frac_high};
                interpolate true;
            }}
            Hazard_Cloud_Emergency
            {{
                type isoSurface;
                isoField {primary_gas};
                isoValue {hazard_mass_frac_emerg};
                interpolate true;
            }}
            breathing_zone
            {{
                type cuttingPlane;
                planeType pointAndNormal;
                pointAndNormalDict
                {{
                    point (0 0 1.5);
                    normal (0 0 1);
                }}
                interpolate true;
            }}
        );
    }}
"""

        lx, ly, lz = leak_loc

        streamlines_str = f"""
    streamlines_leak
    {{
        type streamlines;
        libs ("libfieldFunctionObjects.so");
        writeControl writeTime;
        setFormat vtk;
        U U;
        direction forward;
        lifeTime 100000;
        nTrackingSteps 100000;
        trackLength 10000.0;
        fields ( p U {primary_gas} );
        seedSampleSet
        {{
            type sphereRandom;
            centre ({lx} {ly} {lz});
            radius 0.15;
            nPoints 50;
        }}
    }}
    
    streamlines_hvac
    {{
        type streamlines;
        libs ("libfieldFunctionObjects.so");
        writeControl writeTime;
        setFormat vtk;
        U U;
        direction forward;
        lifeTime 100000;
        nTrackingSteps 100000;
        trackLength 10000.0;
        fields ( p U {primary_gas} );
        seedSampleSet
        {{
            type boundaryRandom;
            patches ( {inlet_name} );
            nPoints 50;
        }}
    }}
"""

        functions_str = f"""
    outletFlow
    {{
        type surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        writeControl timeStep; writeInterval 1;
        log true; writeFields true;
        surfaceFormat vtk; regionType patches; patches ( {outlet_name} );
        operation sum;
        fields ( phi );
    }}
"""

        for spec in species_list:
            functions_str += f"""
    outlet_{spec}
    {{
        type surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        writeControl timeStep; writeInterval 1; log true; writeFields false;
        regionType patches; patches ( {outlet_name} ); operation areaAverage;
        fields ( {spec} );
    }}
"""

        probes_str = ""
        if probe_locations:
            pts_str = " ".join([f"({p[0]} {p[1]} {p[2]})" for p in probe_locations])
            fields_str = " ".join(species_list)

            probes_str = f"""
    gasSensors
    {{
        type probes;
        libs ("libsampling.so");
        writeControl timeStep; writeInterval 1;
        fields ( {fields_str} ); probeLocations ( {pts_str} );
    }}
"""

        max_dt = 0.05
        
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object controlDict; }}

application foamRun;
solver multicomponentFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {sim_duration};
deltaT {delta_t};
writeControl runTime;
writeInterval 5;
purgeWrite 0;
writeFormat ascii;
writePrecision 6;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable true;
adjustTimeStep yes;
maxCo {max_courant};
maxAlphaCo {max_courant};
maxDeltaT {max_dt};
timeStepChangeRatio 1.1;

functions
{{
{functions_str}
{probes_str}
{field_avg_str}
{surfaces_str}
{streamlines_str}
}}
"""

        with open(os.path.join(case_path, "system", "controlDict"), "w") as f:
            f.write(content)

    def _edit_thermo_multi(self, case_path, species_list):
        path = os.path.join(case_path, "constant", "thermophysicalProperties")
        spec_list_str = "\n    ".join(species_list)
        coeffs_blocks = ""

        for spec in species_list:
            c = self._get_coeffs(spec)
            coeffs_blocks += f"    {spec} {{ specie {{ molWeight {c['molWeight']}; }} thermodynamics {{ Tlow 200; Thigh 6000; Tcommon 1000; highCpCoeffs {c['high']}; lowCpCoeffs {c['low']}; }} transport {{ mu 1.8e-05; Pr 0.7; }} }}\n"

        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object thermophysicalProperties; }}

thermoType {{ type heRhoThermo; mixture multicomponentMixture; transport const; thermo janaf; equationOfState perfectGas; specie specie; energy sensibleEnthalpy; }}

defaultSpecie Air;
species ( {spec_list_str} );

{coeffs_blocks}
"""
        with open(path, "w") as f:
            f.write(content)

    def _setup_species_fields_multi(self, case_path, species_list):
        base_header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/\n"""

        for spec in species_list:
            is_air = (spec == "Air")
            val_internal = "1" if is_air else "0"
            val_inlet = "1" if is_air else "0"
            val_outlet = "1" if is_air else "0"

            content = base_header + f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object {spec}; }}

dimensions [0 0 0 0 0 0 0];

internalField uniform {val_internal};

boundaryField
{{
    ".*" {{ type zeroGradient; }}
    {self.vent_inlet} {{ type fixedValue; value uniform {val_inlet}; }}
    {self.vent_outlet} {{ type inletOutlet; inletValue uniform {val_outlet}; value uniform {val_outlet}; }}
}}
"""
            new_file = os.path.join(case_path, "0", spec)
            with open(new_file, "w") as f:
                f.write(content)