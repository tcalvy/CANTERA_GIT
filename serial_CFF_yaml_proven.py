from pathlib import Path
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import yaml
import os
import shutil



with open("cff_config.yaml", "r") as f:
    config = yaml.safe_load(f)

print(config)



folder_name = config['path']['output_folder']

comp_f = config['species']['composition_fuel']
comp_o = config['species']['composition_oxidizer']

width = config['geometry']['width']  # Distance between inlets 

type = config['type']



#%%CONFIG DEBIT INCREMENTATION
if type == 'debit' : 
    p = config['d_properties']['pressure']
    tin_f = config['d_properties']['fuel_temperature']
    tin_o = config['d_properties']['oxidizer_temperature']

    loglevel = 1  # amount of diagnostic output (0 to 5)

    min_d_fuel = config['d_debits']['min_fuel']
    max_d_fuel = config['d_debits']['max_fuel']
    nb_pts_d_fuel = config['d_debits']['nb_points_fuel']

    min_d_ox = config['d_debits']['min_oxidizer']
    max_d_ox = config['d_debits']['max_oxidizer']
    nb_pts_d_ox = config['d_debits']['nb_points_oxidizer']

 

    # Create the gas object used to evaluate all thermodynamic, kinetic, and
    # transport properties.
    gas = ct.Solution('gri30.yaml')
    gas.TP = gas.T, p

    list_mdot_o = np.linspace(min_d_ox, max_d_ox, nb_pts_d_ox)    # kg/m^2/s
    list_mdot_f = np.linspace(min_d_fuel, max_d_fuel, nb_pts_d_fuel)  # kg/m^2/s

    for mdot_o in list_mdot_o :
        for mdot_f in list_mdot_f : 

            # Create an object representing the counterflow flame configuration,
            # which consists of a fuel inlet on the left, the flow in the middle,
            # and the oxidizer inlet on the right.
            f = ct.CounterflowDiffusionFlame(gas, width=width)

            # Set the state of the two inlets
            f.fuel_inlet.mdot = mdot_f
            f.fuel_inlet.X = comp_f
            f.fuel_inlet.T = tin_f

            f.oxidizer_inlet.mdot = mdot_o
            f.oxidizer_inlet.X = comp_o
            f.oxidizer_inlet.T = tin_o


            f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)

            f.solve(loglevel, auto=True)

            # --- Write important info into the first line of the CSV ---
            header = (
                f"mdot_f={mdot_f}, mdot_o={mdot_o}, "
                f"comp_f={comp_f}, comp_o={comp_o}, "
                f"tin_f={tin_f}, tin_o={tin_o}"
            )

            if not os.path.exists(folder_name):
                os.makedirs(folder_name)

            csv_path = (f"{folder_name}/diffusion_flame_DEBIT_{mdot_f:.5f}_{mdot_o:.5f}"
            .replace(".", "p")
            + ".csv"
    )
            # Step 1: let Cantera create the CSV
            f.save(csv_path, basis="mole", overwrite=True)

            # Step 2: prepend your header
            with open(csv_path, "r") as f_read:
                original = f_read.read()

            with open(csv_path, "w") as f_write:
                f_write.write("# " + header + "\n" + original)

            # Step 3: continue processing
            f.show_stats()
            no_rad = f.to_array()

#%%CONFIG TEMPERATURE INCREMENTATION
if type == 'temperature' : 
    p = config['t_properties']['pressure']
    d_fuel = config['t_properties']['debit_fuel']
    d_ox = config['t_properties']['debit_oxidizer']

    loglevel = 1  # amount of diagnostic output (0 to 5)

    min_t_fuel = config['t_temperatures']['min_fuel']
    max_t_fuel = config['t_temperatures']['max_fuel']
    nb_pts_t_fuel = config['t_temperatures']['nb_points_fuel']

    min_t_ox = config['t_temperatures']['min_oxidizer']
    max_t_ox = config['t_temperatures']['max_oxidizer']
    nb_pts_t_ox = config['t_temperatures']['nb_points_oxidizer']

 

    # Create the gas object used to evaluate all thermodynamic, kinetic, and
    # transport properties.
    gas = ct.Solution('gri30.yaml')
 

    list_t_o = np.linspace(min_t_ox, max_t_ox, nb_pts_t_ox)    # K
    list_t_f = np.linspace(min_t_fuel, max_t_fuel, nb_pts_t_fuel)  # K

    for temp_o in list_t_o :
        for temp_f in list_t_f : 

            gas.TP = gas.T, p

            # Create an object representing the counterflow flame configuration,
            # which consists of a fuel inlet on the left, the flow in the middle,
            # and the oxidizer inlet on the right.
            f = ct.CounterflowDiffusionFlame(gas, width=width)

            # Set the state of the two inlets
            f.fuel_inlet.mdot = d_fuel
            f.fuel_inlet.X = comp_f
            f.fuel_inlet.T = temp_f

            f.oxidizer_inlet.mdot = d_ox
            f.oxidizer_inlet.X = comp_o
            f.oxidizer_inlet.T = temp_o


            f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)

            f.solve(loglevel, auto=True)

            # --- Write important info into the first line of the CSV ---
            header = (
                f"mdot_f={d_fuel}, mdot_o={d_ox}, "
                f"comp_f={comp_f}, comp_o={comp_o}, "
                f"tin_f={temp_f}, tin_o={temp_o}"
            )

            if not os.path.exists(folder_name):
                os.makedirs(folder_name)

            csv_path = (f"{folder_name}/diffusion_flame_TEMP_{temp_f:04.4f}_{temp_o:04.4f}"
            .replace(".", "p")
            + ".csv"
    )
            # Step 1: let Cantera create the CSV
            f.save(csv_path, basis="mole", overwrite=True)

            # Step 2: prepend your header
            with open(csv_path, "r") as f_read:
                original = f_read.read()

            with open(csv_path, "w") as f_write:
                f_write.write("# " + header + "\n" + original)

            # Step 3: continue processing
            f.show_stats()
            no_rad = f.to_array()

shutil.copy("cff_config.yaml", f"{folder_name}/cff_config.yaml")