from pathlib import Path
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import yaml
import os

with open("cff_config.yaml", "r") as f:
    config = yaml.safe_load(f)

print(config)

comp_f = config['species']['composition_fuel']
comp_o = config['species']['composition_oxidizer']

p = config['properties']['pressure']
tin_f = config['properties']['fuel_temperature']
tin_o = config['properties']['oxidizer_temperature']

width = config['geometry']['width']  # Distance between inlets 

loglevel = 1  # amount of diagnostic output (0 to 5)

min_fuel = config['debits']['min_fuel']
max_fuel = config['debits']['max_fuel']
nb_pts_fuel = config['debits']['nb_points_fuel']

min_ox = config['debits']['min_oxidizer']
max_ox = config['debits']['max_oxidizer']
nb_pts_ox = config['debits']['nb_points_oxidizer']

folder_name = config['path']['output_folder']

# Create the gas object used to evaluate all thermodynamic, kinetic, and
# transport properties.
gas = ct.Solution('gri30.yaml')
gas.TP = gas.T, p

list_mdot_o = np.linspace(min_ox, max_ox, nb_pts_ox)    # kg/m^2/s
list_mdot_f = np.linspace(min_fuel, max_fuel, nb_pts_fuel)  # kg/m^2/s

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

        csv_path = (f"{folder_name}/diffusion_flame_H2AIR_{mdot_f:.5f}_{mdot_o:.5f}"
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
                
