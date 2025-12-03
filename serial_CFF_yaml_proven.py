from pathlib import Path
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import yaml
import os
import shutil



with open("CRECK_MECANISM/cff_config.yaml", "r") as f:
    config = yaml.safe_load(f)

print(config)



folder_name = config['path']['output_folder']

mechanism_path = config['mechanism']['path']

comp_f = config['species']['composition_fuel']
#comp_o = config['species']['composition_oxidizer']

width = config['geometry']['width']  # Distance between inlets 
fuel_inlet_surface = config['geometry']['fuel_inlet_surface']  # Distance between inlets 
oxi_inlet_surface = config['geometry']['oxi_inlet_surface']

type = config['type']

#%%BURNED GASES COMPUTATION

burned_gas_fuel = config['before_combustion']['fuel']
#phi_bb = config['before_combustion']['phi']
T_bb = config['before_combustion']['T']
P_bb = config['before_combustion']['P']

def compute_burned_gas(fuel, debit_m_fuel, debit_m_ox, T, P) :

    burned_gas = ct.Solution(mechanism_path)

    if fuel =='CH4' :
        MW_CH4 = 16.04e-3
        MW_AIR = 28.965e-3

        mdch4 = debit_m_fuel
        mdair = debit_m_ox
        A_stoechio = 2.0  # Moles d'O2 pour 1 mole de CH4 
        a = 3.76  # Rapport N2/O2 dans l'air
        phi_bb = (mdch4/mdair) / (1*MW_CH4/(A_stoechio*(1+3.76)*MW_AIR))
        mol_CH4_bb = 1.0
        mol_O2_bb = A_stoechio / phi_bb 
        mol_N2_bb = mol_O2_bb * a

        total_moles_bb = mol_CH4_bb + mol_O2_bb + mol_N2_bb

        X_CH4 = mol_CH4_bb / total_moles_bb
        X_O2 = mol_O2_bb / total_moles_bb
        X_N2 = mol_N2_bb / total_moles_bb

        burned_gas.TPX = T, P, {'CH4': mol_CH4_bb, 'O2': mol_O2_bb, 'N2': mol_N2_bb}
        compo={'CH4': mol_CH4_bb, 'O2': mol_O2_bb, 'N2': mol_N2_bb}
        burned_gas.equilibrate('HP')
        return burned_gas, compo, burned_gas.T, phi_bb
        


    else :
        print('ERROR : Fuel for burned gases not recognized')
        exit

#%%PLOT THE BURNED GAS COMPOSITION

def plot_burned_gas_compo(gas, plot_filename) :
    species_names = gas.species_names
    mole_fractions = gas.X

    # Combiner, filtrer les traces (X < 1e-6) et trier
    data = zip(species_names, mole_fractions)
    data = [(name, X) for name, X in data if X > 1e-6] 

    # Trier par fraction molaire décroissante
    sorted_data = sorted(data, key=lambda item: item[1], reverse=True)

    plot_species = [item[0] for item in sorted_data]
    plot_fractions = [item[1] for item in sorted_data]

    # Limiter aux 15 espèces principales pour la clarté
    if len(plot_species) > 15:
        plot_species = plot_species[:15]
        plot_fractions = plot_fractions[:15]

    # --- 5. Tracé des Résultats ---
    plt.figure(figsize=(12, 6))
    plt.bar(plot_species, plot_fractions, color='skyblue')

    plt.xlabel('Specie', fontsize=12)
    plt.ylabel('Molar Fraction', fontsize=12)
    # Affichage de la Température d'équilibre sur le titre
    T_ad = gas.T
    plt.title(f'Composition of burned gases - $\\text{{{burned_gas_fuel}}}/Air$ ($\\phi={phi_bb:.2f}$, $T_{{ad}} = {T_ad:.1f}\text{{ K}}$)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Sauvegarde de la figure
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = f"{folder_name}/{plot_filename}"
    plt.savefig(plot_filename)
    #print(f"Burned gas composition plot saved in  {plot_filename}")



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
    gas = ct.Solution(mechanism_path)
    gas.TP = gas.T, p

    list_mdot_o = np.linspace(min_d_ox, max_d_ox, nb_pts_d_ox)    # kg/m^2/s
    list_mdot_f = np.linspace(min_d_fuel, max_d_fuel, nb_pts_d_fuel)  # kg/m^2/s

    for mdot_o in list_mdot_o :
        for mdot_f in list_mdot_f : 

            

            burned_gas, comp_o, T_adiab, phi_bb = compute_burned_gas(burned_gas_fuel, mdot_f, mdot_o, T_bb, P_bb)
            plot_burned_gas_compo(burned_gas, f"burnt_gas_composition_{mdot_f:.5f}_{mdot_o:.5f}_{phi_bb:.2f}"
            .replace(".", "p")+ ".png")

            if tin_o == 'adiabatic' :
                tin_o = T_adiab


            # Create an object representing the counterflow flame configuration,
            # which consists of a fuel inlet on the left, the flow in the middle,
            # and the oxidizer inlet on the right.
            f = ct.CounterflowDiffusionFlame(gas, width=width)

            # Set the state of the two inlets
            f.fuel_inlet.mdot = mdot_f/fuel_inlet_surface
            f.fuel_inlet.X = comp_f
            f.fuel_inlet.T = tin_f

            f.oxidizer_inlet.mdot = mdot_o/oxi_inlet_surface
            #f.oxidizer_inlet.X = comp_o
            #f.oxidizer_inlet.T = tin_o
            f.oxidizer_inlet.X = burned_gas.X
            f.oxidizer_inlet.T = burned_gas.T
            


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
    gas = ct.Solution(mechanism_path)
 

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