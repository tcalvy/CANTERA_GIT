import numpy as np
import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt
import yaml

##-----------------PLOT PARAM------------------

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams["font.family"] = "Times New Roman"

SMALL_SIZE = 15
MEDIUM_SIZE = 20
BIGGER_SIZE = 25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

##-------------------------------------------


with open("plot_config.yaml", "r") as f:
    config = yaml.safe_load(f)

# ----------------------------------------------------------
# CONFIG
# ----------------------------------------------------------
folder = config['path']['output_folder'] 
type = config['type']

XOH_figsizex = config['XOH_PLOT']['figsizex']
XOH_figsizey = config['XOH_PLOT']['figsizey']
XOH_xlimactiv = str(config['XOH_PLOT']['xlimactiv'])
XOH_xmin = config['XOH_PLOT']['xmin']
XOH_xmax = config['XOH_PLOT']['xmax']
XOH_ylimactiv = str(config['XOH_PLOT']['ylimactiv'])
XOH_ymin = config['XOH_PLOT']['ymin']
XOH_ymax = config['XOH_PLOT']['ymax']

MAX_figsizex = config['MAX_PLOT']['figsizex']
MAX_figsizey = config['MAX_PLOT']['figsizey']
MAX_xlimactiv = str(config['MAX_PLOT']['xlimactiv'])
MAX_xmin = config['MAX_PLOT']['xmin']
MAX_xmax = config['MAX_PLOT']['xmax']
MAX_ylimactiv = str(config['MAX_PLOT']['ylimactiv'])
MAX_ymin = config['MAX_PLOT']['ymin']
MAX_ymax = config['MAX_PLOT']['ymax']


output_folder = folder + "/PLOTS"



#%%------------------------PLOT FOR DEBIT CONFIGURATION------------------------------------------------------
if type == 'debit' :
      # Regex to extract mdot_fuel and mdot_oxy from filename
      pattern = re.compile(r'diffusion_flame_H2AIR_DEBIT_(\d+)p(\d+)_(\d+)p(\d+)\.csv')
      ratios = []
      max_locations = []

      plt.figure(figsize=(XOH_figsizex,XOH_figsizey))
      # ----------------------------------------------------------
      # PROCESS FILES
      # ----------------------------------------------------------
      for filepath in glob.glob(os.path.join(folder, "diffusion_flame_H2AIR_DEBIT_*.csv")):
            filename = os.path.basename(filepath)
            match = pattern.match(filename)
            if not match:
                  continue



            # Convert '0p00127' to 0.00127
            mdot_fuel = float(f"{match.group(1)}.{match.group(2)}")
            mdot_oxy  = float(f"{match.group(3)}.{match.group(4)}")
            
            ratio = mdot_fuel / mdot_oxy

            # Load into DataFrame
            df = pd.read_csv(filepath, comment='#')
            
            # Extract numpy arrays
            grid = df["grid"]
            X_OH = df["X_OH"]

            plt.plot(grid, X_OH, label='$\dot{m}_{fuel} = $'+ str(mdot_fuel))

            # Compute index of maximum
            idx_max = np.argmax(X_OH)
            grid_at_max = grid[idx_max]

            ratios.append(ratio)
            max_locations.append(grid_at_max)

      # Convert to arrays & sort by mass-flow-ratio
      ratios = np.array(ratios)
      max_locations = np.array(max_locations)
      order = np.argsort(ratios)

      ratios = ratios[order]
      max_locations = max_locations[order]


      if not os.path.exists(output_folder):
                  os.makedirs(output_folder)


      plt.xlabel("Grid location ($m$)")
      if XOH_xlimactiv == 'yes' :
            plt.xlim([XOH_xmin,XOH_xmax])
      if XOH_ylimactiv == 'yes' :
            plt.xlim([XOH_ymin,XOH_ymax])
      plt.ylabel("$X_{OH}$")
      plt.legend()
      plt.grid(True)
      plt.tight_layout()
      plt.savefig(os.path.join(output_folder, "XCO_grid_debit.png"), dpi=300)
      plt.show()

      # ----------------------------------------------------------
      # PLOT RESULT
      # ----------------------------------------------------------

      plt.figure(figsize=(MAX_figsizex,MAX_figsizey))
      plt.plot(ratios, max_locations, marker="+")
      if MAX_xlimactiv == 'yes' :
            plt.xlim([MAX_xmin,MAX_xmax])
      if MAX_ylimactiv == 'yes' :
            plt.xlim([MAX_ymin,MAX_ymax])
      plt.xlabel("Mass-flow-rate ratio  mdot_fuel / mdot_oxy")
      plt.ylabel("Grid location of max(X_OH)")
      plt.savefig(os.path.join(output_folder, "max_XCO_ratio_debit.png"), dpi=300)
      plt.grid(True)
      plt.tight_layout()
      plt.show()




#%%------------------------PLOT FOR TEMPERATURE CONFIGURATION--------------------------------------------
if type == 'temperature' :
      # Regex to extract mdot_fuel and mdot_oxy from filename
      pattern = re.compile(r'diffusion_flame_H2AIR_TEMP_(\d+)p(\d+)_(\d+)p(\d+)\.csv')
      ratios = []
      max_locations = []


      plt.figure(figsize=(XOH_figsizex,XOH_figsizey))
      # ----------------------------------------------------------
      # PROCESS FILES
      # ----------------------------------------------------------
      for filepath in glob.glob(os.path.join(folder, "diffusion_flame_H2AIR_TEMP_*.csv")):
            filename = os.path.basename(filepath)
            match = pattern.match(filename)
            if not match:
                  continue
                  



            # Convert '0p00127' to 0.00127
            t_fuel = float(f"{match.group(1)}.{match.group(2)}")
            t_oxi  = float(f"{match.group(3)}.{match.group(4)}")
            ratio = t_fuel / t_oxi

            # Load into DataFrame
            df = pd.read_csv(filepath, comment='#')
            
            # Extract numpy arrays
            grid = df["grid"]
            X_OH = df["X_OH"]

            plt.plot(grid, X_OH, label='$T_{oxidizer} = $'+ str(t_oxi))

            # Compute index of maximum
            idx_max = np.argmax(X_OH)
            grid_at_max = grid[idx_max]

            ratios.append(ratio)
            max_locations.append(grid_at_max)

      # Convert to arrays & sort by mass-flow-ratio
      ratios = np.array(ratios)
      max_locations = np.array(max_locations)
      order = np.argsort(ratios)

      ratios = ratios[order]
      max_locations = max_locations[order]


      if not os.path.exists(output_folder):
                  os.makedirs(output_folder)


      plt.xlabel("Grid location ($m$)")
      if XOH_xlimactiv == 'yes' :
            plt.xlim([XOH_xmin,XOH_xmax])
      if XOH_ylimactiv == 'yes' :
            plt.xlim([XOH_ymin,XOH_ymax])
      plt.ylabel("$X_{OH}$")
      plt.legend()
      plt.grid(True)
      plt.tight_layout()
      plt.savefig(os.path.join(output_folder, "XCO_grid_temp.png"), dpi=300)
      plt.show()

      # ----------------------------------------------------------
      # PLOT RESULT
      # ----------------------------------------------------------

      plt.figure(figsize=(MAX_figsizex,MAX_figsizey))
      plt.plot(ratios, max_locations, marker="+")
      if MAX_xlimactiv == 'yes' :
            plt.xlim([MAX_xmin,MAX_xmax])
      if MAX_ylimactiv == 'yes' :
            plt.xlim([MAX_ymin,MAX_ymax])
      plt.xlabel("Mass-flow-rate ratio  mdot_fuel / mdot_oxy")
      plt.ylabel("Grid location of max(X_OH)")
      plt.savefig(os.path.join(output_folder, "max_XCO_ratio_temperature.png"), dpi=300)
      plt.grid(True)
      plt.tight_layout()
      plt.show()

