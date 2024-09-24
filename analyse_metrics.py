import os
import pandas as pd
import matplotlib.pyplot as plt

# Define the base directory where all experiment folders are located
base_dir = r"C:\Users\Utilisateur\Desktop\wetransfer_data_2024-07-26_1546\twin_exp\Step_1\Output"

# Get the directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Initialize lists to store data for all experiments
all_params = []
all_E_metrics = []
all_L_metrics = []

# Define experiment groups based on jtlag, lenght_scale, and NiterOut
exp_groups = {
    "jtlag": [1, 2, 3, 4, 5],
    "length_scale": [6, 1, 7, 8, 9],
    "NiterOut": [11, 12, 1, 13]
}

# Loop through each experiment folder
for exp_num in range(1, 14):  # Assuming there are 13 experiments
    exp_dir = os.path.join(base_dir, f"{exp_num:02d}_exp")
    
    # Read the rapport.csv file
    rapport_file = os.path.join(exp_dir, "rapport.csv")
    rapport_df = pd.read_csv(rapport_file)
    
    # Extract parameters and metrics from rapport.csv, convert to float for numeric operations
    params = rapport_df.set_index("Parameter")["Value"].apply(pd.to_numeric, errors='coerce').to_dict()
    all_params.append(params)
    
    # Read the E_metricsArray.csv file
    e_metrics_file = os.path.join(exp_dir, "E_metricsArray.csv")
    e_metrics_df = pd.read_csv(e_metrics_file, header=None)
    all_E_metrics.append(e_metrics_df.values.flatten())
    
    # Read the L_metricsArray.csv file
    l_metrics_file = os.path.join(exp_dir, "L_metricsArray.csv")
    l_metrics_df = pd.read_csv(l_metrics_file, header=None)
    all_L_metrics.append(l_metrics_df.values.flatten())

# Plot E_metrics and L_metrics according to jtlag
plt.figure()
for exp_num in exp_groups["jtlag"]:
    jtlag_value = round(float(all_params[exp_num - 1]["jtlag"]), 1)
    plt.plot(all_E_metrics[exp_num - 1], label=f"Exp {exp_num} (jtlag={jtlag_value})")
plt.xlabel("Time Step")
plt.ylabel("E_metrics")
plt.title("Evolution of E_metrics according to jtlag")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "E_metrics_jtlag.png"))
plt.show()

plt.figure()
for exp_num in exp_groups["jtlag"]:
    jtlag_value = round(float(all_params[exp_num - 1]["jtlag"]), 1)
    plt.plot(all_L_metrics[exp_num - 1], label=f"Exp {exp_num} (jtlag={jtlag_value})")
plt.xlabel("Time Step")
plt.ylabel("L_metrics (m)")
plt.title("Evolution of L_metrics according to jtlag")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "L_metrics_jtlag.png"))
plt.show()

# Plot E_metrics according to lenght_scale
plt.figure()
for exp_num in exp_groups["length_scale"]:
    lenght_scale_value = round(float(all_params[exp_num - 1]["lenfil"]), 1)
    plt.plot(all_E_metrics[exp_num - 1], label=f"Exp {exp_num} (lenght_scale={lenght_scale_value})")
plt.xlabel("Time Step")
plt.ylabel("E_metrics")
plt.title("Evolution of E_metrics according to length_scale")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "E_metrics_length_scale.png"))
plt.show()

# Plot L_metrics according to lenght_scale
plt.figure()
for exp_num in exp_groups["length_scale"]:
    lenght_scale_value = round(float(all_params[exp_num - 1]["lenfil"]), 1)
    plt.plot(all_L_metrics[exp_num - 1], label=f"Exp {exp_num} (length_scale={lenght_scale_value})")
plt.xlabel("Time Step")
plt.ylabel("L_metrics (m)")
plt.title("Evolution of L_metrics according to length_scale")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "L_metrics_length_scale.png"))
plt.show()

# Plot E_metrics according to NiterOut
plt.figure()
for exp_num in exp_groups["NiterOut"]:
    NiterOut_value = round(float(all_params[exp_num - 1]["NiterOut"]), 1)
    plt.plot(all_E_metrics[exp_num - 1], label=f"Exp {exp_num} (NiterOut={NiterOut_value})")
plt.xlabel("Time Step")
plt.ylabel("E_metrics")
plt.title("Evolution of E_metrics according to NiterOut")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "E_metrics_NiterOut.png"))
plt.show()

# Plot L_metrics according to NiterOut
plt.figure()
for exp_num in exp_groups["NiterOut"]:
    NiterOut_value = round(float(all_params[exp_num - 1]["NiterOut"]), 1)
    plt.plot(all_L_metrics[exp_num - 1], label=f"Exp {exp_num} (NiterOut={NiterOut_value})")
plt.xlabel("Time Step")
plt.ylabel("L_metrics (m)")
plt.title("Evolution of L_metrics according to NiterOut")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(script_dir, "L_metrics_NiterOut.png"))
plt.show()

