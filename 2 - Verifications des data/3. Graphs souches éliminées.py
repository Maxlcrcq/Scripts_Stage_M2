import matplotlib.pyplot as plt
import pandas as pd

# Read the Excel file
df = pd.read_excel(r"C:\Users\lecle\Desktop\Clean\2 - Verificiations Data\Rapport.xlsx")

# Define the parameter columns
param_cols = ['Longueur Totale', 'GC (%)', '# contigs (>= 1000 bp)']

# Define the Souche column
souce_col = 'Souche'

# Define the specific souche values to highlight in red
highlight_souches = ['A3', 'B3', 'B4', 'B5', 'G3', 'G5', 'II_E1', 'II_F4']

# Create a figure and axis object
fig, ax = plt.subplots(nrows=len(param_cols), ncols=1, figsize=(12, 8*len(param_cols)))

# Loop through each parameter column
for i, param in enumerate(param_cols):
    # Create a scatter plot for each parameter
    ax[i].scatter(df[souce_col], df[param], c=['red' if souce in highlight_souches else 'blue' for souce in df[souce_col]])
    ax[i].set_xlabel(souce_col)
    ax[i].set_ylabel(param)
    ax[i].set_title(param)
    ax[i].tick_params(axis='both', labelsize=8)
    ax[i].tick_params(axis='x', labelrotation=45)

# Show the plot
plt.tight_layout()
plt.show()