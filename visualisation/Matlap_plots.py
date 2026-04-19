# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:55:34 2025
@author: Samantha Crouch
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# ------------------------------------------------------------------
FS = 16                         # uniform font size
XTICK_START, XTICK_END = 0, 10.5  # x-axis range
# ------------------------------------------------------------------

os.chdir(r"YourFolder")

df = pd.read_csv('chi2_full.csv')
df_sorted = df.sort_values(by='Time')
variables = df_sorted['Variable'].unique()

n_cols = 3  # change if needed!
n_rows = int(np.ceil(len(variables) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), sharex=False)
axes = axes.flatten()

colors = plt.cm.viridis(np.linspace(0, 1, len(df_sorted['Time'].unique())))
legend_elements = []

for ax, var in zip(axes, variables):
    df_var = df_sorted[df_sorted['Variable'] == var]

    for idx, (time, exp_data, error) in enumerate(zip(df_var['Time'], df_var['ExpData'], df_var['Std'])):
        color = colors[idx]
        label = f'{time} h'
        if var == variables[0]:
            legend_elements.append(
                plt.Line2D([0], [0], color=color, marker='o', linestyle='None', label=label)
            )
        ax.errorbar(time, exp_data, yerr=error, fmt='o', color=color, capsize=5)

    ax.plot(df_var['Time'], df_var['SimulatedData'], color='red', linestyle='-', label='Fitting curve')
    if var == variables[0]:
        legend_elements.append(plt.Line2D([0], [0], color='red', label='Fitting curve'))

    chi2_sum = df_var['Chi2'].sum()
    ax.set_title(f"{var} | χ² {chi2_sum:.4f}", fontsize=FS)
    ax.set_xlabel('days', fontsize=FS)
    ax.set_ylabel('state', fontsize=FS)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.set_xlim(XTICK_START, XTICK_END)
    ax.set_xticks(np.arange(XTICK_START, int(np.floor(XTICK_END)) + 1, 1))  # 1-unit grid
    ax.tick_params(axis='both', labelsize=FS)
    ax.grid(True)

for ax in axes[len(variables):]:
    fig.delaxes(ax)

fig.subplots_adjust(bottom=0.20, hspace=0.4, wspace=0.3)

fig.legend(
    handles=legend_elements,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.01),
    ncol=4,
    fontsize=FS,
    markerscale=1.5,
    borderpad=0.5,
    frameon=False
)

fig.savefig('gplot1.pdf')
fig.savefig('gplot.png', dpi=600)
fig.savefig('gplot.jpeg', dpi=600)
fig.savefig('gplot.tiff', dpi=600)

plt.show()

#### End