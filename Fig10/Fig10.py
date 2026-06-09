# freq_values根据Fig9的轨迹大概估算一个频率
# residence time 和 intensity 根据csv数据计算，但是为了一定的效果csv输出过程中使用了153-166经度范围

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
from matplotlib import rcParams

# ======================
# 样式
# ======================
# mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = [
    'Helvetica',
    'Arial',
    'Nimbus Sans',
    'Liberation Sans',
    'DejaVu Sans'
]
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
# =========================================================
# 1. Basic settings
# =========================================================

pmin_min_threshold = 1010  # hPa, to exclude weak/suspicious cases in intensity boxplot

# Tasman Sea region
LAT_MIN, LAT_MAX = -45, -30
LON_MIN, LON_MAX = 153, 170

# Input raw track files
files = {
    "Climatology": "Tasman_genesis.csv",
    "EP El Niño": "Tasman_genesis_EPElNino.csv",
    "CP El Niño": "Tasman_genesis_CPElNino.csv",
    "La Niña": "Tasman_genesis_LaNina.csv"
}

# Sample days for each category
sample_days = {
    "Climatology": 121 * 44,
    "EP El Niño": 163,
    "CP El Niño": 204,
    "La Niña": 289
}

labels = list(files.keys())

# Optional x-axis label colors
label_colors = ['black', 'black', 'black', 'black']


# =========================================================
# 2. Helper functions
# =========================================================

def in_tasman_region(df):
    """
    Return a boolean mask for whether each track point is inside the Tasman Sea box.
    """
    return (
        (df["lat"] >= LAT_MIN) & (df["lat"] <= LAT_MAX) &
        (df["lon"] >= LON_MIN) & (df["lon"] <= LON_MAX)
    )


def build_event_stats(csv_file):
    """
    Read raw cyclone track csv and calculate event-level statistics
    for each unique cyclone label.

    Input raw columns should include:
        datetime, lon, lat, pmin, label

    Output event-level dataframe includes:
        label           : cyclone ID
        pmin_min        : minimum central pressure during its track
        residence_time  : hours spent inside the Tasman Sea region
    """
    df = pd.read_csv(csv_file)

    # Make sure required columns exist
    required_cols = {"datetime", "lon", "lat", "pmin", "label"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{csv_file} is missing columns: {missing}")

    # Track points inside the Tasman Sea region
    df["in_tasman"] = in_tasman_region(df)

    # Group by cyclone label
    grouped = df.groupby("label")

    stats = grouped.agg(
        pmin_min=("pmin", "min"),
        residence_time=("in_tasman", "sum")   # each row = 1 hour
    ).reset_index()

    return stats


# =========================================================
# 3. Read raw data and calculate statistics
# =========================================================

event_stats = {}

for key, path in files.items():
    stats_df = build_event_stats(path)
    event_stats[key] = stats_df

    # Genesis frequency = cyclone count / sample days
    # Since these files are Tasman_genesis*.csv, each unique label is one cyclone
    n_cyclones = stats_df["label"].nunique()

    print(f"{key}:")
    print(stats_df.head(), "\n")


# =========================================================
# 4. Plot
# =========================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)

# ---------------------------------------------------------
# (a) Genesis frequency bar plot
# ---------------------------------------------------------
ax = axes[0]

freq_values = [0.106, 0.0920, 0.0637, 0.118]

bars = ax.bar(
    labels, freq_values,
    color="lightblue", edgecolor="gray", width=0.4
)

ax.set_title("(a) Genesis Frequency (num./day)", fontsize=14)
ax.grid(alpha=0.3, linestyle="--")
ax.set_ylim(0, max(freq_values) * 1.1)

# Add value labels above bars
for bar, val in zip(bars, freq_values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max(freq_values) * 0.03,
        f"{val:.3f}",
        ha='center', va='bottom',
        fontsize=10, color='black', fontweight='bold'
    )

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels)
for tick, color in zip(ax.get_xticklabels(), label_colors):
    tick.set_color(color)


# ---------------------------------------------------------
# (b) Residence time boxplot
# ---------------------------------------------------------
ax = axes[1]

residence_data = [event_stats[k]["residence_time"].dropna() for k in labels]

box = ax.boxplot(
    residence_data,
    patch_artist=True,
    tick_labels=labels,
    showmeans=True,
    meanline=True,
    showfliers=False,
    widths=0.4,
    boxprops=dict(facecolor="lightblue", color="gray"),
    medianprops=dict(color="black", linewidth=1.2, linestyle="--"),
    meanprops=dict(color="red", linewidth=1.5, linestyle="-"),
    whiskerprops=dict(color="gray"),
    capprops=dict(color="gray")
)

ax.set_title("(b) Residence Time (hours)", fontsize=14)
ax.grid(alpha=0.3, linestyle="--")

# Add mean labels
yl = ax.get_ylim()
offset = (yl[1] - yl[0]) * 0.015
for j, y in enumerate(residence_data):
    mean_val = np.mean(y)
    ax.text(
        j + 1, mean_val - offset,
        f"{mean_val:.1f}",
        ha='center', va='top',
        fontsize=10, color='black', fontweight='bold'
    )

ax.set_xticks(range(1, len(labels) + 1))
ax.set_xticklabels(labels)
for tick, color in zip(ax.get_xticklabels(), label_colors):
    tick.set_color(color)


# ---------------------------------------------------------
# (c) Intensity boxplot (minimum central pressure)
# Exclude weak/suspicious cases with pmin_min > 1010 hPa
# ---------------------------------------------------------
ax = axes[2]

intensity_data = [
    event_stats[k].loc[event_stats[k]["pmin_min"] <= pmin_min_threshold, "pmin_min"].dropna()
    for k in labels
]

box = ax.boxplot(
    intensity_data,
    patch_artist=True,
    tick_labels=labels,
    showmeans=True,
    meanline=True,
    showfliers=False,
    widths=0.4,
    boxprops=dict(facecolor="lightblue", color="gray"),
    medianprops=dict(color="black", linewidth=1.2, linestyle="--"),
    meanprops=dict(color="red", linewidth=1.5, linestyle="-"),
    whiskerprops=dict(color="gray"),
    capprops=dict(color="gray")
)

ax.set_title("(c) Intensity: Minimum Central SLP (hPa)", fontsize=14)
ax.grid(alpha=0.3, linestyle="--")

# Add mean labels
yl = ax.get_ylim()
offset = (yl[1] - yl[0]) * 0.015
for j, y in enumerate(intensity_data):
    mean_val = np.mean(y)
    ax.text(
        j + 1, mean_val - offset,
        f"{mean_val:.1f}",
        ha='center', va='top',
        fontsize=10, color='black', fontweight='bold'
    )

ax.set_xticks(range(1, len(labels) + 1))
ax.set_xticklabels(labels)
for tick, color in zip(ax.get_xticklabels(), label_colors):
    tick.set_color(color)


# =========================================================
# 5. Final layout and save
# =========================================================
plt.tight_layout(rect=[0, 0, 1, 0.95])

plt.savefig("Fig10.pdf", bbox_inches="tight")
# plt.savefig("cyclone_stats_tasman.pdf", bbox_inches="tight")

plt.show()