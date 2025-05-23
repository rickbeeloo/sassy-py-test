#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd

f = "benchmarks/results.csv"
df = pd.read_csv(f, delimiter=",")

profile = "dna"
k = 1

filtered = df[
    (df["profile"] == profile)
    & (df["alphabet"] == "Dna")
    & (df["k"] == k)
    & (df["match_fraction"] == 1)
].copy()

# print(filtered)

# Create a custom colormap from gray to #fcc007
colors = ["#fcc007", "black"]
custom_cmap = LinearSegmentedColormap.from_list("custom_gray_yellow", colors)

# Plot
plt.figure(figsize=(12, 7))

# Plot each query length with a gradient color and connect points with lines
for i, ql in enumerate(filtered["query_length"].unique()):
    subset = filtered[filtered["query_length"] == ql]
    color = custom_cmap(i / len(filtered["query_length"].unique()))

    # Plot Edlib times with dotted line
    plt.plot(
        subset["text_length"],
        subset["edlib_ms"],
        color=color,
        marker="^",
        markersize=10,
        linestyle=":",
        label=f"Edlib (query_length={ql})",
    )

    # Plot Sassy times with solid line
    plt.plot(
        subset["text_length"],
        subset["sassy_ms"],
        color=color,
        marker="o",
        markersize=10,
        linestyle="-",
        label=f"Sassy (query_length={ql})",
    )

# Log scale for x-axis with original labels
xticks = sorted(filtered["text_length"].unique())
plt.xscale("log")
plt.yscale("log")
plt.xticks(xticks, [str(x) for x in xticks])

plt.xlabel("Text Length (log scale)", fontsize=12)
plt.ylabel("Execution Time (ms, log scale)", fontsize=12)
plt.title(f"Edlib vs Sassy Execution Times (bounded, {profile}, k={k})", fontsize=14)
plt.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.tight_layout()
plt.savefig("results.png", dpi=300)
plt.show()