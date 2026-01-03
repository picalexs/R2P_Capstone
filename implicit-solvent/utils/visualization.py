#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Some helper functions to inspect your circuit layout to make sure it seems sane. You can visualize the selected layout, 
# then go inspect the hardware map online at https://quantum.cloud.ibm.com/ to make sure there are no obvious hardware errors sources on the layout
import seaborn as sns
from utils.gate_map import plot_gate_map
from qiskit.converters import circuit_to_dag

def used_qubits(qc, include={"unitary"}):
    """
    include: {"unitary","measure","reset","barrier","delay"}  (pick any)
    """
    wanted = set(include)
    dag = circuit_to_dag(qc)

    # collect qubits touched by selected ops
    used_qubits = set()
    for node in dag.topological_op_nodes():
        k = node.name
        is_unitary = k not in {"measure", "reset", "barrier", "delay"}
        if ("unitary" in wanted and is_unitary) or (k in wanted):
            used_qubits.update(node.qargs)

    idxs = set()
    for qb in used_qubits:
        loc = qc.find_bit(qb)
        idx = getattr(loc, "index", None) 
        if idx is None:
            idx = loc[1]
        idxs.add(idx)
    return sorted(idxs)

def color_batches(batches):
    """
    Given a list of sets of nodes (batches), assign a unique color to each node
    according to its batch.

    Returns:
        list indexed by node, value is a color name or hex code.
    """
    flat_palette = sns.color_palette("tab20", len(batches)).as_hex()
    color_map = {}
    for color_idx, batch in enumerate(batches):
        for node in batch:
            color_map[node] = flat_palette[color_idx]
    return [color_map[key] for key in sorted(color_map.keys())]

def plot_data(data, baseline_1=0, baseline_2=0, name=None, save=False):
    x_vals, y_vals, job_ids = zip(*data)
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot line + markers
    ax.plot(x_vals, y_vals, color='navy', linewidth=2, marker='o', markersize=5, label="Energy trajectory")
    ax.axhline(baseline_1, color='red', linestyle='--', linewidth=1.5, label="Hartree Fock energy")
    ax.axhline(baseline_2, color='blue', linestyle='--', linewidth=1.5, label="CASCI-derived energy")

    # Force plain formatting
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='plain', axis='y')

    # Annotate each point with its exact value
    for x, y, job_id in zip(x_vals, y_vals, job_ids):
        ax.annotate(f"{y:.8f}, ID: {job_id}", 
                    (x, y), 
                    textcoords="offset points", 
                    xytext=(0, 8),  # vertical offset
                    ha='center', fontsize=8, rotation=25, color='navy')

    # Annotate the Hartree Fock Classical Reference line
    for x,y in zip([0.5], [baseline_1]):
        ax.annotate(f"{y:.5f}", 
                (x, y), 
                textcoords="offset points", 
                xytext=(0, 8),  # vertical offset
                ha='center', fontsize=8, rotation=25, color='red')
        
    # Annotate the CASCI Classical Reference line
    for x,y in zip([0.5], [baseline_2]):
        ax.annotate(f"{y:.5f}", 
                (x, y), 
                textcoords="offset points", 
                xytext=(0, 8),  # vertical offset
                ha='center', fontsize=8, rotation=25, color='blue')

    # Titles, labels, etc 
    ax.set_title(f"SQD/IEF-PCM(cc-pVDZ) - {name}\nEnergy Convergence", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("Update Iterations", fontsize=12)
    ax.set_ylabel("Total Energy (Hartrees)", fontsize=12)
    ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
    ax.legend(frameon=True, loc="best")
    plt.tight_layout()

    if save:
        plt.savefig(f"./results/{name}_energy_convergence.png")
    return fig, ax

if __name__ == "__main__":
    raise SystemExit("This file is a utility module, not meant to be run directly.")
