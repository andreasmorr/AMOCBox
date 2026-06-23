"""
plotting_paper.py  –  AMOCBox 4-panel paper figure.

Reads data exported by amoc3box_export_paper_data.jl:
  data/paper/basin_{co2}ppm.csv
  data/paper/trajectories_{co2}ppm.csv
  data/paper/attractors_{co2}ppm.csv

Output: plots/amocbox_paper.png

Run from the AMOCBox submodule root or the project root:
    python scripts/plotting_paper.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup — find amoc_plot_style in the umbrella repo
# ---------------------------------------------------------------------------

SCRIPT_DIR  = Path(__file__).resolve().parent          # AMOCBox/scripts/
AMOC_BOX    = SCRIPT_DIR.parent                        # AMOCBox/
UMBRELLA    = AMOC_BOX.parent                          # AMOCResilience/
DATA_DIR    = AMOC_BOX / "data" / "paper"
PLOTS_DIR   = AMOC_BOX / "plots"

sys.path.insert(0, str(UMBRELLA))
from amoc_plot_style import (
    COL_ON, COL_OFF, COL_EDGE,
    BASIN_ON_FILL, BASIN_OFF_FILL, TRAJ_COLORS,
    LW_TRAJ, LW_TRAJ_PHASE, LW_ATTRACTOR, LS_ATTRACTOR, ALPHA_ATTRACTOR,
    LW_EQUIL, LS_EQUIL, ALPHA_EQUIL, IC_EDGE_COLOR, IC_EDGE_LW,
    make_paper_figure, add_panel_label, savefig_pdf,
)



# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

SCENARIOS = [
    ("280ppm",  "280 ppm"),
    ("448ppm", "448 ppm"),
]


def load_scenario(scenario: str) -> dict | None:
    """Load all CSV data for one scenario. Returns None if files are missing."""
    basin_path = DATA_DIR / f"basin_{scenario}.csv"
    traj_path  = DATA_DIR / f"trajectories_{scenario}.csv"
    attr_path  = DATA_DIR / f"attractors_{scenario}.csv"

    missing = [p for p in (basin_path, traj_path, attr_path) if not p.exists()]
    if missing:
        print(f"[{scenario}] Missing data files:")
        for p in missing:
            print(f"  {p}")
        return None

    return {
        "basin":      pd.read_csv(basin_path),
        "trajs":      pd.read_csv(traj_path),
        "attractors": pd.read_csv(attr_path),
    }


# ---------------------------------------------------------------------------
# Main plotting
# ---------------------------------------------------------------------------

def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load data for both scenarios
    data = {}
    for scenario, _ in SCENARIOS:
        d = load_scenario(scenario)
        if d is None:
            print(
                f"\nRun amoc3box_export_paper_data.jl first to generate data "
                f"in {DATA_DIR}."
            )
            sys.exit(1)
        data[scenario] = d

    # Build figure
    fig, axes_top, axes_bottom = make_paper_figure()
    for ax in axes_bottom:
        ax.set_aspect("auto")
    axes_bottom[1].sharex(axes_bottom[0])

    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for col, (scenario, title) in enumerate(SCENARIOS):
        d         = data[scenario]
        df_basin  = d["basin"]
        df_trajs  = d["trajs"]
        df_attr   = d["attractors"]

        ax_top    = axes_top[col]
        ax_bot    = axes_bottom[col]

        # ── Attractor equilibrium q values ───────────────────────────────────
        row_on  = df_attr[df_attr["state"] == "on"].iloc[0]
        row_off = df_attr[df_attr["state"] == "off"].iloc[0]
        q_on    = float(row_on["q"])
        q_off   = float(row_off["q"])

        # ── TOP panel: q(t) for 4 trajectories ───────────────────────────────
        traj_ids = sorted(df_trajs["traj_id"].unique())
        legend_handles = []
        legend_labels  = []
        for k, tid in enumerate(traj_ids[:4]):
            color  = TRAJ_COLORS[k % len(TRAJ_COLORS)]
            sub    = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            sub    = sub[sub["time"] <= 3000]
            t_arr  = sub["time"].values
            q_arr  = sub["q"].values
            lbl    = f"IC {k + 1}"
            (line,) = ax_top.plot(t_arr, q_arr, color=color, lw=LW_TRAJ, label=lbl)
            legend_handles.append(line)
            legend_labels.append(lbl)

        # Horizontal lines for equilibria
        ax_top.axhline(q_on,  color=COL_ON,  lw=LW_EQUIL, ls=LS_EQUIL, alpha=ALPHA_EQUIL)
        ax_top.axhline(q_off, color=COL_OFF, lw=LW_EQUIL, ls=LS_EQUIL, alpha=ALPHA_EQUIL)

        ax_top.set_title(title, fontsize=9)
        if col == 0:
            ax_top.set_ylabel("AMOC strength (Sv)")
        else:
            ax_top.tick_params(labelleft=False)
        ax_top.set_xlabel("Time (model years)")
        add_panel_label(ax_top, panel_labels[col])


        # ── BOTTOM panel: phase space ─────────────────────────────────────────
        # Basin fill
        sn_vals  = np.sort(df_basin["S_N"].unique())
        st_vals  = np.sort(df_basin["S_T"].unique())
        SN, ST   = np.meshgrid(sn_vals, st_vals, indexing="ij")
        basin_mat = np.full(SN.shape, np.nan)
        for _, row in df_basin.iterrows():
            ix = np.searchsorted(sn_vals, row["S_N"])
            iy = np.searchsorted(st_vals, row["S_T"])
            if 0 <= ix < len(sn_vals) and 0 <= iy < len(st_vals):
                basin_mat[ix, iy] = row["basin_label"]

        # Determine which label corresponds to the on-attractor
        on_label_val  = None
        off_label_val = None
        # Find which basin label is closest to the on-attractor position
        sn_on  = float(row_on["S_N"])
        st_on  = float(row_on["S_T"])
        sn_off = float(row_off["S_N"])
        st_off = float(row_off["S_T"])

        # Nearest grid cell to on-attractor
        ix_on = np.argmin(np.abs(sn_vals - sn_on))
        iy_on = np.argmin(np.abs(st_vals - st_on))
        on_label_val  = basin_mat[ix_on, iy_on]
        off_label_val = 1.0 if on_label_val == 2.0 else 2.0

        # Custom colormap: on=pale blue, off=pale orange
        cmap_basin = matplotlib.colors.ListedColormap(
            [BASIN_ON_FILL, BASIN_OFF_FILL]
            if on_label_val == 1.0 else
            [BASIN_OFF_FILL, BASIN_ON_FILL]
        )
        bounds = [0.5, 1.5, 2.5]
        norm   = matplotlib.colors.BoundaryNorm(bounds, cmap_basin.N)
        ax_bot.pcolormesh(SN, ST, basin_mat, cmap=cmap_basin, norm=norm,
                          shading="nearest", zorder=0)

        # Trajectories
        for k, tid in enumerate(traj_ids[:4]):
            color = TRAJ_COLORS[k % len(TRAJ_COLORS)]
            sub   = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            ax_bot.plot(sub["S_N"].values, sub["S_T"].values,
                        color=color, alpha=0.8, lw=LW_TRAJ_PHASE)
            ax_bot.scatter(sub["S_N"].values[0], sub["S_T"].values[0],
                           marker="o", s=25, color=color, zorder=3,
                           edgecolors=IC_EDGE_COLOR, linewidths=IC_EDGE_LW)

        # Attractor markers

        # Neighbourhood circles of radius 0.05 psu (= 0.005 model units)
        theta = np.linspace(0, 2 * np.pi, 300)
        r = 0.005
        for center, color in [((sn_on, st_on), COL_ON), ((sn_off, st_off), COL_OFF)]:
            ax_bot.plot(center[0] + r * np.cos(theta),
                        center[1] + r * np.sin(theta),
                        color=color, lw=LW_ATTRACTOR, ls=LS_ATTRACTOR, alpha=ALPHA_ATTRACTOR, zorder=4)

        # Convert raw model salinity units to psu: 0 → 35 psu, step 0.1 → 1 psu
        sal_fmt = mticker.FuncFormatter(lambda x, _: f"{35 + x * 10:.0f}")
        ax_bot.xaxis.set_major_formatter(sal_fmt)
        ax_bot.yaxis.set_major_formatter(sal_fmt)

        ax_bot.set_xlabel("North Atlantic salinity (psu)")
        if col == 0:
            ax_bot.set_ylabel("Tropical salinity (psu)")
        else:
            ax_bot.tick_params(labelleft=False)
        add_panel_label(ax_bot, panel_labels[col + 2])

    x_lo = min(ax.get_xlim()[0] for ax in axes_top)
    x_hi = max(ax.get_xlim()[1] for ax in axes_top)
    for ax in axes_top:
        ax.set_xlim(x_lo, x_hi)

    out_path = PLOTS_DIR / "amocbox_paper.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
