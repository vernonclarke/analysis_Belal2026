'''
functions
'''
import os
import sys
import re

from typing import Dict, List, Tuple

import h5py
import json
import matplotlib.pyplot as plt
from pathlib import Path
import xml.etree.ElementTree as et
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import remfile
import seaborn as sns
from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv
from pynwb import NWBHDF5IO
from tqdm import tqdm


# figure 5 functions
def setup_figure_style():
    """Setup matplotlib parameters to match paper style"""
    plt.rcParams.update({
        'font.size': 8,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.titlesize': 12,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.minor.width': 0.6,
    })

def get_session_id(asset_path: str) -> str:
    """Extract session ID from DANDI asset path."""
    if not asset_path:
        return ""
    bottom_level_path = asset_path.split("/")[1]
    session_id_with_ses_prefix = bottom_level_path.split("_")[1]
    session_id = session_id_with_ses_prefix.split("-")[1]
    return session_id

def get_figure_number(session_id: str):
    """Extract which figure this data corresponds to."""
    return session_id.split("+"+"+")[0]

def get_measurement(session_id: str) -> str:
    """Extract measurement type."""
    if not session_id:
        return ""
    return session_id.split("+"+"+")[1]

def get_cell_type(session_id: str) -> str:
    """Extract cell type."""
    if not session_id:
        return ""
    return session_id.split("+"+"+")[2]

def get_state(session_id: str) -> str:
    """Extract experimental state."""
    if not session_id:
        return ""
    return session_id.split("+"+"+")[3]

def is_f5_achfp(session_id: str) -> bool:
    """Check if data belongs to Figure 5 acetylcholine fluorescent protein experiments."""
    return get_figure_number(session_id) == "F5" and get_measurement(session_id) == "AChFP"

def get_condition_label(session_id: str) -> str:
    """Convert session metadata to condition label."""
    parts = session_id.split("+"+"+")
    if len(parts) < 4:
        return "unknown"
    
    state = parts[3]
    
    if state == "CTRL":
        return "UL control"
    elif state == "PD":
        return "6-OHDA"
    elif state == "OFF":
        return "LID off-state"
    else:
        return "unknown"

# CatalystNeuro original
# def process_acetylcholine_trial(
#     timestamps: np.ndarray,
#     fluorescence: np.ndarray,
#     stim_time: float = 10.0,
# ) -> Dict:
#     """Process a single acetylcholine biosensor trial."""
#     # Calculate stimulus-relative windows
#     baseline_window = (stim_time - 1.0, stim_time)
#     response_window = (stim_time, stim_time + 2.0)

#     # Find indices for time windows
#     time_interval = timestamps[1] - timestamps[0]
#     idx_baseline_start = np.argmin(np.abs(timestamps - baseline_window[0]))
#     idx_baseline_end = np.argmin(np.abs(timestamps - baseline_window[1]))
#     idx_response_start = np.argmin(np.abs(timestamps - response_window[0]))
#     idx_response_end = np.argmin(np.abs(timestamps - response_window[1]))

#     # Calculate F0 (baseline fluorescence)
#     F0 = float(np.mean(fluorescence[idx_baseline_start:idx_baseline_end]))

#     # Calculate dF/F0
#     dF_over_F0 = (fluorescence - F0) / F0 if F0 > 0 else np.zeros_like(fluorescence)

#     # Unnormalized AUC (optional diagnostic)
#     auc = float(np.sum(dF_over_F0[idx_response_start:idx_response_end]) * time_interval)

#     # Peak response
#     peak_response = float(np.max(dF_over_F0[idx_response_start:idx_response_end]))

#     return {
#         "F0": F0,
#         "dF_over_F0": dF_over_F0,
#         "auc": auc,
#         "peak_response": peak_response,
#         "timestamps": timestamps,
#         "fluorescence": fluorescence,
#     }

# added PMT_background default 0
def process_acetylcholine_trial(
    timestamps: np.ndarray,
    fluorescence: np.ndarray,
    stim_time: float = 10.0,
    PMT_background: float = 0.0,
) -> Dict:
    """Process a single acetylcholine biosensor trial with optional PMT background subtraction."""

    # Subtract PMT background first
    fluorescence = fluorescence - PMT_background

    # Calculate stimulus-relative windows
    baseline_window = (stim_time - 1.0, stim_time)
    response_window = (stim_time, stim_time + 2.0)

    # Find indices for time windows
    time_interval = timestamps[1] - timestamps[0]
    idx_baseline_start = np.argmin(np.abs(timestamps - baseline_window[0]))
    idx_baseline_end = np.argmin(np.abs(timestamps - baseline_window[1]))
    idx_response_start = np.argmin(np.abs(timestamps - response_window[0]))
    idx_response_end = np.argmin(np.abs(timestamps - response_window[1]))

    # Baseline fluorescence (F0)
    F0 = float(np.mean(fluorescence[idx_baseline_start:idx_baseline_end]))

    # dF/F0 trace
    dF_over_F0 = (fluorescence - F0) / F0 if F0 > 0 else np.zeros_like(fluorescence)
    
    # Unnormalized AUC in response window
    rs = dF_over_F0[idx_response_start:idx_response_end]
    ts = timestamps[idx_response_start:idx_response_end]
    auc = float(np.trapezoid(rs, ts))
    
    # # previous method
    # # Unnormalized AUC in response window
    # auc = float(np.sum(dF_over_F0[idx_response_start:idx_response_end]) * time_interval)

    # Peak response above baseline (ΔF = Fpeak - F0)
    F = float(np.max(fluorescence[idx_response_start:idx_response_end]) - F0)

    return {
        "F0": F0,
        "F": F,
        "dF_over_F0": dF_over_F0,
        "auc": auc,
        "timestamps": timestamps,
        "fluorescence": fluorescence,
        "PMT_background": PMT_background,
    }


def get_calibration_values(trials_data: List[Dict], PMT_background: float = 0.0) -> Tuple[float, float]:
    """Extract Fmax and Fmin using full-trace calibration averages with background subtraction."""
    ach_trials = [t for t in trials_data if t["treatment"] == "ACh_calibration"]
    ttx_trials = [t for t in trials_data if t["treatment"] == "TTX_calibration"]

    if ach_trials and ttx_trials:
        ach_vals = []
        for t in ach_trials:
            ach_vals.extend((t["fluorescence"] - PMT_background).tolist())
        ttx_vals = []
        for t in ttx_trials:
            ttx_vals.extend((t["fluorescence"] - PMT_background).tolist())
        if ach_vals and ttx_vals:
            return float(np.mean(ach_vals)), float(np.mean(ttx_vals))
    # 2-4 seconds not entire / if more than 1 take highest ACH max
    # Fallback constants (from original script)
    return 493.47, 212.39


def normalize_fluorescence(trial_data: Dict, Fmax: float, Fmin: float, stim_time: float = 10.0) -> Dict:
    """Normalize fluorescence; compute AUC via trapezoidal integration in [stim, stim+2s]."""
    
    FI = Fmax - Fmin
    dF_over_FI = (trial_data["fluorescence"] - trial_data["F0"]) / FI if FI > 0 else np.zeros_like(trial_data["fluorescence"])
    trial_data["dF_over_FI"] = dF_over_FI

    # Response window indices
    idx_start = np.argmin(np.abs(trial_data["timestamps"] - stim_time))
    idx_end = np.argmin(np.abs(trial_data["timestamps"] - (stim_time + 2.0)))
    rs = dF_over_FI[idx_start:idx_end]
    ts = trial_data["timestamps"][idx_start:idx_end]

    if len(rs) > 1:
        trial_data["auc_normalized"] = float(np.trapezoid(rs, ts))
    else:
        trial_data["auc_normalized"] = 0.0

    return trial_data

# Map R's quantile types (1–9) to NumPy methods
_RTYPE2METHOD = {
    1: "inverted_cdf", 2: "averaged_inverted_cdf", 3: "closest_observation",
    4: "interpolated_inverted_cdf", 5: "hazen", 6: "weibull",
    7: "linear", 8: "median_unbiased", 9: "normal_unbiased",
}

# def boxplot_rtype(
#     ax,
#     data,
#     rtype: int = 7,
#     whis=1.5,
#     tick_labels=None,
#     showfliers: bool = True,
#     manage_ticks: bool = True,
#     box_linewidth: float = 1.5,
#     median_overhang: float = 0.05,
#     showpoints: bool = False,
#     jitter_frac: float = 1.0,
#     point_diameter: float = 6.0,
#     point_alpha: float = 0.5,
#     point_color: str = "gray",
#     point_edgecolor: str = "black",
#     positions=None,
#     **bxp_kwargs,
# ):
#     if rtype not in _RTYPE2METHOD:
#         raise ValueError(f"rtype must be 1..9, got {rtype}")
#     method = _RTYPE2METHOD[rtype]

#     if hasattr(data, "__array__") and np.ndim(data) == 1:
#         datasets = [np.asarray(data, float)]
#     else:
#         datasets = [np.asarray(d, float) for d in data]
#     datasets = [d[~np.isnan(d)] for d in datasets if len(d) > 0]

#     stats = []
#     for d in datasets:
#         if d.size == 0:
#             continue
#         q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75], method=method)
#         iqr = q3 - q1
#         if isinstance(whis, str) and whis.lower() in {"minmax", "min-max"}:
#             whislo, whishi = float(np.min(d)), float(np.max(d))
#         elif isinstance(whis, (list, tuple, np.ndarray)) and len(whis) == 2:
#             low_p, high_p = whis
#             whislo = float(np.quantile(d, low_p/100.0, method=method))
#             whishi = float(np.quantile(d, high_p/100.0, method=method))
#         else:
#             w = float(whis)
#             low_b = q1 - w * iqr
#             high_b = q3 + w * iqr
#             d_lo = d[d >= low_b]
#             d_hi = d[d <= high_b]
#             whislo = float(d_lo.min() if d_lo.size else d.min())
#             whishi = float(d_hi.max() if d_hi.size else d.max())
#         fliers = d[(d < whislo) | (d > whishi)]
#         stats.append(dict(med=float(med), q1=float(q1), q3=float(q3),
#                           whislo=float(whislo), whishi=float(whishi),
#                           fliers=fliers))

#     defaults = dict(
#         patch_artist=True,
#         boxprops=dict(edgecolor="black", linewidth=box_linewidth),
#         whiskerprops=dict(color="black", linewidth=box_linewidth),
#         capprops=dict(color="black", linewidth=box_linewidth),
#         medianprops=dict(color="black", linewidth=3*box_linewidth),
#     )
#     defaults.update(bxp_kwargs)

#     if positions is None:
#         positions = list(range(1, len(stats) + 1))

#     artists = ax.bxp(stats, showfliers=showfliers, positions=positions, **defaults)

#     for m in artists["medians"]:
#         x0, x1 = m.get_xdata()
#         m.set_xdata([x0 - median_overhang, x1 + median_overhang])
#         m.set_solid_capstyle("round")
#     for ln in artists["whiskers"] + artists["caps"]:
#         ln.set_solid_capstyle("round")

#     if tick_labels is not None and manage_ticks:
#         ax.set_xticks(positions)
#         ax.set_xticklabels(tick_labels)

#     scatters = []
#     if showpoints:
#         widths = bxp_kwargs.get("widths", 0.6)
#         if isinstance(widths, (list, tuple, np.ndarray)):
#             widths_seq = list(widths)
#         else:
#             widths_seq = [float(widths)] * len(stats)
#         capwidths = bxp_kwargs.get("capwidths", [0.25 * w for w in widths_seq])

#         for xpos, d, w, cw in zip(positions, datasets, widths_seq, capwidths):
#             if d.size == 0:
#                 continue
#             half_box = 0.5 * float(w)
#             half_cap = 0.5 * float(cw)
#             half_span = min(half_box, half_cap) * jitter_frac
#             offsets = (np.random.rand(d.size) - 0.5) * (2.0 * half_span)
#             x = xpos + offsets
#             sc = ax.scatter(
#                 x, d,
#                 s=point_diameter**2,
#                 alpha=point_alpha,
#                 c=point_color,
#                 edgecolors=point_edgecolor,
#                 linewidths=0.5,
#                 zorder=3
#             )
#             scatters.append(sc)

#     return artists, scatters

def boxplot_rtype(
    ax,
    data,
    rtype: int = 7,
    whis=1.5,
    tick_labels=None,
    showfliers: bool = True,
    manage_ticks: bool = True,
    box_linewidth: float = 1.5,
    median_overhang: float = 0.05,
    showpoints: bool = False,
    jitter_frac: float = 1.0,
    jitter_seed: int = 42,
    point_diameter: float = 6.0,
    point_alpha: float = 0.5,
    point_color: str = "gray",
    point_edgecolor: str = "black",
    positions=None,
    paired: bool = False,
    subjects=None,
    **bxp_kwargs,
):
    if rtype not in _RTYPE2METHOD:
        raise ValueError(f"rtype must be 1..9, got {rtype}")
    method = _RTYPE2METHOD[rtype]

    if hasattr(data, "__array__") and np.ndim(data) == 1:
        datasets = [np.asarray(data, float)]
    else:
        datasets = [np.asarray(d, float) for d in data]
    datasets = [d[~np.isnan(d)] for d in datasets if len(d) > 0]

    if paired and subjects is None:
        raise ValueError("If paired=True, subjects must be provided")
    if paired:
        if len(subjects) != len(datasets):
            raise ValueError("subjects must align with data (one array per dataset)")

    stats = []
    for d in datasets:
        if d.size == 0:
            continue
        q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75], method=method)
        iqr = q3 - q1
        if isinstance(whis, str) and whis.lower() in {"minmax", "min-max"}:
            whislo, whishi = float(np.min(d)), float(np.max(d))
        elif isinstance(whis, (list, tuple, np.ndarray)) and len(whis) == 2:
            low_p, high_p = whis
            whislo = float(np.quantile(d, low_p/100.0, method=method))
            whishi = float(np.quantile(d, high_p/100.0, method=method))
        else:
            w = float(whis)
            low_b = q1 - w * iqr
            high_b = q3 + w * iqr
            d_lo = d[d >= low_b]
            d_hi = d[d <= high_b]
            whislo = float(d_lo.min() if d_lo.size else d.min())
            whishi = float(d_hi.max() if d_hi.size else d.max())
        fliers = d[(d < whislo) | (d > whishi)]
        stats.append(dict(med=float(med), q1=float(q1), q3=float(q3),
                          whislo=float(whislo), whishi=float(whishi),
                          fliers=fliers))

    defaults = dict(
        patch_artist=True,
        boxprops=dict(edgecolor="black", linewidth=box_linewidth),
        whiskerprops=dict(color="black", linewidth=box_linewidth),
        capprops=dict(color="black", linewidth=box_linewidth),
        medianprops=dict(color="black", linewidth=3*box_linewidth),
    )
    defaults.update(bxp_kwargs)

    if positions is None:
        positions = list(range(1, len(stats) + 1))

    artists = ax.bxp(stats, showfliers=showfliers, positions=positions, **defaults)

    for m in artists["medians"]:
        x0, x1 = m.get_xdata()
        m.set_xdata([x0 - median_overhang, x1 + median_overhang])
        m.set_solid_capstyle("round")
    for ln in artists["whiskers"] + artists["caps"]:
        ln.set_solid_capstyle("round")

    if tick_labels is not None and manage_ticks:
        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels)

    scatters = []
    subject_coords = {} if paired else None

    if showpoints:
        rng = np.random.default_rng(jitter_seed) if jitter_seed is not None else np.random.default_rng()      
        widths = bxp_kwargs.get("widths", 0.6)
        if isinstance(widths, (list, tuple, np.ndarray)):
            widths_seq = list(widths)
        else:
            widths_seq = [float(widths)] * len(stats)
        capwidths = bxp_kwargs.get("capwidths", [0.25 * w for w in widths_seq])

        if paired:
            all_subjects = np.concatenate(subjects)
            unique_subjects = np.unique(all_subjects)
            half_span = min(0.5 * float(widths_seq[0]), 0.5 * float(capwidths[0])) * jitter_frac
            subject_offsets = {s: (rng.random() - 0.5) * (2.0 * half_span) for s in unique_subjects}
            subject_coords = {s: [] for s in unique_subjects}

        for xpos, d, subs, w, cw in zip(positions, datasets, subjects if paired else [None]*len(datasets), widths_seq, capwidths):
            if d.size == 0:
                continue
            half_box = 0.5 * float(w)
            half_cap = 0.5 * float(cw)
            half_span = min(half_box, half_cap) * jitter_frac

            if paired:
                offsets = np.array([subject_offsets[s] for s in subs])
            else:
                offsets = (rng.random(d.size) - 0.5) * (2.0 * half_span)

            x = xpos + offsets
            ax.scatter(
                x, d,
                s=point_diameter**2,
                alpha=point_alpha,
                c=point_color,
                edgecolors=point_edgecolor,
                linewidths=0.5,
                zorder=3
            )

            if paired:
                for xi, yi, si in zip(x, d, subs):
                    subject_coords[si].append((xi, yi))

    if paired and subject_coords is not None:
        for si, coords in subject_coords.items():
            if len(coords) > 1:
                coords = sorted(coords, key=lambda c: c[0])
                ax.plot([c[0] for c in coords], [c[1] for c in coords],
                        linestyle=":", color="gray", alpha=0.6, linewidth=1)

    return artists

def boxplot_rtype_plotly(
    data,
    rtype: int = 7,
    whis=1.5,
    tick_labels=None,
    showfliers: bool = True,
    manage_ticks: bool = True,
    box_linewidth: float = 1.5,
    median_overhang: float = 0.05,
    showpoints: bool = False,
    jitter_frac: float = 1.0,
    point_diameter: float = 6.0,
    point_alpha: float = 0.5,
    point_color: str = "gray",
    point_edgecolor: str = "black",
    positions=None,
    paired: bool = False,
    subjects=None,
):
    if rtype not in _RTYPE2METHOD:
        raise ValueError(f"rtype must be one of {_RTYPE2METHOD.keys()}, got {rtype}")
    method = _RTYPE2METHOD[rtype]

    if hasattr(data, "__array__") and np.ndim(data) == 1:
        datasets = [np.asarray(data, float)]
    else:
        datasets = [np.asarray(d, float) for d in data]
    datasets = [d[~np.isnan(d)] for d in datasets if len(d) > 0]

    if paired and subjects is None:
        raise ValueError("If paired=True, subjects must be provided")
    if paired and len(subjects) != len(datasets):
        raise ValueError("subjects must align with data (one array per dataset)")

    stats = []
    for d in datasets:
        if d.size == 0:
            continue
        q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75], method=method)
        iqr = q3 - q1
        if isinstance(whis, str) and whis.lower() in {"minmax", "min-max"}:
            whislo, whishi = float(np.min(d)), float(np.max(d))
        elif isinstance(whis, (list, tuple, np.ndarray)) and len(whis) == 2:
            low_p, high_p = whis
            whislo = float(np.quantile(d, low_p / 100.0, method=method))
            whishi = float(np.quantile(d, high_p / 100.0, method=method))
        else:
            w = float(whis)
            low_b = q1 - w * iqr
            high_b = q3 + w * iqr
            d_lo = d[d >= low_b]
            d_hi = d[d <= high_b]
            whislo = float(d_lo.min() if d_lo.size else d.min())
            whishi = float(d_hi.max() if d_hi.size else d.max())
        fliers = d[(d < whislo) | (d > whishi)]
        stats.append(dict(med=float(med), q1=float(q1), q3=float(q3),
                          whislo=float(whislo), whishi=float(whishi),
                          fliers=fliers))

    if positions is None:
        positions = list(range(1, len(stats) + 1))

    fig = go.Figure()

    for xpos, st in zip(positions, stats):
        # Box
        fig.add_shape(type="rect",
                      x0=xpos - 0.3, x1=xpos + 0.3,
                      y0=st["q1"], y1=st["q3"],
                      line=dict(color="black", width=box_linewidth),
                      fillcolor="rgba(255,255,255,0)",
                      layer="above")

        # Median
        fig.add_shape(type="line",
                      x0=xpos - 0.3 - median_overhang,
                      x1=xpos + 0.3 + median_overhang,
                      y0=st["med"], y1=st["med"],
                      line=dict(color="black", width=3*box_linewidth))

        # Whiskers
        fig.add_shape(type="line",
                      x0=xpos, x1=xpos,
                      y0=st["whislo"], y1=st["q1"],
                      line=dict(color="black", width=box_linewidth))
        fig.add_shape(type="line",
                      x0=xpos, x1=xpos,
                      y0=st["q3"], y1=st["whishi"],
                      line=dict(color="black", width=box_linewidth))

        # Caps
        fig.add_shape(type="line",
                      x0=xpos - 0.15, x1=xpos + 0.15,
                      y0=st["whislo"], y1=st["whislo"],
                      line=dict(color="black", width=box_linewidth))
        fig.add_shape(type="line",
                      x0=xpos - 0.15, x1=xpos + 0.15,
                      y0=st["whishi"], y1=st["whishi"],
                      line=dict(color="black", width=box_linewidth))

        # Fliers
        if showfliers and st["fliers"].size > 0:
            fig.add_trace(go.Scatter(
                x=[xpos] * len(st["fliers"]),
                y=st["fliers"],
                mode="markers",
                marker=dict(size=point_diameter,
                            color=point_color,
                            opacity=point_alpha,
                            line=dict(color=point_edgecolor, width=0.5)),
                showlegend=False
            ))

    # Scatter points + paired lines
    subject_coords = {} if paired else None
    if showpoints:
        for xpos, d, subs in zip(positions, datasets, subjects if paired else [None]*len(datasets)):
            half_span = 0.25 * jitter_frac
            if paired:
                offsets = np.linspace(-half_span, half_span, num=len(d))
            else:
                offsets = (np.random.rand(d.size) - 0.5) * (2.0 * half_span)

            x = xpos + offsets
            fig.add_trace(go.Scatter(
                x=x,
                y=d,
                mode="markers",
                marker=dict(size=point_diameter,
                            color=point_color,
                            opacity=point_alpha,
                            line=dict(color=point_edgecolor, width=0.5)),
                showlegend=False
            ))

            if paired:
                for xi, yi, si in zip(x, d, subs):
                    subject_coords.setdefault(si, []).append((xi, yi))

    if paired and subject_coords is not None:
        for coords in subject_coords.values():
            if len(coords) > 1:
                coords = sorted(coords, key=lambda c: c[0])
                fig.add_trace(go.Scatter(
                    x=[c[0] for c in coords],
                    y=[c[1] for c in coords],
                    mode="lines",
                    line=dict(color="gray", dash="dot", width=1),
                    opacity=0.6,
                    showlegend=False
                ))

    # Axis labels
    if tick_labels is not None and manage_ticks:
        fig.update_xaxes(tickmode="array", tickvals=positions, ticktext=tick_labels)

    fig.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(showline=True, linewidth=1, linecolor="black", ticks="outside"),
        yaxis=dict(showline=True, linewidth=1, linecolor="black", ticks="outside")
    )

    return fig


def analyse_fluorescence(trials_data, PMT_background=162, verbose=False, plot2screen=False):
    # define helper to compute mean within window
    def mean_fun(t, y, a, b):
        mask = (t >= a) & (t < b)
        return np.nanmean(y[mask]) if np.any(mask) else np.nan

    # define helper to average traces with common time base
    def mean_traces_fun(trials):
        if not trials:
            return None, None
        t_common = np.linspace(
            max(trial['timestamps'][0] for trial in trials),
            min(trial['timestamps'][-1] for trial in trials),
            num=min(len(trial['timestamps']) for trial in trials)
        )
        ys = [np.interp(t_common, trial['timestamps'], trial['fluorescence']) for trial in trials]
        return t_common, np.mean(ys, axis=0)

    # calculate fmax from ach trials
    ach_vals = [
        mean_fun(trial['timestamps'], trial['fluorescence'] - PMT_background, 2, 4)
        for trial in trials_data if 'ach' in trial['treatment'].lower()
    ]
    Fmax = np.nanmax(ach_vals) if ach_vals else np.nan

    # calculate fmin from ttx trials
    ttx_vals = [
        mean_fun(trial['timestamps'], trial['fluorescence'] - PMT_background, 2, 4)
        for trial in trials_data if 'ttx' in trial['treatment'].lower()
    ]
    Fmin_TTX = np.nanmin(ttx_vals) if ttx_vals else np.nan

    # calculate f0_quin with fallback to last 1 s if needed
    quin_trials = [t for t in trials_data if 'quin' in t['treatment'].lower()]
    if len(quin_trials) >= 2:
        t_q, ymean_q = mean_traces_fun(quin_trials[:2])
        if t_q is not None:
            F0_quin = mean_fun(t_q, ymean_q - PMT_background, 9, 10)
            if np.isnan(F0_quin):
                F0_quin = mean_fun(t_q, ymean_q - PMT_background, t_q[-1] - 1, t_q[-1])
        else:
            F0_quin = np.nan
    else:
        F0_quin = np.nan

    # determine final fmin combining ttx and quin
    if np.isnan(Fmin_TTX) and np.isnan(F0_quin):
        Fmin = np.nan
    elif np.isnan(Fmin_TTX):
        Fmin = F0_quin
    elif np.isnan(F0_quin):
        Fmin = Fmin_TTX
    else:
        Fmin = min(Fmin_TTX, F0_quin)

    # calculate fi from fmax and fmin
    FI = Fmax - Fmin if np.isfinite(Fmax) and np.isfinite(Fmin) else np.nan

    # build results table
    df_meta = pd.DataFrame(trials_data)
    results = []

    # loop through treatments
    for treatment, subset in df_meta.groupby('treatment'):
        single_pulse = subset[subset['stimulation'] == 'single_pulse']
        if len(single_pulse) == 0:
            continue
        avg_n = min(2, len(single_pulse))
        t, ymean = mean_traces_fun(single_pulse.iloc[:avg_n].to_dict('records'))
        if t is None:
            continue

        # extract stim time and compute background corrected signal
        stim_time = float(single_pulse.iloc[0]['stim_time'])
        y_bg = ymean - PMT_background
        dt = t[1] - t[0]

        # define baseline and response windows
        baseline_start = stim_time - 1.0 - dt
        baseline_end = stim_time - dt
        post_end = stim_time + 2.0

        # compute baseline f0 and delta f
        F0 = mean_fun(t, y_bg, baseline_start, baseline_end)
        dF = y_bg - F0

        # adjust response window for sulp
        if 'sulp' in treatment.lower():
            resp_a, resp_b = baseline_start, stim_time + 2.0
        else:
            resp_a, resp_b = baseline_start, stim_time + 1.0

        # compute peak and auc metrics
        mask_resp = (t >= resp_a) & (t <= resp_b)
        peak_dF = float(np.nanmax(dF[mask_resp])) if np.any(mask_resp) else np.nan

        if verbose:
            print(f"{treatment}: F0={F0:.2f}, peak_dF={peak_dF:.2f}")

        # plot traces if enabled
        if plot2screen:
            plt.figure(figsize=(7, 4))
            plt.plot(t, y_bg, color='black', lw=1.2, label=f'{treatment}')
            plt.axvline(stim_time, color='red', linestyle='--', lw=1)
            plt.xlim(baseline_start, post_end)
            plt.xlabel('time (s)')
            plt.ylabel('fluorescence – pmt background')
            plt.title(f'corrected fluorescence – {treatment}')
            plt.tight_layout()
            plt.show()
            plt.figure(figsize=(7, 4))
            plt.plot(t, dF, color='#404040', lw=1.2, label=f'Δf ({treatment})')
            plt.axhline(0, color='black', lw=0.8)
            plt.axvline(stim_time, color='red', linestyle='--', lw=1)
            plt.xlim(baseline_start, post_end)
            plt.xlabel('time (s)')
            plt.ylabel('Δf (fluorescence − f₀)')
            plt.title(f'Δf trace – {treatment}')
            plt.tight_layout()
            plt.show()

        # calculate auc metrics
        AUC = np.trapezoid(dF[mask_resp], t[mask_resp]) if np.any(mask_resp) else np.nan
        AUC_norm_F0 = np.trapezoid((dF[mask_resp] / F0), t[mask_resp]) if (np.any(mask_resp) and F0) else np.nan
        AUC_norm_FI = np.trapezoid((dF[mask_resp] / FI), t[mask_resp]) if (np.any(mask_resp) and np.isfinite(FI) and FI != 0) else np.nan

        # collect metadata and results
        meta = single_pulse.iloc[0]
        results.append({
            'condition': meta['condition'],
            'treatment': treatment,
            'subject_id': meta['subject_id'],
            'session_id': meta['session_id'],
            'series_name': meta['series_name'],
            'file': meta['file'],
            'N': avg_n,
            'F0': F0,
            'peak dF': peak_dF,
            'AUC_norm_F0': AUC_norm_F0,
            'AUC_norm_FI': AUC_norm_FI,
            'Fmax': Fmax,
            'Fmin': Fmin,
            'FI': FI,
            'PMT_background': PMT_background
        })
    return pd.DataFrame(results)


def get_stimulation_info_any(folder):
    voltage_xmls = sorted(folder.glob("*VoltageOutput_001.xml"))
    if not voltage_xmls:
        voltage_xmls = sorted(folder.glob("*VoltageOutput*.xml"))
    if not voltage_xmls:
        return {
            "stimulation_type": "calibration",
            "stimulus_delay_ms": None,
            "pulse_count": 0,
            "pulse_width_ms": None,
            "pulse_spacing_ms": None,
            "has_stimulus": False,
        }

    root = et.parse(voltage_xmls[0]).getroot()
    for waveform in root.findall("Waveform"):
        name = (waveform.findtext("Name") or "").strip().lower()
        enabled = (waveform.findtext("Enabled") or "").strip().lower() == "true"
        if enabled and ("stim" in name or "led" in name):
            pulse_train = waveform.find("WaveformComponent_PulseTrain")
            if pulse_train is None:
                continue
            pulse_count = int(float(pulse_train.findtext("PulseCount", "1")))
            return {
                "stimulation_type": "single_pulse" if pulse_count == 1 else "burst_stimulation",
                "stimulus_delay_ms": int(float(pulse_train.findtext("FirstPulseDelay", "0"))),
                "pulse_count": pulse_count,
                "pulse_width_ms": int(float(pulse_train.findtext("PulseWidth", "0"))),
                "pulse_spacing_ms": int(float(pulse_train.findtext("PulseSpacing", "0"))),
                "has_stimulus": True,
            }

    return {
        "stimulation_type": "calibration",
        "stimulus_delay_ms": None,
        "pulse_count": 0,
        "pulse_width_ms": None,
        "pulse_spacing_ms": None,
        "has_stimulus": False,
    }

def first(row, keys, default=None):
    for k in keys:
        if k in row and row[k] not in [None, ""]:
            return row[k]
    return default

def to_int(x):
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else None

def group_from_row(row):
    raw = str(first(row, ["group", "condition_raw", "treatment", "condition"], "ctrl")).strip().lower()
    if raw in {"ctrl", "ctr", "control", "ul control"}:
        return "ctrl"
    if raw in {"test", "pd", "6-ohda", "lid off", "off-state"}:
        return "test"
    return re.sub(r"[^a-z0-9]+", "", raw) or "ctrl"

def slice_roi_from_row(row):
    slice_id = to_int(first(row, ["slice_id", "slice", "slice_label"]))
    roi_id = to_int(first(row, ["roi_id", "roi", "slice_label"]))
    if roi_id is None:
        xlsx_file = str(first(row, ["xlsx_file"], ""))
        m = re.search(r"roi\s*[-# ]*\s*(\d+)", Path(xlsx_file).stem, flags=re.IGNORECASE)
        if m:
            roi_id = int(m.group(1))
    return slice_id, roi_id

def repetition_from_row(row):
    rep = first(row, ["repetition", "repetition_number"])
    if rep not in [None, ""]:
        return int(rep)
    col = str(first(row, ["excel_col", "col"], "")).strip().upper()
    return {"B": 1, "C": 2, "D": 3}.get(col, 1)

def treatment_from_row(row):
    t = str(first(row, ["treatment"], "")).strip()
    if t:
        return re.sub(r"[^A-Za-z0-9]+", "", t)
    return "ctr" if group_from_row(row) == "ctrl" else "test"

def condition_for_f5(row):
    c = str(first(row, ["condition"], "")).strip().lower()
    if c in {"pd", "6-ohda"}:
        return "PD"
    if c in {"lid off", "off-state", "lid_off"}:
        return "LID off"
    return "UL control"

def slice_label_from_row(row):
    sl = first(row, ["slice_label"])
    if sl:
        return re.sub(r"\s+", "", str(sl))
    s, r = slice_roi_from_row(row)
    if s is not None and r is not None:
        return f"slice{s}ROI{r}"
    return "slice1ROI1"

def row_from_json(json_path, raw_folder):
    raw_folder = Path(raw_folder)
    project = raw_folder.parent.name
    trial = raw_folder.name
    key = f"{project}/{trial}"

    with open(json_path, "r") as f:
        meta = json.load(f)

    if isinstance(meta, dict):
        if key not in meta:
            raise KeyError(f"no metadata key: {key}")
        row = meta[key]
        if isinstance(row, list):
            if len(row) == 0:
                raise KeyError(f"metadata key exists but empty list: {key}")
            row = row[0]
        return row

    if isinstance(meta, list):
        for row in meta:
            p = Path(str(row.get("matched_csv_file", "")))
            if p.parent.name == trial and p.parent.parent.name == project:
                return row
        raise KeyError(f"no metadata row for: {key}")

    raise ValueError("metadata json must be dict or list")

def build_out_file(raw_folder, dest_root, row):
    raw_folder = Path(raw_folder)
    project = raw_folder.parent.name
    trial = raw_folder.name

    group = group_from_row(row)
    animal_id = to_int(first(row, ["animal_id", "animal", "subject_id"]))
    slice_id, roi_id = slice_roi_from_row(row)
    repetition = repetition_from_row(row)

    animal_txt = f"{animal_id:02d}" if animal_id is not None else "xx"
    slice_txt = f"{slice_id:02d}" if slice_id is not None else "xx"
    roi_txt = f"{roi_id:02d}" if roi_id is not None else "xx"
    rep_txt = f"{repetition:03d}" if repetition is not None else "xxx"

    m = re.match(r"^BrightnessOverTime-(\d{8})-(\d{4})-(\d+)$", trial)
    if m:
        date_token = m.group(1)
        run_token = f"{int(m.group(3)):03d}"
        out_name = (
            f"bot_{date_token}_animal{animal_txt}_{group}_"
            f"slice{slice_txt}_roi{roi_txt}_rep{rep_txt}-{run_token}.nwb"
        )
    else:
        out_name = (
            f"bot_animal{animal_txt}_{group}_"
            f"slice{slice_txt}_roi{roi_txt}_rep{rep_txt}.nwb"
        )

    out_dir = Path(dest_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / out_name

def pick_files(raw_folder):
    all_files = [p for p in raw_folder.rglob("*") if p.is_file()]

    csv_files = [p for p in all_files if p.suffix.lower() == ".csv"]
    bot_csv = [p for p in csv_files if p.name.endswith("-botData.csv")]
    csv_file = bot_csv[0] if bot_csv else (csv_files[0] if csv_files else None)

    xml_files = [p for p in all_files if p.suffix.lower() == ".xml"]
    main_xml = [p for p in xml_files if "voltageoutput" not in p.name.lower() and "voltagerecording" not in p.name.lower()]
    xml_file = main_xml[0] if main_xml else None

    voltage_xml = [p for p in xml_files if "voltageoutput" in p.name.lower()]
    voltage_xml_file = voltage_xml[0] if voltage_xml else None

    tifs = [p for p in all_files if p.suffix.lower() in [".tif", ".tiff"]]

    return csv_file, xml_file, voltage_xml_file, tifs

def _balance_tifs(tifs, max_frames=3533):
    """
    Truncate ch2 and ch3 TIF lists to the same length, capped at max_frames.
    Handles mismatched frame counts by discarding trailing frames from the
    longer channel. Sorts each channel by frame number before truncating.
    """
    import re as _re

    def _frame_num(p):
        m = _re.search(r"_(\d{6})\.ome\.tif$", p.name, _re.IGNORECASE)
        return int(m.group(1)) if m else 0

    ch2 = sorted([p for p in tifs if "_Ch2_" in p.name], key=_frame_num)
    ch3 = sorted([p for p in tifs if "_Ch3_" in p.name], key=_frame_num)
    other = [p for p in tifs if "_Ch2_" not in p.name and "_Ch3_" not in p.name]

    n = min(len(ch2), len(ch3), max_frames) if ch2 and ch3 else max_frames
    if len(ch2) != len(ch3):
        print(
            f"  [balance_tifs] ch2={len(ch2)} ch3={len(ch3)} frames → "
            f"truncating both to {n}"
        )
    return ch2[:n] + ch3[:n] + other
    
def build_nwb_index(trial_metadata_path, dest_root):
    with open(trial_metadata_path, "r") as f:
        meta = json.load(f)
    index = {}
    for _key, row in meta.items():
        raw_folder = Path(row["project"]) / row["trial_folder"]
        nwb_filename = build_out_file(raw_folder, dest_root, row).name
        entry = dict(row)
        entry.pop("matched_csv_file", None)
        if "sources" in entry:
            sources = [s for s in entry["sources"] if "paper representation" not in s.lower()]
            entry["sources"] = sources[0] if sources else None
        entry["nwb_file"] = nwb_filename
        index[nwb_filename] = entry
    return index


def compute_dff(f_raw, time, stim_time, pmt_background=0.0, baseline_window=1.0):
    """
    Compute dF/F0.

    Parameters
    ----------
    f_raw : np.ndarray, shape (n_frames,)
        Raw fluorescence trace (mean over ROI).
    time : np.ndarray, shape (n_frames,)
        Time vector in seconds.
    stim_time : float
        Stimulus onset time (s). Baseline is [stim_time - baseline_window, stim_time).
    pmt_background : float
        Constant PMT background to subtract before computing dF/F. Default 0.
    baseline_window : float
        Duration (s) of the pre-stimulus baseline. Default 1.0 s.

    Returns
    -------
    dff : np.ndarray, shape (n_frames,)
    f0  : float   (baseline mean, useful for diagnostics)
    """
    f = f_raw - pmt_background
    baseline_mask = (time >= stim_time - baseline_window) & (time < stim_time)
    if not baseline_mask.any():
        raise ValueError(
            f"No baseline frames found in [{stim_time - baseline_window}, {stim_time}). "
            f"Check stim_time ({stim_time}) and t range ({time.min():.2f}–{time.max():.2f})."
        )
    f0 = float(np.mean(f[baseline_mask]))
    if abs(f0) < 1e-9:
        raise ValueError(f"Baseline f0 ≈ 0 ({f0:.4g}); check ROI or background subtraction.")
    dff = (f - f0) / f0
    return dff
    
def plot_dff(time, dff, stim_time, xlim=None, ylim=None, lwd=1, tick_len=4, axis_color="black",
             vmin=None, vmax=None, dp=2, savepath=None):
    """
    Plot dF/F0 as a line trace (top) and a 1-row heatmap (bottom).
    Parameters
    ----------
    time      : np.ndarray   time vector (s)
    dff       : np.ndarray   dF/F0 trace
    stim_time : float        stimulus onset time (s)
    xlim      : [float, float] or None   x-axis limits; defaults to [floor(t.min()), ceil(t.max())]
    ylim      : [float, float] or None   y-axis limits for trace panel (default None -> auto)
    lwd       : float        line / spine width (default 1)
    tick_len  : float        tick length in points (default 4)
    axis_color: str          colour for spines and tick labels (default "black")
    vmin      : float or None  colorbar minimum (default None -> 0.0)
    vmax      : float or None  colorbar maximum (default None -> data max)
    dp        : int          decimal places for colorbar min/max labels (default 1)
    """
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    plt.close('all')
    
    x0, x1 = xlim if xlim is not None else [np.floor(time.min()), np.ceil(time.max())]
    
    plt.rcParams['savefig.format'] = 'svg'
    fig, ax = plt.subplots(
        2, 1, figsize=(10, 6.2), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2.1, 1.2]},
    )
    # main plot
    ax[0].plot(time, dff, color="black", lw=lwd)
    ax[0].axvline(stim_time, color="firebrick", ls="--", lw=lwd)
    ax[0].axhline(0, color="gray", ls=":", lw=lwd)
    ax[0].set_ylabel("dF/F0")
    ax[0].set_title("dF/F0 trace")
    ax[0].grid(False)
    ax[0].format_coord = lambda x, y: ""
    if ylim is not None:
        ax[0].set_ylim(ylim)
    # heatmap
    _vmin = round(0.0 if vmin is None else vmin, dp)
    _vmax = round(max(float(np.nanmax(dff)), 1e-9) if vmax is None else vmax, dp)
    im = ax[1].imshow(
        dff[np.newaxis, :],
        aspect="auto", origin="lower",
        extent=[time.min(), time.max(), -0.5, 0.5],
        cmap="magma", vmin=_vmin, vmax=_vmax, interpolation="bicubic",
    )
    ax[1].format_coord = lambda x, y: ""
    im.format_cursor_data = lambda data: ""
    ax[1].axvline(stim_time, color="firebrick", ls="--", lw=lwd)
    ax[1].set_xlabel("time (s)")
    ax[1].set_title("dF/F0 heatmap")
    ax[1].set_yticks([])
    ax[1].set_ylabel("")
    ax[1].spines["left"].set_visible(False)
    ax[1].tick_params(axis="y", left=False, right=False, labelleft=False)
    # style
    for a in ax:
        a.set_xlim(x0, x1)
        a.xaxis.set_major_locator(MultipleLocator(1.0))
        a.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        a.minorticks_off()
        a.tick_params(axis="x", which="major", bottom=True, top=False,
                      length=tick_len, width=lwd, direction="out",
                      color=axis_color, labelcolor=axis_color)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_linewidth(lwd)
        a.spines["bottom"].set_color(axis_color)
    ax[0].spines["left"].set_linewidth(lwd)
    ax[0].spines["left"].set_color(axis_color)
    ax[0].tick_params(axis="y", which="major", left=True, right=False,
                      length=tick_len, width=lwd, direction="out",
                      color=axis_color, labelcolor=axis_color)
    ax[1].tick_params(axis="x", length=0)
    cbar = fig.colorbar(im, ax=ax[1], fraction=0.32, pad=0.02)
    cbar.set_label("dF/F0")
    cbar.outline.set_visible(False)
    cbar.set_ticks([_vmin, _vmax])
    cbar.ax.tick_params(length=0)

    if savepath is not None:
        fig.savefig(savepath, format="svg", bbox_inches="tight")

# ---------------------------------------------------------------------------
# Auto-patches applied on import
# ---------------------------------------------------------------------------
from prairie_view_fluorescence_interface import (
    PrairieViewFluorescenceInterface as _PVFInterface,
)
from neuroconv.converters import BrukerTiffSinglePlaneConverter as _BrukerConverter

# Rename whatever channel key is present to "Ch2"
if not hasattr(_PVFInterface, "_orig_load_roi_region_data"):
    _PVFInterface._orig_load_roi_region_data = _PVFInterface._load_roi_region_data

def load_roi_region_data_force_ch2(self):
    self._orig_load_roi_region_data()
    keys = list(self.regions_info.keys())
    if "Ch2" not in keys and len(keys) > 0:
        gaasp = [k for k in keys if "gaasp" in k.lower()]
        src = gaasp[0] if gaasp else keys[0]
        self.regions_info["Ch2"] = self.regions_info.pop(src)

_PVFInterface._load_roi_region_data = load_roi_region_data_force_ch2

# Rename TwoPhotonSeriesGaAsP -> TwoPhotonSeriesCh2 in metadata
if not hasattr(_BrukerConverter, "_orig_get_metadata"):
    _BrukerConverter._orig_get_metadata = _BrukerConverter.get_metadata

def get_metadata_force_grabch(self, *args, **kwargs):
    md = self._orig_get_metadata(*args, **kwargs)
    for s in md.get("Ophys", {}).get("TwoPhotonSeries", []):
        if s.get("name") == "TwoPhotonSeriesGaAsP":
            s["name"] = "TwoPhotonSeriesCh2"
    return md

_BrukerConverter.get_metadata = get_metadata_force_grabch
    
