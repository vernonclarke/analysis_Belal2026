import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go


import json
import re
from pathlib import Path
from pynwb import NWBHDF5IO

__all__ = [
    "compute_dff",
    "plot_dff",
    "boxplot_rtype",
    "boxplot_rtype_plotly",
    "parse_filename_meta",
    "load_nwb_fluorescence",
    "load_dataset_description_index",
    "load_all_nwb",
]


# Map R's quantile types (1-9) to NumPy methods.
_RTYPE2METHOD = {
    1: "inverted_cdf",
    2: "averaged_inverted_cdf",
    3: "closest_observation",
    4: "interpolated_inverted_cdf",
    5: "hazen",
    6: "weibull",
    7: "linear",
    8: "median_unbiased",
    9: "normal_unbiased",
}


def compute_dff(f_raw, time, stim_time, pmt_background=0.0, baseline_window=1.0):
    """Compute dF/F0 using the pre-stimulus baseline window."""
    f_raw = np.asarray(f_raw, dtype=float)
    time = np.asarray(time, dtype=float)
    f = f_raw - pmt_background

    baseline_mask = (time >= stim_time - baseline_window) & (time < stim_time)
    if not baseline_mask.any():
        raise ValueError(
            f"No baseline frames found in [{stim_time - baseline_window}, {stim_time}). "
            f"Check stim_time ({stim_time}) and t range ({time.min():.2f}-{time.max():.2f})."
        )

    f0 = float(np.mean(f[baseline_mask]))
    if abs(f0) < 1e-9:
        raise ValueError(f"Baseline f0 is near zero ({f0:.4g}); check ROI or background subtraction.")

    return (f - f0) / f0


def plot_dff(
    time,
    dff,
    stim_time,
    xlim=None,
    ylim=None,
    lwd=1,
    tick_len=4,
    axis_color="black",
    vmin=None,
    vmax=None,
    dp=2,
    savepath=None,
):
    """Plot dF/F0 as a line trace plus a one-row heatmap."""
    from matplotlib.ticker import FormatStrFormatter, MultipleLocator

    time = np.asarray(time, dtype=float)
    dff = np.asarray(dff, dtype=float)
    x0, x1 = xlim if xlim is not None else [np.floor(time.min()), np.ceil(time.max())]

    plt.close("all")
    plt.rcParams["savefig.format"] = "svg"
    fig, ax = plt.subplots(
        2,
        1,
        figsize=(10, 6.2),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.1, 1.2]},
    )

    ax[0].plot(time, dff, color="black", lw=lwd)
    ax[0].axvline(stim_time, color="firebrick", ls="--", lw=lwd)
    ax[0].axhline(0, color="gray", ls=":", lw=lwd)
    ax[0].set_ylabel("dF/F0")
    ax[0].set_title("dF/F0 trace")
    ax[0].grid(False)
    ax[0].format_coord = lambda x, y: ""
    if ylim is not None:
        ax[0].set_ylim(ylim)

    _vmin = round(0.0 if vmin is None else vmin, dp)
    _vmax = round(max(float(np.nanmax(dff)), 1e-9) if vmax is None else vmax, dp)
    im = ax[1].imshow(
        dff[np.newaxis, :],
        aspect="auto",
        origin="lower",
        extent=[time.min(), time.max(), -0.5, 0.5],
        cmap="magma",
        vmin=_vmin,
        vmax=_vmax,
        interpolation="bicubic",
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

    for a in ax:
        a.set_xlim(x0, x1)
        a.xaxis.set_major_locator(MultipleLocator(1.0))
        a.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        a.minorticks_off()
        a.tick_params(
            axis="x",
            which="major",
            bottom=True,
            top=False,
            length=tick_len,
            width=lwd,
            direction="out",
            color=axis_color,
            labelcolor=axis_color,
        )
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_linewidth(lwd)
        a.spines["bottom"].set_color(axis_color)

    ax[0].spines["left"].set_linewidth(lwd)
    ax[0].spines["left"].set_color(axis_color)
    ax[0].tick_params(
        axis="y",
        which="major",
        left=True,
        right=False,
        length=tick_len,
        width=lwd,
        direction="out",
        color=axis_color,
        labelcolor=axis_color,
    )
    ax[1].tick_params(axis="x", length=0)

    cbar = fig.colorbar(im, ax=ax[1], fraction=0.32, pad=0.02)
    cbar.set_label("dF/F0")
    cbar.outline.set_visible(False)
    cbar.set_ticks([_vmin, _vmax])
    cbar.ax.tick_params(length=0)

    if savepath is not None:
        fig.savefig(savepath, format="svg", bbox_inches="tight")

    plt.show()


def _boxplot_stats(data, rtype, whis):
    if rtype not in _RTYPE2METHOD:
        raise ValueError(f"rtype must be 1..9, got {rtype}")
    method = _RTYPE2METHOD[rtype]

    if hasattr(data, "__array__") and np.ndim(data) == 1:
        datasets = [np.asarray(data, float)]
    else:
        datasets = [np.asarray(d, float) for d in data]
    datasets = [d[~np.isnan(d)] for d in datasets if len(d) > 0]

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
        stats.append(
            {
                "med": float(med),
                "q1": float(q1),
                "q3": float(q3),
                "whislo": float(whislo),
                "whishi": float(whishi),
                "fliers": fliers,
            }
        )

    return datasets, stats

def boxplot_rtype(
    ax,
    data,
    rtype=7,
    whis=1.5,
    tick_labels=None,
    showfliers=True,
    manage_ticks=True,
    box_linewidth=1.5,
    median_overhang=0.05,
    showpoints=False,
    jitter_frac=1.0,
    jitter_seed=42,
    point_diameter=6.0,
    point_alpha=0.5,
    point_color="gray",
    point_edgecolor="black",
    positions=None,
    paired=False,
    subjects=None,
    **bxp_kwargs,
):
    """Matplotlib boxplot using R-compatible quantile definitions."""
    datasets, stats = _boxplot_stats(data=data, rtype=rtype, whis=whis)

    if paired and subjects is None:
        raise ValueError("If paired=True, subjects must be provided")
    if paired and len(subjects) != len(datasets):
        raise ValueError("subjects must align with data (one array per dataset)")

    defaults = {
        "patch_artist": True,
        "boxprops": {"edgecolor": "black", "linewidth": box_linewidth},
        "whiskerprops": {"color": "black", "linewidth": box_linewidth},
        "capprops": {"color": "black", "linewidth": box_linewidth},
        "medianprops": {"color": "black", "linewidth": 3 * box_linewidth},
    }
    defaults.update(bxp_kwargs)

    if positions is None:
        positions = list(range(1, len(stats) + 1))

    artists = ax.bxp(stats, showfliers=showfliers, positions=positions, **defaults)

    for median in artists["medians"]:
        x0, x1 = median.get_xdata()
        median.set_xdata([x0 - median_overhang, x1 + median_overhang])
        median.set_solid_capstyle("round")
    for line in artists["whiskers"] + artists["caps"]:
        line.set_solid_capstyle("round")

    if tick_labels is not None and manage_ticks:
        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels)

    subject_coords = {} if paired else None
    if showpoints:
        rng = np.random.default_rng(jitter_seed) if jitter_seed is not None else np.random.default_rng()
        widths = bxp_kwargs.get("widths", 0.6)
        widths_seq = list(widths) if isinstance(widths, (list, tuple, np.ndarray)) else [float(widths)] * len(stats)
        capwidths = bxp_kwargs.get("capwidths", [0.25 * w for w in widths_seq])

        if paired:
            all_subjects = np.concatenate(subjects)
            unique_subjects = np.unique(all_subjects)
            half_span = min(0.5 * float(widths_seq[0]), 0.5 * float(capwidths[0])) * jitter_frac
            subject_offsets = {s: (rng.random() - 0.5) * (2.0 * half_span) for s in unique_subjects}
            subject_coords = {s: [] for s in unique_subjects}

        subject_iter = subjects if paired else [None] * len(datasets)
        for xpos, d, subs, w, cap_width in zip(positions, datasets, subject_iter, widths_seq, capwidths):
            if d.size == 0:
                continue
            half_span = min(0.5 * float(w), 0.5 * float(cap_width)) * jitter_frac
            offsets = np.array([subject_offsets[s] for s in subs]) if paired else (rng.random(d.size) - 0.5) * (2.0 * half_span)
            x = xpos + offsets
            ax.scatter(
                x,
                d,
                s=point_diameter**2,
                alpha=point_alpha,
                c=point_color,
                edgecolors=point_edgecolor,
                linewidths=0.5,
                zorder=3,
            )

            if paired:
                for xi, yi, subject in zip(x, d, subs):
                    subject_coords[subject].append((xi, yi))

    if paired and subject_coords is not None:
        for coords in subject_coords.values():
            if len(coords) > 1:
                coords = sorted(coords, key=lambda c: c[0])
                ax.plot(
                    [c[0] for c in coords],
                    [c[1] for c in coords],
                    linestyle=":",
                    color="gray",
                    alpha=0.6,
                    linewidth=1,
                )

    return artists

def boxplot_rtype_plotly(
    data,
    rtype=7,
    whis=1.5,
    tick_labels=None,
    showfliers=True,
    manage_ticks=True,
    box_linewidth=1.5,
    median_overhang=0.05,
    showpoints=False,
    jitter_frac=1.0,
    point_diameter=6.0,
    point_alpha=0.5,
    point_color="gray",
    point_edgecolor="black",
    positions=None,
    widths=0.6,
    xrange=None,
    yrange=None,
    title=None,
    xaxis_title=None,
    yaxis_title=None,
    width=None,
    height=None,
    paired=False,
    subjects=None,
):
    """Plotly boxplot using R-compatible quantile definitions."""
    datasets, stats = _boxplot_stats(data=data, rtype=rtype, whis=whis)

    if paired and subjects is None:
        raise ValueError("If paired=True, subjects must be provided")
    if paired and len(subjects) != len(datasets):
        raise ValueError("subjects must align with data (one array per dataset)")
    if positions is None:
        positions = list(range(1, len(stats) + 1))
    xaxis_range = [min(positions) - 0.5, max(positions) + 0.5] if xrange is None else list(xrange)
    widths_seq = list(widths) if isinstance(widths, (list, tuple, np.ndarray)) else [float(widths)] * len(stats)
    capwidths = [0.25 * w for w in widths_seq]

    def _point_color_for_dataset(index, values):
        if isinstance(point_color, (list, tuple, np.ndarray, pd.Series)):
            if len(point_color) == len(datasets):
                dataset_color = point_color[index]
                if isinstance(dataset_color, (list, tuple, np.ndarray, pd.Series)):
                    if len(dataset_color) != len(values):
                        raise ValueError("Each point_color sequence must match its dataset length")
                    return list(dataset_color)
                return dataset_color
            if len(point_color) == len(values):
                return list(point_color)
        return point_color

    fig = go.Figure()
    for xpos, stat, box_width, cap_width in zip(positions, stats, widths_seq, capwidths):
        fig.add_shape(
            type="rect",
            x0=xpos - 0.5 * box_width,
            x1=xpos + 0.5 * box_width,
            y0=stat["q1"],
            y1=stat["q3"],
            line={"color": "gray", "width": box_linewidth},
            fillcolor="rgba(0,0,0,0)",
            layer="above",
        )
        fig.add_shape(
            type="line",
            x0=xpos - 0.5 * box_width - median_overhang,
            x1=xpos + 0.5 * box_width + median_overhang,
            y0=stat["med"],
            y1=stat["med"],
            line={"color": "gray", "width": 3 * box_linewidth},
        )
        fig.add_shape(
            type="line",
            x0=xpos,
            x1=xpos,
            y0=stat["whislo"],
            y1=stat["q1"],
            line={"color": "gray", "width": box_linewidth},
        )
        fig.add_shape(
            type="line",
            x0=xpos,
            x1=xpos,
            y0=stat["q3"],
            y1=stat["whishi"],
            line={"color": "gray", "width": box_linewidth},
        )
        fig.add_shape(
            type="line",
            x0=xpos - 0.5 * cap_width,
            x1=xpos + 0.5 * cap_width,
            y0=stat["whislo"],
            y1=stat["whislo"],
            line={"color": "gray", "width": box_linewidth},
        )
        fig.add_shape(
            type="line",
            x0=xpos - 0.5 * cap_width,
            x1=xpos + 0.5 * cap_width,
            y0=stat["whishi"],
            y1=stat["whishi"],
            line={"color": "gray", "width": box_linewidth},
        )

        if showfliers and stat["fliers"].size > 0:
            fig.add_trace(
                go.Scatter(
                    x=[xpos] * len(stat["fliers"]),
                    y=stat["fliers"],
                    mode="markers",
                    marker={
                        "size": point_diameter,
                        "color": point_color,
                        "opacity": point_alpha,
                        "line": {"color": point_edgecolor, "width": 0.5},
                    },
                    showlegend=False,
                )
            )

    subject_coords = {} if paired else None
    if showpoints:
        subject_iter = subjects if paired else [None] * len(datasets)
        for dataset_index, (xpos, d, subs, w, cap_width) in enumerate(zip(positions, datasets, subject_iter, widths_seq, capwidths)):
            half_span = min(0.5 * float(w), 0.5 * float(cap_width)) * jitter_frac
            offsets = np.linspace(-half_span, half_span, num=len(d)) if paired else (np.random.rand(d.size) - 0.5) * (2.0 * half_span)
            x = xpos + offsets
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=d,
                    mode="markers",
                    marker={
                        "size": point_diameter,
                        "color": _point_color_for_dataset(dataset_index, d),
                        "opacity": point_alpha,
                        "line": {"color": point_edgecolor, "width": 0.5},
                    },
                    showlegend=False,
                )
            )

            if paired:
                for xi, yi, subject in zip(x, d, subs):
                    subject_coords.setdefault(subject, []).append((xi, yi))

    if paired and subject_coords is not None:
        for coords in subject_coords.values():
            if len(coords) > 1:
                coords = sorted(coords, key=lambda c: c[0])
                fig.add_trace(
                    go.Scatter(
                        x=[c[0] for c in coords],
                        y=[c[1] for c in coords],
                        mode="lines",
                        line={"color": "gray", "dash": "dot", "width": 1},
                        opacity=0.6,
                        showlegend=False,
                    )
                )

    if tick_labels is not None and manage_ticks:
        fig.update_xaxes(tickmode="array", tickvals=positions, ticktext=tick_labels)

    y_max_values = []
    for stat in stats:
        y_max_values.extend([stat["whishi"], stat["q3"], stat["med"]])
        if showfliers and stat["fliers"].size > 0:
            y_max_values.extend(stat["fliers"].tolist())
    for dataset in datasets:
        if dataset.size:
            y_max_values.append(float(np.nanmax(dataset)))
    y_max = max(y_max_values) if y_max_values else 1.0
    y_range_top = y_max * 1.05 if y_max > 0 else 1.0
    yaxis_range = [0, y_range_top] if yrange is None else list(yrange)

    fig.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        title={"text": title, "x": 0.5, "xanchor": "center", "xref": "paper"},
        font={"color": "gray"},
        width=width,
        height=height,
        xaxis={
            "showline": True,
            "showgrid": False,
            "zeroline": False,
            "range": xaxis_range,
            "linewidth": 1,
            "linecolor": "gray",
            "ticks": "outside",
            "tickcolor": "gray",
            "tickfont": {"color": "gray"},
            "tickangle": 45,
            "title": {"text": xaxis_title, "font": {"color": "gray"}},
        },
        yaxis={
            "showline": True,
            "showgrid": False,
            "zeroline": False,
            "range": yaxis_range,
            "linewidth": 1,
            "linecolor": "gray",
            "ticks": "outside",
            "tickcolor": "gray",
            "tickfont": {"color": "gray"},
            "title": {"text": yaxis_title, "font": {"color": "gray"}},
        },
    )
    return fig

def parse_filename_meta(nwb_path):
    stem = Path(nwb_path).stem

    m = re.match(
        r"sub-(?P<group>[^-]+)-animal(?P<animal>\d+)-slice(?P<slice>\d+)-roi(?P<roi>\d+)"
        r"_ses-.*?(?P<datetime>\d{14})_ophys$",
        stem,
    )
    if not m:
        return {}

    return {
        "date": m.group("datetime")[:8],
        "time": m.group("datetime")[8:],
        "animal_id": int(m.group("animal")),
        "group": m.group("group"),
        "slice_id": int(m.group("slice")),
        "roi_id": int(m.group("roi")),
        "slice_label": f"slice{int(m.group('slice'))}ROI{int(m.group('roi'))}",
    }

def load_nwb_fluorescence(nwb_path, json_meta=None):
    nwb_path = Path(nwb_path)
    records  = []
    with NWBHDF5IO(str(nwb_path), mode="r", load_namespaces=True) as io:
        nwbfile = io.read()
        meta = {}
        json_meta = json_meta or {}
        meta["nwb_file"]   = nwb_path.name
        meta["key"]        = json_meta.get("dandi_path", nwbfile.identifier)
        meta["nwb_identifier"] = nwbfile.identifier
        meta.update(json_meta)
        meta["session_id"] = nwbfile.session_id
        meta["project"]    = nwbfile.session_id.split("/")[0] if nwbfile.session_id else None
        if nwbfile.subject is not None:
            meta["subject_id"] = nwbfile.subject.subject_id
        meta.update({k: v for k, v in parse_filename_meta(nwb_path).items()
                     if k not in meta})

        trials_df = nwbfile.trials.to_dataframe() if nwbfile.trials is not None else None

        if "ophys" not in nwbfile.processing:
            print(f"[WARN] no ophys in {nwb_path.name}")
            return records

        fluor = nwbfile.processing["ophys"]["Fluorescence"]
        for series_name, series in fluor.roi_response_series.items():
            fluorescence = series.data[:]
            timestamps   = series.get_timestamps()

            stim_time = None
            if trials_df is not None and "roi_series_name" in trials_df.columns:
                match = trials_df[trials_df["roi_series_name"] == series_name]
                if not match.empty:
                    stim_time = float(match.iloc[0]["stimulus_start_time"])
                    for col in ("treatment", "stimulation", "group"):
                        if col in match.columns:
                            meta[col] = match.iloc[0][col]

            records.append({
                **meta,
                "series_name":  series_name,
                "fluorescence": fluorescence,
                "timestamps":   timestamps,
                "stim_time":    stim_time,
                "rate_hz":      getattr(series, "rate", None),
                "n_frames":     len(fluorescence),
            })
    return records

def load_dataset_description_index(nwb_dir):
    desc_path = Path(nwb_dir) / "dataset_description" / "dataset_description.json"
    with open(desc_path, "r") as f:
        desc = json.load(f)

    index = {}
    for figure, entries in desc.get("FigureMappings", {}).items():
        for entry in entries:
            dandi_path = entry["dandi_path"]
            index[Path(nwb_dir) / dandi_path] = {
                "key": dandi_path,
                "figure": figure,
                "original_path": entry.get("original_path"),
                "dandi_path": dandi_path,
            }
    return index


def load_all_nwb(nwb_dir):
    nwb_dir = Path(nwb_dir)
    nwb_index = load_dataset_description_index(nwb_dir)
    nwb_files = sorted(
        path for path in nwb_index
        if path.exists() and path.name.endswith("_ophys.nwb")
    )

    if not nwb_files:
        print(f"No .nwb files found from dataset_description.json under {nwb_dir}")
        return pd.DataFrame()

    all_records = []
    for f in nwb_files:
        print(f"  loading {f.name} …", end=" ")
        recs = load_nwb_fluorescence(f, nwb_index[f])
        print(f"{len(recs)} series")
        all_records.extend(recs)

    df = pd.DataFrame(all_records)

    if not df.empty:
        df = df.sort_values("nwb_file").reset_index(drop=True)
        df["run_id"] = (
            df.groupby(["group", "animal_id", "slice_id", "roi_id", "series_name"])
              .cumcount()
            + 1
        )

    print(f"\n→ {len(df)} total records from {len(nwb_files)} files")
    return df
