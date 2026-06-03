"""
Standalone RNAscope conversion/source-loading workflow.

This module separates:
1. Conversion of Olympus RNAscope folders to a self-contained NWB
2. Loading field data for analysis from either the original folder or the NWB

Stored in NWB per field:
- raw source TIFFs: s_C00x.tif
- raw-image metadata: s_C00x.pty
- exact display reconstruction metadata: s_LUTx.lut
- exported per-channel display TIFFs: ..._C001.tif, ..._C002.tif, ...

Not stored as canonical conversion assets:
- composite display TIFFs such as GB / RB
"""

import base64
import json
import re
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector
from matplotlib.patches import Polygon, Circle
from matplotlib.path import Path as MplPath
from IPython.display import display, Markdown
from PIL import Image
from scipy import ndimage as ndi
from pynwb import NWBHDF5IO, NWBFile
from pynwb.file import Subject
from pynwb.image import GrayscaleImage, Images, RGBImage
from pynwb.misc import DynamicTable


FIELD_RE = re.compile(r"^(?P<slice_id>.+?)[._](?P<hemisphere>UL|L)_60x\.?(?P<field_index>\d+)$")


def read_utf16_text(path: Path) -> str:
    return Path(path).read_bytes().decode("utf-16", errors="ignore")


def _field_id_from_name(field_name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(field_name).strip()).strip("_")


def _normalize_field(field: str) -> str:
    field = str(field)
    for suffix in (".oif.files", ".oif", ".tif.frames"):
        if field.endswith(suffix):
            field = field[: -len(suffix)]
    return field


def _normalize_channel_name(channel: str) -> str:
    channel = str(channel).strip()
    match = re.fullmatch(r"(?:s_)?[Cc](\d{3})", channel)
    if not match:
        raise ValueError(f"Invalid channel name {channel!r}. Expected like 's_C001'")
    return f"s_C{match.group(1)}"


def _channel_sort_key(channel: str) -> int:
    return int(_normalize_channel_name(channel)[-3:])


def _normalize_roi_layout(roi_groups=("GB", "RB"), roi_channels=("s_C002", "s_C003")):
    if isinstance(roi_groups, str):
        roi_groups = (roi_groups,)
    if isinstance(roi_channels, str):
        roi_channels = (roi_channels,)

    roi_groups = tuple(str(g).strip() for g in roi_groups)
    roi_channels = tuple(_normalize_channel_name(ch) for ch in roi_channels)

    if not roi_groups:
        raise ValueError("roi_groups must contain at least one group")
    if len(set(roi_groups)) != len(roi_groups):
        raise ValueError("roi_groups must be unique")
    if len(roi_groups) != len(roi_channels):
        raise ValueError("roi_groups and roi_channels must have the same length")

    return roi_groups, roi_channels


def _normalize_roi_specs(roi_specs=None, roi_groups=("GB", "RB"), roi_channels=("s_C002", "s_C003")):
    if roi_specs is None:
        roi_groups, roi_channels = _normalize_roi_layout(roi_groups, roi_channels)
        return tuple({"group": g, "channel": ch} for g, ch in zip(roi_groups, roi_channels))

    if isinstance(roi_specs, dict):
        roi_specs = (roi_specs,)

    normalized = []
    for spec in roi_specs:
        if isinstance(spec, (tuple, list)) and len(spec) == 2:
            group, channel = spec
        elif isinstance(spec, dict):
            if "group" not in spec or "channel" not in spec:
                raise ValueError("Each roi_spec must contain 'group' and 'channel'")
            group, channel = spec["group"], spec["channel"]
        else:
            raise ValueError("Each roi_spec must be a dict or a (group, channel) pair")

        group = str(group).strip()
        if not group:
            raise ValueError("ROI group names cannot be empty")

        normalized.append(
            {
                "group": group,
                "channel": _normalize_channel_name(channel),
            }
        )

    groups = [spec["group"] for spec in normalized]
    if len(set(groups)) != len(groups):
        raise ValueError("ROI group names must be unique")

    return tuple(normalized)


def _roi_groups_from_specs(roi_specs):
    return tuple(spec["group"] for spec in roi_specs)


def _roi_channels_from_specs(roi_specs):
    return tuple(spec["channel"] for spec in roi_specs)


def _compact_channel_label(channel):
    channel = str(channel).strip()
    match = re.fullmatch(r"(?:s_)?[Cc](\d{3})", channel)
    return f"C{match.group(1)}" if match else channel


def _format_start_title(group, roi_channel):
    return (
        f"{group} ({_compact_channel_label(roi_channel)}) | "
        "draw polygon, double-click close, Backspace undo, q finish"
    )


def _draw_roi_editor(state, group, draw_key, image_used_key, title_text=None):
    img = state["images"][draw_key]
    rois = [np.asarray(v, dtype=float) for v in state["roi_sets"][group]]

    fig, ax = plt.subplots(figsize=(8, 8))
    selector_holder = {"selector": None}
    key_cid = {"id": None}

    def redraw():
        ax.clear()
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(
            title_text
            or f"{group} on {draw_key} | draw polygon, double-click close, Backspace undo, q finish"
        )
        for idx, verts in enumerate(rois, start=1):
            patch = Polygon(
                verts,
                closed=True,
                fill=False,
                edgecolor="white",
                linewidth=1,
                linestyle=(0, (1.2, 2.0)),
            )
            ax.add_patch(patch)
            center = verts.mean(axis=0)
            ax.text(center[0], center[1], str(idx), color="white", ha="center", va="center", fontsize=10)
        fig.canvas.draw_idle()

    def onselect(verts):
        verts = np.asarray(verts, dtype=float)
        if len(verts) < 3:
            return
        if rois and rois[-1].shape == verts.shape and np.allclose(rois[-1], verts):
            return
        rois.append(verts.copy())
        state["roi_sets"][group] = [r.copy() for r in rois]
        state["roi_image_used"][group] = image_used_key
        print(f"{group}: stored {len(rois)} ROI(s)")
        redraw()
        reset_selector()

    def reset_selector():
        if selector_holder["selector"] is not None:
            selector_holder["selector"].disconnect_events()
        selector_holder["selector"] = PolygonSelector(
            ax,
            onselect,
            useblit=True,
            props=dict(color="white", linewidth=1, alpha=0.6),
            handle_props=dict(marker="o", markersize=5, markerfacecolor="cyan", markeredgecolor="black"),
        )

    def cleanup():
        if selector_holder["selector"] is not None:
            selector_holder["selector"].disconnect_events()
            selector_holder["selector"] = None
        if key_cid["id"] is not None:
            fig.canvas.mpl_disconnect(key_cid["id"])
            key_cid["id"] = None

    def on_key(event):
        if event.key == "backspace" and rois:
            rois.pop()
            state["roi_sets"][group] = [r.copy() for r in rois]
            if not rois:
                state["roi_image_used"][group] = None
            print(f"{group}: remaining {len(rois)} ROI(s)")
            redraw()
            reset_selector()
        elif event.key in ("q", "escape", "enter"):
            cleanup()
            plt.close(fig)

    key_cid["id"] = fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("close_event", lambda event: cleanup())

    redraw()
    reset_selector()
    plt.show()

    state["_ui"][group] = {"fig": fig, "draw_key": draw_key, "image_used_key": image_used_key}


def parse_capture_time(oif_text: str, tz: ZoneInfo):
    match = re.search(r"ImageCaputreDate='([^']+)'", oif_text) or re.search(
        r"ImageCaptureDate='([^']+)'", oif_text
    )
    if not match:
        return None

    dt = datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S").replace(tzinfo=tz)
    ms_match = re.search(r"ImageCaputreDate\+MilliSec=(\d+)", oif_text) or re.search(
        r"ImageCaptureDate\+MilliSec=(\d+)", oif_text
    )
    if ms_match:
        dt += timedelta(milliseconds=int(ms_match.group(1)))
    return dt


def parse_side_and_index(oif_name: str):
    oif_name = Path(oif_name).name
    match = re.search(r"(?:_C)?\.(UL|L)_60x\.?(\d+)\.oif$", oif_name)
    if not match:
        raise ValueError(f"Could not parse side/index from OIF filename: {oif_name}")
    return match.group(1), match.group(2).zfill(2)


def parse_channel_dyes(oif_text: str):
    out = {}
    for match in re.finditer(r"\[Channel (\d+) Parameters\](.*?)(?:\n\[|$)", oif_text, re.S):
        channel_number = int(match.group(1))
        block = match.group(2)
        dye_match = re.search(r'DyeName="([^"]+)"', block)
        out[f"s_C{channel_number:03d}"] = dye_match.group(1) if dye_match else None
    return out


def parse_pty_text(pty_text: str):
    parsed = {}
    section = None
    for raw_line in pty_text.splitlines():
        line = raw_line.strip("\r")
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            parsed.setdefault(section, {})
            continue
        if "=" in line and section is not None:
            key, value = line.split("=", 1)
            parsed[section][key] = value.strip().strip('"')
    return parsed


def parse_pty_file(pty_path: Path):
    text = read_utf16_text(pty_path)
    return text, parse_pty_text(text)


def parse_pixel_size_um_from_pty(pty_path: Path):
    _, parsed = parse_pty_file(pty_path)
    return _maybe_float(parsed.get("Image Parameters", {}).get("WidthConvertValue"))


def parse_resolution_cm_from_pty(pty_path: Path):
    pixel_size_um = parse_pixel_size_um_from_pty(pty_path)
    if pixel_size_um is None:
        return None
    return 10000.0 / pixel_size_um


def parse_lut_bytes(lut_bytes: bytes):
    marker = "[ColorLUTData]\r\n".encode("utf-16le")
    start = lut_bytes.find(marker)
    if start < 0:
        raise ValueError("Could not locate [ColorLUTData] marker in LUT bytes")

    header = lut_bytes[:start].decode("utf-16le", errors="ignore").replace("\ufeff", "")
    payload = lut_bytes[start + len(marker): start + len(marker) + 65536 * 4]
    if len(payload) != 65536 * 4:
        raise ValueError(f"Unexpected LUT payload length: {len(payload)}")

    lut = np.frombuffer(payload, dtype=np.uint8).reshape(65536, 4)
    params = {}
    section = None
    for raw_line in header.splitlines():
        line = raw_line.strip("\r")
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            params.setdefault(section, {})
            continue
        if "=" in line and section is not None:
            key, value = line.split("=", 1)
            params[section][key] = value.strip().strip('"')
    return lut, params


def parse_lut_file(lut_path: Path):
    lut_bytes = Path(lut_path).read_bytes()
    lut, params = parse_lut_bytes(lut_bytes)
    return lut_bytes, lut, params


def render_channel(raw_u16, lut, params):
    active = max(
        ("Blue", "Green", "Red"),
        key=lambda color: float(params.get(color, {}).get("Contrast", 0) or 0),
    )
    section = params.get(active, {})
    shadow = float(section.get("Shadow", 0))
    highlight = float(section.get("HiLight", 65535))
    gamma = float(section.get("Gamma", 100))

    idx16 = np.clip(raw_u16.astype(np.float64) * 16.0, 0, 65535)
    z = np.clip((idx16 - shadow) / max(highlight - shadow, 1.0), 0, 1)
    z = z ** (100.0 / max(gamma, 1e-6))
    lut_idx = np.clip(z * 65535.0, 0, 65535).astype(np.uint16)

    bgr = lut[lut_idx]
    return np.stack([bgr[..., 2], bgr[..., 1], bgr[..., 0]], axis=-1).astype(np.uint8)


def discover_exported_channel_tiffs(frames_dir: Path, field_name: str):
    frames_dir = Path(frames_dir)
    if not frames_dir.exists():
        return {}

    out = {}
    pattern = re.compile(rf"^{re.escape(field_name)}_C(\d{{3}})\.tif$")
    for tif_path in sorted(frames_dir.glob("*.tif")):
        match = pattern.match(tif_path.name)
        if not match:
            continue
        channel = f"s_C{match.group(1)}"
        out[channel] = tif_path
    return dict(sorted(out.items(), key=lambda item: _channel_sort_key(item[0])))


def _infer_source_type(source) -> str:
    source = Path(source)
    return "nwb" if source.suffix == ".nwb" else "folder"


def _resolve_source_path(source, source_type=None):
    source = Path(source)
    if source_type is None:
        source_type = _infer_source_type(source)
    if source_type == "folder":
        if source.suffix == ".nwb":
            return source.parent
        return source
    if source_type == "nwb":
        if source.suffix == ".nwb":
            return source
        return source / f"{source.name}.nwb"
    raise ValueError("source_type must be 'folder' or 'nwb'")


def _maybe_float(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except Exception:
        return None


def discover_rnascope_fields_from_folder(session_dir: Path, timezone="America/Chicago"):
    session_dir = Path(session_dir)
    tz = ZoneInfo(timezone) if isinstance(timezone, str) else timezone

    fields = []
    for oif_path in sorted(session_dir.glob("*.oif")):
        field_name = oif_path.stem
        field_id = _field_id_from_name(field_name)
        oif_text = read_utf16_text(oif_path)
        oif_files_dir = oif_path.with_suffix(".oif.files")
        frames_dir = session_dir / f"{field_name}.tif.frames"

        raw_tiffs = {
            tif_path.stem: tif_path
            for tif_path in sorted(oif_files_dir.glob("s_C*.tif"))
        }
        raw_tiffs = {
            _normalize_channel_name(channel): path
            for channel, path in raw_tiffs.items()
        }
        if not raw_tiffs:
            raise ValueError(f"No raw channel TIFFs found in {oif_files_dir}")

        pty_files = {
            _normalize_channel_name(path.stem): path
            for path in sorted(oif_files_dir.glob("s_C*.pty"))
        }
        lut_files = {
            f"s_C{int(match.group(1)):03d}": path
            for path in sorted(oif_files_dir.glob("s_LUT*.lut"))
            for match in [re.search(r"s_LUT(\d+)\.lut$", path.name)]
            if match is not None
        }
        exported_channel_tiffs = discover_exported_channel_tiffs(frames_dir, field_name)

        pty_path = pty_files.get("s_C001")
        pixel_size_um = parse_pixel_size_um_from_pty(pty_path) if pty_path else None
        resolution_cm = parse_resolution_cm_from_pty(pty_path) if pty_path else None
        side, index = parse_side_and_index(oif_path.name)

        fields.append(
            {
                "session_dir": session_dir,
                "session_name": session_dir.name,
                "field_name": field_name,
                "field_id": field_id,
                "side": side,
                "index": index,
                "oif_path": oif_path,
                "oif_text": oif_text,
                "oif_files_dir": oif_files_dir,
                "frames_dir": frames_dir,
                "channel_dyes": parse_channel_dyes(oif_text),
                "raw_tiffs": dict(sorted(raw_tiffs.items(), key=lambda item: _channel_sort_key(item[0]))),
                "pty_files": dict(sorted(pty_files.items(), key=lambda item: _channel_sort_key(item[0]))),
                "lut_files": dict(sorted(lut_files.items(), key=lambda item: _channel_sort_key(item[0]))),
                "exported_channel_tiffs": exported_channel_tiffs,
                "capture_time": parse_capture_time(oif_text, tz),
                "pixel_size_um": pixel_size_um,
                "resolution_cm": resolution_cm,
            }
        )

    if not fields:
        raise ValueError(f"No .oif files found in {session_dir}")
    return fields


def _create_base_nwbfile(
    session_dir: Path,
    fields,
    timezone="America/Chicago",
    metadata=None,
    session_description=None,
    experiment_description=None,
    subject_id=None,
    subject_sex="U",
    subject_species="Mus musculus",
):
    metadata = metadata or {}
    nwb_md = metadata.get("NWBFile", {})
    subject_md = metadata.get("Subject", {})
    custom_md = metadata.get("Custom", {})

    tz = ZoneInfo(timezone) if isinstance(timezone, str) else timezone
    capture_times = [field["capture_time"] for field in fields if field["capture_time"] is not None]
    session_start_time = min(capture_times) if capture_times else datetime.now(tz)

    nwbfile = NWBFile(
        session_description=(
            nwb_md.get("session_description")
            or session_description
            or f"RNAscope imaging session converted from {Path(session_dir).name}"
        ),
        identifier=Path(session_dir).name,
        session_start_time=session_start_time,
        experiment_description=(
            nwb_md.get("experiment_description")
            or experiment_description
            or "RNAscope source conversion with raw TIFFs, PTY metadata, LUT metadata, and exported per-channel display TIFFs."
        ),
        lab=nwb_md.get("lab", "Surmeier Lab"),
        institution=nwb_md.get("institution", "Northwestern University"),
    )
    nwbfile.subject = Subject(
        subject_id=(
            subject_md.get("subject_id")
            or subject_id
            or Path(session_dir).name
        ),
        description=subject_md.get("description"),
        species=subject_md.get("species", subject_species),
        sex=subject_md.get("sex", subject_sex),
        genotype=subject_md.get("genotype"),
        age=subject_md.get("age"),
        strain=subject_md.get("strain"),
    )
    nwbfile.keywords = nwb_md.get(
        "keywords",
        [
            "RNAscope",
            "confocal microscopy",
            "raw image conversion",
            "Olympus FV10i",
        ],
    )
    if nwb_md.get("experimenter") is not None:
        nwbfile.experimenter = list(nwb_md["experimenter"])
    if nwb_md.get("protocol") is not None:
        nwbfile.protocol = nwb_md["protocol"]
    if nwb_md.get("notes") is not None:
        nwbfile.notes = nwb_md["notes"]
    if custom_md:
        nwbfile.add_scratch(
            json.dumps(custom_md, indent=2),
            name="session_metadata_custom",
            description="Custom session metadata imported from session_metadata.json",
        )
    return nwbfile


def _get_or_create_processing_module(nwbfile, name: str, description: str):
    mod = nwbfile.processing.get(name)
    if mod is None:
        mod = nwbfile.create_processing_module(name=name, description=description)
    return mod


def _get_or_create_field_metadata_table(nwbfile):
    mod = _get_or_create_processing_module(
        nwbfile,
        "rnascope_source_metadata",
        "RNAscope source metadata including field, PTY, and LUT information",
    )
    if "field_metadata" in mod.data_interfaces:
        return mod["field_metadata"]

    table = DynamicTable(
        name="field_metadata",
        description="One row per RNAscope imaging field",
    )
    table.add_column("session_name", "Session folder name")
    table.add_column("field_name", "Full field name")
    table.add_column("field_id", "Sanitized field identifier")
    table.add_column("side", "Side token from field name")
    table.add_column("field_index", "Field index token")
    table.add_column("capture_time", "Capture time in ISO format")
    table.add_column("oif_path", "Source OIF path")
    table.add_column("oif_text", "Raw OIF text")
    table.add_column("oif_files_dir", "Source .oif.files path")
    table.add_column("frames_dir", "Source .tif.frames path")
    table.add_column("pixel_size_um", "Pixel size in micrometers")
    table.add_column("resolution_cm", "Pixels per centimeter")
    table.add_column("channel_codes_json", "JSON list of raw channel codes")
    table.add_column("exported_channel_codes_json", "JSON list of exported per-channel TIFF codes")
    table.add_column("channel_dyes_json", "JSON mapping of channel code to dye name")
    table.add_column("raw_container_name", "Acquisition Images container holding raw channel TIFFs")
    table.add_column("exported_container_name", "Acquisition Images container holding exported per-channel TIFFs")
    mod.add(table)
    return table


def _get_or_create_pty_metadata_table(nwbfile):
    mod = _get_or_create_processing_module(
        nwbfile,
        "rnascope_source_metadata",
        "RNAscope source metadata including field, PTY, and LUT information",
    )
    if "channel_pty_metadata" in mod.data_interfaces:
        return mod["channel_pty_metadata"]

    table = DynamicTable(
        name="channel_pty_metadata",
        description="One row per field/channel with saved PTY metadata",
    )
    table.add_column("field_name", "Full field name")
    table.add_column("field_id", "Sanitized field identifier")
    table.add_column("channel", "Channel code such as s_C001")
    table.add_column("pty_file_name", "Original PTY file name")
    table.add_column("pty_text", "Raw PTY UTF-16 text")
    table.add_column("pty_metadata_json", "Parsed PTY metadata as JSON")
    mod.add(table)
    return table


def _get_or_create_lut_metadata_table(nwbfile):
    mod = _get_or_create_processing_module(
        nwbfile,
        "rnascope_source_metadata",
        "RNAscope source metadata including field, PTY, and LUT information",
    )
    if "channel_lut_metadata" in mod.data_interfaces:
        return mod["channel_lut_metadata"]

    table = DynamicTable(
        name="channel_lut_metadata",
        description="One row per field/channel with saved LUT metadata",
    )
    table.add_column("field_name", "Full field name")
    table.add_column("field_id", "Sanitized field identifier")
    table.add_column("channel", "Channel code such as s_C001")
    table.add_column("lut_file_name", "Original LUT file name")
    table.add_column("lut_bytes_b64", "Raw LUT bytes encoded as base64")
    table.add_column("lut_params_json", "Parsed LUT parameters as JSON")
    mod.add(table)
    return table


def add_raw_field_to_nwb(nwbfile, field_info):
    container_name = f"RawField_{field_info['field_id']}"
    if container_name in nwbfile.acquisition:
        return container_name

    images = Images(
        name=container_name,
        description=f"Raw single-channel RNAscope TIFFs for field {field_info['field_name']}",
    )

    for channel, tif_path in field_info["raw_tiffs"].items():
        array = np.asarray(Image.open(tif_path), dtype=np.uint16)
        images.add_image(
            GrayscaleImage(
                name=f"{field_info['field_id']}_{channel}_raw",
                data=array,
                description=f"Raw source TIFF for {field_info['field_name']} channel {channel}",
            )
        )

    nwbfile.add_acquisition(images)
    return container_name


def add_exported_channel_tiffs_to_nwb(nwbfile, field_info):
    if not field_info["exported_channel_tiffs"]:
        return ""

    container_name = f"ExportedField_{field_info['field_id']}"
    if container_name in nwbfile.acquisition:
        return container_name

    images = Images(
        name=container_name,
        description=f"Exported per-channel display TIFFs for field {field_info['field_name']}",
    )

    for channel, tif_path in field_info["exported_channel_tiffs"].items():
        array = np.asarray(Image.open(tif_path), dtype=np.uint8)
        images.add_image(
            RGBImage(
                name=f"{field_info['field_id']}_{channel}_exported",
                data=array,
                description=f"Exported per-channel display TIFF for {field_info['field_name']} channel {channel}",
            )
        )

    nwbfile.add_acquisition(images)
    return container_name


def add_pty_metadata_to_nwb(nwbfile, field_info):
    table = _get_or_create_pty_metadata_table(nwbfile)
    for channel, pty_path in field_info["pty_files"].items():
        pty_text, pty_metadata = parse_pty_file(pty_path)
        table.add_row(
            field_name=field_info["field_name"],
            field_id=field_info["field_id"],
            channel=channel,
            pty_file_name=pty_path.name,
            pty_text=pty_text,
            pty_metadata_json=json.dumps(pty_metadata),
        )


def add_lut_metadata_to_nwb(nwbfile, field_info):
    table = _get_or_create_lut_metadata_table(nwbfile)
    for channel, lut_path in field_info["lut_files"].items():
        lut_bytes, _, lut_params = parse_lut_file(lut_path)
        table.add_row(
            field_name=field_info["field_name"],
            field_id=field_info["field_id"],
            channel=channel,
            lut_file_name=lut_path.name,
            lut_bytes_b64=base64.b64encode(lut_bytes).decode("ascii"),
            lut_params_json=json.dumps(lut_params),
        )


def add_field_metadata_to_nwb(nwbfile, field_info):
    table = _get_or_create_field_metadata_table(nwbfile)
    raw_container_name = f"RawField_{field_info['field_id']}"
    exported_container_name = (
        f"ExportedField_{field_info['field_id']}"
        if field_info["exported_channel_tiffs"]
        else ""
    )
    table.add_row(
        session_name=field_info["session_name"],
        field_name=field_info["field_name"],
        field_id=field_info["field_id"],
        side=field_info["side"],
        field_index=field_info["index"],
        capture_time=field_info["capture_time"].isoformat() if field_info["capture_time"] else "",
        oif_path=str(field_info["oif_path"]),
        oif_text=field_info["oif_text"],
        oif_files_dir=str(field_info["oif_files_dir"]),
        frames_dir=str(field_info["frames_dir"]),
        pixel_size_um=float(field_info["pixel_size_um"]) if field_info["pixel_size_um"] is not None else np.nan,
        resolution_cm=float(field_info["resolution_cm"]) if field_info["resolution_cm"] is not None else np.nan,
        channel_codes_json=json.dumps(list(field_info["raw_tiffs"].keys())),
        exported_channel_codes_json=json.dumps(list(field_info["exported_channel_tiffs"].keys())),
        channel_dyes_json=json.dumps(field_info["channel_dyes"]),
        raw_container_name=raw_container_name,
        exported_container_name=exported_container_name,
    )


def convert_rnascope_session_to_nwb(
    session_dir: Path,
    nwb_path: Path | None = None,
    timezone: str = "America/Chicago",
    metadata: dict | None = None,
    session_description: str | None = None,
    experiment_description: str | None = None,
    subject_id: str | None = None,
    subject_sex: str = "U",
    subject_species: str = "Mus musculus",
    include_exported_channel_tiffs: bool = True,
    overwrite: bool = True,
):
    session_dir = Path(session_dir)
    nwb_path = Path(nwb_path) if nwb_path is not None else session_dir / f"{session_dir.name}.nwb"
    fields = discover_rnascope_fields_from_folder(session_dir, timezone=timezone)

    if nwb_path.exists() and not overwrite:
        raise FileExistsError(f"NWB already exists and overwrite=False: {nwb_path}")

    nwbfile = _create_base_nwbfile(
        session_dir=session_dir,
        fields=fields,
        timezone=timezone,
        metadata=metadata,
        session_description=session_description,
        experiment_description=experiment_description,
        subject_id=subject_id,
        subject_sex=subject_sex,
        subject_species=subject_species,
    )

    for field_info in fields:
        add_raw_field_to_nwb(nwbfile, field_info)
        if include_exported_channel_tiffs:
            add_exported_channel_tiffs_to_nwb(nwbfile, field_info)
        add_pty_metadata_to_nwb(nwbfile, field_info)
        add_lut_metadata_to_nwb(nwbfile, field_info)
        add_field_metadata_to_nwb(nwbfile, field_info)

    with NWBHDF5IO(str(nwb_path), "w") as io:
        io.write(nwbfile)

    return nwb_path, fields


def _load_images_container(container):
    out = {}
    for image in container.images.values():
        name = getattr(image, "name", "") or ""
        match = re.search(r"(s_C\d{3})", name)
        if match is None:
            continue
        out[match.group(1)] = np.asarray(image.data)
    return dict(sorted(out.items(), key=lambda item: _channel_sort_key(item[0])))


def _build_field_payload(
    session_name,
    session_path,
    field_name,
    field_id,
    raw_channels,
    pty_text_by_channel,
    pty_metadata_by_channel,
    lut_bytes_by_channel,
    lut_params_by_channel,
    exported_channels,
    display_mode,
    metadata=None,
):
    if display_mode == "rendered_from_raw":
        display_channels = {
            channel: render_channel(
                raw_channels[channel],
                parse_lut_bytes(lut_bytes_by_channel[channel])[0],
                lut_params_by_channel[channel],
            )
            for channel in raw_channels
            if channel in lut_bytes_by_channel and channel in lut_params_by_channel
        }
    elif display_mode == "exported_channel_tiffs":
        if not exported_channels:
            raise ValueError(f"No exported per-channel TIFFs are available for field {field_name}")
        display_channels = dict(exported_channels)
    else:
        raise ValueError(
            "display_mode must be 'rendered_from_raw' or 'exported_channel_tiffs'"
        )

    return {
        "session_name": session_name,
        "session_path": str(session_path),
        "field_name": field_name,
        "field_id": field_id,
        "raw_channels": raw_channels,
        "pty_text_by_channel": pty_text_by_channel,
        "pty_metadata_by_channel": pty_metadata_by_channel,
        "lut_bytes_by_channel": lut_bytes_by_channel,
        "lut_params_by_channel": lut_params_by_channel,
        "exported_channels": exported_channels,
        "display_mode": display_mode,
        "display_channels": display_channels,
        "metadata": metadata or {},
    }


def load_field_from_folder(session_dir: Path, field: str, display_mode="rendered_from_raw"):
    session_dir = Path(session_dir)
    field = _normalize_field(field)
    fields = discover_rnascope_fields_from_folder(session_dir)
    matches = [field_info for field_info in fields if field_info["field_name"] == field]
    if not matches:
        raise KeyError(f"Field {field!r} not found in {session_dir}")
    field_info = matches[0]

    raw_channels = {
        channel: np.asarray(Image.open(path), dtype=np.uint16)
        for channel, path in field_info["raw_tiffs"].items()
    }
    exported_channels = {
        channel: np.asarray(Image.open(path), dtype=np.uint8)
        for channel, path in field_info["exported_channel_tiffs"].items()
    }

    pty_text_by_channel = {}
    pty_metadata_by_channel = {}
    for channel, path in field_info["pty_files"].items():
        text, metadata = parse_pty_file(path)
        pty_text_by_channel[channel] = text
        pty_metadata_by_channel[channel] = metadata

    lut_bytes_by_channel = {}
    lut_params_by_channel = {}
    for channel, path in field_info["lut_files"].items():
        lut_bytes, _, lut_params = parse_lut_file(path)
        lut_bytes_by_channel[channel] = lut_bytes
        lut_params_by_channel[channel] = lut_params

    return _build_field_payload(
        session_name=field_info["session_name"],
        session_path=session_dir,
        field_name=field_info["field_name"],
        field_id=field_info["field_id"],
        raw_channels=raw_channels,
        pty_text_by_channel=pty_text_by_channel,
        pty_metadata_by_channel=pty_metadata_by_channel,
        lut_bytes_by_channel=lut_bytes_by_channel,
        lut_params_by_channel=lut_params_by_channel,
        exported_channels=exported_channels,
        display_mode=display_mode,
        metadata={
            "side": field_info["side"],
            "field_index": field_info["index"],
            "capture_time": field_info["capture_time"].isoformat() if field_info["capture_time"] else "",
            "pixel_size_um": field_info["pixel_size_um"],
            "resolution_cm": field_info["resolution_cm"],
            "channel_dyes": field_info["channel_dyes"],
        },
    )


def load_field_from_nwb(nwb_path: Path, field: str, display_mode="rendered_from_raw"):
    nwb_path = Path(nwb_path)
    field = _normalize_field(field)

    with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
        nwbfile = io.read()
        mod = nwbfile.processing.get("rnascope_source_metadata")
        if mod is None:
            raise KeyError("NWB is missing processing module `rnascope_source_metadata`")

        field_table = mod["field_metadata"].to_dataframe().reset_index(drop=True)
        row_matches = field_table[field_table["field_name"].astype(str) == field]
        if row_matches.empty:
            raise KeyError(f"Field {field!r} not found in {nwb_path}")
        row = row_matches.iloc[0]

        raw_container_name = str(row["raw_container_name"])
        if raw_container_name not in nwbfile.acquisition:
            raise KeyError(f"Missing raw acquisition container {raw_container_name!r}")
        raw_channels = {
            channel: np.asarray(array, dtype=np.uint16)
            for channel, array in _load_images_container(nwbfile.acquisition[raw_container_name]).items()
        }

        exported_channels = {}
        exported_container_name = str(row.get("exported_container_name", "") or "")
        if exported_container_name:
            if exported_container_name not in nwbfile.acquisition:
                raise KeyError(f"Missing exported acquisition container {exported_container_name!r}")
            exported_channels = {
                channel: np.asarray(array, dtype=np.uint8)
                for channel, array in _load_images_container(nwbfile.acquisition[exported_container_name]).items()
            }

        pty_table = mod["channel_pty_metadata"].to_dataframe().reset_index(drop=True)
        pty_rows = pty_table[pty_table["field_name"].astype(str) == field]
        pty_text_by_channel = {}
        pty_metadata_by_channel = {}
        for _, pty_row in pty_rows.iterrows():
            channel = str(pty_row["channel"])
            pty_text_by_channel[channel] = str(pty_row["pty_text"])
            pty_metadata_by_channel[channel] = json.loads(pty_row["pty_metadata_json"])

        lut_table = mod["channel_lut_metadata"].to_dataframe().reset_index(drop=True)
        lut_rows = lut_table[lut_table["field_name"].astype(str) == field]
        lut_bytes_by_channel = {}
        lut_params_by_channel = {}
        for _, lut_row in lut_rows.iterrows():
            channel = str(lut_row["channel"])
            lut_bytes_by_channel[channel] = base64.b64decode(str(lut_row["lut_bytes_b64"]))
            lut_params_by_channel[channel] = json.loads(lut_row["lut_params_json"])

        return _build_field_payload(
            session_name=str(row["session_name"]),
            session_path=nwb_path.parent,
            field_name=str(row["field_name"]),
            field_id=str(row["field_id"]),
            raw_channels=raw_channels,
            pty_text_by_channel=pty_text_by_channel,
            pty_metadata_by_channel=pty_metadata_by_channel,
            lut_bytes_by_channel=lut_bytes_by_channel,
            lut_params_by_channel=lut_params_by_channel,
            exported_channels=exported_channels,
            display_mode=display_mode,
            metadata={
                "side": str(row["side"]),
                "field_index": str(row["field_index"]),
                "capture_time": str(row["capture_time"]),
                "pixel_size_um": _maybe_float(row["pixel_size_um"]),
                "resolution_cm": _maybe_float(row["resolution_cm"]),
                "channel_codes": json.loads(row["channel_codes_json"]),
                "exported_channel_codes": json.loads(row["exported_channel_codes_json"]),
                "channel_dyes": json.loads(row["channel_dyes_json"]),
            },
        )


def get_rnascope_field(source, field: str, source_type=None, display_mode="rendered_from_raw"):
    if source_type is None:
        source_type = _infer_source_type(source)
    source = _resolve_source_path(source, source_type=source_type)
    if source_type == "folder":
        return load_field_from_folder(source, field=field, display_mode=display_mode)
    if source_type == "nwb":
        return load_field_from_nwb(source, field=field, display_mode=display_mode)
    raise ValueError("source_type must be 'folder' or 'nwb'")


def _display_channel_key(channel):
    return f"display_channel::{_normalize_channel_name(channel)}"


def _display_group_key(group):
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(group).strip()).strip("_").lower()
    return f"display_group::{slug}"


def _compose_display_channels(display_channels, roi_channel, count_channel):
    roi_channel = _normalize_channel_name(roi_channel)
    count_channel = _normalize_channel_name(count_channel)

    if roi_channel not in display_channels:
        raise KeyError(f"Missing display image for ROI channel {roi_channel}")
    if count_channel not in display_channels:
        raise KeyError(f"Missing display image for count channel {count_channel}")

    out = display_channels[roi_channel].astype(np.uint16).copy()
    if count_channel != roi_channel:
        out += display_channels[count_channel].astype(np.uint16)
    return np.clip(out, 0, 255).astype(np.uint8)


def build_analysis_state_from_field_data(
    field_data,
    roi_specs=None,
    roi_groups=("GB", "RB"),
    roi_channels=("s_C002", "s_C003"),
    count_channel="s_C004",
    blind=False,
):
    roi_specs = _normalize_roi_specs(
        roi_specs=roi_specs,
        roi_groups=roi_groups,
        roi_channels=roi_channels,
    )
    roi_groups = _roi_groups_from_specs(roi_specs)
    roi_channels = _roi_channels_from_specs(roi_specs)
    count_channel = _normalize_channel_name(count_channel)

    display_channels = field_data.get("display_channels", {})
    if count_channel not in field_data["raw_channels"]:
        raise KeyError(f"Missing raw count channel {count_channel} in field data")
    images = {}
    for channel, display_array in display_channels.items():
        images[_display_channel_key(channel)] = np.asarray(display_array)

    for spec in roi_specs:
        group = spec["group"]
        roi_channel = spec["channel"]
        images[_display_group_key(group)] = _compose_display_channels(
            display_channels=display_channels,
            roi_channel=roi_channel,
            count_channel=count_channel,
        )

    state = {
        "root": Path(field_data["session_path"]),
        "field": field_data["field_name"],
        "oif": None,
        "frames": None,
        "images": images,
        "offsets": {},
        "roi_specs": roi_specs,
        "roi_groups": roi_groups,
        "roi_channels": roi_channels,
        "count_channel": count_channel,
        "count_image": np.asarray(field_data["raw_channels"][count_channel], dtype=np.float32),
        "display_mode": field_data["display_mode"],
        "blind": bool(blind),
        "verify_keys": {group: _display_group_key(group) for group in roi_groups},
        "roi_sets": {group: [] for group in roi_groups},
        "roi_image_used": {group: None for group in roi_groups},
        "_ui": {},
        "field_data": field_data,
    }
    return state


def RNAscopeAnalysisStartFromFieldData(
    field_data,
    roi_specs=None,
    roi_groups=("GB", "RB"),
    roi_channels=("s_C002", "s_C003"),
    count_channel="s_C004",
    blind=False,
):
    backend = plt.get_backend().lower()
    if "ipympl" not in backend and "widget" not in backend:
        raise RuntimeError(f"Run `%matplotlib widget` first. Current backend: {backend}")

    state = build_analysis_state_from_field_data(
        field_data=field_data,
        roi_specs=roi_specs,
        roi_groups=roi_groups,
        roi_channels=roi_channels,
        count_channel=count_channel,
        blind=blind,
    )

    for spec in state["roi_specs"]:
        group = spec["group"]
        draw_key = _display_channel_key(spec["channel"]) if state["blind"] else _display_group_key(group)
        image_used_label = _compact_channel_label(spec["channel"])
        title_text = _format_start_title(group, spec["channel"])
        _draw_roi_editor(
            state,
            group,
            draw_key,
            image_used_label,
            title_text=title_text,
        )

    print("Draw ROIs in all widget figures.")
    print("Press q, escape, or enter inside each figure when finished.")
    print("Then run: results, state = RNAscopeAnalysisFinish(state)")
    return state

def RNAscopeAnalysisStartFromSource(
    source,
    field,
    source_type=None,
    display_mode="exported_channel_tiffs",
    roi_specs=None,
    roi_groups=("GB", "RB"),
    roi_channels=("s_C002", "s_C003"),
    count_channel="s_C004",
    blind=False,
):
    if source_type is None:
        source_type = _infer_source_type(source)
    field_data = get_rnascope_field(
        source,
        field=field,
        source_type=source_type,
        display_mode=display_mode,
    )
    return RNAscopeAnalysisStartFromFieldData(
        field_data=field_data,
        roi_specs=roi_specs,
        roi_groups=roi_groups,
        roi_channels=roi_channels,
        count_channel=count_channel,
        blind=blind,
    )


def _count_image_for_group(state):
    count_channel = _normalize_channel_name(state.get("count_channel", "s_C004"))
    if state.get("count_image") is not None:
        return np.asarray(state["count_image"], dtype=np.float32), count_channel
    raise KeyError("State is missing in-memory count_image")


def _peak_mask_legacy(dog, threshold_percentile=99, peak_footprint=4):
    positive = dog[dog > 0]
    threshold = np.percentile(positive, threshold_percentile) if positive.size else np.inf
    maxima = ndi.maximum_filter(dog, size=peak_footprint, mode="nearest")
    return (dog == maxima) & (dog >= threshold)


def _flood_fill_mask(image, start_r, start_c, threshold):
    from collections import deque

    h, w = image.shape
    mask = np.zeros((h, w), dtype=bool)
    queue = deque([(start_r, start_c)])
    mask[start_r, start_c] = True

    while queue:
        r, c = queue.popleft()
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w and not mask[nr, nc] and image[nr, nc] >= threshold:
                    mask[nr, nc] = True
                    queue.append((nr, nc))

    return mask


def _peak_mask_fiji(img, maxima_tolerance=200):
    maxima = ndi.maximum_filter(img, size=3, mode="nearest")
    candidates = np.argwhere((img == maxima) & (img > 0))
    peak_mask = np.zeros(img.shape, dtype=bool)

    if len(candidates) == 0:
        return peak_mask

    heights = img[candidates[:, 0], candidates[:, 1]]
    order = np.argsort(heights)[::-1]
    candidates = candidates[order]
    heights = heights[order]
    claimed = np.zeros(img.shape, dtype=bool)

    for i in range(len(candidates)):
        r, c = candidates[i]
        if claimed[r, c]:
            continue
        peak_mask[r, c] = True
        region = _flood_fill_mask(img, r, c, heights[i] - maxima_tolerance)
        claimed |= region

    return peak_mask


def _detect_centers(
    state,
    sigma_small=1.0,
    sigma_large=2.8,
    detection_method="DoG",
    dog_mode="tolerance",
    threshold_percentile=99,
    peak_footprint=4,
    maxima_tolerance=200,
):
    raw, src_label = _count_image_for_group(state)
    method = detection_method.strip().lower()

    if method == "dog":
        g1 = ndi.gaussian_filter(raw, sigma_small)
        g2 = ndi.gaussian_filter(raw, sigma_large)
        dog = g1 - g2
        dog[dog < 0] = 0

        mode = dog_mode.strip().lower()
        if mode == "legacy":
            peak_mask = _peak_mask_legacy(
                dog,
                threshold_percentile=threshold_percentile,
                peak_footprint=peak_footprint,
            )
        elif mode == "tolerance":
            peak_mask = _peak_mask_fiji(
                dog,
                maxima_tolerance=maxima_tolerance,
            )
        else:
            raise ValueError("dog_mode must be 'legacy' or 'tolerance'")
    elif method == "fiji":
        peak_mask = _peak_mask_fiji(
            raw.astype(np.float32),
            maxima_tolerance=maxima_tolerance,
        )
    else:
        raise ValueError("detection_method must be 'DoG' or 'fiji'")

    labels, n_peaks = ndi.label(peak_mask)
    centers = (
        np.array(ndi.center_of_mass(peak_mask, labels, np.arange(1, n_peaks + 1)))
        if n_peaks
        else np.empty((0, 2))
    )

    height, width = raw.shape
    rr = np.clip(np.rint(centers[:, 0]).astype(int), 0, height - 1) if n_peaks else np.array([], dtype=int)
    cc = np.clip(np.rint(centers[:, 1]).astype(int), 0, width - 1) if n_peaks else np.array([], dtype=int)

    return src_label, centers, rr, cc, height, width


def _centers_in_group_rois(state, group, centers, rr, cc, height, width):
    if len(centers) == 0 or not state["roi_sets"].get(group):
        return np.empty((0, 2))

    center_xy = np.c_[cc, rr]
    keep = np.zeros(len(centers), dtype=bool)

    for verts in state["roi_sets"][group]:
        path = MplPath(np.asarray(verts, dtype=float))
        keep |= path.contains_points(center_xy)

    return centers[keep]


def _roi_channel_for_group(state, group):
    roi_specs = state.get("roi_specs")
    if not roi_specs:
        return None

    for spec in _normalize_roi_specs(roi_specs=roi_specs):
        if spec["group"] == group:
            return spec["channel"]
    return None


def _image_used_label_for_group(state, group):
    channel = _roi_channel_for_group(state, group)
    if channel is not None:
        return _compact_channel_label(channel)

    image_used = str(state.get("roi_image_used", {}).get(group, "")).strip()
    match = re.search(r"[Cc](\d{3})", image_used)
    return f"C{match.group(1)}" if match else image_used


def _format_finish_title(group, roi_channel, count_channel):
    roi_label = _compact_channel_label(roi_channel)
    count_label = _compact_channel_label(count_channel)
    if roi_label == count_label:
        return f"{group} ({roi_label})"
    return f"{group} ({roi_label} + {count_label})"


def measure_rois(
    state,
    sigma_small=1.0,
    sigma_large=2.8,
    detection_method="DoG",
    dog_mode="tolerance",
    threshold_percentile=99,
    peak_footprint=4,
    maxima_tolerance=200,
):
    _, centers, rr, cc, height, width = _detect_centers(
        state,
        sigma_small=sigma_small,
        sigma_large=sigma_large,
        detection_method=detection_method,
        dog_mode=dog_mode,
        threshold_percentile=threshold_percentile,
        peak_footprint=peak_footprint,
        maxima_tolerance=maxima_tolerance,
    )

    center_xy = np.c_[cc, rr] if len(centers) else np.empty((0, 2))
    yy, xx = np.mgrid[:height, :width]
    all_xy = np.c_[xx.ravel(), yy.ravel()]
    groups = tuple(state.get("roi_groups", tuple(state["roi_sets"].keys())))

    rows = []
    for group in groups:
        for idx, verts in enumerate(state["roi_sets"].get(group, []), start=1):
            path = MplPath(np.asarray(verts, dtype=float))
            dot_count = int(path.contains_points(center_xy).sum()) if len(centers) else 0
            pixel_count = int(path.contains_points(all_xy).sum())
            rows.append(
                {
                    "field": state["field"],
                    "group": group,
                    "image used": _image_used_label_for_group(state, group),
                    "roi": idx,
                    "pixel count": pixel_count,
                    "count": dot_count,
                }
            )

    return pd.DataFrame(rows)


def RNAscopeAnalysisFinish(
    state,
    sigma_small=1.0,
    sigma_large=2.8,
    detection_method="DoG",
    dog_mode="tolerance",
    threshold_percentile=99,
    peak_footprint=4,
    maxima_tolerance=170,
    show_verify=True,
    show_detected=False,
    verify_circle_radius=2,
    verify_circle_color="white",
    verify_image="group",
):
    for ui in state.get("_ui", {}).values():
        fig = ui.get("fig")
        if fig is not None and plt.fignum_exists(fig.number):
            plt.close(fig)

    groups = tuple(state.get("roi_groups", tuple(state["roi_sets"].keys())))
    results = measure_rois(
        state,
        sigma_small=sigma_small,
        sigma_large=sigma_large,
        detection_method=detection_method,
        dog_mode=dog_mode,
        threshold_percentile=threshold_percentile,
        peak_footprint=peak_footprint,
        maxima_tolerance=maxima_tolerance,
    )

    if show_verify:
        detected_by_group = {group: np.empty((0, 2)) for group in groups}
        if show_detected:
            _, centers, rr, cc, height, width = _detect_centers(
                state,
                sigma_small=sigma_small,
                sigma_large=sigma_large,
                detection_method=detection_method,
                dog_mode=dog_mode,
                threshold_percentile=threshold_percentile,
                peak_footprint=peak_footprint,
                maxima_tolerance=maxima_tolerance,
            )
            for group in groups:
                detected_by_group[group] = _centers_in_group_rois(
                    state, group, centers, rr, cc, height, width
                )

        count_channel = _normalize_channel_name(state.get("count_channel", "s_C004"))
        fig, axes = plt.subplots(1, len(groups), figsize=(6 * len(groups), 6))
        axes = np.atleast_1d(axes)

        for idx, group in enumerate(groups):
            if verify_image == "count":
                key = _display_channel_key(count_channel)
            elif verify_image == "group":
                key = state["verify_keys"][group]
            else:
                raise ValueError("verify_image must be 'group' or 'count'")

            axes[idx].imshow(state["images"][key])
            axes[idx].axis("off")
            roi_channel = _roi_channel_for_group(state, group)
            axes[idx].set_title(_format_finish_title(group, roi_channel, count_channel))

            for verts in state["roi_sets"].get(group, []):
                verts = np.asarray(verts, dtype=float)
                patch = Polygon(
                    verts,
                    closed=True,
                    fill=False,
                    edgecolor="white",
                    linewidth=1,
                    linestyle=(0, (1.2, 2.0)),
                )
                axes[idx].add_patch(patch)

            if show_detected:
                for center in detected_by_group[group]:
                    axes[idx].add_patch(
                        Circle(
                            (center[1], center[0]),
                            radius=verify_circle_radius,
                            fill=False,
                            edgecolor=verify_circle_color,
                            linewidth=0.5,
                        )
                    )

        plt.tight_layout()
        plt.show()

    if results.empty:
        print("[info] no ROIs selected; returning empty results")

    return results, state


def save_rnascope_field_analysis(state, results, analysis_params=None, out_dir=None):
    if hasattr(results, "to_dict"):
        rows = results.to_dict(orient="records")
    else:
        rows = list(results)

    row_lookup = {
        (str(row["group"]), int(row["roi"])): row
        for row in rows
    }

    payload = {
        "schema_version": 1,
        "session": Path(state["root"]).name,
        "field": state["field"],
        "display_mode": state.get("display_mode"),
        "count_channel": state.get("count_channel", "s_C004"),
        "roi_specs": [dict(spec) for spec in state.get("roi_specs", ())],
        "blind": bool(state.get("blind", False)),
        "analysis_timestamp": datetime.now().isoformat(timespec="seconds"),
        "analysis_params": dict(analysis_params or {}),
        "groups": {group: [] for group in state.get("roi_groups", ())},
    }

    for group in state.get("roi_groups", ()):
        for roi_index, verts in enumerate(state["roi_sets"].get(group, []), start=1):
            row = row_lookup.get((group, roi_index), {})
            payload["groups"][group].append(
                {
                    "roi_index": roi_index,
                    "vertices_xy": np.asarray(verts, dtype=float).tolist(),
                    "pixel_count": int(row.get("pixel count", 0)),
                    "dot_count": int(row.get("count", 0)),
                    "image_used": row.get("image used", state.get("roi_image_used", {}).get(group)),
                }
            )

    out_dir = Path(out_dir) if out_dir is not None else Path(state["root"]) / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{state['field']}.roi_analysis.json"
    out_path.write_text(json.dumps(payload, indent=2))
    return out_path


def load_rnascope_field_analysis(path):
    return json.loads(Path(path).read_text())


def load_session_analysis_jsons(session_root):
    session_root = Path(session_root)
    analysis_dir = session_root / "analysis"
    json_paths = sorted(analysis_dir.glob("*.roi_analysis.json"))

    rows = []
    for json_path in json_paths:
        payload = load_rnascope_field_analysis(json_path)

        for group, rois in payload.get("groups", {}).items():
            for roi in rois:
                rows.append(
                    {
                        "field": str(payload["field"]),
                        "group": str(group),
                        "image used": str(roi.get("image_used", "")),
                        "roi": int(roi["roi_index"]),
                        "pixel count": int(roi.get("pixel_count", 0)),
                        "count": int(roi.get("dot_count", 0)),
                        "count_channel": str(payload.get("count_channel", "")),
                        "display_mode": str(payload.get("display_mode", "")),
                        "blind": bool(payload.get("blind", False)),
                        "analysis_timestamp": str(payload.get("analysis_timestamp", "")),
                        "analysis_json": json_path.name,
                    }
                )

    return pd.DataFrame(rows)


def apply_saved_analysis_to_state(state, analysis_payload):
    if str(analysis_payload.get("field", "")) != str(state["field"]):
        raise ValueError(
            f"Analysis payload field {analysis_payload.get('field')!r} does not match state field {state['field']!r}"
        )

    if analysis_payload.get("count_channel") is not None:
        state["count_channel"] = _normalize_channel_name(analysis_payload["count_channel"])
        state["count_image"] = np.asarray(
            state["field_data"]["raw_channels"][state["count_channel"]],
            dtype=np.float32,
        )

    if analysis_payload.get("roi_specs"):
        roi_specs = _normalize_roi_specs(analysis_payload["roi_specs"])
        state["roi_specs"] = roi_specs
        state["roi_groups"] = _roi_groups_from_specs(roi_specs)
        state["roi_channels"] = _roi_channels_from_specs(roi_specs)

    if analysis_payload.get("blind") is not None:
        state["blind"] = bool(analysis_payload["blind"])

    state["roi_sets"] = {}
    state["roi_image_used"] = {}
    for group in state.get("roi_groups", ()):
        group_rois = analysis_payload.get("groups", {}).get(group, [])
        state["roi_sets"][group] = [
            np.asarray(roi["vertices_xy"], dtype=float)
            for roi in group_rois
        ]
        state["roi_image_used"][group] = (
            group_rois[0].get("image_used") if group_rois else None
        )

    state["analysis_payload"] = analysis_payload
    return state


def reconstruct_state_from_saved_analysis(source, analysis_path, source_type=None, display_mode=None):
    analysis_payload = load_rnascope_field_analysis(analysis_path)
    field = analysis_payload["field"]
    display_mode = display_mode or analysis_payload.get("display_mode", "rendered_from_raw")

    field_data = get_rnascope_field(
        source,
        field=field,
        source_type=source_type,
        display_mode=display_mode,
    )

    state = build_analysis_state_from_field_data(
        field_data=field_data,
        roi_specs=analysis_payload.get("roi_specs"),
        count_channel=analysis_payload.get("count_channel", "s_C004"),
        blind=bool(analysis_payload.get("blind", False)),
    )
    state = apply_saved_analysis_to_state(state, analysis_payload)
    return state, analysis_payload

def apply_session_metadata_to_nwb(nwb_path, metadata):
    nwb_path = Path(nwb_path)

    with NWBHDF5IO(str(nwb_path), "r+", load_namespaces=True) as io:
        nwbfile = io.read()

        # Top-level NWBFile metadata
        for key, value in metadata.get("NWBFile", {}).items():
            if value is None:
                continue
            nwbfile.fields[key] = list(value) if key in ("experimenter", "keywords") else value

        # Subject metadata
        subject_md = metadata.get("Subject", {})
        if subject_md:
            if nwbfile.subject is None:
                nwbfile.subject = Subject(
                    subject_id=subject_md.get("subject_id") or nwbfile.identifier,
                    species=subject_md.get("species", "Mus musculus"),
                    sex=subject_md.get("sex", "U"),
                )

            for key, value in subject_md.items():
                if value is None:
                    continue
                nwbfile.subject.fields[key] = value

        # Custom metadata
        custom_md = metadata.get("Custom", {})
        if custom_md:
            custom_json = json.dumps(custom_md, indent=2)
            if "session_metadata_custom" in nwbfile.scratch:
                del nwbfile.scratch["session_metadata_custom"]
            nwbfile.add_scratch(
                custom_json,
                name="session_metadata_custom",
                description="Custom session metadata imported from session metadata",
            )

        nwbfile.set_modified()
        io.write(nwbfile)

def pretty(value):
    if value is None:
        return "_null_"
    if isinstance(value, (list, tuple)):
        if not value:
            return "_[]_"
        return "\n".join(f"- {v}" for v in value)
    if isinstance(value, dict):
        if not value:
            return "_{}_"
        lines = []
        for k, v in value.items():
            lines.append(f"**{k}:** {pretty(v)}")
        return "\n".join(lines)
    return str(value)


def show_block(title, mapping):
    display(Markdown(f"## {title}"))
    if not mapping:
        display(Markdown("_No entries_"))
        return
    parts = []
    for k, v in mapping.items():
        parts.append(f"### {k}\n{pretty(v)}")
    display(Markdown("\n\n".join(parts)))


def sort_experimenter_style(df):
    df = df.copy()

    if "group" in df.columns and "cell_type" not in df.columns:
        df = df.rename(columns={"group": "cell_type"})
    if "roi" in df.columns and "replicate" not in df.columns:
        df = df.rename(columns={"roi": "replicate"})

    pattern = re.compile(
        r"^(?P<slice_id>.+)\.(?P<hemisphere>UL|L)_60x\.?(?P<field_index>\d+)$"
    )

    parsed = df["field"].astype(str).str.extract(pattern)
    df["slice_id"] = parsed["slice_id"]
    df["hemisphere"] = parsed["hemisphere"]
    df["field_index"] = parsed["field_index"].astype(int)
    df["condition"] = df["hemisphere"].map({"UL": "Intact", "L": "Lesioned"})

    df["cell_type"] = pd.Categorical(
        df["cell_type"],
        categories=["NDNF+", "TH+"],
        ordered=True,
    )
    df["condition"] = pd.Categorical(
        df["condition"],
        categories=["Intact", "Lesioned"],
        ordered=True,
    )
    df["hemisphere"] = pd.Categorical(
        df["hemisphere"],
        categories=["UL", "L"],
        ordered=True,
    )
    df["replicate"] = df["replicate"].astype(int)

    return df.sort_values(
        [
            "session",
            "cell_type",
            "condition",
            "slice_id",
            "hemisphere",
            "field_index",
            "field",
            "replicate",
        ],
        kind="stable",
    ).reset_index(drop=True)

def parse_field_metadata(field_name):
    match = FIELD_RE.match(str(field_name))
    if match is None:
        raise ValueError(f"Could not parse field name: {field_name}")

    hemisphere = match.group("hemisphere")
    return {
        "slice_id": match.group("slice_id"),
        "hemisphere": hemisphere,
        "field_index": int(match.group("field_index")),
        "condition": "Intact" if hemisphere == "UL" else "Lesioned",
    }

def counts_from_roi_jsons(json_paths, session_group=None):
    rows = []

    if json_paths is None:
        json_paths = []

    for json_path in sorted(map(Path, json_paths)):
        payload = json.loads(json_path.read_text())
        field_name = str(payload["field"])
        field_meta = parse_field_metadata(field_name)
        session = str(payload.get("session") or session_group or json_path.parents[1].name)

        for group, rois in payload.get("groups", {}).items():
            for roi in rois:
                rows.append(
                    {
                        "condition": field_meta["condition"],
                        "cell_type": str(group),
                        "slice_id": field_meta["slice_id"],
                        "field": field_name,
                        "hemisphere": field_meta["hemisphere"],
                        "field_index": field_meta["field_index"],
                        "replicate": int(roi["roi_index"]),
                        "count": int(roi.get("dot_count", 0)),
                        "session": session,
                        "session_group": session_group or session,
                        "count_channel": str(payload.get("count_channel", "")),
                        "analysis_json": json_path.name,
                    }
                )

    columns = [
        "condition",
        "cell_type",
        "slice_id",
        "field",
        "hemisphere",
        "field_index",
        "replicate",
        "count",
        "session",
        "session_group",
        "count_channel",
        "analysis_json",
    ]

    if not rows:
        return pd.DataFrame(columns=columns)

    return (
        pd.DataFrame(rows, columns=columns)
        .sort_values(
            ["session", "condition", "cell_type", "hemisphere", "field_index", "field", "replicate"],
            kind="stable",
        )
        .reset_index(drop=True)
    )
