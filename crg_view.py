#!/usr/bin/env python3
"""
crg_view.py - lightweight reader and 3D viewer for ASAM OpenCRG (.crg) files.

Dependencies: numpy, matplotlib  (pip install numpy matplotlib)

Supported: all four data formats (LRFI, LDFI, KRBI, KDBI and friends, i.e.
long/compact, single/double, ASCII/binary), long sections defined by index or
by explicit v position, reference line given by heading (phi), by explicit
x/y channels, or straight, plus slope and banking (channel or constant).
Not supported: $ROAD_CRG_FILE includes and $ROAD_CRG_MODS / _OPTS modifiers
(they are ignored, so the grid is shown as stored in the file).

The reconstruction follows the reference C loader of the OpenCRG repository
(https://github.com/asam-ev/OpenCRG, c-api/baselib/src/crgLoader.c).

Usage:
    python crg_view.py road.crg
    python crg_view.py road.crg --zscale 10 --max-points 400
    python crg_view.py road.crg --mode uv --save road.png
    python crg_view.py road.crg --export road.ply          # point cloud for CloudCompare
    python crg_view.py road.crg --export road.ply --mesh   # triangulated surface
"""

import argparse
import re
import sys

import numpy as np

_NUMERIC_CHARS = set("0123456789+-.eEdD ")


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def _split_header(raw):
    """Return (header lines, byte offset of the data section)."""
    pos = 0
    lines = []
    while pos < len(raw):
        end = raw.find(b"\n", pos)
        end = len(raw) if end < 0 else end
        line = raw[pos:end].rstrip(b"\r").decode("latin-1")
        pos = end + 1
        if line.startswith("$$$$"):
            return lines, pos
        lines.append(line)
    return lines, None


def _strip_comment(line):
    return line.split("!", 1)[0].strip()


def _parse_header(lines):
    """Return (road parameters dict, data format string, channel name list)."""
    params, channels, fmt = {}, [], None
    section = None
    for line in lines:
        if line.startswith("$") and not line.startswith("$!"):
            tag = _strip_comment(line)[1:].strip().upper()
            section = tag or None
            if section == "ROAD_CRG_FILE":
                raise NotImplementedError("$ROAD_CRG_FILE references are not supported")
            continue
        if section == "CT" or line.startswith(("*", "%", "$!")):
            continue
        text = _strip_comment(line)
        if not text:
            continue
        if section == "ROAD_CRG" and "=" in text:
            key, val = text.split("=", 1)
            params[key.strip().lower()] = float(val.replace("d", "e").replace("D", "e"))
        elif section == "KD_DEFINITION":
            if text.startswith("#:"):
                fmt = text[2:6].upper()
            elif text.upper().startswith("D:"):
                channels.append(text[2:].split(",", 1)[0].strip().lower())
    if fmt is None:
        raise ValueError("no data format ('#:LRFI' etc.) in $KD_DEFINITION")
    if not channels:
        raise ValueError("no data channels ('D:...') in $KD_DEFINITION")
    return params, fmt, channels


def _ascii_value(field):
    if set(field) - _NUMERIC_CHARS:
        return 0.0 if field == "**unused**" else np.nan
    field = field.strip().replace("d", "e").replace("D", "e")
    return float(field) if field else 0.0


def _read_records(data, fmt, nch):
    """Decode the data section into an (n_records, nch) float array."""
    long_fmt = "L" in fmt
    double = "D" in fmt
    if "F" in fmt:  # ASCII, fixed width fields
        width = 20 if double else 10
        text = data.decode("latin-1")
        if long_fmt:  # each record spans ceil(nch*width/80) lines
            lines_per_rec = -(-nch * width // 80)
            lines = [ln.rstrip("\r") for ln in text.split("\n")]
            lines = [ln for ln in lines if ln.strip()]
            fields = []
            for r in range(len(lines) // lines_per_rec):
                rec = []
                for ln in lines[r * lines_per_rec:(r + 1) * lines_per_rec]:
                    rec += [ln[i:i + width] for i in range(0, len(ln) - width + 1, width)]
                fields.append(rec[:nch])
            fields = [f for f in fields if len(f) == nch]
        else:  # compact: one continuous stream of fields
            stream = text.replace("\r", "").replace("\n", "")
            nrec = len(stream) // (width * nch)
            fields = [[stream[(r * nch + c) * width:(r * nch + c + 1) * width]
                       for c in range(nch)] for r in range(nrec)]
        return np.array([[_ascii_value(f) for f in rec] for rec in fields], dtype=float)

    # binary, big endian IEEE
    size = 8 if double else 4
    rec_bytes = nch * size
    if long_fmt and rec_bytes % 80:
        rec_bytes = (rec_bytes // 80 + 1) * 80
    nrec = len(data) // rec_bytes
    dtype = np.dtype(">f8" if double else ">f4")
    buf = np.frombuffer(data[:nrec * rec_bytes], dtype=np.uint8).reshape(nrec, rec_bytes)
    return buf[:, :nch * size].copy().view(dtype).astype(float)


def _forward_backward(start, inc, end=None):
    """Integrate increments like the C loader (Euler, optionally end-weighted)."""
    n = len(inc)
    out = np.empty(n)
    out[0] = start
    out[1:] = start + np.cumsum(inc[1:])
    if end is None or n < 2:
        return out
    back = np.empty(n)
    back[-1] = end
    back[:-1] = end - np.cumsum(inc[1:][::-1])[::-1]
    for i in range(n - 1):  # weighted forward integration towards end point
        f = (i + 1) / (n - 1)
        out[i + 1] = (1 - f) * (out[i] + inc[i + 1]) + f * back[i + 1]
    return out


def read_crg(path):
    """Read an OpenCRG file and return a dict with grid and 3D coordinates.

    Keys: u (nu), v (nv), z (nu x nv grid elevation as stored), x, y,
    zref, bank (nu each, reference line), X, Y, Z (nu x nv, road surface).
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    lines, offset = _split_header(raw)
    p, fmt, channels = _parse_header(lines)
    if offset is None:
        raise ValueError("no data section ('$$$$...') found")
    rec = _read_records(raw[offset:], fmt, len(channels))

    # expected record count from the header (binary may carry padding bytes)
    if "reference_line_end_u" in p and "reference_line_increment" in p:
        nexp = int((p["reference_line_end_u"] - p.get("reference_line_start_u", 0.0))
                   / p["reference_line_increment"] + 0.5) + 1
        rec = rec[:nexp]
    if len(rec) < 2:
        raise ValueError("less than two data records found")

    # --- classify channels
    idx = {}
    sec_pos, sec_idx, sec_num = [], [], []
    for i, name in enumerate(channels):
        if name.startswith("long section"):
            m = re.search(r"v\s*=\s*([-+0-9.eEdD]+)", name)
            if m:
                sec_pos.append(float(m.group(1).replace("d", "e").replace("D", "e")))
                sec_idx.append(i)
            else:
                m = re.search(r"(\d+)", name)
                sec_num.append(int(m.group(1)) if m else len(sec_num) + 1)
                sec_idx.append(i)
        elif name.startswith("reference line"):
            word = name[len("reference line"):].strip().split()[0]
            idx[word] = i
    if sec_pos and sec_num:
        raise ValueError("mixed long section definitions (by position and by index)")
    if len(sec_idx) < 2:
        raise ValueError("at least two long sections required")

    if sec_pos:
        order = np.argsort(sec_pos)
        v = np.asarray(sec_pos)[order]
    else:
        order = np.argsort(sec_num)
        n = len(sec_num)
        dv = p.get("long_section_v_increment")
        if dv is None and "long_section_v_right" in p and "long_section_v_left" in p:
            dv = (p["long_section_v_left"] - p["long_section_v_right"]) / (n - 1)
        dv = 0.01 if dv is None else dv
        v = p.get("long_section_v_right", -0.5 * (n - 1) * dv) + dv * np.arange(n)
    z = rec[:, np.asarray(sec_idx)[order]]

    # NaN handling like the C default (grid_nan_mode = keep last): copy the
    # nearest valid value along each cross section
    if np.isnan(z).any():
        for row in z:
            ok = np.flatnonzero(~np.isnan(row))
            if ok.size:
                row[:] = row[ok[np.clip(np.searchsorted(ok, np.arange(row.size), side="right") - 1, 0, None)]]
                row[:ok[0]] = row[ok[0]]

    nu = len(rec)
    # --- reference line (x, y, phi) and u
    if "x" in idx and "y" in idx:
        x, y = rec[:, idx["x"]], rec[:, idx["y"]]
        ds = np.hypot(np.diff(x), np.diff(y))
        u = p.get("reference_line_start_u", 0.0) + np.concatenate(([0.0], np.cumsum(ds)))
    else:
        if "u" in idx:
            u = rec[:, idx["u"]]
            du = np.diff(u).max()
        else:
            du = p.get("reference_line_increment", 0.01)
            u = p.get("reference_line_start_u", 0.0) + du * np.arange(nu)
        x0 = p.get("reference_line_start_x", 0.0)
        y0 = p.get("reference_line_start_y", 0.0)
        phi0 = p.get("reference_line_start_phi", 0.0)
        if "phi" in idx:
            phi = rec[:, idx["phi"]].copy()
            phi[0] = phi0
            has_end = "reference_line_end_x" in p and "reference_line_end_y" in p
            x = _forward_backward(x0, du * np.cos(phi), p.get("reference_line_end_x") if has_end else None)
            y = _forward_backward(y0, du * np.sin(phi), p.get("reference_line_end_y") if has_end else None)
        else:
            x = x0 + (u - u[0]) * np.cos(phi0)
            y = y0 + (u - u[0]) * np.sin(phi0)

    # --- reference line elevation from slope, banking
    du_mean = (u[-1] - u[0]) / (nu - 1)
    if "slope" in idx:
        slope = rec[:, idx["slope"]]
    else:
        slope = np.full(nu, p.get("reference_line_start_s", 0.0))
    zref = _forward_backward(p.get("reference_line_start_z", 0.0), slope * du_mean,
                             p.get("reference_line_end_z"))
    if "banking" in idx:
        bank = rec[:, idx["banking"]]
    else:
        bank = np.full(nu, p.get("reference_line_start_b", 0.0))

    # --- road surface in (x, y, z); v is positive to the left of the reference line.
    # Like the C library, offset each node along the normal of the chord between
    # its neighbours, mitred against the normal of the following segment.
    nx, ny = -np.gradient(y), np.gradient(x)
    nlen = np.hypot(nx, ny)
    nx, ny = nx / nlen, ny / nlen
    sx, sy = -np.diff(y), np.diff(x)
    slen = np.hypot(sx, sy)
    sx, sy = np.append(sx / slen, sx[-1] / slen[-1]), np.append(sy / slen, sy[-1] / slen[-1])
    dot = nx * sx + ny * sy
    dot = np.where(np.abs(dot) > 1e-10, dot, 1.0)
    nx, ny = nx / dot, ny / dot
    X = x[:, None] + v[None, :] * nx[:, None]
    Y = y[:, None] + v[None, :] * ny[:, None]
    Z = z + zref[:, None] + bank[:, None] * v[None, :]

    return dict(u=u, v=v, z=z, x=x, y=y, zref=zref, bank=bank,
                X=X, Y=Y, Z=Z, params=p, format=fmt)


# --------------------------------------------------------------------------
# export
# --------------------------------------------------------------------------

def export_cloud(crg, path, step=1, mesh=False):
    """Write the road surface as a point cloud (or mesh) for CloudCompare etc.

    .ply : binary PLY, double x/y/z plus u and v as scalar fields; with
           mesh=True also triangle faces
    .xyz, .txt, .csv : ASCII columns x y z u v
    step decimates the grid (every step-th point in u and v).
    Returns the number of points written.
    """
    sl = (slice(None, None, step), slice(None, None, step))
    X, Y, Z = crg["X"][sl], crg["Y"][sl], crg["Z"][sl]
    U, V = np.meshgrid(crg["u"][::step], crg["v"][::step], indexing="ij")
    ok = np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
    n = int(ok.sum())
    ext = path.lower().rsplit(".", 1)[-1]

    if ext in ("xyz", "txt", "csv"):
        if mesh:
            raise ValueError("--mesh needs a .ply file")
        sep = "," if ext == "csv" else " "
        cols = np.column_stack([X[ok], Y[ok], Z[ok], U[ok], V[ok]])
        np.savetxt(path, cols, fmt="%.6f", delimiter=sep,
                   header=sep.join("xyzuv"), comments="//" if ext != "csv" else "")
        return n
    if ext != "ply":
        raise ValueError(f"unsupported export format '.{ext}' (use .ply, .xyz, .txt or .csv)")

    verts = np.empty(n, dtype=[("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("u", "<f4"), ("v", "<f4")])
    verts["x"], verts["y"], verts["z"], verts["u"], verts["v"] = X[ok], Y[ok], Z[ok], U[ok], V[ok]

    faces = None
    if mesh:  # two triangles per grid cell whose four corners are valid
        vid = np.full(ok.shape, -1, dtype=np.int64)
        vid[ok] = np.arange(n)
        a, b, c, d = vid[:-1, :-1], vid[1:, :-1], vid[1:, 1:], vid[:-1, 1:]
        cell = (a >= 0) & (b >= 0) & (c >= 0) & (d >= 0)
        tri = np.concatenate([np.stack([a[cell], b[cell], c[cell]], -1),
                              np.stack([a[cell], c[cell], d[cell]], -1)])
        faces = np.empty(len(tri), dtype=[("n", "u1"), ("i", "<i4", (3,))])
        faces["n"], faces["i"] = 3, tri

    header = ["ply", "format binary_little_endian 1.0",
              "comment exported from OpenCRG by crg_view.py",
              f"element vertex {n}",
              "property double x", "property double y", "property double z",
              "property float u", "property float v"]
    if faces is not None:
        header += [f"element face {len(faces)}", "property list uchar int vertex_indices"]
    header.append("end_header")
    with open(path, "wb") as fh:
        fh.write(("\n".join(header) + "\n").encode("ascii"))
        verts.tofile(fh)
        if faces is not None:
            faces.tofile(fh)
    return n


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def plot_crg(crg, title="", mode="xyz", zscale=1.0, max_points=300, cmap="viridis", save=None):
    import matplotlib
    if save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # decimate so that the mesh stays interactive
    su = max(1, int(np.ceil(len(crg["u"]) / max_points)))
    sv = max(1, int(np.ceil(len(crg["v"]) / max_points)))
    sl = (slice(None, None, su), slice(None, None, sv))

    if mode == "uv":  # elevation grid as stored, unrolled along u
        U, V = np.meshgrid(crg["u"], crg["v"], indexing="ij")
        X, Y, Z = U[sl], V[sl], crg["z"][sl]
        labels = ("u [m]", "v [m]", "z [m]")
    else:
        X, Y, Z = crg["X"][sl], crg["Y"][sl], crg["Z"][sl]
        labels = ("x [m]", "y [m]", "z [m]")

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(projection="3d")
    surf = ax.plot_surface(X, Y, Z * zscale, cmap=cmap, linewidth=0,
                           antialiased=False, rcount=Z.shape[0], ccount=Z.shape[1])
    if mode != "uv":
        ax.plot(crg["x"][::su], crg["y"][::su], crg["zref"][::su] * zscale, "r-", lw=1, label="reference line")
        ax.legend(loc="upper left")
    cb = fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1)
    cb.set_label("z [m]" + (f" (x{zscale:g})" if zscale != 1 else ""))
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2] + (f" x{zscale:g}" if zscale != 1 else ""))

    # true proportions (z exaggerated by zscale only), but keep each axis at
    # least 2 % of the largest one so that thin or flat roads stay visible
    ext = np.array([np.ptp(X), np.ptp(Y), np.nanmax(Z) * zscale - np.nanmin(Z) * zscale])
    ax.set_box_aspect(np.maximum(ext, 0.02 * ext.max()) if ext.max() > 0 else None)
    ax.set_title(f"{title}  [{crg['format']}, {len(crg['u'])} x {len(crg['v'])} grid"
                 + (f", shown every {su}x{sv}" if su > 1 or sv > 1 else "") + "]")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150)
        print(f"saved {save}")
    else:
        plt.show()


def main(argv=None):
    ap = argparse.ArgumentParser(description="3D viewer for ASAM OpenCRG files")
    ap.add_argument("file", help="OpenCRG file (.crg)")
    ap.add_argument("--mode", choices=("xyz", "uv"), default="xyz",
                    help="xyz: road in world coordinates (default); uv: unrolled elevation grid")
    ap.add_argument("--zscale", type=float, default=1.0, help="vertical exaggeration factor")
    ap.add_argument("--max-points", type=int, default=300,
                    help="max mesh points per direction (decimates large grids)")
    ap.add_argument("--cmap", default="viridis", help="matplotlib colormap")
    ap.add_argument("--save", metavar="PNG", help="write an image instead of opening a window")
    ap.add_argument("--export", metavar="FILE",
                    help="export the surface as a point cloud (.ply, .xyz, .txt, .csv) and exit")
    ap.add_argument("--step", type=int, default=1,
                    help="export every STEP-th grid point in u and v (default 1 = full resolution)")
    ap.add_argument("--mesh", action="store_true", help="with --export to .ply: also write triangle faces")
    ap.add_argument("--show", action="store_true", help="with --export: also open the 3D plot")
    args = ap.parse_args(argv)

    crg = read_crg(args.file)
    u, v, Z = crg["u"], crg["v"], crg["Z"]
    print(f"{args.file}: format {crg['format']}, {len(u)} x {len(v)} points, "
          f"u = {u[0]:.3f} .. {u[-1]:.3f} m, v = {v[0]:.3f} .. {v[-1]:.3f} m, "
          f"z = {np.nanmin(Z):.4f} .. {np.nanmax(Z):.4f} m")
    if args.export:
        n = export_cloud(crg, args.export, step=max(1, args.step), mesh=args.mesh)
        print(f"wrote {n} points to {args.export}")
        if not args.show:
            return 0
    plot_crg(crg, title=args.file.replace("\\", "/").split("/")[-1], mode=args.mode,
             zscale=args.zscale, max_points=args.max_points, cmap=args.cmap, save=args.save)


if __name__ == "__main__":
    sys.exit(main())
