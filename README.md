# OpenCRG Viewer

A lightweight 3D viewer for [ASAM OpenCRG](https://github.com/asam-ev/OpenCRG) road surface files (`.crg`).

## Python (`crg_view.py`)

A single file that needs only `numpy` and `matplotlib`. It includes a pure Python reader for the OpenCRG format, so neither the C nor the MATLAB library has to be installed.

```bash
pip install numpy matplotlib
python crg_view.py path/to/road.crg                  # 3D road surface in x/y/z
python crg_view.py road.crg --zscale 10              # 10x vertical exaggeration
python crg_view.py road.crg --mode uv                # unrolled elevation grid (u/v/z)
python crg_view.py road.crg --max-points 600         # finer mesh (default 300 per axis)
python crg_view.py road.crg --save road.png          # write an image, no window
```

### Export to CloudCompare

```bash
python crg_view.py road.crg --export road.ply                  # binary PLY point cloud (recommended)
python crg_view.py road.crg --export road.ply --mesh           # triangulated surface instead of points
python crg_view.py road.crg --export road.xyz                  # ASCII: x y z u v
python crg_view.py road.crg --export road.ply --step 5         # every 5th point in u and v
```

- **What's written:** every grid node at full resolution, unless you set `--step`. Points are in world coordinates, with x, y and z stored as doubles.
- **Scalar fields:** `u` (distance along the road) and `v` (lateral offset) are exported as scalar fields that CloudCompare can colour by.
- **Large coordinates:** CloudCompare will offer a global shift for georeferenced files; accept it.
- **File size:** a full-resolution PLY takes 32 bytes per point. For example, `country_road.crg` (57k x 370 points) comes to about 670 MB, so use `--step` for a quick look.

You can also use it as a library:

```python
from crg_view import read_crg
crg = read_crg("road.crg")
crg["X"], crg["Y"], crg["Z"]   # (nu x nv) road surface in world coordinates
crg["u"], crg["v"], crg["z"]   # grid axes and elevation grid as stored
```

What the reader supports:

- **Data formats:** `LRFI`, `LDFI`, `KRFI`, `KDFI`, `KRBI`, `KDBI`, `LRBI` and `LDBI`, which covers every ASCII and binary format at single or double precision.
- **Long sections:** defined either by index (`long section N`) or by position (`long section at v = ...`).
- **Reference line:** a straight line from the start heading, integrated from the `phi` channel, or taken from explicit `x`/`y` channels.
- **Elevation:** slope and banking, from a channel or from the constant start values in the header.
- **Missing values:** `*missing*` and NaN values are filled like the C library default (`grid_nan_mode = keep last`).

What it does not do:

- `$ROAD_CRG_FILE` include files are rejected.
- `$ROAD_CRG_MODS` and `$ROAD_CRG_OPTS` are ignored, so the grid is shown exactly as stored.

Large grids are decimated for display only.

The reader was checked against the official OpenCRG C library (`crgEvaluv2xy` / `crgEvaluv2z`) on every sample file in `crg-txt/` and `crg-bin/` of the OpenCRG repository. At the grid nodes, x, y and z agree to about 1e-9 m. The one exception is `handmade_straight_minimalist.crg`, where the header gives no `long_section_v_right`. There the reader centres the grid on the reference line, as the MATLAB `crg_read` does, while the C library does not.

## MATLAB (`crg_view.m`)

This is a thin wrapper around the toolbox shipped in the OpenCRG repository:

```matlab
data = crg_view('road.crg', 'C:\path\to\OpenCRG');            % add OpenCRG/matlab to the path, read, plot
data = crg_view('road.crg', 'ZScale', 10);                      % 10x vertical exaggeration
data = crg_view('road.crg', 'Export', 'road.ply');              % point cloud for CloudCompare
data = crg_view('road.crg', 'Export', 'road.xyz', 'Step', 5);   % ASCII, every 5th point
```

The plot does not use the toolbox's own surface helper. That helper fixes the aspect ratio but not the camera, so the road drifts out of the axes when you rotate it. `crg_view.m` instead plots around a local origin with `axis vis3d`, so the road stays in the axes while you rotate and zoom. Double-click the plot to restore the initial view.
