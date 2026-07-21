import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from pathlib import Path
from typing import Optional, Sequence, Tuple, List
from matplotlib.colors import ListedColormap, BoundaryNorm
import os

def dist2(p: Tuple[float, float]) -> float:
    return p[0]**2 + p[1]**2

def Coordinates_computer(
    xCyC: Tuple[float, float],
    no_faces: int,
    innr: Optional[float] = None,
    outr: Optional[float] = None,
    rotation: Optional[float] = None,
) -> List[Tuple[float, float]]:
    """
    Return vertices of a regular n-gon centered at xCyC.
    Ordering starts from the vertex closest to the origin (0,0), then CCW.

    Provide either `outr` (circumradius) or `innr` (inradius).
    """
    cx, cy = xCyC

    if no_faces < 3:
        raise ValueError("no_faces must be >= 3")

    if innr is None and outr is None:
        raise ValueError("Provide either innr (inradius) or outr (outer/circumradius).")

    if outr is not None and innr is not None:
        # Optional: verify consistency (not required, but helpful)
        expected_outr = innr / math.cos(math.pi / no_faces)
        # tolerate small numeric differences
        if abs(outr - expected_outr) > 1e-9:
            raise ValueError("innr and outr are inconsistent for the given no_faces.")

    # Use circumradius R
    if outr is not None:
        R = outr
    else:
        R = innr / math.cos(math.pi / no_faces)  # convert inradius -> circumradius

    if R <= 0:
        raise ValueError("Radius must be positive.")

    theta0 = 0.0 if rotation is None else float(rotation)

    # Generate CCW vertices (starting at theta0)
    verts = []
    for k in range(no_faces):
        ang = theta0 + 2 * math.pi * k / no_faces
        x = cx + R * math.cos(ang)
        y = cy + R * math.sin(ang)
        verts.append((x, y))

    # Start from the vertex closest to the origin; tie-break by smaller atan2
    min_idx = min(
        range(no_faces),
        key=lambda i: (dist2(verts[i]), math.atan2(verts[i][1], verts[i][0]))
    )

    # Rotate list so the closest vertex comes first
    ordered = verts[min_idx:] + verts[:min_idx]
    return ordered

def plot_flux2(
    vertices_list: Sequence[np.ndarray],
    values: Optional[Sequence[float]] = None,
    cmap: str = 'viridis',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_mesh: bool = True,
    mesh_color: str = 'k',
    mesh_width: float = 0.5,
    show_centroids: bool = False,
    title: Optional[str] = r'$\phi_{1},\ dx\ =\ 5.00[cm]$',
    save_as: Optional[str] = None,
    dpi: int = 300,
    colorbar_label: str = 'value'
) -> Tuple[plt.Figure, plt.Axes, List[Path]]:
    """
    Plot polygonal mesh with optional scalar coloring.

    Parameters
    ----------
    vertices_list : list of (Ni,2) float arrays
        Each item is the (x,y) vertices of one polygon (open or closed; it is closed).
    values : (M,) array-like or None
        One scalar per polygon for coloring. If None, all faces share a single fill color.
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : float or None
        Color scale limits. If None, inferred from `values`.
    show_mesh : bool
        If False, hides polygon edges (mesh outline removal).
    mesh_color : str
        Edge color when show_mesh=True.
    mesh_width : float
        Edge line width when show_mesh=True.
    show_centroids : bool
        If True, plot centroids of polygons.
    title : str or None
        Figure title. If None, no title.
    save_as : str or None
        Base filename (with or without extension). If provided, saves both
        PNG and SVG to '<base>.png' and '<base>.svg'.
    dpi : int
        PNG export resolution (ignored for SVG).
    colorbar_label : str
        Label for the colorbar when `values` is not None.

    Returns
    -------
    fig, ax, saved_paths : (Figure, Axes, list[Path])
        Matplotlib figure/axes and a list of saved file paths (empty if not saved).
    """
    # Build patches
    patches = [Polygon(np.asarray(verts), closed=True) for verts in vertices_list]

    # Edge appearance
    if show_mesh:
        pc = PatchCollection(patches, edgecolor=mesh_color, linewidth=mesh_width)
    else:
        pc = PatchCollection(patches, edgecolor='none', linewidth=0.0)

    # Face coloring
    if values is not None:
        values = np.asarray(values)
        if vmin is None:
            vmin = np.nanmin(values)
        if vmax is None:
            vmax = np.nanmax(values)
        pc.set_cmap(cmap)
        pc.set_norm(Normalize(vmin=vmin, vmax=vmax))
        pc.set_array(values)  # color per polygon
    else:
        pc.set_facecolor('#cccccc')  # single fill if no values

    # Figure
    fig, ax = plt.subplots()
    ax.add_collection(pc)

    # Auto-zoom to all vertices
    all_xy = np.vstack(vertices_list)
    ax.set_xlim(all_xy[:, 0].min(), all_xy[:, 0].max())
    ax.set_ylim(all_xy[:, 1].min(), all_xy[:, 1].max())
    ax.set_aspect('equal', adjustable='box')

    # Colorbar
    if values is not None:
        cbar = plt.colorbar(pc, ax=ax)
        cbar.set_label(colorbar_label)

    # Centroids
    if show_centroids:
        cents = [np.asarray(verts).mean(axis=0) for verts in vertices_list]
        ax.scatter([c[0] for c in cents], [c[1] for c in cents], s=10, zorder=3)

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()

    # Save to disk (both PNG and SVG) if requested
    saved_paths: List[Path] = []
    if save_as is not None:
        base = Path(save_as)
        # If user passed an extension, strip it to ensure both formats are created
        base = base.with_suffix('')
        png_path = base.with_suffix('.png')
        svg_path = base.with_suffix('.svg')

        fig.savefig(png_path, dpi=dpi, bbox_inches='tight')
        fig.savefig(svg_path, bbox_inches='tight')
        saved_paths.extend([png_path, svg_path])

    return fig, ax, saved_paths


## For octagons
#no_faces=8
## Box
#xy_arr = np.array([(5.0, 5.0),(122.5, 5.0),(245.0, 5.0), (245.0, 122.5), (245.0, 245.0), (122.5, 245.0), (5.0, 245.0), (5.0, 122.5)], dtype=float)[..., None]
#tilt = math.pi/8
#inradius=35.0

# For squares
no_faces=4
# Box
xy_arr = np.array([(0.0, 0.0), (250.0, 0.0), (250.0, 250.0), (0.0, 250.0)], dtype=float)[..., None]
tilt = 0 #math.pi/4 #
inradius=28.28 #20.0 #


coordinates=np.zeros((no_faces,2,0), dtype=float)
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(45.0, 45.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(125.0, 45.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(205.0, 45.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(45.0, 125.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(125.0, 125.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(205.0, 125.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(45.0, 205.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(125.0, 205.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)

xy_arr = np.array(Coordinates_computer(xCyC=(205.0, 205.0), no_faces=no_faces, innr=inradius, rotation=tilt))[..., None]
coordinates=np.concatenate([coordinates, xy_arr], axis=-1)


vertices_list=[]
vert2=[coordinates[:,:,i] for i in range(coordinates.shape[-1])]
vertices_list=vert2

one=np.ones((np.size(coordinates,axis=-1)))
one[0]=0

prog=np.arange(0, np.size(coordinates,axis=-1), 1)
prog[1:]=prog[1:]+20


#blue   = "#8fdaf8"
#orange = "#F0CAB6"

blue   = "#A6CAEC"
orange = "#F2AA84"
cmap = ListedColormap([blue, orange], name="blue_orange_01")


os.chdir(r'C:\Users\user\Desktop\Codes\Official_FVM_Code\Deterministic_FVM_Code\Deterministic_FVM_Code')
#plot_flux2(vertices_list, values=one, show_mesh= True, cmap='bwr', mesh_width = 0.75, show_centroids=False, save_as='Reactor_crossview_view_Octa',  title='Reactor Cross-sectional View',  dpi = 600)
plot_flux2(vertices_list, values=prog, show_mesh= True, cmap=cmap, mesh_width = 0.75, show_centroids=False, save_as='Reactor_crossview_view_Octa',  title='Reactor Cross-sectional View',  dpi = 600)

for i in range(1, coordinates.shape[-1]):
    for j in range(coordinates.shape[0]):
        print(coordinates[j,0,i], coordinates[j,1,i])
    print(' ')
    print(' ')
    print(' ')

plt.show()