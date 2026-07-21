import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from pathlib import Path
from typing import Optional, Sequence, Tuple, List

def plot_flux(vertices_list, values=None, cmap='viridis', vmin=None, vmax=None,
                          edgecolor='k', linewidth=0.5, show_centroids=False):
    """
    vertices_list : list of (Ni,2) float arrays
        Each item is the (x,y) vertices of one polygon (closed or open; it will be closed).
    values : (M,) array-like or None
        One scalar per polygon for coloring. If None, all faces are same color.
    """
    patches = [Polygon(verts, closed=True) for verts in vertices_list]
    pc = PatchCollection(patches, edgecolor=edgecolor, linewidth=linewidth)

    if values is not None:
        values = np.asarray(values)
        if vmin is None: vmin = np.nanmin(values)
        if vmax is None: vmax = np.nanmax(values)
        pc.set_cmap(cmap)
        pc.set_norm(Normalize(vmin=vmin, vmax=vmax))
        pc.set_array(values)             # color per polygon
    else:
        pc.set_facecolor('#cccccc')      # single fill if no values

    fig, ax = plt.subplots()
    ax.add_collection(pc)

    # auto-zoom to all vertices
    all_xy = np.vstack(vertices_list)
    ax.set_xlim(all_xy[:,0].min(), all_xy[:,0].max())
    ax.set_ylim(all_xy[:,1].min(), all_xy[:,1].max())
    ax.set_aspect('equal', adjustable='box')

    if values is not None:
        cbar = plt.colorbar(pc, ax=ax)
        cbar.set_label('value')

    if show_centroids:
        # geometric centroid of each polygon (simple polygon, no holes)
        cents = [verts.mean(axis=0) for verts in vertices_list]
        ax.scatter([c[0] for c in cents], [c[1] for c in cents], s=10, zorder=3)

    plt.title(r'$\phi_{1},\ dx\ =\ 5.00[cm]$')
    plt.tight_layout

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


# Main

import os
import pandas as pd
os.chdir(r'C:\Users\user\Desktop\Codes\Official_FVM_Code\Deterministic_FVM_Code\Deterministic_FVM_Code')

xy_ctrd=pd.read_csv('Coordinates_InCntrs.out', sep='\s+', header=None)
xy_ctrd.columns=['x', 'y']
xy_vert=pd.read_csv('Coordinates_Vert.out', sep='\s+', header=None)
xy_vert.columns=['x_1','y_1','x_2','y_2','x_3','y_3']

Phi = pd.read_csv('Eigfcns.out', sep='\s+', skiprows=2, header=None)
Phi.columns=['Flux']

Phi['x_c']=np.concatenate((xy_ctrd['x'], xy_ctrd['x']), axis=None)
Phi['y_c']=np.concatenate((xy_ctrd['y'], xy_ctrd['y']), axis=None)
Phi['x_1']=np.concatenate((xy_vert['x_1'], xy_vert['x_1']), axis=None)
Phi['y_1']=np.concatenate((xy_vert['y_1'], xy_vert['y_1']), axis=None)
Phi['x_2']=np.concatenate((xy_vert['x_2'], xy_vert['x_2']), axis=None)
Phi['y_2']=np.concatenate((xy_vert['y_2'], xy_vert['y_2']), axis=None)
Phi['x_3']=np.concatenate((xy_vert['x_3'], xy_vert['x_3']), axis=None)
Phi['y_3']=np.concatenate((xy_vert['y_3'], xy_vert['y_3']), axis=None)


Univ=pd.read_csv('Universes.out', sep='\s+', header=None)
#Phi['uni']=np.concatenate((Univ, Univ), axis=None)
Phi['uni']=Univ

vertices_list=[]
for i in range(np.size(xy_ctrd['x'])):
    vert=np.array([[Phi['x_1'][i], Phi['y_1'][i]], [Phi['x_2'][i], Phi['y_2'][i]], [Phi['x_3'][i], Phi['y_3'][i]]])
    vertices_list.append(vert)

progr=np.arange(1, np.size(xy_ctrd['x'])+1)
Phi['ID']=np.concatenate((progr, progr), axis=None)

#plot_flux(vertices_list, values=Phi["uni"][0:np.size(xy_ctrd['x'])], cmap="plasma", show_centroids=False)

#plot_flux(vertices_list, values=Phi["Flux"][0:np.size(xy_ctrd['x'])], cmap="plasma", show_centroids=False)
#plt.savefig("Flux_g1_Trng.png")
#plot_flux(vertices_list, values=Phi["Flux"][np.size(xy_ctrd['x'])+1:], cmap="plasma", show_centroids=False)
#plt.savefig("Flux_g2_Trng.png")
#print(Univ)

plot_flux2(vertices_list, values=Phi["uni"][0:np.size(xy_ctrd['x'])], cmap="plasma", show_mesh=True, mesh_width=0.1, show_centroids=False, title=r'$universes$', save_as="Flux_g1_Trng", dpi=600)
plot_flux2(vertices_list, values=Phi["Flux"][0:np.size(xy_ctrd['x'])], cmap="plasma", show_mesh=True, mesh_width=0.1, show_centroids=False, title=r'$\phi_{1},\ dx\ =\ 1.25[cm]$', save_as="Flux_g1_Trng", dpi=400)
plot_flux2(vertices_list, values=Phi["Flux"][np.size(xy_ctrd['x'])+1:], cmap="plasma", show_mesh=True, mesh_width=0.1, show_centroids=False, title=r'$\phi_{1},\ dx\ =\ 1.25[cm]$', save_as="Flux_g2_Trng", dpi=400)



plt.show()