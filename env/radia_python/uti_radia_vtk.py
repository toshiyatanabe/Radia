"""Utility to display Radia 3D geometry using pyvista instead of ObjDrwOpenGL.

Usage:
    import radia as rad
    from uti_radia_vtk import ObjDrwPyVista

    ObjDrwPyVista(my_radia_obj)
"""

import atexit
import gc
import numpy as np
import pyvista as pv

# Ensure all VTK/pyvista objects are destroyed before Python tears down modules.
# Without this, Python 3.14 + pyvista produces spurious AttributeError /
# ImportError messages in __del__ during interpreter shutdown.
atexit.register(gc.collect)


def ObjDrwPyVista(obj, title='Radia 3D Viewer', show_edges=True):
    """Display a Radia object in a pyvista window.

    :param obj: Radia object index (return value of rad.Obj* calls)
    :param title: window title
    :param show_edges: whether to draw edge lines
    """
    import radia as rad
    data = rad.ObjDrwVTK(obj)

    plotter = pv.Plotter(title=title)

    # --- polygons (faces) ---
    poly = data.get('polygons', {})
    verts = poly.get('vertices', [])
    lengths = poly.get('lengths', [])
    colors = poly.get('colors', [])

    if verts and lengths:
        pts = np.array(verts, dtype=float).reshape(-1, 3)
        # build VTK face connectivity array: [n, i0, i1, ..., n, i0, ...]
        faces = []
        idx = 0
        for n in lengths:
            faces.append(n)
            faces.extend(range(idx, idx + n))
            idx += n
        mesh = pv.PolyData(pts, faces=np.array(faces))

        if colors:
            # colors per face: flat list [r,g,b, r,g,b, ...]
            col = np.array(colors, dtype=float).reshape(-1, 3)
            if col.max() <= 1.0:
                col = (col * 255).astype(np.uint8)
            mesh.cell_data['colors'] = col
            plotter.add_mesh(mesh, scalars='colors', rgb=True,
                             show_edges=show_edges, edge_color='black')
        else:
            plotter.add_mesh(mesh, color='lightblue',
                             show_edges=show_edges, edge_color='black')

    # --- lines (edges/coils) ---
    lines_data = data.get('lines', {})
    lverts = lines_data.get('vertices', [])
    llengths = lines_data.get('lengths', [])
    lcolors = lines_data.get('colors', [])

    if lverts and llengths:
        lpts = np.array(lverts, dtype=float).reshape(-1, 3)
        conn = []
        idx = 0
        for n in llengths:
            conn.append(n)
            conn.extend(range(idx, idx + n))
            idx += n
        lmesh = pv.PolyData(lpts, lines=np.array(conn))
        if lcolors:
            lcol = np.array(lcolors, dtype=float).reshape(-1, 3)
            if lcol.max() <= 1.0:
                lcol = (lcol * 255).astype(np.uint8)
            lmesh.cell_data['colors'] = lcol
            plotter.add_mesh(lmesh, scalars='colors', rgb=True, line_width=2)
        else:
            plotter.add_mesh(lmesh, color='black', line_width=2)

    plotter.show_axes()
    plotter.show()
    plotter.close()
    del plotter
    gc.collect()
