"""Display conventions for the 24-La / 16-Ni sqrt2 bilayer cells.

Only origin/lattice-basis changes and periodic images are used. No atom is
displaced relative to the periodic structure; original atom indices survive.
"""
import itertools
import numpy as np


def is_sqrt2_lno(atoms):
    symbols = list(atoms.get_chemical_symbols())
    return (symbols.count('La') == 24 and symbols.count('Ni') == 16
            and symbols.count('O') >= 56 and set(symbols) == {'La', 'Ni', 'O'})


def prepare_lno_display(atoms):
    from ni_octahedra import cluster_layers
    if not is_sqrt2_lno(atoms):
        raise ValueError('LNO display requires La24 Ni16 O(56+delta)')
    aligned = atoms.copy()
    aligned.set_pbc(True)
    old_cell = np.array(atoms.cell)
    if np.linalg.det(old_cell) <= 0:
        raise ValueError('LNO display requires a right-handed input cell')
    normal = np.cross(old_cell[0], old_cell[1])
    normal /= np.linalg.norm(normal)
    ni = [i for i, s in enumerate(atoms.symbols) if s == 'Ni']
    frac = atoms.get_scaled_positions(wrap=True)
    layers = cluster_layers(frac, ni, old_cell[2] @ normal, .6)
    if len(layers) != 4 or any(len(l['Ni_indices_1based']) != 4 for l in layers):
        raise ValueError('LNO display requires four Ni layers with four independent Ni each')
    centers = np.array([l['fractional_c'] for l in layers])
    gaps = np.diff(np.r_[centers, centers[0]+1])
    # Preserve an existing complete-bilayer origin. Otherwise place the new
    # boundary halfway inside the first wide gap between bilayers.
    origin_z = 0.
    if gaps[-1] < .8*gaps.max():
        k = int(np.flatnonzero(gaps >= .9*gaps.max())[0])
        origin_z = float((centers[k]+gaps[k]/2) % 1)
    frac[:, 2] = (frac[:, 2]-origin_z) % 1.
    layers = cluster_layers(frac, ni, old_cell[2] @ normal, .6)
    bottom = [i-1 for i in layers[0]['Ni_indices_1based']]
    # The three-image edge has independent rows at 0 and 1/2, with 1 a copy
    # of 0. The perpendicular edge has rows near 1/4 and 3/4.
    folded = frac[bottom, :2] % .5
    scores = np.mean(np.minimum(folded, .5-folded), axis=0)
    x_axis = int(np.argmin(scores))
    if scores[x_axis] > .08 or scores[1-x_axis] < .17:
        raise ValueError('Bottom Ni 2x3 registry is ambiguous; inspect origin/axes before using LNO display')
    basis = np.eye(3, dtype=int) if x_axis == 0 else np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    aligned.positions -= origin_z*old_cell[2]
    aligned.set_cell(basis @ old_cell, scale_atoms=False)
    aligned.wrap()
    info = {
        'mode': 'lno_boundary_complete',
        'x_parallel_input_axis': 'ab'[x_axis],
        'y_toward_input_axis': '+b' if x_axis == 0 else '-a',
        'basis_rows_in_input_basis': basis.tolist(),
        'origin_input_fractional': [0., 0., origin_z],
        'input_cell_A': old_cell.tolist(),
        'fractional_display_bounds': [[-.05, 1.05], [-.05, 1.05], [-.005, 1.005]],
        'note': 'x follows the bottom Ni three-image edge; y is perpendicular in the layer plane. Periodic image copies are display-only.',
    }
    return aligned, info


def oxygen_markers(atoms, result):
    """Identify excess O only when the complete La3Ni2O7 framework exists."""
    symbols = list(atoms.get_chemical_symbols())
    oxygen = {i+1 for i, s in enumerate(symbols) if s == 'O'}
    rows = result['octahedra']
    used = {b['O_index_1based'] for row in rows for b in row['neighbors']}
    ni_count = symbols.count('Ni')
    host_oxygen = 7*ni_count//2
    complete = (ni_count in (8, 16) and symbols.count('La') == 3*ni_count//2
                and set(symbols) == {'La', 'Ni', 'O'} and len(rows) == ni_count
                and all(row['orientation'] is not None for row in rows)
                and len(used) == host_oxygen and used <= oxygen)
    if not complete:
        return [], {'status': 'unresolved', 'reason': 'A complete La3Ni2O7 NiO6 framework was not established; unused O is not automatically interstitial.'}
    frac = atoms.get_scaled_positions(wrap=True)
    markers = [{'O_index_1based': i, 'fractional_position': frac[i-1].tolist()}
               for i in sorted(oxygen-used)]
    return markers, {'status': 'identified', 'count': len(markers),
                     'method': f'Excess oxygen outside the complete {host_oxygen}-O framework of {ni_count} six-coordinate Ni; no atom-order assumption.'}


def _panel_images(records, plane, panel, bounds):
    hidden, horizontal = (0, 1) if plane == 'yz' else (1, 0)
    start = panel['start_fraction']
    width = panel['width_fraction']
    images = []
    for row in records:
        f = np.asarray(row['fractional_position'], dtype=float)
        if (f[hidden]-start) % 1. >= width:
            continue
        # Exactly one representative along depth: the third boundary row
        # is a periodic image, not a third independent slice.
        depth_shift = int(np.ceil(start-f[hidden]-1e-10))
        for h, z in itertools.product((-1, 0, 1), repeat=2):
            shift = np.zeros(3, dtype=int)
            shift[hidden], shift[horizontal], shift[2] = depth_shift, h, z
            position = f+shift
            if not all(bounds[k][0]-1e-10 <= position[k] <= bounds[k][1]+1e-10 for k in (horizontal, 2)):
                continue
            image = dict(row)
            image.update(fractional_position=position.tolist(), display_image_shift=shift.tolist(),
                         periodic_copy=bool(h or z))
            images.append(image)
    return images


def periodic_panel_rows(result, plane, panel):
    return _panel_images(result['octahedra'], plane, panel,
                         result['display']['fractional_display_bounds'])


def periodic_panel_oxygen(result, plane, panel):
    return _panel_images(result.get('interstitial_oxygen', []), plane, panel,
                         result['display']['fractional_display_bounds'])
