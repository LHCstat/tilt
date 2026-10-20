#!/usr/bin/env python3
"""NiO6 geometry from VASP structures. No VASP execution. See README.md for conventions."""
import argparse
import csv
import hashlib
import itertools
import io
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read
from ase.neighborlist import neighbor_list
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

VERSION = '1.0.1'


def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    if np.any(n < 1e-12):
        raise ValueError('Zero-length vector / coincident atoms')
    return v / n


def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1))))


def sign_relation(a, b, threshold):
    if a is None or b is None:
        return 'unresolved'
    if min(abs(a), abs(b)) < threshold:
        return 'near_zero'
    return 'same' if a * b > 0 else 'opposite'


def fit_orientation(vectors, azimuth):
    """Vectors expressed in the right-handed x,y,z layer frame; column active R."""
    phi = np.radians(azimuth)
    basis = np.array([[np.cos(phi), np.sin(phi), 0],
                      [-np.sin(phi), np.cos(phi), 0], [0, 0, 1.]])
    ideal = np.array([basis[0], -basis[0], basis[1], -basis[1], basis[2], -basis[2]])
    directions = unit(vectors)
    rows, cols = linear_sum_assignment(-(ideal @ directions.T))
    target = directions[cols]
    u, _, vt = np.linalg.svd(ideal.T @ target)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    row_rotation = u @ correction @ vt
    residual = np.degrees(np.arccos(np.clip(np.sum((ideal @ row_rotation) * target, axis=1), -1, 1)))
    rv = np.degrees(Rotation.from_matrix(row_rotation.T).as_rotvec())
    apical = unit(vectors[cols[4]] - vectors[cols[5]])
    return {
        'rotation_matrix': row_rotation.T.tolist(),
        'omega_frame_deg': rv.tolist(), 'omega_pc_deg': (basis @ rv).tolist(),
        'rotation_magnitude_deg': float(np.linalg.norm(rv)),
        'fit_rms_deg': float(np.sqrt(np.mean(residual**2))),
        'fit_max_deg': float(residual.max()),
        'assignment': cols.tolist(),
        'apical_axis_frame': apical.tolist(),
        'apical_tilt_deg': angle(apical, [0, 0, 1]),
        # Screen horizontal +y, vertical +z: top apical atom right -> positive clockwise.
        'projected_cw_yz_deg': float(np.degrees(np.arctan2(apical[1], apical[2]))),
        # Screen horizontal +x, vertical +z.
        'projected_cw_xz_deg': float(np.degrees(np.arctan2(apical[0], apical[2]))),
        'trans_angles_deg': [angle(vectors[cols[k]], vectors[cols[k+1]]) for k in (0, 2, 4)],
    }


def cluster_layers(frac, ni_indices, height, tolerance):
    """Periodic clustering in fractional c; labels are origin dependent, physical gaps are not."""
    parent = {i: i for i in ni_indices}
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in itertools.combinations(ni_indices, 2):
        d = frac[i, 2] - frac[j, 2]
        if abs(d - np.rint(d)) * height < tolerance:
            parent[root(j)] = root(i)
    groups = {}
    for i in ni_indices:
        groups.setdefault(root(i), []).append(i)
    layers = []
    for ids in groups.values():
        z = frac[ids, 2] % 1
        center = float(np.angle(np.mean(np.exp(2j * np.pi * z))) / (2*np.pi) % 1)
        if abs(center-1) < 1e-10:
            center = 0.
        delta = z-center
        delta -= np.rint(delta)
        layers.append({'Ni_indices_1based': [i+1 for i in ids], 'fractional_c': center,
                       'normal_spread_A': float(np.ptp(delta)*height)})
    layers.sort(key=lambda r: r['fractional_c'])
    for k, layer in enumerate(layers):
        layer['layer'] = k+1
    return layers


def analyze(atoms, azimuth=None, cutoff=2.6, search_radius=4.2, layer_tol=.6,
            sign_threshold=.5, excluded_oxygen=()):
    if not np.isfinite([cutoff, search_radius, layer_tol, sign_threshold]).all() or cutoff <= 0 or search_radius <= cutoff or layer_tol <= 0 or sign_threshold <= 0:
        raise ValueError('Require 0 < cutoff < search-radius; positive layer-tol/sign-threshold')
    atoms = atoms.copy()
    atoms.set_pbc(True)
    atoms.wrap()
    cell = atoms.cell.array
    if np.linalg.det(cell) <= 1e-8:
        raise ValueError('Require a nondegenerate right-handed cell; reorder lattice axes explicitly')
    x = unit(cell[0]); z = unit(np.cross(cell[0], cell[1])); y = np.cross(z, x)
    frame = np.array([x, y, z])
    height = float(cell[2] @ z)
    symbols = np.array(atoms.get_chemical_symbols())
    nis = np.flatnonzero(symbols == 'Ni').tolist()
    if not nis or 'O' not in symbols:
        raise ValueError('Input must contain Ni and O with explicit element symbols')
    excluded = set(excluded_oxygen)
    if any(i < 1 or i > len(atoms) or symbols[i-1] != 'O' for i in excluded):
        raise ValueError('--exclude-o must contain valid 1-based oxygen indices')
    ii, jj, dd, vv, ss = neighbor_list('ijdDS', atoms, search_radius, self_interaction=False)
    raw, azimuth_samples = {}, []
    for ni in nis:
        ix = np.flatnonzero((ii == ni) & (symbols[jj] == 'O') & ~np.isin(jj+1, list(excluded)))
        ix = ix[np.argsort(dd[ix], kind='stable')]
        selected = ix[dd[ix] <= cutoff]
        raw[ni] = (ix, selected)
        if len(selected) == 6:
            v = vv[selected] @ frame.T
            u = unit(v)
            for q in u[np.abs(u[:, 2]) < .5]:
                azimuth_samples.append(np.arctan2(q[1], q[0]))
    strength = None
    automatic = azimuth is None
    if automatic:
        if not azimuth_samples:
            raise ValueError('No complete near-upright NiO6 to infer axes; supply --azimuth and inspect coordination')
        m = np.mean(np.exp(4j*np.array(azimuth_samples)))
        strength = float(abs(m))
        if strength < .7:
            raise ValueError('In-plane axes ambiguous (fourfold coherence < 0.7); supply an explicit --azimuth')
        azimuth = float((np.degrees(np.angle(m))/4) % 90)
        if abs(azimuth-90) < 1e-8:
            azimuth = 0.
    if not np.isfinite(azimuth):
        raise ValueError('azimuth must be finite')
    frac = atoms.get_scaled_positions()
    layers = cluster_layers(frac, nis, height, layer_tol)
    layer_for = {i-1: r['layer'] for r in layers for i in r['Ni_indices_1based']}
    rows, bond_records = [], {}
    for ni in nis:
        ix, selected = raw[ni]
        r = {'Ni_index_1based': ni+1, 'fractional_position': frac[ni].tolist(),
             'layer': layer_for[ni], 'coordination_within_cutoff': len(selected),
             'sixth_O_distance_A': float(dd[ix[5]]) if len(ix)>5 else None,
             'seventh_O_distance_A': float(dd[ix[6]]) if len(ix)>6 else None,
             'flags': [], 'orientation': None, 'neighbors': []}
        for k in selected:
            r['neighbors'].append({'O_index_1based': int(jj[k]+1), 'image_shift': ss[k].tolist(),
                'distance_A': float(dd[k]), 'vector_cart_A': vv[k].tolist()})
        if len(selected) != 6:
            r['status'] = 'invalid_coordination'
            r['flags'].append('Not exactly six O within cutoff; no octahedral orientation assigned')
        else:
            lengths = dd[selected]
            r['orientation'] = fit_orientation(vv[selected] @ frame.T, azimuth)
            r['bond_mean_A'] = float(lengths.mean())
            r['bond_distortion'] = float(np.mean((lengths/lengths.mean()-1)**2))
            if len(ix)>6 and dd[ix[6]]-dd[ix[5]] < .15:
                r['flags'].append('Sixth/seventh oxygen separation < 0.15 A; coordination ambiguous')
            if cutoff-lengths.max() < .1:
                r['flags'].append('Longest bond within 0.1 A of cutoff; check cutoff sensitivity')
            if r['orientation']['fit_rms_deg'] > 5:
                r['flags'].append('Rigid-fit RMS > 5 deg: strong non-rigid distortion or wrong reference')
            if r['orientation']['fit_max_deg'] > 15 or r['orientation']['rotation_magnitude_deg'] > 35:
                r['flags'].append('Large rotation/misfit: nearest-axis assignment may be ambiguous')
            r['status'] = 'warning' if r['flags'] else 'ok'
            for k in selected:
                bond_records.setdefault(int(jj[k]), []).append((ni, ss[k], vv[k]))
        rows.append(r)
    by_ni = {r['Ni_index_1based']-1: r for r in rows}
    # Actual shared oxygen topology, with separate periodic images retained.
    bridges, seen = [], set()
    for oxygen, entries in bond_records.items():
        for (i, si, vi), (j, sj, vj) in itertools.combinations(entries, 2):
            shift = si-sj
            if i > j or (i == j and tuple(shift) < (0, 0, 0)):
                i, j, si, sj, vi, vj, shift = j, i, sj, si, vj, vi, -shift
            key = (i, j, *shift, oxygen, *si)
            if key in seen or (i == j and np.all(shift == 0)):
                continue
            seen.add(key)
            disp = atoms.positions[j]+shift@cell-atoms.positions[i]
            bridges.append({'Ni_i': i+1, 'Ni_j': j+1, 'Ni_j_shift': shift.tolist(),
                'O_index': oxygen+1, 'O_image_from_i': si.tolist(),
                'Ni_O_Ni_deg': angle(vi, vj), 'Ni_distance_A': float(np.linalg.norm(disp)),
                'kind': 'interlayer_shared_O' if abs(disp@z) > layer_tol else 'intralayer_shared_O'})
    pairs = []
    # Across disconnected RP slabs there need not be a unique 3D nearest registry.
    # Preserve every tied nearest lateral candidate instead of asserting a global +/- symbol.
    for la in layers if len(layers)>1 else []:
        lb = layers[la['layer'] % len(layers)]
        shift_c = int(lb['fractional_c'] <= la['fractional_c'])
        for i1 in la['Ni_indices_1based']:
            i = i1-1; candidates = []
            for j1 in lb['Ni_indices_1based']:
                j = j1-1
                # Place atoms near their own circular layer centers first.
                iz = int(np.rint(la['fractional_c']-frac[i, 2]))
                jz = int(np.rint(lb['fractional_c']-frac[j, 2])) + shift_c
                shift0 = np.array([0, 0, jz-iz])
                disp0 = atoms.positions[j] + shift0@cell - atoms.positions[i]
                planar = disp0-(disp0@z)*z
                mic, distance = find_mic(planar, cell, pbc=[True, True, False])
                abshift = np.rint((mic-planar)@np.linalg.inv(cell)).astype(int)
                # Enumerate equivalent nearest image registries, not only the MIC tie chosen by ASE.
                for a, b in itertools.product(range(-1, 2), repeat=2):
                    shift = shift0+abshift+np.array([a,b,0])
                    disp = atoms.positions[j]+shift@cell-atoms.positions[i]
                    lateral = np.linalg.norm(disp-(disp@z)*z)
                    candidates.append((lateral, j, shift, disp))
            nearest = min(t[0] for t in candidates)
            chosen = [t for t in candidates if t[0] <= nearest + .02]
            for lateral, j, shift, disp in chosen:
                oi, oj = by_ni[i]['orientation'], by_ni[j]['orientation']
                shared = any(b['Ni_i']==i+1 and b['Ni_j']==j+1 and b['Ni_j_shift']==shift.tolist()
                             or b['Ni_i']==j+1 and b['Ni_j']==i+1 and b['Ni_j_shift']==(-shift).tolist()
                             for b in bridges)
                p = {'Ni_i': i+1, 'Ni_j': j+1, 'Ni_j_shift': shift.tolist(),
                     'layer_i': la['layer'], 'layer_j': lb['layer'],
                     'normal_gap_A': float(disp@z), 'lateral_offset_A': float(lateral),
                     'registry_candidate_count': len(chosen),
                     'kind': 'shared_O' if shared else 'unbridged_candidate',
                     'projected_yz_relation': sign_relation(oi['projected_cw_yz_deg'] if oi else None,
                         oj['projected_cw_yz_deg'] if oj else None, sign_threshold)}
                for axis, k in zip('xyz', range(3)):
                    p['omega_'+axis+'_relation'] = sign_relation(oi['omega_frame_deg'][k] if oi else None,
                         oj['omega_frame_deg'][k] if oj else None, sign_threshold)
                pairs.append(p)
    valid = [r['orientation'] for r in rows if r['orientation']]
    for la in layers:
        members = [r for r in rows if r['layer']==la['layer'] and r['orientation']]
        values = [r['orientation']['projected_cw_yz_deg'] for r in members]
        la.update({'valid_octahedra':len(members), 'cw_yz_mean_deg':float(np.mean(values)) if values else None,
                   'cw_yz_rms_deg':float(np.sqrt(np.mean(np.square(values)))) if values else None})
    return {'program_version':VERSION, 'atom_count':len(atoms),
            'composition':{s:int(sum(symbols==s)) for s in sorted(set(symbols))},
            'cell_A':cell.tolist(), 'frame_rows_cartesian':frame.tolist(),
            'parameters':{'cutoff_A':cutoff,'search_radius_A':search_radius,'layer_tol_A':layer_tol,
                          'sign_threshold_deg':sign_threshold,'excluded_O_1based':sorted(excluded)},
            'reference':{'azimuth_deg':azimuth, 'mode':'auto_fourfold' if automatic else 'explicit',
                         'fourfold_coherence':strength,
                         'warning':'Auto azimuth absorbs uniform rotation about z; use fixed azimuth for comparing structures.' if automatic else None},
            'summary':{'Ni_count':len(rows),'valid_NiO6_count':len(valid),
                       'warning_count':sum(r['status']!='ok' for r in rows),
                       'apical_tilt_rms_deg':float(np.sqrt(np.mean([o['apical_tilt_deg']**2 for o in valid]))) if valid else None,
                       'omega_frame_rms_deg':np.sqrt(np.mean(np.square([o['omega_frame_deg'] for o in valid]),axis=0)).tolist() if valid else None},
            'glazer_symbol':None,
            'scope_note':'Layer analysis assumes ab is the layer plane and c stacks the layers; reorder axes before using other settings.',
            'glazer_note':'No automatic Glazer/space-group assignment: RP stacking, axis registration, equal amplitudes and translational phase relations require validation.',
            'layers':layers, 'octahedra':rows, 'bridges':bridges, 'interlayer_pairs':pairs}


def write_csv(path, rows, fallback_header):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else fallback_header)
        w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()})


def plot_projection(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.cm import ScalarMappable
    frame=np.array(result['frame_rows_cartesian']);cell=np.array(result['cell_A'])
    fig, ax = plt.subplots(figsize=(6,9))
    values=[r['orientation']['projected_cw_yz_deg'] for r in result['octahedra'] if r['orientation']]
    limit=max(1.,max(map(abs,values),default=1.))
    norm=TwoSlopeNorm(0,-limit,limit);cmap=plt.get_cmap('coolwarm')
    for r in result['octahedra']:
        p=np.array(r['fractional_position'])@cell@frame.T
        ori=r['orientation']
        color=cmap(norm(ori['projected_cw_yz_deg'])) if ori else 'gray'
        for bond in r['neighbors']:
            v=np.array(bond['vector_cart_A'])@frame.T
            ax.plot([p[1],p[1]+v[1]],[p[2],p[2]+v[2]],c=color,lw=1.3,alpha=.65)
        ax.scatter(p[1],p[2],c=[color],s=28,zorder=3)
        ax.annotate(str(r['Ni_index_1based']),(p[1],p[2]),xytext=(3,3),textcoords='offset points',fontsize=7)
    ax.set(xlabel='y in layer frame (A)',ylabel='z = normal to ab (A)',
           title='Ni-O bonds projected along +x\nCW-positive apical-axis tilt; labels = POSCAR Ni indices')
    ax.set_aspect('equal');fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,label='Projected CW tilt (deg)',shrink=.65)
    fig.tight_layout();fig.savefig(output/'tilt_yz.png',dpi=180);fig.savefig(output/'tilt_yz.svg');plt.close(fig)


def export(result, output, plot=True):
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    (output/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n', encoding='utf-8')
    flat=[]
    for r in result['octahedra']:
        o=r['orientation'] or {}
        d={k:r[k] for k in ['Ni_index_1based','layer','fractional_position','status','coordination_within_cutoff','sixth_O_distance_A','seventh_O_distance_A']}
        for k in ['apical_tilt_deg','projected_cw_yz_deg','projected_cw_xz_deg','rotation_magnitude_deg','fit_rms_deg','fit_max_deg']:
            d[k]=o.get(k)
        for k,a in enumerate('xyz'):
            d['omega_'+a+'_deg']=o.get('omega_frame_deg',[None]*3)[k]
            d['omega_pc_'+a+'_deg']=o.get('omega_pc_deg',[None]*3)[k]
        d['flags']=r['flags'];flat.append(d)
    write_csv(output/'octahedra.csv',flat,['Ni_index_1based'])
    write_csv(output/'layers.csv',result['layers'],['layer'])
    write_csv(output/'bridges.csv',result['bridges'],['Ni_i','Ni_j','Ni_O_Ni_deg'])
    write_csv(output/'interlayer_pairs.csv',result['interlayer_pairs'],['Ni_i','Ni_j'])
    if plot:
        plot_projection(result,output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('poscar',type=Path)
    p.add_argument('--out',type=Path,default=Path('tilt_result'),help='New output directory; existing paths are never overwritten')
    p.add_argument('--azimuth',type=float,help='Reference in-plane Ni-O azimuth from frame x, degrees (Amam conventional cell: commonly 45; inspect)')
    p.add_argument('--cutoff',type=float,default=2.6)
    p.add_argument('--search-radius',type=float,default=4.2)
    p.add_argument('--layer-tol',type=float,default=.6)
    p.add_argument('--sign-threshold',type=float,default=.5)
    p.add_argument('--exclude-o',default='',help='Comma-separated 1-based oxygen indices explicitly excluded from NiO6 identification')
    p.add_argument('--no-plot',action='store_true')
    args=p.parse_args()
    if args.out.exists():
        p.error('Output directory already exists; choose a new --out')
    try:
        data=args.poscar.read_bytes()
        # ASE's VASP reader supports Direct/Cartesian, selective dynamics and scaling.
        atoms=read(io.StringIO(data.decode('utf-8-sig')),format='vasp')
        result=analyze(atoms,azimuth=args.azimuth,cutoff=args.cutoff,search_radius=args.search_radius,
                       layer_tol=args.layer_tol,sign_threshold=args.sign_threshold,
                       excluded_oxygen=[int(x) for x in args.exclude_o.split(',') if x.strip()])
        result['source']={'path':str(args.poscar.resolve()),'sha256':hashlib.sha256(data).hexdigest(),
                          'geometry_only':'POSCAR cannot establish electronic/ionic convergence, pressure, U or magnetism.'}
        export(result,args.out,not args.no_plot)
    except (ValueError,OSError) as e:
        p.error(str(e))
    print(json.dumps({'output':str(args.out.resolve()),'reference':result['reference'],'summary':result['summary']},ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
