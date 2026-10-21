import copy
import unittest
import numpy as np
from ase import Atoms


def host():
    ni = []
    for z in (.155, .345, .655, .845):
        xy = [(x, y) for x in (.25, .75) for y in (0., .5)] if z < .5 else [
            (x, y) for x in (0., .5) for y in (.25, .75)]
        ni.extend([[x, y, z] for x, y in xy])
    rng = np.random.default_rng(11)
    return Atoms('La24Ni16O56', scaled_positions=np.vstack([
        rng.random((24, 3)), ni, rng.random((56, 3))]),
        cell=[[8, 0, 0], [-.08, 8, 0], [0, 0, 20]], pbc=True)


class DisplayTests(unittest.TestCase):
    def test_aligns_bottom_three_atom_edge_and_preserves_geometry(self):
        from lno_display import prepare_lno_display
        atoms = host()
        original = atoms.copy()
        aligned, info = prepare_lno_display(atoms)
        self.assertEqual(info['x_parallel_input_axis'], 'b')
        np.testing.assert_allclose(atoms.positions, original.positions)
        old_d = atoms.get_all_distances(mic=True)
        np.testing.assert_allclose(aligned.get_all_distances(mic=True), old_d, atol=1e-10)
        self.assertGreater(np.linalg.det(aligned.cell), 0)
        f = aligned.get_scaled_positions()[24:28]
        self.assertTrue(np.all(np.minimum(f[:, 0] % .5, .5-f[:, 0] % .5) < 1e-8))
        np.testing.assert_allclose(np.sort(f[:, 1]), [.25, .25, .75, .75], atol=1e-10)

    def test_repairs_bilayer_cut_by_z_boundary_and_swapped_axes(self):
        from lno_display import prepare_lno_display
        atoms = host()
        atoms.positions += .25*atoms.cell[2]
        atoms.wrap()
        atoms.set_cell([atoms.cell[1], -atoms.cell[0], atoms.cell[2]], scale_atoms=False)
        aligned, info = prepare_lno_display(atoms)
        levels = sorted(set(np.round(aligned.get_scaled_positions()[24:40, 2], 3)))
        np.testing.assert_allclose(levels, [.155, .345, .655, .845], atol=1e-8)
        self.assertNotEqual(info['origin_input_fractional'][2], 0)

    def test_boundary_images_keep_two_slices_and_complete_grid(self):
        from lno_display import prepare_lno_display, periodic_panel_rows
        from ni_octahedra import projection_slices
        aligned, info = prepare_lno_display(host())
        rows = [{'Ni_index_1based': i+25, 'fractional_position': f.tolist()}
                for i, f in enumerate(aligned.get_scaled_positions()[24:40])]
        cell = np.array(aligned.cell)
        x = cell[0]/np.linalg.norm(cell[0]); z = np.cross(cell[0], cell[1]); z /= np.linalg.norm(z)
        result = {'octahedra': rows, 'cell_A': cell.tolist(),
                  'frame_rows_cartesian': [x.tolist(), np.cross(z, x).tolist(), z.tolist()],
                  'display': info}
        before = copy.deepcopy(result)
        for plane in ('yz', 'xz'):
            plan = projection_slices(result, plane)
            self.assertEqual(len(plan['slices']), 2)
            allrows = []
            for panel in plan['slices']:
                drawn = periodic_panel_rows(result, plane, panel)
                self.assertEqual(len(drawn), 10)
                allrows.extend(drawn)
            self.assertEqual(len(set(r['Ni_index_1based'] for r in allrows)), 16)
            self.assertEqual(sum(r['periodic_copy'] for r in allrows), 4)
        self.assertEqual(result, before)

    def test_interstitial_detection_requires_complete_host(self):
        from lno_display import oxygen_markers
        atoms = host(); atoms += Atoms('O', positions=[[2, 2, 10]])
        rows = [{'orientation': {}, 'neighbors': [{'O_index_1based': i} for i in range(41, 97)]} for _ in range(16)]
        markers, info = oxygen_markers(atoms, {'octahedra': rows})
        self.assertEqual([m['O_index_1based'] for m in markers], [97])
        self.assertEqual(info['status'], 'identified')
        rows[0]['orientation'] = None
        markers, info = oxygen_markers(atoms, {'octahedra': rows})
        self.assertEqual(markers, [])
        self.assertEqual(info['status'], 'unresolved')

    def test_unit_cell_interstitial_count_uses_28_oxygen_framework(self):
        from lno_display import oxygen_markers
        atoms = Atoms('La12Ni8O30', positions=np.random.default_rng(5).random((50, 3))*8, cell=[8,8,20], pbc=True)
        rows = [{'orientation': {}, 'neighbors': [{'O_index_1based': i} for i in range(21,49)]} for _ in range(8)]
        markers, info = oxygen_markers(atoms, {'octahedra': rows})
        self.assertEqual([m['O_index_1based'] for m in markers], [49,50])
        self.assertEqual(info['status'], 'identified')


if __name__ == '__main__':
    unittest.main()
