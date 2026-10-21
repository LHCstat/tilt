import copy
import json
from pathlib import Path
import tempfile
import unittest

from ni_octahedra import analyze, export
from test_ni_octahedra import isolated


def supercell_result():
    result = analyze(isolated((8, -3, 0)), azimuth=0)
    template = result['octahedra'][0]
    result['cell_A'] = [[8, 0, 0], [-.1, 8, 0], [0, 0, 20]]
    result['octahedra'] = []
    for layer, coords in enumerate([
        [(.001, .25), (.501, .25), (.499, .75), (.999, .75)],
        [(.25, .001), (.25, .501), (.75, .499), (.75, .999)],
    ], 1):
        for x, y in coords:
            row = copy.deepcopy(template)
            row.update(Ni_index_1based=len(result['octahedra'])+1,
                       fractional_position=[x, y, .2*layer], layer=layer)
            result['octahedra'].append(row)
    return result


class SliceTests(unittest.TestCase):
    def test_periodic_slabs_separate_both_views_without_mutation(self):
        from ni_octahedra import projection_slices
        result = supercell_result()
        before = copy.deepcopy(result)
        for plane in ('yz', 'xz'):
            plan = projection_slices(result, plane)
            groups = [s['Ni_indices_1based'] for s in plan['slices']]
            self.assertEqual(sorted(sum(groups, [])), list(range(1, 9)))
            self.assertEqual([len(g) for g in groups], [4, 4])
            self.assertGreater(len(plan['overlapping_pairs_before']), 0)
            self.assertEqual(plan['overlapping_pairs_after'], [])
            # Near-zero and near-one depths stay together across the seam.
            pair = [1, 4] if plane == 'yz' else [5, 8]
            self.assertTrue(any(set(pair) <= set(g) for g in groups))
        self.assertEqual(result, before)

    def test_single_octahedron_needs_no_split(self):
        from ni_octahedra import projection_slices
        result = analyze(isolated(), azimuth=0)
        for plane in ('yz', 'xz'):
            self.assertIsNone(projection_slices(result, plane))

    def test_two_slabs_report_unresolved_threefold_overlap(self):
        from ni_octahedra import projection_slices
        result = supercell_result()
        result['octahedra'] = result['octahedra'][:3]
        for row, depth in zip(result['octahedra'], [.1, .4, .7]):
            row['fractional_position'] = [depth, .25, .2]
        plan = projection_slices(result, 'yz')
        self.assertEqual(len(plan['overlapping_pairs_before']), 3)
        self.assertEqual(len(plan['overlapping_pairs_after']), 1)
        self.assertIsNotNone(plan['warning'])

    def test_export_draws_each_ni_once_in_each_slice_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'slices'
            export(supercell_result(), output)
            metadata = json.loads((output/'projection_slices.json').read_text())
            for plane in ('yz', 'xz'):
                self.assertEqual(len(metadata[plane]['slices']), 2)
                svg = (output/f'tilt_{plane}_slices.svg').read_text(encoding='utf-8')
                for ni in range(1, 9):
                    self.assertEqual(svg.count(f'id="octahedron-{ni}"'), 1)
                self.assertNotIn('stroke-dasharray', svg)
                self.assertTrue((output/f'tilt_{plane}.png').exists())
                self.assertGreater((output/f'tilt_{plane}_slices.png').stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
