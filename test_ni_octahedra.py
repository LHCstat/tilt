import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read, write
from scipy.spatial.transform import Rotation
from ni_octahedra import analyze, fit_orientation, sign_relation

BONDS=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1.]])
def isolated(rv=(0,0,0), lengths=None):
    rot=Rotation.from_rotvec(np.radians(rv)).as_matrix()
    vectors=BONDS*(np.array(lengths)[:,None] if lengths is not None else 2.)
    return Atoms('NiO6',positions=np.vstack([[5,5,5],vectors@rot.T+[5,5,5]]),cell=[10,10,10],pbc=True)
def ori(atoms,azimuth=0):
    return analyze(atoms,azimuth=azimuth)['octahedra'][0]['orientation']

class GeometryTests(unittest.TestCase):
    def test_exact_known_rotations(self):
        for rv in [(0,0,0),(8,0,0),(0,-9,0),(0,0,12),(5,-7,3)]:
            o=ori(isolated(rv))
            np.testing.assert_allclose(o['omega_frame_deg'],rv,atol=1e-10)
            self.assertLess(o['fit_rms_deg'],2e-6)
    def test_projection_sign_is_not_rotation_sign(self):
        o=ori(isolated((8,0,0)))
        self.assertAlmostEqual(o['projected_cw_yz_deg'],-8)
        self.assertAlmostEqual(o['apical_tilt_deg'],8)
    def test_bond_stretch_not_rotation(self):
        o=ori(isolated((4,3,2),[1.8,1.9,2,2.1,2.2,2.3]))
        np.testing.assert_allclose(o['omega_frame_deg'],[4,3,2],atol=1e-10)
    def test_nonrigid_distortion_reported(self):
        a=isolated();a.positions[1]+=[0,.3,.1]
        self.assertGreater(ori(a)['fit_rms_deg'],1)
    def test_five_and_seven_coordination_rejected(self):
        a=isolated();del a[-1]
        self.assertIsNone(ori(a))
        a=isolated();a+=Atoms('O',positions=[[6,6,6]])
        self.assertIsNone(ori(a))
    def test_explicit_interstitial_exclusion(self):
        a=isolated();a+=Atoms('O',positions=[[6,6,6]])
        result=analyze(a,azimuth=0,excluded_oxygen=[8])
        self.assertEqual(result['summary']['valid_NiO6_count'],1)
        with self.assertRaises(ValueError):analyze(a,azimuth=0,excluded_oxygen=[1])
    def test_periodic_oxygen_images_not_deduplicated(self):
        a=Atoms('NiO3',scaled_positions=[[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5]],cell=[4,4,4],pbc=True)
        r=analyze(a,azimuth=0)
        self.assertEqual(r['summary']['valid_NiO6_count'],1)
        self.assertEqual(len(r['octahedra'][0]['neighbors']),6)
        self.assertEqual(len(r['bridges']),3)
        np.testing.assert_allclose([b['Ni_O_Ni_deg'] for b in r['bridges']],180)
    def test_translation_and_atom_permutation(self):
        a=isolated((4,3,2));a.positions+=[5.1,4.9,5.3];a.wrap()
        a=a[[6,3,1,0,5,2,4]]
        np.testing.assert_allclose(ori(a)['omega_frame_deg'],[4,3,2],atol=1e-10)
    def test_global_rotation_covariance(self):
        a=isolated((4,3,2));q=Rotation.from_rotvec([.3,-.2,.1]).as_matrix()
        a.positions=a.positions@q.T;a.set_cell(a.cell.array@q.T)
        np.testing.assert_allclose(ori(a)['omega_frame_deg'],[4,3,2],atol=1e-10)
    def test_azimuth_and_auto_gauge(self):
        a=isolated((0,0,45))
        np.testing.assert_allclose(ori(a,45)['omega_frame_deg'],0,atol=1e-10)
        r=analyze(a)
        self.assertAlmostEqual(r['reference']['azimuth_deg'],45)
        self.assertIsNotNone(r['reference']['warning'])
    def test_direct_cartesian_selective_and_scale(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=isolated((4,3,2))
            for direct in [False,True]:
                p=Path(tmp)/'POSCAR';write(p,a,format='vasp',direct=direct)
                text=p.read_text().splitlines();text[1]='-1000.0'
                index=7;text.insert(index,'Selective dynamics')
                text[index+2:]=[line+' T F T' for line in text[index+2:]]
                p.write_text('\n'.join(text)+'\n')
                np.testing.assert_allclose(ori(read(p,format='vasp'))['omega_frame_deg'],[4,3,2],atol=1e-9)
    def test_layer_pair_and_sign_relation(self):
        a=isolated((8,0,0));b=isolated((-8,0,0));b.positions+=[0,0,5]
        a+=b;a.set_cell([10,10,20])
        r=analyze(a,azimuth=0)
        self.assertEqual(len(r['layers']),2)
        self.assertTrue(all(p['projected_yz_relation']=='opposite' for p in r['interlayer_pairs']))
        self.assertEqual(sign_relation(.1,-8,.5),'near_zero')
        self.assertEqual(sign_relation(None,8,.5),'unresolved')
    def test_skew_cell_periodic_invariance(self):
        a=isolated((4,3,2));a.set_cell([[10,0,0],[2,10,0],[1,1,10]])
        a.positions += a.cell[0]+a.cell[2]
        a.wrap()
        np.testing.assert_allclose(ori(a)['omega_frame_deg'],[4,3,2],atol=1e-10)
    def test_periodic_layer_cluster_at_boundary(self):
        a=isolated();a.positions += [0,0,4.95];a.wrap()
        b=isolated();b.positions += [4,0,-4.95];b.wrap()
        a+=b
        r=analyze(a,azimuth=0)
        self.assertEqual(len(r['layers']),1)
        self.assertLess(r['layers'][0]['normal_spread_A'],.11)
    def test_invalid_reference(self):
        with self.assertRaises(ValueError):analyze(isolated(),azimuth=float('nan'))
        with self.assertRaises(ValueError):analyze(isolated(),cutoff=4.3)

if __name__=='__main__':unittest.main(verbosity=2)
