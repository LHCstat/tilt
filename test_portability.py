import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parent

class PortabilityTests(unittest.TestCase):
    def test_unicode_spaces_plot_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/'中文 folder';folder.mkdir()
            source=folder/'POSCAR 测试'
            lines=(ROOT/'examples/POSCAR_rx8').read_text(encoding='utf-8').splitlines()
            lines[0]='人工八面体测试'
            source.write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
            initial=hashlib.sha256(source.read_bytes()).hexdigest()
            out=folder/'结果 result'
            cmd=[sys.executable,'-X','utf8',str(ROOT/'ni_octahedra.py'),str(source),'--azimuth','0','--out',str(out)]
            run=subprocess.run(cmd,capture_output=True,encoding='utf-8')
            self.assertEqual(run.returncode,0,run.stderr)
            data=json.loads((out/'analysis.json').read_text(encoding='utf-8'))
            self.assertAlmostEqual(data['octahedra'][0]['orientation']['omega_frame_deg'][0],8,places=8)
            self.assertAlmostEqual(data['octahedra'][0]['orientation']['projected_cw_yz_deg'],-8,places=8)
            for file in ['analysis.json','octahedra.csv','layers.csv','bridges.csv','interlayer_pairs.csv','tilt_yz.png','tilt_yz.svg','tilt_xz.png','tilt_xz.svg']:
                self.assertGreater((out/file).stat().st_size,0)
            saved=(out/'analysis.json').read_bytes()
            again=subprocess.run(cmd,capture_output=True,encoding='utf-8')
            self.assertNotEqual(again.returncode,0)
            self.assertEqual(saved,(out/'analysis.json').read_bytes())
            self.assertEqual(initial,hashlib.sha256(source.read_bytes()).hexdigest())
    def test_projection_labels_show_only_projected_angle(self):
        from ni_octahedra import analyze, export
        from test_ni_octahedra import isolated
        with tempfile.TemporaryDirectory() as tmp:
            for name, rotation, yz, xz in [
                ('rx', (8,0,0), '-8.00', '0.00'),
                ('ry', (0,-9,0), '0.00', '-9.00'),
            ]:
                out=Path(tmp)/name
                export(analyze(isolated(rotation),azimuth=0),out)
                for plane, expected in [('yz',yz),('xz',xz)]:
                    self.assertTrue((out/f'tilt_{plane}.svg').exists())
                    svg=(out/f'tilt_{plane}.svg').read_text(encoding='utf-8')
                    self.assertTrue(f'<!-- {expected}° -->' in svg, f'{plane}: missing compact angle')
                    self.assertTrue('id="octahedron-1"' in svg, f'{plane}: missing octahedral outline')
                    self.assertNotIn('alpha_',svg)
                    self.assertNotIn('theta=',svg)
                    self.assertNotIn('stroke-dasharray',svg)
    def test_launcher_demo_and_error_propagation(self):
        with tempfile.TemporaryDirectory() as t:
            args=[sys.executable,'-X','utf8',str(ROOT/'run_local.py'),'--demo','--out',str(Path(t)/'demo')]
            run=subprocess.run(args,capture_output=True,encoding='utf-8')
            self.assertEqual(run.returncode,0,run.stderr)
            run=subprocess.run(args,capture_output=True,encoding='utf-8')
            self.assertNotEqual(run.returncode,0)

if __name__=='__main__':unittest.main(verbosity=2)
