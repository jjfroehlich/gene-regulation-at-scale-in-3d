"""First-run and recovery tests; run with python -m unittest discover -s tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('pipeline',ROOT/'animation/scripts/run_dark_pipeline.py')
pipeline=importlib.util.module_from_spec(spec);spec.loader.exec_module(pipeline)

class BootstrapTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.out=self.root/'outputs/animation';self.rollback=self.out/'rollback/before_dark'
        for name,value in [('ROOT',self.root),('OUT',self.out),('ROLLBACK',self.rollback)]:
            guard=patch.object(pipeline,name,value);guard.start();self.addCleanup(guard.stop)
        self.args=SimpleNamespace(blender_exe='blender')

    def canonical(self):
        folder=self.root/'outputs/canonical';folder.mkdir(parents=True,exist_ok=True)
        (folder/'gene_expression_surface_style.blend').write_bytes(b'canonical')
        (folder/'gene_expression_surface_scene_report.json').write_text('{}')

    def baseline_build(self,args,extra):
        self.assertIn('--build-baseline',extra)
        folder=Path(extra[extra.index('--output-dir')+1]);folder.mkdir(parents=True,exist_ok=True)
        (folder/'flythrough_animation.blend').write_bytes(b'original route')
        (folder/'flythrough_animation_report.json').write_text(json.dumps({'duration_seconds':66,'frame_end':1585,'camera_motion_continuity':{'failures':[]}}))

    def test_retained_animation_geometry_matches_declared_checksum(self):
        folder=ROOT/"animation/assets"
        metadata=json.loads((folder/"baseline.json").read_text())
        self.assertEqual(hashlib.sha256((folder/"compact_mrna.blend").read_bytes()).hexdigest(),metadata["geometry_sha256"])
        self.assertEqual(len(metadata["objects"]),4)

    def test_fresh_workspace_builds_canonical_and_original_route(self):
        with patch.object(pipeline,'run',side_effect=lambda command:self.canonical()) as canonical, patch.object(pipeline,'blender',side_effect=self.baseline_build) as build:
            pipeline.preserve(self.args)
            canonical.assert_called_once();build.assert_called_once()
        self.assertEqual((self.rollback/'flythrough_animation.blend').read_bytes(),b'original route')

    def test_existing_pair_is_preserved_without_building(self):
        self.rollback.mkdir(parents=True)
        for name in ('flythrough_animation.blend','flythrough_animation_report.json'):(self.rollback/name).write_bytes(b'keep')
        with patch.object(pipeline,'run') as canonical,patch.object(pipeline,'blender') as build:
            pipeline.preserve(self.args);canonical.assert_not_called();build.assert_not_called()
        self.assertEqual((self.rollback/'flythrough_animation.blend').read_bytes(),b'keep')

    def test_partial_baseline_rebuilds_without_copying_treated_delivery(self):
        self.canonical();self.rollback.mkdir(parents=True)
        (self.rollback/'flythrough_animation_report.json').write_text('{}')
        (self.out/'flythrough_animation.blend').write_bytes(b'dark film must not become original baseline')
        with patch.object(pipeline,'run') as canonical,patch.object(pipeline,'blender',side_effect=self.baseline_build):
            pipeline.preserve(self.args);canonical.assert_not_called()
        self.assertEqual((self.rollback/'flythrough_animation.blend').read_bytes(),b'original route')
        self.assertEqual((self.out/'flythrough_animation.blend').read_bytes(),b'dark film must not become original baseline')

    def test_failed_build_does_not_publish_a_baseline(self):
        self.canonical()
        with patch.object(pipeline,'blender',side_effect=RuntimeError('failed')):
            with self.assertRaises(RuntimeError):pipeline.preserve(self.args)
        self.assertFalse((self.rollback/'flythrough_animation.blend').exists())

if __name__=='__main__':unittest.main()
