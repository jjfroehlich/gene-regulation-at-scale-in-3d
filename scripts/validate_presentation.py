"""Validate saved still/animation scenes in background Blender without saving them.

Run: blender --background --factory-startup --python-exit-code 1
             --python scripts/validate_presentation.py -- --kind all
"""
import argparse
import json
import sys
from pathlib import Path
import bpy
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'animation/scripts')]


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def still():
    import build_gene_expression_surface_scene as builder
    builder.configure_renderer()
    directory=ROOT/'outputs/canonical'
    report=load(directory/'gene_expression_surface_scene_report.json')
    bpy.ops.wm.open_mainfile(filepath=str(directory/'gene_expression_surface_style.blend'))
    scene=bpy.context.scene
    assert scene.camera.name=='Camera_canonical_full_overview'
    builder.validate_canonical_report(report)
    scales=[]
    for name,row in report['closeup_rendering'].items():
        camera=bpy.data.objects[row['camera']]
        # Close-up overlays are temporary; their measured lengths are recorded at render time.
        length=row['scale_bar']['measured_mm']
        expected=row['scale_bar']['length_nm']*report['units']['nm_to_mm']
        assert abs(length-expected)<1e-5,(name,length,expected)
        fraction=length/camera.data.ortho_scale
        assert abs(fraction-row['scale_bar']['fraction_of_width'])<1e-5
        assert .15-1e-5<=fraction<=.25+1e-5,(name,fraction)
        assert row['scale_bar']['camera_plane']
        scales.append({'subject':name,'length_mm':length,'width_fraction':fraction})
    return {'camera':scene.camera.name,'canonical_checks_passed':True,'closeup_scale_bars':scales}


def animation(directory):
    import build_flythrough_animation as builder
    from dark_presentation import snapshot
    report=load(directory/'flythrough_animation_report.json')
    source=load(directory/'staging/baseline_snapshot.json') if (directory/'staging/baseline_snapshot.json').exists() else load(directory/'baseline_snapshot.json')
    bpy.ops.wm.open_mainfile(filepath=str(directory/'flythrough_animation.blend'))
    scene=bpy.context.scene;camera=scene.camera
    assert scene.frame_current==1 and camera.name==report['camera']
    assert [scene.frame_start,scene.frame_end]==[report['frame_start'],report['frame_end']]
    assert scene.frame_end==round(report['duration_seconds']*report['fps'])+1
    assert [scene.render.resolution_x,scene.render.resolution_y]==report['resolution']
    assert scene.render.engine=='CYCLES' and scene.cycles.samples==256
    assert scene.cycles.use_denoising and not scene.render.use_motion_blur
    assert snapshot()==source['molecules'],'Molecular geometry or transforms changed'
    labels=report['individual_asset_label_coverage']
    assert labels['passed'] and labels['actual_count']==labels['expected_count']
    fonts=[o for o in bpy.data.objects if o.type=='FONT']
    hold_start=round((report['duration_seconds']-report['actin_hold_seconds'])*report['fps'])+1
    hold=None;previous=-1.;matrix_error=0.;fov_error=0.
    assert len(source['camera_frames'])==scene.frame_end
    for frame,row in enumerate(source['camera_frames'],1):
        scene.frame_set(frame);bpy.context.view_layer.update()
        matrix=[list(v) for v in camera.matrix_world]
        matrix_error=max(matrix_error,max(abs(a-b) for ra,rb in zip(matrix,row['matrix']) for a,b in zip(ra,rb)))
        fov_error=max(fov_error,abs(camera.data.lens/camera.data.sensor_width-row['lens']/row['sensor']))
        assert row['labels']=={o.name:not o.hide_render for o in fonts},('Label visibility',frame)
        assert bpy.data.objects['Animation_highlight_DNA_3954bp_path'].hide_render
        progress=camera.constraints['Follow Path'].offset_factor
        assert progress>=previous,('Path reverses',frame)
        previous=progress
        assert camera.data.dof.aperture_fstop>=1.19999
        if frame==1 or frame>=hold_start:assert not camera.data.dof.use_dof
        if frame==hold_start:hold=matrix
        if frame>=hold_start:assert matrix==hold,('Final hold moves',frame)
    assert matrix_error<1e-6 and fov_error<1e-6
    canonical=load(ROOT/'outputs/canonical/gene_expression_surface_scene_report.json')
    actin=Vector(next(a['location_mm'] for a in canonical['pdb_assets'] if a['name']=='Actin protein'))
    pointer=report['overview_labels'][0]
    pointer_check=builder.camera_pointer_projection_validation(bpy.data.objects[pointer['pointer']],camera,actin,1,int(pointer['end_frame']))
    assert pointer_check['passed']
    return {'frames':scene.frame_end,'camera_matrix_max_error':matrix_error,'fov_ratio_max_error':fov_error,'geometry_preserved':True,'labels_passed':True,'path_progress_passed':True,'hold_frames':scene.frame_end-hold_start+1,'pointer':pointer_check}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind',choices=['all','still','animation'],default='all')
    parser.add_argument('--animation-dir',type=Path,default=ROOT/'outputs/animation')
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/validation/presentation.json')
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    result={}
    if args.kind in ('all','still'):result['still']=still()
    if args.kind in ('all','animation'):result['animation']=animation(args.animation_dir.resolve())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('PRESENTATION_VALIDATION_PASSED '+str(args.output))

if __name__=='__main__':main()
