"""Dark presentation of the saved flythrough, preserving its evaluated camera route."""
import array
import hashlib
import json
import math
import sys
from pathlib import Path
import bpy
from mathutils import Vector, Quaternion
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import build_gene_expression_surface_scene as canonical
canonical.configure_renderer()


def curves(action):
    result=list(getattr(action,'fcurves',[]))
    for layer in getattr(action,'layers',[]):
        for strip in layer.strips:
            for slot in action.slots:
                bag=strip.channelbag(slot)
                if bag: result.extend(bag.fcurves)
    return result


def snapshot():
    data={}
    for o in bpy.data.objects:
        if o.type not in {'MESH','CURVE','EMPTY'} or any(c.name in {'Beauty','Labels','Scale bars','Detail context','Dark presentation'} or c.name.startswith('Animation') for c in o.users_collection):continue
        h=hashlib.sha256()
        if o.type=='MESH':
            for seq,attr,n,kind in [(o.data.vertices,'co',3,'f'),(o.data.loops,'vertex_index',1,'i'),(o.data.polygons,'loop_total',1,'i')]:
                buf=array.array(kind,[0])*(len(seq)*n);seq.foreach_get(attr,buf);h.update(buf.tobytes())
        elif o.type=='CURVE':
            h.update(repr([(s.type,[tuple(p.co) for p in s.points],[(tuple(p.co),tuple(p.handle_left),tuple(p.handle_right)) for p in s.bezier_points]) for s in o.data.splines]).encode())
        data[o.name]={'matrix':[list(v) for v in o.matrix_world],'hash':h.hexdigest()}
    return data


def smooth(x):
    x=max(0.,min(1.,x));return x*x*(3-2*x)


def presentation_seconds(story):
    """Keep 66 seconds; smoothly redistribute two seconds from the ribosome pass."""
    if story <= 46: return story + 2*smooth(story/46)
    if story < 54: return story + 2*(1-smooth((story-46)/8))
    return story


def source_seconds(presentation):
    lo,hi=0.,66.
    for _ in range(52):
        mid=(lo+hi)/2
        if presentation_seconds(mid)<presentation: lo=mid
        else: hi=mid
    return (lo+hi)/2


def look(o,p):o.rotation_euler=(p-o.location).to_track_quat('-Z','Y').to_euler()


def configure(profile,width,height):
    settings=canonical.renderer.configure_canonical_beauty_render()
    s=bpy.context.scene
    assert s.render.engine=='CYCLES','Cycles is required for the approved treatment'
    s.render.resolution_x=width;s.render.resolution_y=height;s.render.resolution_percentage=100
    s.cycles.samples={'smoke':16,'review':64,'final':256}[profile]
    s.cycles.adaptive_threshold=.015 if profile=='final' else .035
    s.cycles.use_denoising=True;s.cycles.seed=73;s.cycles.use_animated_seed=False
    s.render.use_motion_blur=False;s.render.use_persistent_data=True
    s.render.image_settings.file_format='PNG';s.render.image_settings.color_mode='RGB';s.render.image_settings.color_depth='8'
    s.view_settings.exposure=.25;s.render.use_compositing=False;s.render.use_sequencer=False
    return settings['cycles']['device_backend']


def build(args):
    args.output_dir=args.output_dir.resolve()
    import build_flythrough_animation as original
    bpy.ops.wm.open_mainfile(filepath=str(args.baseline_blend))
    s=bpy.context.scene;cam=s.camera;target=bpy.data.objects['Animation_Flythrough_Target']
    old=json.loads(Path(args.baseline_report).read_text())
    source_fps=s.render.fps/s.render.fps_base
    source_duration=(s.frame_end-1)/source_fps
    duration=args.duration_seconds;fps=args.fps;end=round(duration*fps)+1
    ratio=duration*fps/(source_duration*source_fps)
    s.frame_set(round(27*source_fps)+1);bpy.context.view_layer.update()
    ref_width=(cam.matrix_world.translation-target.matrix_world.translation).length*cam.data.sensor_width/cam.data.lens
    baseline=[]
    fonts=[o for o in bpy.data.objects if o.type=='FONT']
    for f in range(1,end+1):
        sf=1+source_seconds((f-1)/(duration*fps)*66)*source_fps
        s.frame_set(math.floor(sf),subframe=sf-math.floor(sf));bpy.context.view_layer.update()
        baseline.append({'source_frame':sf,'offset_factor':cam.constraints['Follow Path'].offset_factor,'target_local':list(target.location),'matrix':[list(v) for v in cam.matrix_world],'lens':cam.data.lens,'sensor':cam.data.sensor_width,'target':list(target.matrix_world.translation),'labels':{o.name:not o.hide_render for o in fonts}})
    s.frame_set(1);geo=snapshot()
    # Sample original channels at the same monotonic time map as the camera audit.
    # Baking at output frames preserves the evaluated route even across smooth retiming.
    for action in bpy.data.actions:
        for curve in curves(action):
            values=[curve.evaluate(row['source_frame']) for row in baseline]
            if action==cam.animation_data.action and curve.data_path.endswith('offset_factor'):
                values=[row['offset_factor'] for row in baseline]
            if action==target.animation_data.action and curve.data_path=='location':
                values=[row['target_local'][curve.array_index] for row in baseline]
            if action==cam.data.animation_data.action and curve.data_path=='lens':
                values=[row['lens'] for row in baseline]
            curve.keyframe_points.clear()
            curve.keyframe_points.add(end)
            curve.keyframe_points.foreach_set('co',[v for f,value in enumerate(values,1) for v in (f,value)])
            for k in curve.keyframe_points:k.interpolation='LINEAR'
            curve.update()
    s.render.fps=fps;s.render.fps_base=1;s.frame_start=1;s.frame_end=end
    # Confirm timing alone before touching the optical gate or presentation.
    matrix_error=0.;label_errors=[]
    for f,row in enumerate(baseline,1):
        s.frame_set(f);bpy.context.view_layer.update()
        matrix_error=max(matrix_error,max(abs(a-b) for ra,rb in zip(row['matrix'],cam.matrix_world) for a,b in zip(ra,rb)))
        if row['labels']!={o.name:not o.hide_render for o in fonts}:label_errors.append(f)
    assert matrix_error<1e-6,(matrix_error,'Retiming changed the camera route')
    assert not label_errors,label_errors
    for o in bpy.data.objects:
        if o.type=='LIGHT' or any(c.name=='Animation Atmosphere' for c in o.users_collection):o.hide_render=True
    blue=bpy.data.objects.get('Animation_highlight_DNA_3954bp_path')
    if blue:blue.hide_render=True;blue.hide_viewport=True
    copied={}
    for o in bpy.data.objects:
        if o.type!='MESH' or 'surface' not in o.name:continue
        for slot in o.material_slots:
            src=slot.material
            if src is None:continue
            if src.name not in copied:
                mat=src.copy();mat.animation_data_clear()
                if mat.node_tree:mat.node_tree.animation_data_clear()
                bsdf=mat.node_tree.nodes.get('Principled BSDF') if mat.use_nodes else None
                if bsdf:
                    if 'protein' in o.name:bsdf.inputs['Base Color'].default_value=(.25,.34,.46,1)
                    elif not o.name.startswith('DNA B-form'):bsdf.inputs['Base Color'].default_value=(.74,.42,.10,1)
                    bsdf.inputs['Roughness'].default_value=.44;bsdf.inputs['Metallic'].default_value=0;bsdf.inputs['Emission Strength'].default_value=0
                copied[src.name]=mat
            slot.material=copied[src.name]
    s.world=s.world.copy();nodes=s.world.node_tree.nodes;links=s.world.node_tree.links;nodes.clear()
    out=nodes.new('ShaderNodeOutputWorld');ambient=nodes.new('ShaderNodeBackground');ambient.inputs['Color'].default_value=(1,1,1,1);ambient.inputs['Strength'].default_value=.09
    bg=nodes.new('ShaderNodeBackground');bg.inputs['Color'].default_value=(.008,.012,.018,1)
    ray=nodes.new('ShaderNodeLightPath');mix=nodes.new('ShaderNodeMixShader')
    links.new(ray.outputs['Is Camera Ray'],mix.inputs[0]);links.new(ambient.outputs[0],mix.inputs[1]);links.new(bg.outputs[0],mix.inputs[2]);links.new(mix.outputs[0],out.inputs['Surface'])
    collection=bpy.data.collections.new('Dark presentation');s.collection.children.link(collection)
    lights=[]
    for name,energy,size,color,offset in [('Key',2800,.42,(1,.93,.82),(-.85,.85,1.3)),('Fill',650,1.6,(.8,.89,1),(1,-.1,1.2)),('Rim',2600,.3,(.76,.85,1),(.65,.8,-.65))]:
        data=bpy.data.lights.new('Dark '+name,'AREA');data.shape='DISK';data.color=color
        o=bpy.data.objects.new('Dark '+name,data);collection.objects.link(o);lights.append((o,energy,size,Vector(offset)))
    windows=[(11.7,13.2,14.7,'p53'),(16.5,18,19.5,'Nucleosome'),(25,27,29.5,'RNA polymerase II'),(47,50.5,54,'Ribosome')]
    centers={}
    for _,_,_,name in windows:
        objects=[o for o in bpy.data.objects if name.lower() in o.name.lower() and o.type=='MESH' and 'surface protein' in o.name]
        points=[o.matrix_world@Vector(c) for o in objects for c in o.bound_box]
        assert points,name
        centers[name]=Vector(tuple((min(p[i] for p in points)+max(p[i] for p in points))/2 for i in range(3)))
    cam.data.dof.focus_object=None;cam.data.sensor_width=baseline[0]['sensor']*10
    # Replace only optical channels; leave camera, path and look-target actions intact.
    for action in [cam.data.animation_data.action] if cam.data.animation_data else []:
        for curve in curves(action):
            if curve.data_path=='lens':
                for k in curve.keyframe_points:k.co.y*=10;k.handle_left.y*=10;k.handle_right.y*=10
    focus_rows=[]
    for f,row in enumerate(baseline,1):
        s.frame_set(f);bpy.context.view_layer.update()
        story=(row['source_frame']-1)/source_fps
        q=cam.matrix_world.to_quaternion();anchor=target.matrix_world.translation.copy()
        width=(cam.matrix_world.translation-anchor).length*row['sensor']/row['lens']
        span=max(8.,min(240.,20*width/ref_width))
        weight=0.;focus=anchor.copy();sweep=0.
        for a,c,b,name in windows:
            if a<=story<=b:
                weight=smooth((story-a)/(c-a)) if story<=c else 1-smooth((story-c)/(b-c))
                focus=anchor.lerp(centers[name],weight)
                # Return to neutral after each passage; no abrupt light reset.
                sweep=math.radians(10)*math.sin(2*math.pi*(story-a)/(b-a))*weight
        dof=8<=story<55
        cam.data.dof.use_dof=dof;cam.data.dof.aperture_fstop=8-6.8*weight
        cam.data.dof.focus_distance=max(.1,(cam.matrix_world.translation-focus).dot(q@Vector((0,0,1))))
        for prop in ('use_dof','aperture_fstop','focus_distance'):cam.data.dof.keyframe_insert(data_path=prop)
        for o,energy,size,offset in lights:
            delta=Quaternion((0,1,0),sweep)@offset if o==lights[0][0] else offset
            o.location=anchor+q@(delta*span);look(o,anchor);o.data.energy=energy*(span/20)**2;o.data.size=size*span
            for prop in ('location','rotation_euler'):o.keyframe_insert(data_path=prop)
            for prop in ('energy','size'):o.data.keyframe_insert(data_path=prop)
        focus_rows.append({'frame':f,'story_seconds':story,'focus_weight':weight,'light_span':span,'dof':dof,'sweep_degrees':math.degrees(sweep),'aperture_fstop':cam.data.dof.aperture_fstop})
    fov_error=0.
    for f,row in enumerate(baseline,1):
        s.frame_set(f);bpy.context.view_layer.update()
        matrix_error=max(matrix_error,max(abs(a-b) for ra,rb in zip(row['matrix'],cam.matrix_world) for a,b in zip(ra,rb)))
        fov_error=max(fov_error,abs(cam.data.lens/cam.data.sensor_width-row['lens']/row['sensor']))
        assert blue is None or blue.hide_render
    assert matrix_error<1e-6 and fov_error<1e-6,(matrix_error,fov_error)
    s.frame_set(1);assert snapshot()==geo,'Molecular geometry changed'
    backend=configure('final',args.resolution_x,args.resolution_y)
    # Existing continuity validator adapts its step thresholds to the shorter duration.
    continuity=original.camera_motion_continuity(cam,target,1,end,duration,fps,moving_end_frame=round(62/66*duration*fps)+1)
    assert not continuity['failures'],continuity
    report=dict(old)
    # Preserve report field meanings while mapping frame and timeline coordinates.
    def retime(v):
        if isinstance(v,list):return [retime(x) for x in v]
        if isinstance(v,dict):
            return {k:(1+presentation_seconds((x-1)/source_fps)/66*duration*fps if k in {'frame','start_frame','end_frame'} and isinstance(x,(int,float)) else presentation_seconds(x)/66*duration if k=='timeline_seconds' and isinstance(x,(int,float)) else retime(x)) for k,x in v.items()}
        return v
    report=retime(report)
    pointer_row=report['overview_labels'][0]
    canonical_report=json.loads((ROOT/'outputs/canonical/gene_expression_surface_scene_report.json').read_text())
    actin=Vector(next(a['location_mm'] for a in canonical_report['pdb_assets'] if a['name']=='Actin protein'))
    pointer_check=original.camera_pointer_projection_validation(bpy.data.objects[pointer_row['pointer']],cam,actin,1,int(pointer_row['end_frame']))
    assert pointer_check['passed'],pointer_check
    report['overview_actin_pointer_validation']=pointer_check
    report.update(duration_seconds=duration,fps=fps,frame_start=1,frame_end=end,actin_hold_seconds=4*duration/66,rendered=False,rendered_frames=0,resolution=[args.resolution_x,args.resolution_y],camera_motion_continuity=continuity)
    report['dark_treatment']={'backend':backend,'camera_matrix_max_error':matrix_error,'fov_ratio_max_error':fov_error,'label_timing_errors':label_errors,'geometry_preserved':True,'molecular_objects':len(geo),'source_duration':source_duration,'speed_multiplier':source_duration/duration,'dna_centerline_hidden':True,'frame_settings':focus_rows,'world_strength':.09,'roughness':.44,'exposure':.25,'optical_gate_scale':10,'timing_map':{'kind':'smooth monotonic redistribution','source_ribosome_seconds':[46,54],'delivery_ribosome_seconds':[48*duration/66,54*duration/66],'final_hold_seconds':4*duration/66,'route_sampling':'Original animation channels evaluated at mapped source times, baked at every output frame'}}
    report['lighting_objects']=[o.name for o,_,_,_ in lights];report['atmosphere_density']=0;report['compositor_nodes']=[]
    report['render_profile']='final';report['smoke_test']=False
    args.output_dir.mkdir(parents=True,exist_ok=True)
    output=args.output_dir/'flythrough_animation.blend'
    report['experiment_blend']=str(output);report['mp4']=str(args.output_dir/'flythrough_animation_1080p.mp4');report['frames_dir']=str(args.output_dir/'frames')
    s.frame_set(1);s.render.filepath=str(args.output_dir/'frames'/'frame_')
    bpy.context.preferences.filepaths.save_version=0;bpy.ops.wm.save_as_mainfile(filepath=str(output))
    (args.output_dir/'baseline_snapshot.json').write_text(json.dumps({'molecules':geo,'camera_frames':baseline,'source_report':old},indent=2))
    (args.output_dir/'flythrough_animation_report.json').write_text(json.dumps(report,indent=2))
    print('DARK_BUILD_VALIDATED '+str(output),flush=True)


def render(args):
    bpy.ops.wm.open_mainfile(filepath=str(args.render_staged))
    s=bpy.context.scene
    configure(args.render_profile,args.resolution_x,args.resolution_y)
    directory=args.frames_dir.resolve();directory.mkdir(parents=True,exist_ok=True)
    report=json.loads((args.render_staged.parent/'flythrough_animation_report.json').read_text())
    if args.finalize_staged:
        for o in bpy.data.objects:
            if (o.type=='LIGHT' and o.hide_render) or any(c.name=='Animation Atmosphere' for c in o.users_collection):o.hide_viewport=True
        s.frame_set(1);s.render.filepath=str(directory/'frame_')
        bpy.context.preferences.filepaths.save_version=0
        bpy.ops.wm.save_as_mainfile(filepath=str(args.render_staged.resolve()))
        assert s.camera.name==report['camera'] and s.frame_current==1
        print('FINAL_SCENE_VERIFIED',flush=True)
        return
    frames=range(s.frame_start,s.frame_end+1)
    if args.checkpoints:frames=sorted({int(round(row['frame'])) for row in report['storyboard']})
    if args.sample_check:frames=[round(presentation_seconds(t)/66*report['duration_seconds']*s.render.fps)+1 for t in (13.2,18,27,50.5)]
    for frame in frames:
        path=directory/f'frame_{frame:04d}.png'
        if path.exists() and path.stat().st_size>1000:
            with path.open('rb') as stream:
                stream.seek(-12,2)
                if stream.read()==b'\x00\x00\x00\x00IEND\xaeB`\x82':continue
        s.frame_set(frame);s.render.filepath=str(path);bpy.ops.render.render(write_still=True)
    print('DARK_RENDER_COMPLETE '+str(directory),flush=True)
