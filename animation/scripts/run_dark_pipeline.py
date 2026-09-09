"""Staged, resumable delivery for the existing flythrough runner (standard library only)."""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
BASE=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'animation'
STAGE=OUT/'staging'
ROLLBACK=OUT/'rollback'/'before_dark'
BUILDER=Path(__file__).with_name('build_flythrough_animation.py')
DARK=Path(__file__).with_name('dark_presentation.py')


def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def run(command):
    print('RUN '+str(command[0]),flush=True)
    subprocess.run(list(map(str,command)),cwd=ROOT,check=True)


def protect():
    return {str(p):digest(p) for p in (ROOT/'outputs/canonical/gene_expression_surface_style.blend',ROOT/'outputs/canonical/gene_expression_surface_scene_report.json')}


def preserve(args):
    """Keep the original route when present; bootstrap it without reusing a treated film."""
    names=('flythrough_animation.blend','flythrough_animation_report.json')
    if all((ROLLBACK/name).is_file() for name in names):
        return
    canonical_blend=ROOT/'outputs/canonical/gene_expression_surface_style.blend'
    canonical_report=ROOT/'outputs/canonical/gene_expression_surface_scene_report.json'
    if not canonical_blend.is_file() or not canonical_report.is_file():
        run([args.blender_exe,'--background','--factory-startup','--python-exit-code','1',
             '--python',ROOT/'scripts/build_gene_expression_surface_scene.py'])
    build_dir=ROLLBACK.parent/'baseline_build'
    blender(args,['--build-baseline','--source-blend',canonical_blend,'--source-report',canonical_report,
                  '--output-dir',build_dir,'--duration-seconds',66,'--fps',24,'--skip-video-render'])
    report=json.loads((build_dir/names[1]).read_text())
    assert report['duration_seconds']==66 and report['frame_end']==1585
    assert not report.get('dark_treatment'),'Baseline must be the original untreated route'
    assert not report['camera_motion_continuity']['failures']
    ROLLBACK.mkdir(parents=True,exist_ok=True)
    for name in names:
        source=build_dir/name;pending=ROLLBACK/(name+'.pending')
        shutil.copy2(source,pending)
        assert digest(source)==digest(pending)
        pending.replace(ROLLBACK/name)
    print('ORIGINAL_CAMERA_BASELINE_READY '+str(ROLLBACK),flush=True)


def signature(duration,fps,width=1920,height=1080):
    return hashlib.sha256(json.dumps({'baseline':digest(ROLLBACK/'flythrough_animation.blend'),'baseline_report':digest(ROLLBACK/'flythrough_animation_report.json'),'builder':digest(BUILDER),'treatment':digest(DARK),'duration':float(duration),'fps':int(fps),'resolution':[int(width),int(height)]},sort_keys=True).encode()).hexdigest()


def blender(args,extra):
    run([args.blender_exe,'--background','--factory-startup','--python-exit-code','1','--python',BUILDER,'--',*extra])


def prepare(args):
    preserve(args);sig=signature(args.duration_seconds,args.fps,args.resolution_x,args.resolution_y)
    stamp=STAGE/'generation.json'
    old=json.loads(stamp.read_text()) if stamp.exists() else {}
    if old.get('signature')!=sig:
        if STAGE.exists():
            archive=OUT/('staging_previous_'+old.get('signature',digest(DARK))[:12])
            counter=1
            while archive.exists():archive=OUT/('staging_previous_'+sig[:12]+'_'+str(counter));counter+=1
            # Both endpoints are direct children of the animation output directory.
            assert STAGE.resolve().parent==OUT.resolve() and archive.resolve().parent==OUT.resolve()
            STAGE.rename(archive)
        STAGE.mkdir(parents=True,exist_ok=True)
        blender(args,['--baseline-blend',ROLLBACK/'flythrough_animation.blend','--baseline-report',ROLLBACK/'flythrough_animation_report.json','--duration-seconds',args.duration_seconds,'--fps',args.fps,'--resolution-x',args.resolution_x,'--resolution-y',args.resolution_y,'--output-dir',STAGE,'--skip-video-render'])
        stamp.write_text(json.dumps({'signature':sig,'protected':protect()},indent=2))
    return sig


def render(args,profile,checkpoints=False,sample=False):
    width,height=(640,360) if profile=='smoke' else (1280,720) if profile=='review' else (args.resolution_x,args.resolution_y)
    if sample:width,height=args.resolution_x,args.resolution_y
    name='checkpoints' if checkpoints else 'sample_'+profile if sample else 'frames_review' if profile=='review' else 'frames'
    extra=['--render-staged',STAGE/'flythrough_animation.blend','--render-profile',profile,'--resolution-x',width,'--resolution-y',height,'--frames-dir',STAGE/name]
    if checkpoints:extra+=['--checkpoints']
    if sample:extra+=['--sample-check']
    blender(args,extra)


def probe(args,p):
    data=json.loads(subprocess.check_output([args.ffprobe_exe,'-v','error','-count_frames','-show_streams','-of','json',str(p)]))
    s=next(v for v in data['streams'] if v['codec_type']=='video')
    return {k:s.get(k) for k in ('width','height','duration','nb_read_frames','r_frame_rate','codec_name')}


def encode(args,profile):
    report=json.loads((STAGE/'flythrough_animation_report.json').read_text())
    directory=STAGE/('frames_review' if profile=='review' else 'frames')
    video=STAGE/('flythrough_animation_review.mp4' if profile=='review' else 'flythrough_animation_1080p.mp4')
    count=report['frame_end']-report['frame_start']+1
    assert len(list(directory.glob('frame_*.png')))==count
    run([args.ffmpeg_exe,'-hide_banner','-loglevel','error','-y','-framerate',args.fps,'-start_number',1,'-i',directory/'frame_%04d.png','-frames:v',count,'-c:v','libx264','-pix_fmt','yuv420p','-crf',19,'-movflags','+faststart',video])
    info=probe(args,video)
    assert int(info['nb_read_frames'])==count
    assert info['r_frame_rate']==str(args.fps)+'/1'
    assert abs(float(info['duration'])-count/args.fps)<.01
    assert [info['width'],info['height']]==([1280,720] if profile=='review' else [args.resolution_x,args.resolution_y])
    (STAGE/(profile+'_video_validation.json')).write_text(json.dumps(info,indent=2))
    return info


def gallery():
    report=json.loads((STAGE/'flythrough_animation_report.json').read_text())
    (STAGE/'index.html').write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Dark molecular flythrough</title><style>body{background:#17191c;color:#eee;font:17px/1.6 system-ui;max-width:1400px;margin:40px auto;padding:24px}video,img{width:100%}a{color:#dcc89e}</style><h1>Dark molecular flythrough</h1><p>__DURATION__-second journey with a shorter ribosome pass. Original camera route, matching DNA segment colors, gentle focus shifts and light sweeps.</p><button id="restart" style="padding:12px 20px;margin-bottom:16px">Play from start</button><video controls playsinline preload="metadata" src="flythrough_animation_review.mp4"></video><p><a href="flythrough_animation_1080p.mp4">1080p final</a> · <a href="flythrough_animation_report.json">Validation and rendering report</a></p><img src="review_contact_sheet.png" alt="Storyboard contact sheet"><script>document.querySelector("#restart").onclick=()=>{const v=document.querySelector("video");v.currentTime=0;v.play();};</script></html>'''.replace('__DURATION__',f"{report['duration_seconds']:g}").replace('src="flythrough_animation_review.mp4"','src="flythrough_animation_1080p.mp4"' if (STAGE/'flythrough_animation_1080p.mp4').exists() else 'src="flythrough_animation_review.mp4"'),encoding='utf-8')


def deliver(args):
    report=json.loads((STAGE/'flythrough_animation_report.json').read_text())
    validations={p:json.loads((STAGE/(p+'_video_validation.json')).read_text()) for p in ('review','final')}
    run([args.ffmpeg_exe,'-hide_banner','-loglevel','error','-y','-i',STAGE/'flythrough_animation_review.mp4','-vf',f"fps=20/{report['duration_seconds']},scale=320:180,tile=5x4",'-frames:v',1,STAGE/'review_contact_sheet.png'])
    if not args.skip_readme_gif:
        run([args.ffmpeg_exe,'-hide_banner','-loglevel','error','-y','-i',STAGE/'flythrough_animation_1080p.mp4','-vf','fps=24,scale=280:-2:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=64[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5',STAGE/'flythrough-preview.gif'])
    report.update(rendered=True,rendered_frames=report['frame_end'],encoded_video_validation=validations,resolution=[args.resolution_x,args.resolution_y],render_profile='final',experiment_blend=str(OUT/'flythrough_animation.blend'),mp4=str(OUT/'flythrough_animation_1080p.mp4'),frames_dir=str(STAGE/'frames'))
    report['render_settings']={'review_samples':64,'final_samples':256,'review_threshold':.035,'final_threshold':.015,'seed':73,'denoising':True,'motion_blur':False}
    report['delivery']={'rollback':str(ROLLBACK),'staging':str(STAGE),'gif_speed_multiplier':1,'endpoint_convention':f"Inclusive final endpoint retained: nominal {report['duration_seconds']} seconds plus one frame at {args.fps} fps"}
    (STAGE/'flythrough_animation_report.json').write_text(json.dumps(report,indent=2));gallery()
    blender(args,['--render-staged',STAGE/'flythrough_animation.blend','--render-profile','final','--resolution-x',args.resolution_x,'--resolution-y',args.resolution_y,'--frames-dir',STAGE/'frames','--finalize-staged'])
    expected=json.loads((STAGE/'generation.json').read_text())['protected']
    assert protect()==expected,'Canonical scene or report changed'
    validation={'videos':validations,'duration_seconds':report['duration_seconds'],'frame_count':report['frame_end'],'camera_matrix_max_error':report['dark_treatment']['camera_matrix_max_error'],'fov_ratio_max_error':report['dark_treatment']['fov_ratio_max_error'],'geometry_preserved':report['dark_treatment']['geometry_preserved'],'dna_centerline_hidden':report['dark_treatment']['dna_centerline_hidden'],'timing_map':report['dark_treatment']['timing_map'],'canonical_artifacts_unchanged':True,'scene_opening_frame':1,'scene_camera':report['camera']}
    if not args.skip_readme_gif:
        gif_info=probe(args,STAGE/'flythrough-preview.gif')
        assert abs(float(gif_info['duration'])-report['frame_end']/args.fps)<.05,gif_info
        validation['gif']=gif_info
    (STAGE/'delivery_validation.json').write_text(json.dumps(validation,indent=2))
    names=['delivery_validation.json','flythrough_animation.blend','flythrough_animation_report.json','flythrough_animation_review.mp4','flythrough_animation_1080p.mp4','review_contact_sheet.png','index.html']
    for name in names:
        source=STAGE/name;tmp=OUT/(name+'.pending');shutil.copy2(source,tmp);assert digest(source)==digest(tmp);tmp.replace(OUT/name)
    if not args.skip_readme_gif:
        dest=ROOT/'docs/images/flythrough-preview.gif';tmp=dest.with_suffix('.pending.gif');shutil.copy2(STAGE/'flythrough-preview.gif',tmp);tmp.replace(dest)
    print('DARK_FLYTHROUGH_DELIVERED '+str(OUT),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--blender-exe',default=r'C:\Program Files\Blender Foundation\Blender 5.1\blender.exe')
    p.add_argument('--ffmpeg-exe',default='ffmpeg');p.add_argument('--ffprobe-exe',default='ffprobe')
    p.add_argument('--duration-seconds',type=float,default=66);p.add_argument('--fps',type=int,default=24)
    p.add_argument('--resolution-x',type=int,default=1920);p.add_argument('--resolution-y',type=int,default=1080)
    for flag in ('smoke-test','review-render','keep-frames','skip-video-render','skip-readme-gif'):p.add_argument('--'+flag,action='store_true')
    args=p.parse_args()
    assert not(args.smoke_test and args.review_render)
    prepare(args)
    if args.skip_video_render:return
    if args.smoke_test:render(args,'smoke',checkpoints=True);return
    render(args,'smoke',checkpoints=True)
    render(args,'review');encode(args,'review')
    if args.review_render:gallery();return
    render(args,'final',sample=True)
    render(args,'final');encode(args,'final');deliver(args)

if __name__=='__main__':main()
