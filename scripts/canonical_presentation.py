"""Camera-space annotations and additional close-ups for canonical still renders.

Molecular meshes and transforms are read-only here; context geometry is copied.
"""
from __future__ import annotations

import math
import re

import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree


DETAILS = {
    "p53_dna": ("p53 tetramer + DNA", "Structure: 3TS8"),
    "nucleosome_loop": ("Nucleosome + wrapped DNA", "Structure: 1AOI"),
    "polymerase_gene_end": ("RNA polymerase II at the gene end", "Structure: 2E2I · DNA context and nascent RNA 3′ end"),
    "ribosome_trna": ("Ribosome + reference tRNA", "Structures: 1J5E / 1JJ2 · Detached reference tRNA: 4TNA"),
    "actin_product": ("Actin protein", "Structure: 1J6Z · 375 amino acids"),
    "cas9_dna": ("Cas9 + guide / target DNA", "Structure: 4UN3"),
}


def overlay_material(r, source):
    name = "presentation_" + source.name
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)
        mat.use_nodes = True
        mat.node_tree.nodes.clear()
        out = mat.node_tree.nodes.new("ShaderNodeOutputMaterial")
        emission = mat.node_tree.nodes.new("ShaderNodeEmission")
        emission.inputs["Color"].default_value = source.diffuse_color
        mat.node_tree.links.new(emission.outputs[0], out.inputs["Surface"])
    return mat


def text(r, camera, collections, materials, name, body, x, y, pixels=26, align="LEFT", detail=False, source=None, material="black"):
    obj = source or r.base.create_text(name, body, (0, 0, 0), 1, materials["black"], collections["Labels"])
    obj.data.body = "H"  # One cap-height reference gives every title/label a consistent type size.
    obj.data.size = 1
    obj.data.offset = 0.012  # Modest stroke weight survives downsampling to normal display size.
    obj.data.align_x = align
    obj.data.align_y = "CENTER"
    obj.rotation_euler = camera.rotation_euler
    obj.location = r._camera_overlay_point(camera, x, y, -1.0)
    obj.hide_render = False
    obj["presentation_overlay"] = True
    if detail:
        obj["canonical_detail_title"] = True
    obj.data.materials.clear()
    obj.data.materials.append(overlay_material(r, materials[material]))
    bpy.context.view_layer.update()
    box = r._projected_objects_bounds(camera, [obj])
    actual = (box[3] - box[2]) * bpy.context.scene.render.resolution_y
    ratio = pixels / max(actual, 1e-8)
    obj.data.size *= ratio
    obj.data.offset *= ratio
    obj.data.body = body
    bpy.context.view_layer.update()
    box = r._projected_objects_bounds(camera, [obj])
    obj.location += camera.matrix_world.to_3x3() @ Vector((0, (y - (box[2]+box[3])/2) * r.camera_view_dimensions(camera)[1], 0))
    bpy.context.view_layer.update()
    return obj


def line(r, camera, collections, materials, name, points, detail=False, material="label_grey", thickness=1.2):
    width, _ = r.camera_view_dimensions(camera)
    obj = r.base.create_curve(name, [tuple(r._camera_overlay_point(camera, x, y, -1.02)) for x, y in points],
                              width / bpy.context.scene.render.resolution_x * thickness / 2,
                              overlay_material(r, materials[material]), collections["Labels"], resolution=1)
    obj["presentation_overlay"] = True
    obj["canonical_annotation_leader"] = True
    if detail:
        obj["canonical_detail_title"] = True
    return obj


def point(camera, r, world):
    local = camera.matrix_world.inverted() @ Vector(world)
    w, h = r.camera_view_dimensions(camera)
    return local.x / w + .5, local.y / h + .5


def place_overview_labels(r, camera_name, collections, materials):
    camera = bpy.data.objects[camera_name]
    camera.data.dof.use_dof = False
    r._remove_annotation_guides()
    for obj in list(bpy.data.objects):
        if obj.get("presentation_overlay") and not obj.get("pdb_id") and obj.name not in r.PRIMARY_CALLOUTS and obj.name != r.COMPACT_CALLOUT["object"]:
            bpy.data.objects.remove(obj, do_unlink=True)
    # Original physical comparison bars are retained, but the overview uses camera-plane copies.
    for obj in bpy.data.objects:
        if obj.name.startswith("label_scale_") or any(c.name == "Scale bars" for c in obj.users_collection):
            obj.hide_render = True
    labels = [o for o in bpy.data.objects if o.type == "FONT" and o.get("pdb_id")]
    molecular_boxes = {o.name: r._projected_object_bounds(camera, o) for o in labels}
    compact = [o for o in bpy.data.objects if o.type == "MESH" and o.name.startswith("actin mRNA compact")]
    compact_bounds = r._projected_objects_bounds(camera, compact)
    molecular_boxes["compact"] = compact_bounds
    # Exclude nucleic-acid centerlines from candidate label rectangles as well as protein silhouettes.
    path_points = []
    paths = {}
    for name in ("Canonical_DNA_source_path", "Canonical_mRNA_source_path"):
        obj = bpy.data.objects[name]
        coords = [point(camera, r, obj.matrix_world @ p.co.xyz) for s in obj.data.splines for p in s.points]
        paths[name] = coords
        path_points.extend(coords)
    occupied = []
    all_rows = []
    def record_fixed(obj):
        box = r._projected_objects_bounds(camera, [obj])
        occupied.append(box)
        return box
    width, _ = r.camera_view_dimensions(camera)
    for i, (label, length) in enumerate([("DNA 100 bp", 13.6), ("RNA 100 nt", 12.0), ("Protein 33 aa", 4.752), ("10 nm", 4.0)]):
        y = .145 - .035*i
        x1, x2 = .035, .035 + length/width
        line(r, camera, collections, materials, "overview_scale_"+str(i), [(x1,y),(x2,y)], material="scale_grey", thickness=6)
        record_fixed(text(r, camera, collections, materials, "overview_scale_label_"+str(i), label, .035+13.6/width+.013, y, 26, material="scale_grey"))
        occupied.append((x1,x2,y-.005,y+.005))

    def place(obj, body, anchor_box, pixels, preferred=None):
        text(r, camera, collections, materials, obj.name, body, .5, .5, pixels, "CENTER", source=obj)
        box = r._projected_objects_bounds(camera, [obj])
        w,h = box[1]-box[0],box[3]-box[2]
        ax,ay = (anchor_box[0]+anchor_box[1])/2,(anchor_box[2]+anchor_box[3])/2
        if obj.name in {"label_Ribosome large subunit", "label_Ribosome small subunit"}:
            ribosome_boxes=[molecular_boxes[n] for n in ("label_Ribosome large subunit", "label_Ribosome small subunit")]
            left=min(b[0] for b in ribosome_boxes); top=max(b[3] for b in ribosome_boxes)
            preferred=(left-.012-w/2,top+(.025 if obj.name=="label_Ribosome large subunit" else -.025))
        candidates = [] if preferred is None else [preferred]
        for gap in (.009,.024,.045,.070,.10):
            candidates.extend([(ax,anchor_box[3]+gap+h/2),(anchor_box[1]+gap+w/2,ay),
                               (anchor_box[0]-gap-w/2,ay),(ax,anchor_box[2]-gap-h/2)])
            for dx in (-.07,-.035,.035,.07):
                candidates.extend([(ax+dx,anchor_box[3]+gap+h/2),(ax+dx,anchor_box[2]-gap-h/2)])
        valid=[]
        for x,y in candidates:
            b=r._label_box((x,y),w,h)
            if b[0]<.025 or b[1]>.812 or b[2]<.052 or b[3]>.976: continue
            if any(r._boxes_overlap(b,other,pad=.003) for other in occupied): continue
            if any(r._boxes_overlap(b,other,pad=.003) for other in molecular_boxes.values()): continue
            if any(b[0]-.003<=px<=b[1]+.003 and b[2]-.004<=py<=b[3]+.004 for px,py in path_points): continue
            edge=(max(b[0],min(ax,b[1])),max(b[2],min(ay,b[3])))
            distance=math.hypot(edge[0]-ax,edge[1]-ay)
            valid.append((distance,(x,y),b,edge))
        if not valid: raise RuntimeError("No collision-free overview label position: "+body)
        distance,(x,y),placed,edge=next((c for c in valid if c[1]==preferred),min(valid,key=lambda c:c[0]))
        obj.location=r._camera_overlay_point(camera,x,y)
        occupied.append(placed)
        leader=line(r,camera,collections,materials,"leader_"+obj.name,[(ax,ay),edge])
        row={"object":obj.name,"text":body,"anchor_position":[ax,ay],"view_position":[x,y],
             "estimated_view_box":[w,h],"molecule_projected_bounds":list(anchor_box),"leader":leader.name,
             "displacement":distance,"placement_space":"camera_plane","projected_height_px":pixels}
        all_rows.append(row)
        return row

    # Keep the three general molecular categories in a dedicated right-hand column.
    # Parallel lines show the overlapping vertical extents of DNA and extended mRNA.
    primary_rows=[]
    for name,spec in r.PRIMARY_CALLOUTS.items():
        obj=bpy.data.objects.get(name)
        if "DNA" in name:
            ys=[p[1] for p in paths["Canonical_DNA_source_path"]]
            span=(min(ys),max(ys)); line_x=.835; y=sum(span)/2
        elif "mRNA" in name:
            ys=[p[1] for p in paths["Canonical_mRNA_source_path"]]
            span=(min(ys),max(ys)); line_x=.825; y=span[0]+.69*(span[1]-span[0])
        else:
            b=molecular_boxes["label_Actin protein"]
            y=(b[2]+b[3])/2; span=(y-.016,y+.016); line_x=.825
        obj=text(r,camera,collections,materials,name,spec["text"],.855,y,36,source=obj)
        b=record_fixed(obj)
        group=line(r,camera,collections,materials,"overview_group_line_"+name,
                   [(line_x,span[0]),(line_x,span[1])],material="scale_grey",thickness=3)
        row={"object":name,"text":spec["text"],"view_position":[.855,y],
             "span":list(span),"line":group.name,"line_x":line_x,"leader":None,
             "estimated_view_box":[b[1]-b[0],b[3]-b[2]],"placement_space":"camera_plane",
             "projected_height_px":36}
        primary_rows.append(row)
        all_rows.append(row)
    compact_obj=bpy.data.objects[r.COMPACT_CALLOUT["object"]]
    compact_row=place(compact_obj,"Compact mRNA conformation",compact_bounds,26)
    compact_row.update(overlaps_molecule=False,offset=list(r.COMPACT_CALLOUT["offset"]))
    rows=[]
    for obj in sorted(labels,key=lambda o:-molecular_boxes[o.name][3]):
        body=re.sub(r"\s*\([A-Za-z0-9]{4}\)","",obj.data.body).strip()
        if obj.name=="label_Standalone tRNA": body="tRNA"
        if obj.name=="label_Actin protein": body="Actin"
        if obj.name=="label_Ribosome large subunit": body="Ribosome (large subunit)"
        if obj.name=="label_Ribosome small subunit": body="(30S subunit)"
        rows.append(place(obj,body,molecular_boxes[obj.name],26))
    bpy.context.view_layer.update()
    visible_text = [o for o in bpy.data.objects if o.type == "FONT" and not o.hide_render]
    actual_boxes = {o.name: r._projected_objects_bounds(camera, [o]) for o in visible_text}
    collisions = [(a, b) for i, a in enumerate(actual_boxes) for b in list(actual_boxes)[i+1:]
                  if r._boxes_overlap(actual_boxes[a], actual_boxes[b], pad=.001)]
    overlaps = [(a, b) for a, box in actual_boxes.items() for b, other in molecular_boxes.items()
                if r._boxes_overlap(box, other, pad=.001)]
    clipped = [name for name, b in actual_boxes.items() if b[0]<.01 or b[1]>.99 or b[2]<.005 or b[3]>.99]
    for row in all_rows:
        b = actual_boxes[row["object"]]
        row["projected_height_px"] = (b[3]-b[2]) * bpy.context.scene.render.resolution_y
    if actual_boxes["label_Ribosome large subunit"][2] <= actual_boxes["label_Ribosome small subunit"][3]:
        raise RuntimeError("Ribosome large-subunit label must appear above the 30S label")
    if collisions or overlaps or clipped:
        raise RuntimeError(f"Overview annotation validation: collisions={collisions}, overlaps={overlaps}, clipped={clipped}")
    return {"policy":"right-side molecular categories and vertical extent lines; packed molecule labels; grey lower-left comparison scales",
            "camera":camera_name,"placed_label_count":len(all_rows),"primary_callouts":primary_rows,
            "compact_callout":compact_row,"rows":rows,"collision_pairs":collisions,"molecule_overlap_pairs":overlaps,
            "measured_text_boxes":{name:list(b) for name,b in actual_boxes.items()},
            "maximum_displacement":max(row["displacement"] for row in rows),
            "maximum_displacement_object":max(rows,key=lambda row:row["displacement"])["object"]}


def add_polymerase_dna_context(r, collections):
    source=bpy.data.objects.get("DNA B-form direct fused PyMOL-like surface")
    path=bpy.data.objects["Canonical_DNA_source_path"]
    points=[path.matrix_world @ p.co.xyz - Vector((0,0,.16)) for s in path.data.splines for p in s.points]
    tail=[points[-1]]
    length=0
    for i in range(len(points)-2,-1,-1):
        length+=(points[i+1]-points[i]).length
        if length>8: break
        tail.append(points[i])
    kd=KDTree(len(tail))
    for i,p in enumerate(tail): kd.insert(p,i)
    kd.balance()
    selected={v.index for v in source.data.vertices if kd.find(source.matrix_world @ v.co)[2]<.8}
    faces=[p for p in source.data.polygons if all(v in selected for v in p.vertices)]
    used=sorted({v for p in faces for v in p.vertices}); remap={v:i for i,v in enumerate(used)}
    mesh=bpy.data.meshes.new("detail_polymerase_dna_context_mesh")
    mesh.from_pydata([tuple(source.matrix_world @ source.data.vertices[i].co) for i in used],[],[[remap[v] for v in p.vertices] for p in faces])
    for mat in source.data.materials: mesh.materials.append(mat)
    for dst,src in zip(mesh.polygons,faces): dst.material_index=src.material_index; dst.use_smooth=src.use_smooth
    obj=bpy.data.objects.new("detail_polymerase_dna_context",mesh)
    collections["Detail context"].objects.link(obj)
    obj.hide_render=True
    obj.hide_viewport=True
    obj["presentation_context_copy"]=True
    if not faces: raise RuntimeError("Polymerase DNA context is empty")


def render_closeup(r,key,shared_camera,collections,materials,units):
    for obj in list(bpy.data.objects):
        if obj.get("canonical_detail_title"): bpy.data.objects.remove(obj,do_unlink=True)
    r.apply_focus_visibility(key)
    if key=="polymerase_gene_end": bpy.data.objects["detail_polymerase_dna_context"].hide_render=False
    camera=shared_camera.copy();camera.data=shared_camera.data.copy()
    camera.name=shared_camera.name+"_closeup"
    bpy.context.scene.collection.objects.link(camera)
    camera.data.dof.use_dof=False
    bpy.context.scene.camera=camera
    bounds=r.camera_space_renderable_bounds(camera.name)
    aspect=bpy.context.scene.render.resolution_x/bpy.context.scene.render.resolution_y
    required=max(bounds["width"]/.75,bounds["height"]*aspect/.65)
    # The 1-2-5 series has gaps; expand framing only when needed to satisfy the scale-bar interval.
    candidates=[]
    for power in range(-3,6):
        for mantissa in (1,2,5):
            nm=mantissa*10.0**power; mm=nm*units["nm_to_mm"]
            scale=max(required,mm/.25)
            if .15-1e-8<=mm/scale<=.25+1e-8: candidates.append((scale,nm,mm))
    scale,nm,mm=min(candidates)
    center=Vector(((bounds["min_x"]+bounds["max_x"])/2,(bounds["min_y"]+bounds["max_y"])/2,0))
    camera.location+=camera.matrix_world.to_3x3() @ center
    camera.data.ortho_scale=scale
    bpy.context.view_layer.update()
    environment=r.add_canonical_beauty_environment(camera.name,collections)
    title,caption=DETAILS[key]
    text(r,camera,collections,materials,"closeup_title_"+key,title,.07,.92,52,detail=True)
    text(r,camera,collections,materials,"closeup_caption_"+key,caption,.07,.855,26,detail=True)
    bar=line(r,camera,collections,materials,"closeup_scale_"+key,[(.07,.08),(.07+mm/scale,.08)],detail=True,material="scale_grey",thickness=6)
    text(r,camera,collections,materials,"closeup_scale_label_"+key,f"{nm:g} nm",.07,.115,28,detail=True,material="scale_grey")
    output=r.DETAIL_PREVIEWS[key].with_name(r.DETAIL_PREVIEWS[key].stem+"_closeup.png")
    bpy.context.scene.render.filepath=str(output)
    bpy.ops.render.render(write_still=True)
    endpoints=[bar.matrix_world @ p.co.xyz for p in bar.data.splines[0].points]
    measured=(endpoints[1]-endpoints[0]).length
    if not math.isclose(measured,mm,rel_tol=1e-5): raise RuntimeError("Incorrect close-up scale bar")
    subject=[o for o in bpy.data.objects if o.type in {"MESH","CURVE"} and not o.hide_render and not o.get("presentation_overlay") and not o.get("canonical_beauty_object")]
    b=r._projected_objects_bounds(camera,subject)
    if b[0]<.125-1e-5 or b[1]>.875+1e-5 or b[2]<.175-1e-5 or b[3]>.825+1e-5: raise RuntimeError("Close-up subject outside reserved frame")
    return {"camera":camera.name,"output":str(output),"title":title,"caption":caption,
            "subject_projected_bounds":list(b),"ortho_scale_mm":scale,"environment":environment,
            "scale_bar":{"length_nm":nm,"length_mm":mm,"measured_mm":measured,"fraction_of_width":mm/scale,"camera_plane":True}}


def refresh_documentation(r,closeups):
    """Use Blender's existing image I/O; no Pillow dependency in the build."""
    mapping={"full_overview":"overview","p53_dna":"p53-dna","nucleosome_loop":"nucleosome-loop",
             "polymerase_gene_end":"transcription-end","ribosome_trna":"translation","actin_product":"actin","cas9_dna":"cas9-dna"}
    settings=bpy.context.scene.render.image_settings
    old=(settings.file_format,settings.color_mode,settings.quality)
    try:
        settings.file_format="JPEG";settings.color_mode="RGB";settings.quality=95
        for key,name in mapping.items():
            source=r.DETAIL_PREVIEWS[key] if key=="full_overview" else closeups[key]["output"]
            img=bpy.data.images.load(str(source),check_existing=False)
            try:
                # Force the source pixels to load before changing output format or path.
                # Reassigning filepath_raw can otherwise reload an old destination image.
                _ = img.pixels[0]
                img.file_format="JPEG"
                img.save(filepath=str(r.ROOT/"docs/images"/(name+".jpg")), quality=95, save_copy=True)
            finally: bpy.data.images.remove(img)
    finally:
        settings.file_format,settings.color_mode,settings.quality=old
