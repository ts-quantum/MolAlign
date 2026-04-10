# created with MolAlign 1.1 by (C) 2026 Dr. Tobias Schulz
import bpy
import os

# --- Settings ---
path_to_glb = "/Users/user/python/MolAlign/examples/ex3/bld/" # adapt path to *glb files!
extension = ".glb"

# 1. Empty current Scene (recommended)
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()

# Find and sort all Files
files = sorted([f for f in os.listdir(path_to_glb) if f.endswith(extension)])

for i, filename in enumerate(files):
    filepath = os.path.join(path_to_glb, filename)
    
    # 2. IMPORT (imported objects are automatically selected)
    bpy.ops.import_scene.gltf(filepath=filepath)
    imported_objs = bpy.context.selected_objects
    
    current_frame = i + 1 # Blender starts with Frame 1

    for obj in imported_objs:
        # --- Set KEYFRAMES ---
        
        # A) Predecessor Frame: invisible
        if current_frame > 1:
            obj.hide_viewport = True
            obj.hide_render = True
            obj.keyframe_insert(data_path="hide_viewport", frame=current_frame - 1)
            obj.keyframe_insert(data_path="hide_render", frame=current_frame - 1)
        
        # B) Current Frame: visible
        obj.hide_viewport = False
        obj.hide_render = False
        obj.keyframe_insert(data_path="hide_viewport", frame=current_frame)
        obj.keyframe_insert(data_path="hide_render", frame=current_frame)
        
        # C) Successor Frame: invisible
        obj.hide_viewport = True
        obj.hide_render = True
        obj.keyframe_insert(data_path="hide_viewport", frame=current_frame + 1)
        obj.keyframe_insert(data_path="hide_render", frame=current_frame + 1)

# Adapt Timeline
bpy.context.scene.frame_end = len(files)
print(f"Finished! {len(files)} Frames processed.")
