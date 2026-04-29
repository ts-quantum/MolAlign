# created with MolAlign 1.3 (C) 2026 by Dr. Tobias Schulz
# ==============================================================================
# USER GUIDE for MolVista Blender Animation
# ==============================================================================
# A) Run This Script in Blender
# 1. GLOBAL VISUAL CONTROL: 
#    This script links all imported meshes to "MASTER" materials in your template.
#    Edit these materials in the 'Material Properties' tab to update ALL frames:
#    - 'MASTER_Molecule'  -> Controls atoms and bonds (mol_***)
#
# 2. RETAINING COLORS (CPK):
#    'MASTER_Molecule' use 'Vertex Colors' (Color Attributes).
#    In the Shader Editor, ensure a 'Color Attribute' node is connected to the 
#    'Base Color' and 'Emission Color' of the Principled BSDF.
#
# 3. POSITIONING:
#    Select the 'TRAJECTORY_CONTROL' (Empty) to move, rotate, or scale the 
#    entire animation sequence simultaneously over your scene.
# ==============================================================================
 
import bpy, re, os

# --- Helper: Natural Sort ---
def natural_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', text)]

# 1. Protection & Controller Setup
protected_keywords = ["Camera", "Plane", "Cylinder", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY"]

if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.collection.objects.link(cntrl)
    cntrl.empty_display_type = 'ARROWS'
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. THE DUMMY KILLER: Hide all template dummies from the start
for obj in bpy.data.objects:
    if "DUMMY" in obj.name.upper():
        obj.scale = (0, 0, 0)
        obj.hide_render = True

# 3. Collect and Sort all imported Meshes
all_meshes = [o for o in bpy.data.objects if o.type == 'MESH' and not any(p.upper() in o.name.upper() for p in protected_keywords)]
all_meshes.sort(key=lambda o: natural_key(o.name))

# 4. Master Material Link & Animation
master_mat = bpy.data.materials.get("MASTER_Molecule")

for i, obj in enumerate(all_meshes):
    target_frame = i + 1 
    
    # Parenting & Transform-Fix
    obj.parent = cntrl
    obj.matrix_parent_inverse = cntrl.matrix_world.inverted()
    
    if master_mat:
        obj.data.materials.clear()
        obj.data.materials.append(master_mat)

    # --- Animation (Visibility via Scale) ---
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=target_frame - 1)
    
    # Show Geometry
    obj.scale = (1, 1, 1) if len(obj.data.vertices) > 0 else (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=target_frame)
    
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=target_frame + 1)
    
    # Constant interpolation
    if obj.animation_data and obj.animation_data.action:
        act = obj.animation_data.action
        if hasattr(act, "fcurves"):
            for fc in act.fcurves:
                for kp in fc.keyframe_points:
                    kp.interpolation = 'CONSTANT'

# 5. Scene Finalize
if all_meshes:
    bpy.context.scene.frame_end = len(all_meshes)
    bpy.context.scene.frame_set(1)
    print(f"One-File Sync complete: {len(all_meshes)} Frames ready. Dummies hidden.")
