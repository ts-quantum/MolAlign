# created with MolAlign 1.2 by (C) 2026 Dr. Tobias Schulz
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
# 
# B) CLEANUP: Search 'Camera Node' and other empty objects in Outliner -> Select all (A) 
#    -> Right Click -> 'Delete Hierarchy'. This removes Cameras but keeps Meshes.
# ==============================================================================
 
import bpy, os, re

path_to_glb = "/Users/user/python/MolAlign/examples/ex1/bld-multi"
extension = ".glb"

# 1. Initial Cleanup - Protects "MASTER" oder "DUMMY" etc.
protected_keywords = ["Camera", "Plane", "Cylinder", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY"]

# Collect elements to delete
to_delete = []
for obj in bpy.data.objects:
    is_protected = any(p.upper() in obj.name.upper() for p in protected_keywords)
    if not is_protected:
        to_delete.append(obj)

for obj in to_delete:
    try:
        bpy.data.objects.remove(obj, do_unlink=True)
    except:
        pass

# 2. Controller Setup
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.scene.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 3. File Discovery
files = sorted([f for f in os.listdir(path_to_glb) if f.endswith(extension)])

# 4. Import Loop
for i, filename in enumerate(files):
    filepath = os.path.join(path_to_glb, filename)
    bpy.ops.import_scene.gltf(filepath=filepath)
    
    new_objs = [o for o in bpy.context.selected_objects if o.type == 'MESH']
    current_frame = i + 1 

    for obj in new_objs:
        # --- POSITIONING ---
        old_matrix = obj.matrix_world.copy()
        obj.parent = cntrl
        obj.matrix_world = old_matrix

        # --- MASTER MATERIAL LINK ---
        # Wir suchen das Master-Material
        mat = bpy.data.materials.get("MASTER_Molecule")
        if mat:
            obj.data.materials.clear()
            obj.data.materials.append(mat)
            obj.color = mat.diffuse_color

        # --- ANIMATION (Visibility via Scale) ---
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame - 1)
        
        # Show Geometry
        obj.scale = (1, 1, 1) if len(obj.data.vertices) > 1 else (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame)
        
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame + 1)
        
        # Constant Interpolation
        if obj.animation_data and obj.animation_data.action:
            action = obj.animation_data.action
            if hasattr(action, "fcurves"):
                for fc in action.fcurves:
                    for kp in fc.keyframe_points:
                        kp.interpolation = 'CONSTANT'
    
    # Cleanup leere Nodes aus dem GLTF Import (Nodes/Empties)
    for o in bpy.context.selected_objects:
        if o.type != 'MESH':
            # Nur löschen, wenn es nicht unser geschützter Controller ist
            if not any(p.upper() in o.name.upper() for p in protected_keywords):
                bpy.data.objects.remove(o, do_unlink=True)

# 5. Final Scene Sync
bpy.context.scene.frame_start = 1
bpy.context.scene.frame_end = len(files)
bpy.context.scene.frame_set(1)
print(f"Sync complete. Master-Dummies preserved. {len(files)} frames processed.")
