# MolAlign 
import sys, os, re
import platform, subprocess, psutil, click
from tqdm import tqdm
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QColor
import pyvista as pv
from collections import defaultdict
from pyscf import data, lib
from pyscf.data import elements
import numpy as np
from pathlib import Path
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import math
from PySide6.QtCore import QThread, Signal
from PySide6.QtCore import QCoreApplication, QEventLoop

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

class MoleculeData:
    def __init__(self, *, name=None, atom_points=None, atom_types=None, energies=None):
        self.name = name
        self.atom_points = atom_points
        self.atom_types = atom_types
        self.energies = energies or [] 

    @classmethod
    def from_xyz(cls, filepath):
        base = os.path.basename(filepath)
        atom_types = []
        atom_points = []
        energies = []
        
        with open(filepath, 'r') as f:
            lines = f.readlines()

        cursor = 0
        while cursor < len(lines):
            line = lines[cursor].strip()
            if not line:
                cursor += 1
                continue
                
            try:
                num_atoms = int(line)
                # --- extract energy from second line
                comment = lines[cursor + 1].strip()
                energy_match = re.search(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", comment)
                energy_val = float(energy_match.group(0)) if energy_match else None
                energies.append(energy_val)

                # atom_types and coordinates
                block = lines[cursor + 2 : cursor + 2 + num_atoms]
                types = []
                coords = []
                for entry in block:
                    parts = entry.split()
                    if len(parts) < 4: continue
                    symbol = parts[0]
                    try: 
                        at_num = data.elements.charge(symbol) # map symbol to atomic number
                    except KeyError:
                        at_num = 0 
                    types.append(at_num) 
                    coords.append([float(x) for x in parts[1:4]])
                
                atom_types.append(np.array(types))
                atom_points.append(np.array(coords))
                cursor += num_atoms + 2
            except (ValueError, IndexError):
                break

        return cls(name=base, atom_types=atom_types, atom_points=atom_points, energies=energies)

def get_optimal_cores():
    # 1. macOS (Apple Silicon Check)
    if platform.system() == "Darwin":
        try:
            # try to find performance cores
            output = subprocess.check_output(['sysctl', '-n', 'hw.perflevel0.physicalcpu'], 
                                           stderr=subprocess.DEVNULL) # Fehler im Hintergrund halten
            return max(1, int(output.strip()))
        except Exception:
            # Fallback
            pass

    # 2. Fallback for Win and Linux
    cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
    return max(1, cores - 1)

def reverse(data_):
        data_.atom_types = data_.atom_types[::-1]
        data_.atom_points = data_.atom_points[::-1]
        data_.energies = data_.energies[::-1]
        return data_

#### Alignment #####
def get_euler_angles(R):
        """
        Extracts Euler angles (XYZ) from a 3x3 rotation matrix R.
        Returns angles in degrees.
        """
        # Calculation of Beta (Y-Axis)
        sy = math.sqrt(R[0,0] * R[0,0] + R[1,0] * R[1,0])
        singular = sy < 1e-6

        if not singular:
            x = math.atan2(R[2,1], R[2,2])
            y = math.atan2(-R[2,0], sy)
            z = math.atan2(R[1,0], R[0,0])
        else:
            # Special Case (Gimbal Lock)
            x = math.atan2(-R[1,2], R[1,1])
            y = math.atan2(-R[2,0], sy)
            z = 0

        # Convert to "°" "Grad"
        return [math.degrees(x), math.degrees(y), math.degrees(z)]

def get_min_rmsd_kabsch(coords_A, coords_B):
        """
        Does temporary Kabsch Alignment and returns min.
        RMSD (Root Mean Square Deviation)
        """
        # 1. Center both point sets
        centroid_A = np.mean(coords_A, axis=0)
        centroid_B = np.mean(coords_B, axis=0)
        A_centered = coords_A - centroid_A
        B_centered = coords_B - centroid_B

        # 2. calculate Kabsch-Rotation 
        H = A_centered.T @ B_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Special Case, prevent mirroring
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # 3. Transform and calculate RMSD
        B_rotated = B_centered @ R
        rmsd = np.sqrt(np.mean(np.sum((A_centered - B_rotated)**2, axis=1)))
        return rmsd

def align_structures(ref_coords, target_coords):
    """Aligns target_coords with ref_coords"""
    # 1.  Centering (Translation to the center of mass).
    centroid_ref = np.mean(ref_coords, axis=0)
    centroid_target = np.mean(target_coords, axis=0)
    
    ref_centered = ref_coords - centroid_ref
    target_centered = target_coords - centroid_target

    # 2. Covariance matrix calculation
    h = target_centered.T @ ref_centered

    # 3. Singular Value Decomposition (SVD)
    u, s, vt = np.linalg.svd(h)
    
    # 4. Calculate Rotation matrix
    r = vt.T @ u.T

    # Special Case: correct for mirroring (Determinant must be 1)
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T

    # 5. Calculate transformed coordinates
    aligned_coords = (target_centered @ r.T) + centroid_ref
    
    return aligned_coords, r, centroid_ref,  centroid_target

def transform_trajectory(trajectory_coords, r_matrix, ref_centroid, target_centroid):
    """Applies Transformation on whole trajectory"""
    transformed_traj = []
    for coords in trajectory_coords:
        # 1. Shift to origin (relative to the original target center of gravity).
        centered = np.array(coords) - target_centroid
        # 2. Rotate and shift to the new reference center of gravity
        new_coords = (centered @ r_matrix.T) + ref_centroid
        transformed_traj.append(new_coords.tolist())
    return transformed_traj

def find_best_flip_strategy(data0, data1):
        """
        Checks all 4 possible start/end combinations to find the best 
        chronological alignment, including automatic atom mapping.
        """
        # Define the 4 possible connection points
        combinations = [
            ("normal",    data0.atom_points[-1], data1.atom_points[0]),
            ("flip_1",    data0.atom_points[-1], data1.atom_points[-1]),
            ("flip_0",    data0.atom_points[0],  data1.atom_points[0]),
            ("flip_both", data0.atom_points[0],  data1.atom_points[-1])
        ]
        
        best_rmsd = float('inf')
        best_case = "normal"
        best_mapping = None

        for case, c0, c1 in combinations:
            # Perform atom mapping for the current geometry pair
            # This handles different atom orderings automatically
            mapping = find_mapping(c0, c1, data0.atom_types[0], data1.atom_types[0])
            
            # Calculate RMSD after a temporary Kabsch alignment to verify the fit
            _, R, cent0, cent1 = align_structures(c0, c1[mapping])
            
            # Transform the target coordinates to the reference frame for RMSD check
            aligned_c1 = (c1[mapping] - cent1) @ R + cent0
            rmsd = np.sqrt(np.mean(np.sum((c0 - aligned_c1)**2, axis=1)))
            
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_case = case
                best_mapping = mapping

        print(f"Best alignment strategy: {best_case} (RMSD: {best_rmsd:.4f} Å)", "info")
        return best_case, best_mapping

def find_mapping(coords_ref, coords_target, types_ref, types_target):
        """
        Finds the optimal assignment (permutation) between two point clouds.
        Based on internal distance fingerprints (rotation-invariant).

        Args:
            coords_ref (np.array): (N, 3) coordinates of the reference framework.
            coords_target (np.array): (N, 3) coordinates of the framework to be mapped.
            types_ref (list/np.array): Atom types (elements) of the reference framework.
            types_target (list/np.array): Atom types of the target framework.
            
        Returns:
            np.array: Index vector to be applied to coords_target.
        """
        # 1. Calculate distance matrices (distances between every atom pair within each system)
        dists_ref = cdist(coords_ref, coords_ref)
        dists_target = cdist(coords_target, coords_target)
        
        # 2. Create fingerprints: Sorted distances for each atom
        # This makes the signature invariant to both indexing and rotation
        sig_ref = np.sort(dists_ref, axis=1)
        sig_target = np.sort(dists_target, axis=1)
        
        # 3. Calculate cost matrix: How similar are the fingerprints?
        # Compare signature of atom i (Ref) with atom j (Target)
        cost_matrix = cdist(sig_ref, sig_target, metric='euclidean')
        
        # 4. Implement type constraints (element validation)
        # If elements do not match, set the cost extremely high to penalize the pairing
        for i in range(len(types_ref)):
            for j in range(len(types_target)):
                if str(types_ref[i]) != str(types_target[j]):
                    cost_matrix[i, j] += 1e6  # Prevents incorrect mapping across different elements

        # 5. Apply the Hungarian Algorithm (Kuhn-Munkres)
        # Finds the pairing that minimizes the total cost (structural differences)
        _, col_ind = linear_sum_assignment(cost_matrix)
        
        return col_ind

def chain_alignment(all_datasets):
    combined_data = all_datasets[0]
    alignment_rmsds = []
    log = []
    # check if first segment needs to be flipped (outside main loop)
    best_case, mapping = find_best_flip_strategy(all_datasets[0], all_datasets[1])
    
    # mapping is the result from find_mapping
    # create a reference array [0, 1, 2, ..., N-1]
    identity_mapping = np.arange(len(mapping))
    # Check if the found mapping deviates from the original order
    log.append("pre-check with first and second segment")
    if not np.array_equal(mapping, identity_mapping):
        # Count how many atoms actually switched positions
        swapped_count = np.sum(mapping != identity_mapping)
        log.append(f"Atom re-ordering detected: {swapped_count} atoms rearranged to match reference.")
    else:
        log.append("Atom order is already consistent between segments.")

    match best_case:
        case "flip_1":
            all_datasets[1] = reverse(all_datasets[1])
            log.append("second segment reversed")
        case "flip_0":
            all_datasets[0] = reverse(all_datasets[0])
            log.append("first segment reversed")
        case "flip_both":
            all_datasets[0] = reverse(all_datasets[0])
            all_datasets[1] = reverse(all_datasets[1])
            log.append("first and second segment reversed")

    p_bar = tqdm(range(1, len(all_datasets)),desc="Aligning trajectories")
    log.append("Loop over all segments")
    for i in p_bar:
        log.append(f"Step {i}")
        current_name = getattr(all_datasets[i], 'name', f"Segment_{i}")
        p_bar.set_description(f"Processing {os.path.basename(current_name)}")

        next_data = all_datasets[i]
        
        # Referenz: Letzter Frame der BISHERIGEN Kette
        ref_frame = combined_data.atom_points[-1]
        
        # Define the 2 possible connection points
        combinations = [
            ("normal",    ref_frame, next_data.atom_points[0]),
            ("flip_1",    ref_frame, next_data.atom_points[-1])
        ]
        best_rmsd = float('inf')
        best_case = "normal"
        best_mapping = None

        for case, c0, c1 in combinations:
            # Perform atom mapping for the current geometry pair
            # This handles different atom orderings automatically
            mapping = find_mapping(c0, c1, ref_frame, next_data.atom_types[0])
            
            # Calculate RMSD after a temporary Kabsch alignment to verify the fit
            _, R, cent0, cent1 = align_structures(c0, c1[mapping])
            
            # Transform the target coordinates to the reference frame for RMSD check
            aligned_c1 = (c1[mapping] - cent1) @ R + cent0
            rmsd = np.sqrt(np.mean(np.sum((c0 - aligned_c1)**2, axis=1)))
            
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_case = case
                best_mapping = mapping

        log.append(f"Best alignment strategy: {best_case} (RMSD: {best_rmsd:.4f} Å)")
        
        identity_mapping = np.arange(len(best_mapping))
        # Check if the found mapping deviates from the original order
        if not np.array_equal(best_mapping, identity_mapping):
            # Count how many atoms actually switched positions
            swapped_count = np.sum(best_mapping != identity_mapping)
            log.append(f"Atom re-ordering detected: {swapped_count} atoms rearranged to match reference.")
        else:
            log.append("Atom order is already consistent between segments.")

        match best_case:
            case "flip_1":
                next_data = reverse(next_data)

        # PERMANENTLY re-index next_data to match the reference order
        next_data.atom_points = [p[best_mapping] for p in next_data.atom_points]
        next_data.atom_types = [t[best_mapping] for t in next_data.atom_types] # If it's a list of lists
                
        target_frame = next_data.atom_points[0] # after potential flip, always the first point

        # Jetzt das (ggf. geflippte) Teil alignen und mergen
  
        _, r, centroid_ref,  centroid_target = align_structures(ref_frame, target_frame)

        angles = get_euler_angles(r)
        log.append(f"Alignment step {i}, dataset {next_data.name}")
        log.append(f"Kabsch Alignment rotation: X={angles[0]:.2f}°, Y={angles[1]:.2f}°, Z={angles[2]:.2f}°")
        log.append(f"Kabsch ref. translation x={centroid_ref[0]:.2f} y={centroid_ref[1]:.2f} z={centroid_ref[2]:.2f}")
        log.append(f"Kabsch target translation x={centroid_target[0]:.2f} y={centroid_target[1]:.2f} z={centroid_target[2]:.2f}")

        aligned_coords = transform_trajectory(next_data.atom_points, r, centroid_ref,  centroid_target)
        
        step_rmsd = get_min_rmsd_kabsch(ref_frame, aligned_coords[0])
        alignment_rmsds.append(step_rmsd)

        energy_offset = combined_data.energies[-1] - next_data.energies[0]
        shifted_energies = [e + energy_offset for e in next_data.energies]

        combined_points = list(combined_data.atom_points) + list(aligned_coords[1:])
        combined_types = list(combined_data.atom_types) + list(next_data.atom_types[1:])
        combined_energies = list(combined_data.energies) + list(shifted_energies[1:])

        combined_data = MoleculeData(
            atom_points=combined_points,
            atom_types=combined_types,
            energies=combined_energies
        )

    return combined_data, alignment_rmsds, log

#### export XYZ ####
def export_xyz(data_,path):
    n_steps = len(data_.energies)
    n_atoms = len(data_.atom_types[0])
    with open(path, 'w') as f:
        for j in range(n_steps):
            f.write(f"{n_atoms}\n")
            energy = data_.energies[j] if data_.energies[j] is not None else 0.0
            f.write(f"Energy: {energy:16.10f}\n")
            for i in range(n_atoms):
                atom_coord = data_.atom_points[j][i]
                at_num = data_.atom_types[j][i]
                symbol = elements.ELEMENTS[at_num] 
                f.write(f"{symbol:2} {atom_coord[0]:12.8f} {atom_coord[1]:12.8f} {atom_coord[2]:12.8f}\n")
    log = f"xyz data {data_.name} written as: {os.path.basename(path)}"
    return log

def create_split_template_orca(xyz_path, ver_no):
    """
    Creates a pre-configured ORCA Python script (_split.py) in the same directory
    as the exported XYZ file for further batch processing.
    """
    if not xyz_path:
        return

    work_dir = os.path.dirname(xyz_path)
    xyz_filename = os.path.basename(xyz_path)
    # Name des Hilfsskripts basierend auf der XYZ-Datei
    script_name = os.path.join(work_dir, f"{os.path.splitext(xyz_filename)[0]}_split.py")

    # Das Template als String (mit Platzhaltern)
    template_content = f"""# --- CONFIGURATION: ADJUST BEFORE RUNNING ---
##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script will create a series of ORCA Input files from 
##### a given IRC Trajectory along with a "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
### Edit method, basis set, etc. here
charge =''  # e.g. '0'	
mult ='' 	# e.g. '1'		
basis =''   # e.g. 'def2-SVP'
method =''  # e.g. 'HF', 'B3LYP', 'PBE0' or other XC
mem = '2000'
nproc = '8'
### adapt to your specific environment
ORCA_EXE = "/usr/local/orca_6_1_0/orca"
ORCA_2MKL = "/usr/local/orca_6_1_0/orca_2mkl"
# pre-configured by MolAlign
inp_file = '{xyz_filename}.xyz' 	

import os, sys

def run_split():
    with open(inp_file, 'r') as f:
        lines = f.readlines()
    try:
        num_atoms = int(lines[0].strip())
    except:
        print("Error: Invalid XYZ format.")
        return
    block_size = num_atoms + 2
    steps = len(lines) // block_size
    base = os.path.splitext(inp_file)[0]
    wrapper_name = f"run_{{base}}_batch.sh"
    with open(wrapper_name, 'w') as w:
        w.write("#!/bin/bash\\n\\n")
        for i in range(steps):
            label = f"{{base}}_{{i:03d}}"
            inp = f"{{label}}.inp"
            out = f"{{label}}.out"
            xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
            xyz_body = "".join(xyz_list).rstrip()
            body = f\"\"\"! {{method}} {{basis}}

%maxcore {{mem}}
%pal nprocs {{nproc}} end

* xyz {{charge}} {{mult}}
{{xyz_body}}
*
\"\"\"
            #write input file
            with open(inp, 'w') as f_inp:
                f_inp.writelines(body)
			# Commands for shell script
            w.write(f"{{ORCA_EXE}} {{inp}} > {{out}} 2>&1 && \\\\\\n")
            w.write(f"{{ORCA_2MKL}} {{label}} -molden && \\\\\\n")
            w.write(f"mv {{label}}.molden.input {{label}}.molden && \\\\\\n")
            w.write(f"echo 'Frame {{i:03d}} finished.'\\n\\n")
    os.chmod(wrapper_name, 0o755)
    print(f"Done. Created {{steps}} inputs and shell script: {{wrapper_name}}")	

if __name__ == "__main__": 
    missing = [name for name, val in [("Charge", charge), ("Mult", mult), ("Basis", basis), ("Method", method)] if not str(val)]
    if missing:
        print(f"Error: Input missing: {{', '.join(missing)}}")
        sys.exit()
    run_split()
"""
    with open(script_name, 'w', encoding='utf-8') as f:
        f.write(template_content)
    return script_name

def create_split_template_nw(xyz_path, ver_no):
    """
    Creates a pre-configured NWChem Python script (_split.py) in the same directory
    as the exported XYZ file for further batch processing.
    """
    if not xyz_path:
        return

    work_dir = os.path.dirname(xyz_path)
    xyz_filename = os.path.basename(xyz_path)
    # Name des Hilfsskripts basierend auf der XYZ-Datei
    script_name = os.path.join(work_dir, f"{os.path.splitext(xyz_filename)[0]}_split.py")

    # Das Template als String (mit Platzhaltern)
    template_content = f"""# --- CONFIGURATION: ADJUST BEFORE RUNNING ---
##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script will create a series of NWChem Input files from 
##### a given IRC Trajectory along with a "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
### Edit method, basis set, etc. here
charge = ''	# e.g. '0'	
mult = '' # e.g. '1'		
basis = '' # e.g. '6-31G'
method_type = '' # or 'scf' or 'dft'
dft_method = '' # e.g. 'xc b3lyp'
method_details =f\"\"\"
 {{dft_method}}
 mult {{mult}}   
\"\"\"
### directories will be created by the "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
work_dir = 'calc'
scr_dir = 'scr'
### adapt to your specific environment
nw_cmd = 'mpiexec -np 8 nwchem'
# pre-configured by MolAlign
inp_file = '{xyz_filename}.xyz'	

import os, sys 

def run_split():
    with open(inp_file, 'r') as f:
        lines = f.readlines()
    try:
        num_atoms = int(lines[0].strip())
    except:
        print("Error: Invalid XYZ format.")
        return
    
    block_size = num_atoms + 2
    steps = len(lines) // block_size
    base = os.path.splitext(inp_file)[0]
    
    wrapper_name = f"run_{{base}}_batch.sh"
    with open(wrapper_name, 'w') as w:
        w.write("#!/bin/bash\\n\\n")
        w.write(f"mkdir -p ./{{work_dir}}\\n\\n")
        w.write(f"mkdir -p ./{{scr_dir}}\\n\\n")
        for i in range(steps):
            xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
            xyz_body = "".join(xyz_list).rstrip()
            body = f\"\"\"Title {{base}}_{{i:03d}}
echo
start {{base}}_{{i:03d}} 

memory 10000 mb

permanent_dir ./{{work_dir}}
scratch_dir   ./{{scr_dir}}

charge {{charge}}
geometry
{{xyz_body}} 
end

basis
 * library {{basis}}
end

{{method_type}}
{{method_details}}
end

task {{method_type}} energy ignore

property
 moldenfile
 molden_norm none
end

task {{method_type}} property ignore 
\"\"\"
            inp = f"{{base}}_{{i:03d}}.nw"
            out = f"{{base}}_{{i:03d}}.out"
            with open(inp, 'w') as f_inp:
                f_inp.writelines(body)
            
            w.write(f"{{nw_cmd}} {{inp}} > {{out}} 2>&1 && \\\\\\n" )
            w.write(f"mv ./{{work_dir}}/{{base}}_{{i:03d}}.molden ./ 2>/dev/null && \\\\\\n")
            w.write(f"echo 'Frame {{i:03d}} finished.'\\n\\n")
    os.chmod(wrapper_name, 0o755)
    print(f"Done. Created {{steps}} inputs and shell script: {{wrapper_name}}")	

if __name__ == "__main__":
    missing = [name for name, val in [("Charge", charge), ("Basis", basis), ("Method", method_type)] if not str(val)]
    if missing:
        print(f"Error: Input missing: {{', '.join(missing)}}")
        sys.exit()
    run_split()

"""
    with open(script_name, 'w', encoding='utf-8') as f:
        f.write(template_content)
    return script_name

def create_split_template_psi4(xyz_path, ver_no):
    """
    Creates a pre-configured PSI4-Python script (_split.py) in the same directory
    as the exported XYZ file for further batch processing.
    """
    if not xyz_path:
        return

    work_dir = os.path.dirname(xyz_path)
    xyz_filename = os.path.basename(xyz_path)
    # Name des Hilfsskripts basierend auf der XYZ-Datei
    script_name = os.path.join(work_dir, f"{os.path.splitext(xyz_filename)[0]}_split.py")

    # Das Template als String (mit Platzhaltern)
    template_content = f"""# --- CONFIGURATION: ADJUST BEFORE RUNNING ---
##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script to run a Psi-4 single point calculation
##### for each point on a given IRC Trajectory and to produce a 
##### .molden and .fchk file for each point
##### must be run inside a Psi-Conda environment
### Edit method, basis set, etc. here
charge = '' # e.g. '0'
mult = '' # e.g. '1'
basis = '' # e.g. 'cc-pVDZ'
ref = ''  # e.g. 'rhf', 'rks', 'uhf', 'uks', ...
method = '' # e.g. 'scf', 'pbe0', 'b3lyp', ...
mem = '1000mb'
nproc = 8
#pre-configured by MolAlign
xyz_input = '{xyz_filename}.xyz' 

import os, sys
import psi4
psi4.set_memory(mem)
psi4.core.set_num_threads(nproc)

def run_split():
	with open(xyz_input, 'r') as f:
		lines = f.readlines()

	try:
		num_atoms = int(lines[0].strip())
	except: 
		print("Error: Invalid XYZ format.")
		return
	
	block_size = num_atoms + 2
	steps = len(lines) // block_size
	base = os.path.splitext(xyz_input)[0]

	psi4.set_options({{
	     'basis': basis,
	     'scf_type': 'df',
	     'reference': ref,  
	     'g_convergence': 'gau_tight',
	     'guess' : 'read'
	}})
	initial = "".join(lines[2 : block_size])
	mol= psi4.geometry(f\"\"\"
		{{charge}} {{mult}}
		{{initial}}
		units angstrom
		no_reorient
		no_com
		\"\"\")
	
	last_wfn = None

	for i in range(steps):
		xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
		xyz_body = "".join(xyz_list)
		temp_mol = psi4.geometry(f"units angstrom\\n{{xyz_body}}")
		mol.set_geometry(temp_mol.geometry())
		#energy calculation
		if last_wfn is not None:
			energy, wfn = psi4.energy(method, molecule=mol, return_wfn=True, restart_wfn=last_wfn)
		else:
			energy, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
		last_wfn = wfn
		
		psi4.fchk(wfn, f'{{base}}_{{i:03d}}.fchk')
		psi4.molden(wfn, f'{{base}}_{{i:03d}}.molden')
	    
		psi4.core.clean()
	
if __name__ == "__main__":
    missing = [name for name, val in [("Charge", charge), ("Mult", mult), ("Basis", basis), ("Ref", ref), ("Method", method)] if not str(val)]
    if missing:
        print(f"Error: Input missing: {{', '.join(missing)}}")
        sys.exit()
    run_split()
"""
    with open(script_name, 'w', encoding='utf-8') as f:
        f.write(template_content)
    return script_name

#### export POV ####
def export_pov_header(ver_no, length, filename="test.inc", object_name="name"):

    with open(filename, 'w') as f:
        f.write(f"""\
// created with MolAlign {ver_no} (C)2026 Dr. Tobias Schulz
// IRC Trajectory as Array of Molecules
// #include "{object_name}.inc" into povray
// ---- use "object{{{object_name}[i]}}" in code
//declare molecule object array
#declare {object_name} = array[{length+1}];
//
// ---- Atom and Bond Section
//transparency
#declare trans_bd = 0;
#declare trans_atom = 0;
//atom radius
#declare atom_rad_h = 0.24;
#declare atom_rad_2 = 0.35;
#declare atom_rad_3 = 0.42;
#declare atom_rad_def = 0.5;
#declare bond_rad = 0.08;
// predefined finishes:
#declare Fin_Glassy   = finish {{ phong 0.9 specular 0.8 reflection 0.1 roughness 0.001 }}
#declare Fin_Metallic = finish {{ phong 0.5 metallic 0.7 brilliance 2.0 diffuse 0.3 }}
#declare Fin_Matte    = finish {{ phong 0.0 ambient 0.1 diffuse 0.8 }}
// Define Bond Finishes
#declare Fin_Bd_Std = finish {{ phong 0.2 ambient 0.2 }}
// Select bond Finish
#declare BdFinish = Fin_Bd_Std;
// Definde Atom finishes
#declare Fin_Atom_Std = finish {{ phong 0.6 specular 0.4 ambient 0.2 }}
// select Atom Finish
#declare AtomFinish = Fin_Atom_Std;
///////////////////////////////////
        """) 

def export_pov_mol(points, atom_types,cov_radii=None, default_radius=None, 
                   cpk_colors=None,filename="test.inc", object_name="name", idx=0):
    #AtomsGroup
    with open(filename, 'a') as f:
        f.write(f"//Begin of Section {idx}  \n")
        f.write(f"#declare AtomsGroup_{idx} = union {{\n")
        for i, pos in enumerate(points):
            val = int(atom_types[i])
            # get colors from Color Dictionary
            color_name = cpk_colors.get(val, "magenta")
            # Conversion of Names to RGB for POV-Ray 
            rgb = {"white": "<1,1,1>", "gray": "<.3,.3,.3>", "blue": "<0,0,1>", 
                "red": "<1,0,0>", "orange": "<1, 0.55, 0>", "yellow": "<1,1,0>", 
                "brown": "<1,0.65,0>", "darkred": "<0.5,0,0>", "green": "<0, 0.82,0>"}.get(color_name, "<1,0,1>")

            atomic_num = int(atom_types[i])
            match atomic_num:
                case 1: # Hydrogen
                    rad_var = "atom_rad_h"
                case _ if 3 <= atomic_num <= 10: # 2nd period (He-Ne)
                    rad_var = "atom_rad_2"
                case _ if 11 <= atomic_num <= 18: # 3rd period (Na-Ar)
                    rad_var = "atom_rad_3"
                case _: # every other element
                    rad_var = "atom_rad_def"

            f.write(f"  sphere {{ <{pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}>, {rad_var}\n")
            f.write(f"    pigment {{ color rgb {rgb} filter trans_atom }}\n")
            f.write("    finish { AtomFinish }\n")
            f.write("  }\n")
        f.write("}\n")
    #BondsGroup
    with open(filename, 'a') as f:
        f.write(f"#declare BondsGroup_{idx} = union {{\n")
        for i in range(len(points)):
            type_i = int(atom_types[i])
            rad_i = cov_radii.get(type_i, default_radius)
            for j in range(i + 1, len(points)):
                type_j = int(atom_types[j])
                rad_j = cov_radii.get(type_j, default_radius) 
                bd_threshold = rad_i + rad_j + 0.6
                dist = np.linalg.norm(points[i] - points[j])
                # Threshold for Bonds in Angstrom
                if 0.6 < dist < bd_threshold:
                    p1 = points[i]
                    p2 = points[j]
                    f.write(f"  cylinder {{ <{p1[0]:.4f}, {p1[1]:.4f}, {p1[2]:.4f}>, "
                            f"<{p2[0]:.4f}, {p2[1]:.4f}, {p2[2]:.4f}>, bond_rad\n")
                    f.write("    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }\n")
                    f.write("    finish { BdFinish }\n")
                    f.write("  }\n")
        f.write("}\n")
# combine all objects
    with open(filename, 'a') as f:
        f.write(f"""\
// Combine Atoms, Bonds and Mesh 
#declare {object_name}[{idx}] = union {{
    object {{AtomsGroup_{idx}}}
    object {{BondsGroup_{idx}}}
}}
// End of Section {idx}
        """)

#### export Blender multi ####
class ExportWorker(QThread):
    progress = Signal(int)  # Signal for progressbar (0-100)   
    finished = Signal(bool, str) # Signal, after completion

    def __init__(self, tasks, base_name):
        super().__init__()
        self.tasks = tasks
        self.base_name = base_name
        self._is_running = True
        self.executor = None

    def stop(self):
        self._is_running = False
        if self.executor:
            # Stops all tasks that haven't started yet immediately
            self.executor.shutdown(wait=False, cancel_futures=True)

    def run(self):
        length = len(self.tasks)
        completed = 0

        # Execute in parallel
        n_proc = get_optimal_cores()
        workers = max(1, n_proc-1)
        with ProcessPoolExecutor(max_workers = workers) as self.executor:
            futures = [self.executor.submit(export_single_frame, t) for t in self.tasks]
            for _ in as_completed(futures):
                if not self._is_running:
                    self.finished.emit(False, self.base_name)
                    return
                
                completed += 1
                self.progress.emit(int((completed / length) * 100))
        
        self.finished.emit(True, self.base_name)

def get_radius_by_group(atomic_number):
    # Definition of Atom Radii
    groups = {
        (1, 1): 0.2,    # Hydrogen
        (2, 2): 0.2,    # Helium
        (3, 10): 0.35,   # 2. Period (Li to Ne)
        (11, 18): 0.45,  # 3. Period (Na to Ar)
        (19, 36): 0.55,  # 4. Period (K to Kr)
        (37, 54): 0.65,  # 5. Period (Rb to Xe)
        (55, 86): 0.75, # 6. Period (Cs-Rn, incl. Pt, Au)
    }
    
    for (start, end), radius in groups.items():
        if start <= atomic_number <= end:
            return radius
    return 0.3  # Default value

def draw_mol_bld(atom_points, atom_types, cpk_colors=None, cov_radii=None, default_radius=None):
    all_parts = []
    # create PolyData-Object from all points
    atoms_poly = pv.PolyData(atom_points)
    # add atom-tpye as scalar for coloring
    atoms_poly.point_data["colors"] = atom_types
    #sphere as template
    #sphere_source = pv.Sphere(radius=0.3, theta_resolution=20, phi_resolution=20)
    # color mapping:
    # glyph object contains original atom IDs, Lookup Table (LUT) can be used
    # loop over type  
    u_types = np.unique(atom_types)
    for atom_type in u_types:
        color = cpk_colors[atom_type]
        mask = atom_types == atom_type
        if np.any(mask):
            sub_atoms = atoms_poly.extract_points(mask)
            r=get_radius_by_group(atom_type)
            sphere_source = pv.Sphere(radius=r, theta_resolution=20, phi_resolution=20)   
            glyphs = sub_atoms.glyph(geom=sphere_source, scale=False, orient=False)
            
            rgb = (np.array(pv.Color(color).float_rgb) * 255).astype(np.uint8)
            colors_array = np.tile(rgb, (glyphs.n_points, 1))
            glyphs.point_data["RGB"] = colors_array
            all_parts.append(glyphs)
    
    # --- Bonds as single net
    lines = []
    for i in range(len(atom_points)):
        type_i = int(atom_types[i])
        rad_i = cov_radii.get(type_i, default_radius)
        for j in range(i + 1, len(atom_points)):
            type_j = int(atom_types[j])
            rad_j = cov_radii.get(type_j, default_radius) 
            dist = np.linalg.norm(atom_points[i] - atom_points[j])
            bd_threshold = rad_i + rad_j + 0.6
            if 0.6 < dist < bd_threshold:
                # only indices of linked points are saved
                lines.append([2, i, j]) # 2: line conects two points      
    tubes = None
    if lines:
        # Create PolyData-Object for lines
        bonds_poly = pv.PolyData(atom_points)
        bonds_poly.lines = np.hstack(lines)
        # convert lines into tubes
        tubes = bonds_poly.tube(radius=0.06)
        # bond color
        bond_rgb = (np.array(pv.Color("lightgray").float_rgb) * 255).astype(np.uint8)
        tubes.point_data["RGB"] = np.tile(bond_rgb, (tubes.n_points, 1))
        all_parts.append(tubes)
    # merge all meshes
    if not all_parts:
        return None 
    combined = all_parts[0].merge(all_parts[1:])
    return combined

def export_single_frame(args):
    i, atom_points, atom_types, cpk, radii, def_rad, base_name = args
    
    # create geometry mesh
    mesh = draw_mol_bld(atom_points, atom_types, cpk, radii, def_rad)
    
    # export (needs plotter)
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(mesh, scalars="RGB", rgb=True)

    file_path = f"{base_name}_{i:03d}.glb"

    pl.export_gltf(file_path)
    pl.close()
    return file_path

def generate_blender_script_multi(base_name, ver_no):
    """
    Generates a Blender Python script to batch-import multiple GLB files
    and sequence them in the timeline (One frame per file).
    """

    current_path = os.getcwd()

    script_path = f"{base_name}_animate.py"
    
    blender_script = f"""# created with MolAlign {ver_no} by (C) 2026 Dr. Tobias Schulz
import bpy
import os

# --- Settings ---
path_to_glb = "{current_path}" # adapt path to *glb files!
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
print(f"Finished! {{len(files)}} Frames processed.")
"""
    with open(script_path, "w") as f:
        f.write(blender_script)

def on_export_finished(success, base_name):
        log = []
        if success:
            generate_blender_script_multi(base_name, ver_no)
            print(f"Blender multi file export done")
        else:
            print(f"Export cancelled", "warning")
        return log

pbar = None
def update_progress(val):
    global pbar
    if pbar: 
        pbar.n = val
        pbar.refresh()

#### export Blender one ####
class OneFileExportWorker(QThread): # One File
    finished = Signal(bool, str) # Variables: Success, Path

    def __init__(self, data, path, name, cpk, radii, def_rad):
        super().__init__()
        self.data = data
        self.path = path
        self.name = name
        self.cpk = cpk
        self.radii = radii
        self.def_rad = def_rad

    def run(self):
        try:
            pl = pv.Plotter(off_screen=True)
            length = len(self.data.atom_points)
            
            for i in range(length):
                pts = np.array(self.data.atom_points[i])
            
                types = self.data.atom_types[i]
                # generate mesh
                mesh = draw_mol_bld(pts, types, self.cpk, self.radii, self.def_rad)
                
                # Generate a unique name for the Blender script
                pl.add_mesh(mesh, name=f"{self.name}_{i:03d}", scalars="RGB", rgb=True)
                
            pl.export_gltf(f"{self.path}_one.glb")
            pl.close()
            self.finished.emit(True, self.path)
        except Exception as e:
            self.finished.emit(False, "")

def on_one_file_finished(success, path):
    if success:
        generate_blender_script(path, ver_no) 
        print(f"Success: GLB export and script created for {os.path.basename(path)}")
    else:
        print("Error: GLB export failed.")

def generate_blender_script(path, ver_no):
    """
    Generates a companion Blender Python script for the exported GLB file.
    This script sets up a frame-by-frame animation by toggling object visibility.
    """
    script_path = f"{path}_one_file.py"
        
    blender_script = f"""import bpy
# created with MolAlign {ver_no} (C)2026 Dr. Tobias Schulz
# 1. Identify all objects starting with "mol_" (representing trajectory frames)
steps = [obj for obj in bpy.data.objects if "mol_" in obj.name]
steps.sort(key=lambda x: x.name)

if not steps:
    print("Error: No 'mol_' objects found in the scene!")
else:
    # 2. Detach objects from any hierarchies (Unparent)
    # This ensures each frame can be transformed independently if needed
    bpy.ops.object.select_all(action='DESELECT')
    for obj in steps:
        obj.select_set(True)
    
    # Clear parents while maintaining the current world transformation
    bpy.ops.object.parent_clear(type='CLEAR_KEEP_TRANSFORM')
    bpy.ops.object.select_all(action='DESELECT')

    # 3. Setup Scene Timeline
    bpy.context.scene.frame_start = 1
    bpy.context.scene.frame_end = len(steps)
    
    # 4. Create Frame-by-Frame Visibility Animation
    for i, obj in enumerate(steps):
        if obj.animation_data:
            obj.animation_data_clear()
            
        target_frame = i + 1
        
        # Initial State: Hidden (Scale set to 0)
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=0)
        
        # Visible State: Active only at its specific frame
        obj.scale = (1, 1, 1)
        obj.keyframe_insert(data_path="scale", frame=target_frame)
        
        # Hide immediately after its frame
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=target_frame + 1)

        # Force CONSTANT interpolation to prevent smooth scaling effects
        if obj.animation_data and obj.animation_data.action:
            for fcurve in obj.animation_data.action.fcurves:
                for kp in fcurve.keyframe_points:
                    kp.interpolation = 'CONSTANT'

    bpy.context.scene.frame_set(1)
    print(f"Animation Setup Complete: {{len(steps)}} frames sequenced.")
"""

    with open(script_path, "w") as f:
        f.write(blender_script)
      
### LOG File ###
def write_batch_log(output_path, log_entries):
    log_name = os.path.splitext(output_path)[0] + ".log"
    with open(log_name, "w") as f:
        f.write(f"MolAlign {ver_no} (C) 2026 Dr. Tobias Schulz\n")
        f.write("Batch Processing Log\n")
        f.write("============================\n\n")
        for entry in log_entries:
            f.write(f"{entry}\n")
    print(f"[INFO] Log file created: {log_name}")

# Atom Colors
cpk_colors = defaultdict(lambda: "magenta")
cpk_colors.update({
            1: "white",  #H
            5: "pink",  #B
            6: "gray",   #C
            7: "blue",   #N
            8: "red",    #O
            9: "orange",  #F
            14: "darkgrey", #Si
            12: "darkgreen", #Mg
            15: "brown",  #P
            16: "yellow", #S
            17: "green",  #Cl
            26: "darkorange", #Fe
            24: "darkcyan",    # Cr (Chrom)
            27: "royalblue",   # Co (Cobalt)
            28: "silver",      # Ni (Nickel)
            29: "chocolate",   # Cu (Kupfer)
            40: "cadetblue",   # Zr (Zirconium)
            44: "teal",        # Ru (Ruthenium)
            45: "deeppink",    # Rh (Rhodium) - oft zur Unterscheidung kräftig
            78: "lightgrey",    # Pt (Platin)
            35: "darkred",   #Br
            53: "darkviolet" # I
        })
        # Bond Parameters
cov_radii = {
        1: 0.31,   # H
        5: 0.82,   # B
        6: 0.76,   # C
        7: 0.71,   # N
        8: 0.66,   # O
        9: 0.57,   # F
        14: 1.11,  # Si
        15: 1.06,  # P
        16: 1.05,  # S
        17: 1.02,  # Cl
        35: 1.20,  # Br
        53: 1.39,  # I
        24: 1.39,  # Cr
        27: 1.26,  # Co
        28: 1.21,  # Ni
        29: 1.32,  # Cu
        40: 1.48,  # Zr
        44: 1.26,  # Ru
        45: 1.35,  # Rh
        78: 1.28   # Pt
        }
        # Standard Radius for unknown elements
default_radius = 1.0

# Version:
ver_no = 1.1

@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument('files', nargs=-1) # 'e.g. file1.xyz file2.xyz ...'

@click.option('--pov', is_flag=True, help='Export POV-Ray .inc files')
@click.option('--bld', is_flag=True, help='Export Blender .glb files (sequence)')
@click.option('--bld-one', is_flag=True, help='Export single Blender .glb file')
@click.option('--xyz', is_flag=True, help='Export combined .xyz trajectory')
@click.option('--log', is_flag=True, help='Provide log File')

@click.option('--obj_name', '-n', default='molecule',
              help='object name, default "molecule" ') 
@click.option('--fname', '-f', default='output',
              help='output file name or base name') 
@click.option('--split', default='none', type=click.Choice(['none','orca', 'nw', 'psi4']), 
              help='Generate split script for specific QM package (for --xyz only)')

@click.option('--rev', '-r', multiple=True, type=int, default=None,
              help='Indices of files to reverse, e.g. -r 0 -r 2 will reverse first and third file')

def main(files, pov, bld, bld_one, xyz, log, obj_name, fname, rev, split):
    # Data Import
    all_datasets = []
    log_data = []
    log_data.append(f"Input files: {files}")
    if not files:
        msg="no input found"
        click.echo(msg)
        log_data.append(msg)
        #write_batch_log(f"{fname}.log", log_data)
        return
    for i, file_path in enumerate(files):
        # 1. load file
        data = MoleculeData.from_xyz(file_path)
        # 2. check if index is in reverse 
        if i in rev:
            data = reverse(data)
            log_ = f"  -> Reversed trajectory for: {file_path}"
            click.echo(log_) 
            log_data.append(log_)
        all_datasets.append(data)
    
    # Alignment
    if len(all_datasets) == 1:  # only one input file given
            print(f"[INFO] Single file mode for: {files[0]}")
            combined_data = all_datasets[0]
            # No alignment/mapping needed, jump straight to export
    else:  # multiple inputs
        combined_data, alignment_rmsds, log_ = chain_alignment(all_datasets)
        log_data.append("Alignment:")
        log_data.extend(log_)

        if alignment_rmsds:
            avg_rmsd = sum(alignment_rmsds) / len(alignment_rmsds)
            
            msg1 = f"\n[DONE] Trajektory with {len(combined_data.energies)} Frames created."
            msg2 = f"[INFO] Average RMSD: {avg_rmsd:.4f} Å (Max: {max(alignment_rmsds):.4f} Å)"
            print(msg1)
            print(msg2)
            log_data.append(msg1)
            log_data.append(msg2)
            if max(alignment_rmsds) > 0.1:
                msg3 = "[WARN] High RMSDs detected, Check Geometry"
                print(msg3)
                log_data.append(msg3)
        
    # Export
    if xyz:
        log_ = export_xyz(combined_data, f"{fname}.xyz")
        click.echo(f"Trajectory saved as {fname}.xyz")
        log_data.append(log_)
        if split == 'orca':
            try:
                script_name = create_split_template_orca(fname, ver_no)
                msg = f"Batch template created: {os.path.basename(script_name)}"
                click.echo(msg)
                log_data.append(msg)
            except Exception as e:
                msg = f"Failed to create template: {str(e)}", "error"
                click.echo(msg)
                log_data.append(msg)
            
        elif split == 'nw':
            try:
                script_name = create_split_template_nw(fname, ver_no)
                msg = f"Batch template created: {os.path.basename(script_name)}"
                click.echo(msg)
                log_data.append(msg)
            except Exception as e:
                msg = f"Failed to create template: {str(e)}", "error"
                click.echo(msg)
                log_data.append(msg)
        elif split == "psi4":
            try:
                script_name = create_split_template_psi4(fname, ver_no)
                msg = f"Batch template created: {os.path.basename(script_name)}"
                click.echo(msg)
                log_data.append(msg)
            except Exception as e:
                msg = f"Failed to create template: {str(e)}", "error"
                click.echo(msg)
                log_data.append(msg)
        
    if pov:
        length = int(len(combined_data.energies))
        export_pov_header(ver_no, length,f"{fname}.inc",obj_name)
        p_bar = tqdm(range(length), desc="Exporting POV-Ray frames")
        for i in p_bar:
            export_pov_mol(np.array(combined_data.atom_points[i]),combined_data.atom_types[i],cov_radii,
                           default_radius,cpk_colors, f"{fname}.inc" , obj_name, i+1)
            if i % 10 == 0: # update every 10 frames
                p_bar.set_description(f"Writing Frame {i}/{length}")
        log_data.append(f"POV-Ray sequence {fname}.inc written")
        click.echo("[*] POV-Ray Sequence complete.")
        
    if bld:
        global pbar
        app = QCoreApplication.instance() or QCoreApplication([])

        base_name = fname

        raw_points = [np.array(p) for p in combined_data.atom_points]
        raw_types = list(combined_data.atom_types)
        pure_cpk = dict(cpk_colors) 
        pure_radii = dict(cov_radii)
        def_rad = default_radius

        # 2. Create the list of tasks for EVERY frame
        # This is the list the executor will iterate over later
        tasks = [
            (i, raw_points[i], raw_types[i], pure_cpk, pure_radii, def_rad, base_name)
            for i in range(len(raw_points))
        ]

        click.echo("Bld Export startet... ")
        #initialize tqdm
        pbar = tqdm(total=100, desc="Exporting Blender Frames", unit="%")
    
        worker = ExportWorker(tasks, base_name)
        worker.progress.connect(update_progress)

        # local Event-Loop, wait for 'finished' 
        loop = QEventLoop()
        worker.finished.connect(lambda success, b: [
            on_export_finished(success, b), 
            loop.quit(),
            pbar.close()
        ])
        
        # Start
        worker.start()
        loop.exec() # block until finished!

        log_data.append(f"Blender multi file export done to {fname}_*.glb")
        log_data.append(f"Configure Script {fname}_animate.py created ")

    if bld_one:
        click.echo("Bld One File Export startet... ")
        app = QCoreApplication.instance() or QCoreApplication([])

        one_worker = OneFileExportWorker(combined_data, fname, obj_name, cpk_colors,
                cov_radii, default_radius)
        
        loop = QEventLoop()
        one_worker.finished.connect(lambda success, path: [
            on_one_file_finished(success, path), 
            loop.quit()
        ])
        
        # Start
        one_worker.start()
        loop.exec() # block until finished!

        log_data.append(f"Blender One File export complete, {os.path.basename(fname)}.glb")

    if log:
        write_batch_log(f"{fname}.log", log_data)

if __name__ == '__main__':
    # filter for -B Flag directly from sys.argv, 
    # before it reaches the Argument-Parser 
    if "-B" in sys.argv:
        sys.argv.remove("-B")

    try:
        n_threads = get_optimal_cores()
        lib.num_threads(n_threads)
        print(f"Auto-Config: PySCF uses {n_threads} Threads.")
    except Exception as e:
        print(f"Could not set Threads automatically: {e}")

    main()