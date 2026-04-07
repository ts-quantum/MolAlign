##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script to run a Psi-4 single point calculation
##### for each point on a given IRC Trajectory and to produce a 
##### .molden and .fchk file for each point
##### must be run inside a Psi-Conda environment
##### User section:
### method details (examples, adapt to your case)
charge = 0
mult = 1
basis = 'cc-pVDZ'
ref = 'rhf'  # or 'rks', 'uhf', 'uks', ...
method = 'scf' # for dft: e.g. 'pbe0', 'b3lyp', ...
mem = '1000mb'
nproc = 8
#pre-configured by MolAlign
xyz_input = {xyz_filename} 

import os
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
	mol= psi4.geometry(f"""
		{{charge}} {{mult}}
		{{initial}}
		units angstrom
		no_reorient
		no_com
		""")
	
	last_wfn = None

	for i in range(steps):
		xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
		xyz_body = "".join(xyz_list)
		temp_mol = psi4.geometry(f"units angstrom\n{{xyz_body}}")
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
    run_split()
