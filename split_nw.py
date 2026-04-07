##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script will create a series of NWChem Input files from 
##### a given IRC Trajectory along with a "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
##### User section:
### method details (examples, adapt to your case)
charge = 0				
basis = '6-31G'
method_type = 'dft' # or 'scf' for HF
method_details ="""
 xc b3lyp
 mult 1
"""
### directories will be created by the "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
work_dir = 'calc'
scr_dir = 'scr'
### adapt to your specific environment
nw_cmd = 'mpiexec -np 8 nwchem'
# pre-configured by MolAlign
inp_file = {xyz_filename}	

import os

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
		w.write("#!/bin/bash\n\n")
		w.write(f"mkdir -p ./{{work_dir}}\n")
		w.write(f"mkdir -p ./{{scr_dir}}\n\n")
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
			w.write(f"{{nw_cmd}} {{inp}} > {{out}} 2>&1 && \\\n" )
			w.write(f"mv ./{{work_dir}}/{{base}}_{{i:03d}}.molden ./ 2>/dev/null && \\\n")
			w.write(f"echo 'Frame {{i:03d}} finished.'\n\n")
	os.chmod(wrapper_name, 0o755)
	print(f"Done. Created {{steps}} inputs and shell script: {{wrapper_name}}")	

if __name__ == "__main__":

    run_split()
