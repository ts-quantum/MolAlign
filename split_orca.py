##### created by MolAlign {ver_no} (C) 2026 by Dr. Tobias Schulz
##### This Python Script will create a series of ORCA Input files from 
##### a given IRC Trajectory along with a "run_{os.path.splitext(xyz_filename)[0]}_batch.sh" script
##### User section:
### method details (examples, adapt to your case)
charge = 0	
mult = 1			
basis = '6-31G'
method = 'HF' # for DFT, e.g. 'B3LYP', 'PBE0' or other XC
mem = 2048
nproc = 8
### adapt to your specific environment
ORCA_EXE = "/usr/local/orca_6_1_0/orca"
ORCA_2MKL = "/usr/local/orca_6_1_0/orca_2mkl"
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
		for i in range(steps):
			label = f"{{base}}_{{i:03d}}"
			inp = f"{{label}}.inp"
			out = f"{{label}}.out"
			xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
			xyz_body = "".join(xyz_list).rstrip()
			body = f"""! {method} {basis}

%maxcore {{mem}}
%pal nprocs {{nproc}} end

* xyz {{charge}} {{mult}}
{{xyz_body}}
*
"""
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

    run_split()
