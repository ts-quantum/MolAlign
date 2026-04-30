# --- CONFIGURATION: ADJUST BEFORE RUNNING ---
<<<<<<< HEAD
##### created by MolAlign 1.2 (C) 2026 by Dr. Tobias Schulz
=======
##### created by MolAlign 1.3 (C) 2026 by Dr. Tobias Schulz
>>>>>>> work
##### This Python Script will create a series of ORCA Input files from 
##### a given IRC Trajectory along with a "run_ts1_batch.sh" script
### Edit method, basis set, etc. here
charge ='0'  # e.g. '0'	
mult ='3' 	# e.g. '1'		
basis ='def2-SVP'   # e.g. 'def2-SVP'
method ='PBE0 D3BJ UKS'  # e.g. 'HF', 'B3LYP', 'PBE0' or other XC
mem = '2000'
nproc = '8'
### adapt to your specific environment
ORCA_EXE = "/usr/local/orca_6_1_0/orca"
ORCA_2MKL = "/usr/local/orca_6_1_0/orca_2mkl"
# pre-configured by MolAlign
inp_file = 'ts1.xyz' 	

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
    wrapper_name = f"run_{base}_batch.sh"
    with open(wrapper_name, 'w') as w:
        w.write("#!/bin/bash\n\n")
        for i in range(steps):
            label = f"{base}_{i:03d}"
            inp = f"{label}.inp"
            out = f"{label}.out"
            xyz_list = lines[i*block_size + 2 : (i+1)*block_size]
            xyz_body = "".join(xyz_list).rstrip()
            body = f"""! {method} {basis}

%maxcore {mem}
%pal nprocs {nproc} end

* xyz {charge} {mult}
{xyz_body}
*
"""
            #write input file
            with open(inp, 'w') as f_inp:
                f_inp.writelines(body)
			# Commands for shell script
            w.write(f"{ORCA_EXE} {inp} > {out} 2>&1 && \\\n")
            w.write(f"{ORCA_2MKL} {label} -molden && \\\n")
            w.write(f"mv {label}.molden.input {label}.molden && \\\n")
            w.write(f"echo 'Frame {i:03d} finished.'\n\n")
    os.chmod(wrapper_name, 0o755)
    print(f"Done. Created {steps} inputs and shell script: {wrapper_name}")	

if __name__ == "__main__": 
    missing = [name for name, val in [("Charge", charge), ("Mult", mult), ("Basis", basis), ("Method", method)] if not str(val)]
    if missing:
        print(f"Error: Input missing: {', '.join(missing)}")
        sys.exit()
    run_split()
