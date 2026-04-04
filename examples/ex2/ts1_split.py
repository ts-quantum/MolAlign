import os

# --- CONFIGURATION: ADJUST BEFORE RUNNING ---
# created with MolAlign 1.0 (C)2026 Dr. Tobias Schulz
XYZ_INPUT = "ts1.xyz"  # The trajectory to split
CHARGE = 0
MULT = 1
ORCA_EXE = "/usr/local/orca_6_1_0/orca"
ORCA_2MKL = "/usr/local/orca_6_1_0/orca_2mkl"

# Custom Header (Edit method, basis set, etc. here)
ORCA_HEADER = """! PBE0 D3BJ UKS def2-SVP
%maxcore 1200
%pal nprocs 8 end
"""

def run_split():
    if not os.path.exists(XYZ_INPUT):
        print(f"Error: {XYZ_INPUT} not found.")
        return

    with open(XYZ_INPUT, 'r') as f:
        lines = f.readlines()

    try:
        num_atoms = int(lines[0].strip())
    except:
        print("Error: Invalid XYZ format.")
        return

    block_size = num_atoms + 2
    steps = len(lines) // block_size
    base = os.path.splitext(XYZ_INPUT)[0]
    wrapper_name = f"run_{base}_batch.sh"

    with open(wrapper_name, 'w') as w:
        w.write("#!/bin/bash\n\n")
        for i in range(steps):
            label = f"{base}_{i:03d}"
            inp = f"{label}.inp"
            out = f"{label}.out"
            
            # Write single point input
            with open(inp, 'w') as f_inp:
                f_inp.write(ORCA_HEADER)
                f_inp.write(f"* xyz {CHARGE} {MULT}\n")
                f_inp.writelines(lines[i*block_size + 2 : (i+1)*block_size])
                f_inp.write("*\n")
            
            # Commands for shell script
            w.write(f"{ORCA_EXE} {inp} > {out} 2>&1 && \\\n")
            w.write(f"{ORCA_2MKL} {label} -molden && \\\n")
            w.write(f"mv {label}.molden.input {label}.molden && \\\n")
            w.write(f"echo 'Frame {i:03d} finished.'\n\n")

    os.chmod(wrapper_name, 0o755)
    print(f"Done. Created {steps} inputs and shell script: {wrapper_name}")

if __name__ == "__main__":
    run_split()
