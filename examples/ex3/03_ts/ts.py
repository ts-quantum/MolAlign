import psi4
import os

# --- Resources ---
psi4.set_memory('4 GB')
psi4.core.set_num_threads(8)

with open('ts.xyz', 'r') as f:
    xyz_content = f.read()

#skip line one and two
xyz_body = "\n".join(xyz_content.strip().split("\n")[2:])

mol= psi4.geometry(f"""
0 1
{xyz_body}
units angstrom
no_reorient
no_com
""")

xc = 'b3lyp'

# --- 2. Options for TS Search ---
psi4.set_options({
    'basis': 'cc-pVDZ',
    'scf_type': 'df',
    'reference': 'rks',
    'opt_type': 'ts',            # Switch to Transition State optimization
    'g_convergence': 'gau_tight',
    'full_hess_every': 0         # Calculate Hessian at the beginning (crucial for TS)
})

print(">>> Starting TS Optimization...")
# Transition states require a Hessian to find the correct curvature
energy, wfn = psi4.optimize(xc, molecule=mol, return_wfn=True)
psi4.molden(wfn, 'ts.molden')

# Save the optimized TS geometry
wfn.molecule().save_xyz_file('ts_optimized.xyz', True)
# Save the wavefunction
wfn.to_file('opt.npy')

# --- 3. Frequency Analysis ---
print(">>> Starting Frequency Analysis...")
wfn_rel = psi4.core.Wavefunction.from_file(f'opt.npy')
wfn_rel.to_file('18.npy')
psi4.set_options({
    'normal_modes_write' : True,
    'guess':'read'
    })
energy, wfn = psi4.frequency(xc, molecule=mol, return_wfn=True)

print(">>> TS Workflow finished. Check the output for imaginary frequencies.")

