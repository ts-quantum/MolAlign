import psi4
from ase import io
from ase.calculators.psi4 import Psi4
from ase.optimize import FIRE
from ase.vibrations import Vibrations
import numpy as np

# --- 1. Load TS and Setup Calculator ---
# Ensure 'ts_optimized.xyz' is in the same folder
atoms = io.read('ts_optimized.xyz')
calc = Psi4(method='b3lyp', basis='cc-pVDZ', num_threads=8, memory='4GB', reference='rks')
atoms.calc = calc

# --- 2. Identify the Reaction Mode (Imaginary Frequency) ---
# This creates a folder 'ts_vib' with displacement data
vib = Vibrations(atoms, name='ts_vib')
vib.run()
# Mode 0 is the one with the largest imaginary frequency (the TS mode)
mode = vib.get_mode(0) 

# --- 3. Run "Manual" IRC Forward ---
print(">>> Starting Forward Path...")
atoms_fwd = atoms.copy()
# Move 0.1 Angstrom along the vibration vector
atoms_fwd.set_positions(atoms.get_positions() + 0.1 * mode)
atoms_fwd.calc = calc
opt_fwd = FIRE(atoms_fwd, trajectory='irc_forward.traj')
opt_fwd.run(fmax=0.05, steps=50)

# --- 4. Run "Manual" IRC Backward ---
print(">>> Starting Backward Path...")
atoms_bwd = atoms.copy()
# Move 0.1 Angstrom in the opposite direction
atoms_bwd.set_positions(atoms.get_positions() - 0.1 * mode)
atoms_bwd.calc = calc
opt_bwd = FIRE(atoms_bwd, trajectory='irc_backward.traj')
opt_bwd.run(fmax=0.05, steps=50)

# --- 5. Final Export for MolAlign ---
fwd_path = io.read('irc_forward.traj', index=':')
bwd_path = io.read('irc_backward.traj', index=':')
#full_path = bwd_path[::-1] + [atoms] + fwd_path
#io.write('full_irc_path.xyz', full_path)

io.write('ts.ircf.xyz', fwd_path)
io.write('ts.ircb.xyz', bwd_path)

#print(">>> IRC path saved to 'full_irc_path.xyz'")

