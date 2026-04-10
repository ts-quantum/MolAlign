MolAlign: CLI Batch Processing Manual

MolAlign is the command-line powerhouse of the suite. It is designed to stitch multiple 
reaction segments (e.g., separate IRC calculations) into one seamless, physically 
consistent trajectory.

1. Core Concept: The "Domino" Alignment
Unlike simple concatenation, MolAlign uses a Chained Alignment strategy:
    1. Reference: It starts with the first file as a fixed reference.
    2. Mapping & Flip: For every subsequent file, it automatically detects if the atoms are ordered differently or if the trajectory needs to be reversed (End-to-Start matching).
    3. Transformation: It applies the optimal rotation and translation (Kabsch Algorithm) to align the new segment perfectly to the existing chain.
    4. Energy Stitching: It calculates energy offsets at connection points to create a continuous, jump-free reaction profile.

2. Command Line Syntax
Usage: main.py [file1.xyz] [file2.xyz] ... [options]

Options:
  --pov                        Export POV-Ray .inc files
  --bld                        Export Blender .glb files (sequence)
  --bld-one                    Export single Blender .glb file
  --xyz                        Export combined .xyz trajectory
  --log                        Provide log File
  -n, --obj_name TEXT          object name, default "molecule"
  -f, --fname TEXT             output file name or base name
  --split [none|orca|nw|psi4]  Generate split script for specific QM package
                               (for --xyz only)
  -r, --rev INTEGER            Indices of files to reverse, e.g. -r 0 -r 2
                               will reverse first and third file
  --help                       Show this message and exit.

Key Arguments:
    • files: List of XYZ files in chronological order (e.g., edu_to_ts.xyz ts_to_prod.xyz).
    • --output / -o: Filename for the merged trajectory (Default: combined_trj.xyz).
    • --reverse [index]: Manually force-reverse a specific segment (rarely needed due to Auto-Flip).

3. Advanced Export Options
MolAlign can trigger high-end visualization exports directly from the console:
    • --pov: Generates a POV-Ray scene file for raytracing.
    • --blender-multi: Exports each frame as an individual .glb file and creates an import_and_animate.py 
    script for Blender.
    • --blender-one: Exports the entire trajectory as a single .glb file and a setup_anim.py script.

4. The "Split & Run" Workflow
After a successful alignment, MolAlign offers to create a Quantum Chemistry Bridge:
    1. Split Template: A pre-configured _split.py script is generated.
    2. Customization: Open _split.py to adjust the ORCA_HEADER (Method, Basis Set, etc.).
    3. Execution: Run the script to generate individual .inp files and a .sh batch wrapper for your cluster.

5. Understanding the Output
During execution, MolAlign provides real-time feedback:
    • Atom re-ordering detected: Indicates that the mapping algorithm corrected inconsistent atom indices.
    • Average RMSD: Shows the quality of the fit between segments.
    • Warning: High RMSD: If RMSD > 0.1 Å, check if the segments actually share the same chemical "anchor" 
    or if a wrong file was provided.

Troubleshooting: Blender Integration
If you are using the Blender export scripts:
    1. Import: In Blender, import the .glb file(s).
    2. Scripting: Switch to the Scripting Tab in Blender.
    3. Open & Run: Open the generated .py script (e.g., import_and_animate.py) and click Run Script.
    4. Result: Your reaction path is now fully sequenced on the Blender timeline.

6. Examples:

Example 1: 2-Chlorocyclohexane Ring Flip (NWChem Workflow) B3LYP/6-31G, [1]
    Pathway Description:
    This example constructs the complete reaction pathway for the cyclohexane ring flip using NWChem IRC output files.
    The required input geometries and NWChem files are provided in the /nwchem directory.
    Command Line Execution:
    To merge the segments and generate the visualization files, run:
    python3 ../../main.py 2.irc.fxyz 2.irc.bxyz 4.irc.fxyz 4.irc.bxyz 6.irc.fxyz 6.irc.bxyz --xyz --pov -n cyhex --fname Cyhex --log
    Parameters Explained:
    --xyz: Generates the combined trajectory file.
    --pov: Creates the POV-Ray include file.
    -n cyhex: Defines the prefix for the molecule array in POV-Ray (accessible via cyhex[i]).
    --fname Cyhex: Sets the base filename for all output files.
    --log: Generates a detailed processing log.
    Results:
    Cyhex.xyz: The full, combined reaction trajectory.
    Cyhex.inc: A POV-Ray include file containing all structures. You can call specific frames in your .pov file using the cyhex array.
    Cyhex.log: Log file documenting the alignment and merging process.
    Visualization: See test1.pov for a ready-to-render setup using the generated include file.
[1] E. Aprà, D. Mejía-Rodríguez, et al., "NWChem: Recent and Ongoing Developments", J. Chem. Theory Comput., 
19, 7077–7096 (2023). doi: 10.1021/acs.jctc.3c00421. 

Example 2: Addition of Methyl Carbene to Propene (Splitting for BatchMol) PBE0/def2-SVP, [2]
    Pathway Description:
    This example focuses on the radical attack of methyl carbene on propene.
    The required IRC trajectory is provided as ts1.irc_IRC_Full_trj.xyz.
    Command Line Execution:
    To process the trajectory and prepare it for further quantum chemical calculations:

    python3 ../../main.py ts1.irc_IRC_Full_trj.xyz --fname ts1 --log --xyz --split orca

    Workflow & Integration:
    Splitting: The --split orca flag generates a dedicated Python script: ts1_split.py.
    Input Generation: Run ts1_split.py to automatically create individual ORCA input files for every single 
    point along the trajectory.
    Batch Processing: These inputs are used to obtain molden files for each step, which can then be 
    processed with BatchMol.

    Result:
    Allows for high-resolution animation of electronic properties, such as tracking the radical center or spin density 
    throughout the entire reaction path.
[2] F. Neese, "Software update: the ORCA program system — Version 6.0", Wiley Interdiscip. Rev.: Comput. Mol. Sci., 
15, e70019 (2025). doi: 10.1002/wcms.70019.

Example 3: 1-5-H-Shift (Splitting for BatchMol) B3LYP/cc-pVDZ
    Pathway Description:
    This example demonstrates the hydrogen migration between the terminal methyl group and the oxygen 
    atom of but-2-en-1-one. It includes input and output files for TS optimization and IRC calculations 
    performed with PSI4.
    MolAlign was used to generate:
    irc_split.py: A script to create .molden and .fchk files for BatchMol 
    (enabling ESP or MO visualization along the pathway).
    Blender Assets: Individual .glb files for each step of the pathway, 
    including an irc_animate.py script for automated processing in Blender.
    POV-Ray Files: An .inc file containing a molecule array over the reaction path, 
    along with the corresponding input/ini files and the final rendered MP4 video.
[3] D. G. A. Smith, L. A. Burns, et al., "Psi4 1.4: Open-source software for high-throughput quantum chemistry", 
J. Chem. Phys., 152, 184108 (2020). doi: 10.1063/5.0006002.

## Installation

1. Clone the repository
    git clone https://github.com
    cd MolAlign

2. Install dependencies
    pip install -r requirements.txt

    Requirements
        Python 3.x
        PySide6
        PyScf
        ....

## Usage

    python3 main.py

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.