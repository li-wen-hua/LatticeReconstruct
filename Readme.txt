Readme
2025.6.10


Program Description
-------------------
This program is designed to construct a reference configuration suitable for Wigner-Seitz (WS) defect analysis in crystalline structures containing dislocation lines or even dislocation networks. It “repairs” point defects while preserving the dislocation structure, enabling more accurate identification and statistical analysis of point defects.



Dependencies and Pre-run Configuration
--------------------------------------
- Programming Language: Fortran90
  The program consists of multiple `.f90` source files coordinated by shell scripts.
- Compiler: gfortran (any compiler supporting Fortran90 is sufficient)
- Compilation: No manual compilation is required. The shell script contains the `gfortran *.f90` command. Make sure `gfortran` is available in your system path and can compile `.f90` files.
- Parallelism: Supports single-node parallelism (OpenMP only).
- Model Limitation: Currently applicable only to BCC crystal structures with periodic boundary conditions. Non-periodic models are untested.
- A collection of auxiliary scripts designed to assist this program are provided in LatticeReconstructionTools repository (https://github.com/li-wen-hua/LatticeReconstructTools/tree/main) . These tools help with generating and modify input files, performing Wigner-Seitz (WS) analysis, and filtering true point defects. Each module operates independently and can be used as needed.



Program Structure
-----------------
1. Files in Main Directory:
   - `repair.sh`: Main execution script
   - `jobrun.sh`, `repair.pbs`: Example job run/submission scripts
   - `repair.params`: Parameter configuration file
   - `Figure1_ParallelPerformance.tif`, `Figure2_IdentificationResultsOfExample.tif`: Parrallel performance and identification result of a small simulation box.


2. Subdirectories:
   - `datafiles-theINI`: Stores original input files
   - `datafiles-input`: Intermediate input files generated during the repair process
   - `datafiles-output`: Output files from each repair step and the final result
   - `example`: Example input files

3. Modules:
   - `vacfill`: Fills vacancies
   - `csp`: Identifies and extracts defective atoms (similar to CNA), updates `type0.data` for each iteration
   - `subtract`: Subtracts defective atoms from the original structure



Input Files (Stored in `datafiles-theINI`)
------------------------------------------
1. `dump.0`: Original LAMMPS dump format structure file (with defects and dislocations)
2. `type0.data`: List of atoms with unrecognized structures (CNA/DXA) in LAMMPS data format
3. `disline.input`: Dislocation data in CrystalAnalysis format from DXA
ps: An example structure and input files are provided in "example" directory
ps: Tools for generating and modifying input files via Ovito Python scripts are provided in LatticeReconstructionTools repository.



Output Files (Stored in `datafiles-output`)
-------------------------------------------
1. `filledbox.data`: Reference lattice in LAMMPS data format
2. `filledbox.dump`: Reference lattice in LAMMPS dump format



Parameter File: `repair.params`
-------------------------------
- Parameter descriptions are provided as comments in the file.
- Typically, only the first two and the last parameter need to be modified (lattice constant, number of atom types, number of parallel threads).
- The third parameter defines the "dislocation region," which defaults to twice the lattice constant. Adjust this as needed to balance identification accuracy and computational cost.



Running Instructions
--------------------
Step 1: Prepare input files and place them in `datafiles-theINI/`
Step 2: Set parameters in `repair.params`
Step 3: Execute the program
- Run in foreground: `bash repair.sh`
- Run in background: `bash repair.sh > test.out &` or submit job script
Step 4: Check results
- Output includes:
  - `filledbox.data`: LAMMPS data format
  - `filledbox.dump`: LAMMPS dump format
  - `nohup.out` or `myrepair.log`: Program log file
- Use `grep @@@@ logfile` to quickly check iteration results and termination status




Performance & OpenMP Parallelization
------------------------------------
- Parallel regions are mainly in `vacfill/3-FindingMirrorAtoms`
- Runtime roughly scales with the square or cube of structure size and linearly with defect density
- Small model (~100k atoms): ~10–30 seconds with parallelization
- Large model (~4.4 million atoms): ~8–18 minutes with parallelization
- Parallel efficiency is optimal around 8 cores; larger models benefit from more threads
- Parrallel performance and identification result of a small box (~100k atoms) are shown in Figure1_ParallelPerformance.tif and Figure2_IdentificationResultsOfExample.tif


Author’s Note
-------------
This program is still under active development. The author is not a professional programmer, and the current implementation remains limited in optimization. The method of “reference lattice reconstruction” represents an initial exploration and is yet to be fully developed to handle more complex systems such as grain boundaries or polycrystals. There is significant potential for improvement and expansion. Feedback and contributions are welcome from users interested in further development.


