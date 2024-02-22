# QRef
`QRef` is a plugin for the crystallographic software suite `Phenix` enabling what in the literature commonly is referred to as "quantum refinement" (QR) in both real and reciprocal space, utilising the free software Orca as the quantum chemistry engine. This first version of `QRef` using `Phenix` was implemented under the paradigm "correctness and completeness first, performance and adherence to coding standards later".

## Theory
Refinement of macromolecules in both real and and reciprocal space relies on previous knowledge (i.e. a Bayesian prior) for the structure. This is usually encoded as a (pseudo-energy) penalty term, $E_{restraints}(\mathbf{R})$, giving rise to a target function for the refinement with the general appearance
```math
E_{total}\left(\mathbf{R}\right) = E_{exp}\left(\mathbf{R}\right) + wE_{restraints}\left(\mathbf{R}\right),
```
where $\mathbf{R}$ is the coordinate set for the current model and $w$ is a weight factor.
$E_{restraints}(\mathbf{R})$ can in turn be broken down into its components:
```math
E_{restraints}(\mathbf{R}) = E_{chem}(\mathbf{R}) + E_{SS}(\mathbf{R}) + E_{NCS}(\mathbf{R}) + \dots
```
In this QR implementation $E_{chem}(\mathbf{R})$ (which traditionally is a molecular mechanics force field) is replaced using a subtractive QM/MM scheme (using hydrogen link-atoms) according to
```math
E_{chem}(\mathbf{R}) = w_{QM}\sum_{i} \left(E_{QM1, i}(\mathbf{R_{syst1, i}}(\mathbf{R})) - E_{MM1, i}(\mathbf{R})\right) + E_{MM}(\mathbf{R}),
```
where index 1 in turn indicates small, but interesting, parts of the structure. Additionally another scaling factor, $w_{QM}$, is needed due to the fact that crystallographic MM force fields are of a statistical nature, whereas $E_{QM1, i}$ represents physical energies. $\mathbf{R_{syst1, i}}$ in turn is the coordinate set for the i:th region of QM atoms, where the hydrogen link-atoms are placed at
```math
\overline{r}_{H_L} = \overline{r}_X + g_{bond}\left(\overline{r}_{C_L} - \overline{r}_X\right).
```
The gradient for the chemical restraints is then obtained as
```math
\nabla E_{restraints}(\mathbf{R}) = w_{QM} \sum_{i} \left( \nabla E_{QM1, i}(\mathbf{R_{syst1, i}(R)}) \cdot J(\mathbf{R_{syst1,i}}; \mathbf{R}) - \nabla E_{MM1, i}(\mathbf{R}) \right) + \nabla E_{MM}(\mathbf{R})
```
where $J(\mathbf{R_{syst1,i}}; \mathbf{R})$ is the Jacobian between $\mathbf{R_{syst1,i}}$ and $\mathbf{R}$. While this obviosuly is a matrix of size $3N_{syst1, i} \times 3N$, for a single junction in one QM system this becomes a 6x6 matrix with the general shape
```math
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
(1 - g_{bond}) & 0 & 0 & g_{bond} & 0 & 0 \\
0 & (1 - g_{bond}) & 0 & 0 & g_{bond} & 0 \\
0 & 0 & (1 - g_{bond}) & 0 & 0 & g_{bond}
\end{pmatrix}.
```

## Installation
### Modules
The directory `modules` should be placed under `<path to Phenix>`; `qref` will thus be a new directory under `modules`, whereas the user should manually overwrite `energies.py` in `modules/cctbx_project/cctbx/geometry_restraints` and `model.py` in `modules/cctbx_project/mmtbx/model`, respectively, with the version of the file corresponding to their installation of `Phenix`.

There is a commented out guard clause in `energies.py`:<br>

    # if not os.path.exists('qm.lock') and (os.path.exists('xyz_reciprocal.lock') or os.path.exists('xyz.lock')):

This is the recommended way to use the quantum restraints, as they are not always needed. In order to make this work one has to edit the file `<path to Phenix>/modules/phenix/phenix/refinement/xyz_reciprocal_space.py` and `import os` as well as surround the call to `mmtbx.refinement.minimization.lbfgs(...)` in the method `run_lbfgs` in the class `run_all` with<br>

    with open('xyz_reciprocal.lock', 'w'):
        pass

and

    os.remove('xyz_reciprocal.lock')

Likewise the file `<path to Phenix>/modules/phenix/phenix/refinement/macro_cycle_real_space.py` should be edited in a similar manner, i.e. with an added `import os` as well as surrounding the calls to `self.minimization_no_ncs()` and `self.minimization_ncs()` in the method `refine_xyz` in the class `run` with<br>

    with open('xyz.lock', 'w'):
      pass

and

    os.remove('xyz.lock')

This implementation of `QRef` has been verified to work with `Phenix` version 1.20.1-4487.

### Scripts
The scripts in the folder `scripts` should be placed somewhere accessible by `$PATH`. The shebang in the scripts might need to be updated to point towards wherever `cctbx.python` is located.

### Orca
`Orca` can be found at https://orcaforum.kofo.mpg.de/ - follow their guide for installation.

## Usage
The general procedure to set up a quantum refinement job consists of

1. Build a model.
    - The model in the QM regions needs to make chemical sense. This for example means that the QM regions should be protonated as well as being complete.
    - The model outside of the QM region, bar the carbon link atom, can be incomplete.
    - `phenix.ready_set add_h_to_water=True` can be useful for this purpose.

2. Prepare restraint files for unknown residues and ligands. The script `qref_prep.py` will tell you if there are any missing restraint files.
    - This can be achieved using `phenix.ready_set` and `phenix.elbow`.

3. Prepare `syst1` files; these files define the QM regions.
    - If there is only one QM region the default is to look for a file named `syst1` by the software. For multiple QM regions the recommended, and default, naming scheme is `syst11`, `syst12`, etc.
    - Which atoms to include in the QM regions is defined using the serial number from the PDB file describing the entire model.
        - While setting `sort_atoms = False` in the input to `Phenix` should ensure that the ordering in the input model is preserved, we have encountered instances where this is not adhered to. Thus it is recommended to use `iotbx.pdb.sort_atoms` (supplied with `Phenix`) which will give you a new PDB file with the suffix `_sorted.pdb` where the atoms are sorted in the same order as that which `Phenix` uses internally. It is recommended to use the `_sorted.pdb` file as the input model for refinement, as well as the reference when defining the `syst1` files.
    - The `syst1` files allows for multiple atoms or intervals of atoms to be specified on a single line, where `,` or `blank` works as delimiters; `-` is used to indicate an interval.
    - `#` and `!` can be used to include comments in the `syst1` files.
    - The second occurence of an atom in the `syst1` files will indicate that this is a link atom, i.e. it will be replaced by a hydrogen at the appropriate position in the QM calculation.
    - Examples are included in the `examples` folder.

4. Run `qref_prep.py <model>_sorted.pdb` to generate `qref.dat` (a file containing settings for the QR interface), as well as PDB files describing the QM regions.
    - The `junctfactor` file needs to be present in the same directory as where `qref_prep.py` is run. The `junctfactor` file contain ideal QM distances for the $C_L - H_L$ bonds for the link-atoms. If another `junctfactor` file is to be used this can be specified with the `-j` or `--junctfactor` option.
    - The theory used for the ideal $C_L - H_L$ QM distance can be changed with the `-l` or `--ltype` option. Default is type 12. The options are:
        - 6: B3LYP/6-31G*
        - 7: BP(RI)/6-31G*
        - 8: BP(RI)/SVP
        - 9: BP(RI)/def2-SV(P)
        - 10: PBE(RI)/def2-SVP
        - 11: B3LYP(RI)/def2-SV(P)
        - 12: TPSS(RI)/def2-SV(P)
        - 13: B97-D(RI)/def2-SV(P)

        There needs to be a parametrisation in the `junctfactor` file for the bond one intends to cleave; it is recommended that the user inspects the `junctfactor` file to verify that there is support to cleave the intended bond type. In the case parametrisation is lacking another selection for the QM system (and in particular where the link between QM and MM occurs) needs to be made or appropriate parametrisation added to the `junctfactor` file.
    - Ideally only the input model is needed as an argument for `qref_prep.py`. If there was a need to prepare restraint files for novel residues or ligands in point 2 above, `qref_prep.py` needs to be made aware of these. This can be achieved with the `-c` or `--cif` option.
    - The output from `qref_prep.py` should be \# $\left(1+2 n_{syst1}\right)$ files as well as selection strings:
        - `qref.dat`, which contains the settings for the QR interface. This file can be changed manually and it is a good idea to inspect that the value for `orca_binary` is the correct path for the actual Orca binary file (`qref_prep.py` tries to locate this file automatically but may sometimes fail).
        - `qm_i_c.pdb`, which is the model used to calculate $E_{MM1, i}$.
        - `qm_i_h.pdb`, which is the model used to calculate $E_{QM1, i}$.

        The output PDB files can, and probably should, be used to inspect that the QM selection is proper.
        - Two selection strings are printed on the screen, one for reciprocal space and one for real space. They are intended to be used in regards to which selection of the model to refine when crafting the input to either `phenix.refine` or `phenix.real_space_refine`, see point 6 below.

    - Harmonic bond length restraints can be added through the `-rb` or `--restraint_bond` option, using the syntax `i atom1_serial atom2_serial desired_distance_in_Å force_constant`. Experience has shown that the force constant needs to be $\geq$ 10 to achieve adherence to the restraint.

    - All available options for `qref_prep.py` can be seen through the `-h` or `--help` option.

5. The next step is to prepare the input files for Orca. Examples can be found in the `examples` folder.
    - The input files should be named `qm_i.inp`, where `i` refers to the i:th QM region; counting starts at 1.
    - The level of theory should match the junction type; using the default (type 12) the corresponding input to `Orca` then becomes `! TPSS D4 DEF2-SV(P)`.
    - Energy and gradient needs to be written to disk. This is achieved through the keyword `! ENGRAD`.
    - To read coordinates from a PDB file (`qm_i_h.pdb`) use `*pdbfile <charge> <multiplicity> qm_i_h.pdb`, where again `i` refers to the i:th QM region.
    - Custom settings can also be supplied (see the `examples` folder).

6. At this point it is possible to start a quantum refinement. It is however recommended to first create an empty file named `qm.lock` (this disables QR) through for example the `touch` command, then start a refinement so that an `.eff` file is obtained (this file contains all the `Phenix` refinement settings for the current experimental data and model). A good idea is to rename this file to `phenix.params` or similar, then edit this file and make sure the following options are set:
    - For reciprocal space refinement (`phenix.refine`):
        - `refinement.pdb_interpretation.restraints_library.cdl = False`
        - `refinement.pdb_interpretation.restraints_library.mcl = False`
        - `refinement.pdb_interpretation.restraints_library.cis_pro_eh99 = False`
        - `refinement.pdb_interpretation.secondary_structure.enabled = False`
        - `refinement.pdb_interpretation.sort_atoms = False`
        - `refinement.refine.strategy = *individual_sites individual_sites_real_space rigid_body *individual_adp group_adp tls occupancies group_anomalous`
        - `refinement.refine.sites.individual = <reciprocal selection string>`
        - `refinement.hydrogens.refine = *individual riding Auto`
        - `refinement.main.nqh_flips = False`
    - For real space refinement (`phenix.real_space_refine`):
        - `refinement.run = *minimization_global rigid_body local_grid_search morphing simulated_annealing adp occupancy nqh_flips`
        - `pdb_interpretation.restraints_library.cdl = False`
        - `pdb_interpretation.restraints_library.mdl = False`
        - `pdb_interpretation.restraints_library.cis_pro_eh99 = False`
        - `pdb_interpretation.sort_atoms = False`
        - `pdb_interpretation.secondary_structure = False`
        - `pdb_interpretation.reference_coordinate_restraints.enabled = True`
        - `pdb_interpretation.reference_coordinate_restraints.selection = <real space selection string>`
        - `pdb_interpretation.reference_coordinate_restraints.sigma = 0.01`
    - Other options can be set as needed.

7. To run the quantum refinement job, make sure that the `qm.lock` file has been deleted, then execute either `phenix.refine phenix.params` or `phenix.real_space_refine phenix.params`. If there is a need to restart the job with different settings for `Phenix`, make sure to delete the file `settings.pickle`.

## Citation