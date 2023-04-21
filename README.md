# AFICS
AFICS is a Python 3 tool for the **A**nalysis of **F**irst **I**on **C**oordination **S**phere structural features over the course of a molecular dynamics (MD) simulation. AFICS processes the atomic coordinates of an MD trajectory to calculate the radial distribution function (RDF), coordination number (CN), angle distribution function (ADF) of elements surrounding a central ion. AFICS is also able compute the distortion (measured as root mean square deviation, RMSD) of the MD coordination geometry from ideal, uniform polyhedra over time.

### Dependencies
* Python 3.0 or after
* Libraries: numpy, csv, itertools

## Summary
Two files are of primary interest for running AFICS: one, the user needs to provide a single file containing the atomic coordinates of an MD simulation trajectory in XYZ format; and two, the user needs to provide an input script that contains information on the simulation and what to compute. The latter may be done by modifying the example `run.py` script provided here. There are three locations in the `run.py` script where users will specify parameters. The *outputFile* variable specifies the prefix to use for output files. *rmsd.Trajectory* stores inputs for the following: name of the central ion (*IonID*), the elements to include in analyzing the coordination geometry and RDF/ADF (*elements*), the dimensions of the simulation box (*boxsize*), the width for the RDF histogram (*binSize*), the number of frame iterations to measure the RMSD to before updating the reference polyhedra ordering orientation (*framesForRMSD*), and the frame indices of the trajectory to begin (*startFrame*) and end (*endFrame*) at. The *traj.getAtoms* is used to specify the path location for the MD trajectory file. Executing the `run.py` script (e.g. `python run.py`) yield output files containing the structural and geometric information in .csv format.

## Examples
Demo files for using AFICS are provided under `example`. Two examples are provided: first, an MD simulation an of the Cr(III) ion within water; and second, an MD simulation of Eu(III) ion chelated by ethylenediaminetetraacetic acid (EDTA) within water. The MD trajectory files are provided and may be analyzed by executing the accompanying `run.py` scripts. While non-coordinated solvent and counterion molecules have been removed from these demo files in order to reduce the file sizes, removing non-first sphere molecules is not required to use AFICS.

## Citation
The journal article describing the AFICS toolkit is "Analysis of the First Ion Coordination Sphere: A Toolkit to Analyze the Coordination Sphere of Ions" published in the *Journal of Chemical Information and Modeling* and can be found at the following DOI: [10.1021/acs.jcim.3c00294](https://doi.org/10.1021/acs.jcim.3c00294)

