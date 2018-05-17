# IRCalign

___________________________
Version = 0.3.0

Dennis Svatunek

TU Wien

dennis.svatunek@gmail.com
___________________________

Align and combine two multiple coordinate \*xyz files, e.g. IRCs, to the same coordinate system for better visualization.

Additional features are:
* removal of duplicate structures at the intersection between combined sequences
* automatic detection and reversion of the direction of structure sequences to allow for best alignment
* Using center of mass rather than centroid as origin

## How is this done?

All structures are centered first. Using the Kabsch algorithm the rotation needed for aligning the second set of structures to the first set is calculated using the last structure in the first set and the first structure in the second set. Then all structures in the second set are rotated accordingly.

When desired, a partial procrustes analysis using the Kabsch algorithm is performed between the first and last structures of both sets of structures to determine the correct orientations of the structure sequences. This allows for inversion of one or both structures to ensure correct overlap between those two sets. This set is performed after centering the structures but before calculating the optimal rotation.

Additionally, removal of duplicates at the intersection of both sets is possible. This is sueful when combining forward and reverse IRC which both start at the transition state.

## Usage

IRCalign.py \[options] IRC_1 IRC_2 

IRC_1 is the first multiple structure \*.xyz file, it will be the reference and structures will not be rotated

IRC_2 is the second multiple structure \*.xyz file, structures will be aligned to IRC_1

Options:

-n Name   allows to specify the output name, otherwise a name will be created from the input file names
-o        requests automatical orientation of the sets of structures by reversing the order if needed
-d        removes duplicate at the intersection between sets
-m        uses center of mass rather than centroid for centering structures

## Limitations

Currently only \*.xyz files are allowed.
All structures must have the same atoms in the same order.

## To-Do

Ignore hydrogens during alignment.
Ignore custom atoms during alignment.
Allow for more than one alignment at a time.

## Contact

Dennis Svatunek

TU Wien

dennis.svatunek@gmail.com

https://github.com/dsvatunek/IRCalign

_____________________________

License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
