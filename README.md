# Computational Biology Assignment 4
## Written by Ray Loerke

This program is a BioGRID MITAB File Processor.
A Network file (plaintext) is generated containing the filtered protein-protein interaction network.
The name and extension for theis file is determined by the user.
Protein-protein interactions are tab seperated and use their Entrez gene names.
Summary statistics are also listed at the top of the Network file.

The filtering criteria are as follows, each interaction added into the network must:
* Consist of 2 human (taxid: 9606) proteins
* Be a physical interaction based on its MI code
* Be experimentally detected based on its MI code
* Not be a duplicate of another interaction
* Not be a self-loop

See the example files for the formatting of these file types. 
Generally, a leading '<' means a sequence follows and a leading '#' means a comment follows.
