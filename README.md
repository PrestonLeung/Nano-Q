# ONT-Tool
Python script for conservatively cleaning ONT reads from bam files and estimate variant frequencies.

### Brief Description
This python tool takes in a **bam file** performs a conservative cleaning procedure and then uses a hierachical clustering method based on hamming distance to identify potential variants. For each cluster formed using hierachical clustering, a consensus sequence is produced to represent that cluster and frequency of occurrence is measured based on reads per cluster over total reads extracted from the alignment.

Process of algorithm:            

1) For SNPS, it will check the base quality scores against user specified cutoff and substitute it with reference sequence used for the alignment in the bamfile if it is below the cutoff.

2) It will remove all indels and substitutes it to the reference sequence used for the alignment in the bamfile.

3) Stop codons check after cleaning is also performed. If a stop codon is identified the entire triplet of nucleotide is replaced with the reference used for the alignment. *Note - coding region is controlled by option '-c'*.

4) The reads are then trimmed to the read with the shortest length to keep output read lengths identical.

5) Pairwise hamming distance are then considered for all reads, which are then used in a distance matrix for hierarchical clustering.

### Requirements
Python version:
- Python3 (3.7)

These packages are required to run the tool:
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [numpy and scipy](https://www.scipy.org/install.html)
- [Biopython](https://biopython.org/wiki/Download)
- [matplotlib](https://matplotlib.org/users/installing.html)

*Note - Using [Anaconda](https://www.anaconda.com/distribution/) to create an environment (just do a simple Google search on "Create Python3 environment in Anaconda") for Python3.7 as well as installing all the above python packages might be easier to do (thats what I did!).*

### Tool Arguments
 
| Arguments | Description |
| --------- | ----------- |
|-b BAMFILE | Filename of bam file. |
|-c CODE_START | Start codon position in the reference sequence |
| -l READ_LENGTH | Length cut off for read size |
|-nr NUM_REF | Number of references used in the alignment. |
|-q QUAL_THRESHOLD | Base quality score cut off. |
|-j JUMP | Increase this to make larger read intervals. |
| -ht HD_THRESHOLD | Hamming distance threshold to call clusters. [Default = 234] |
|-mc MINREAD_CLUSTER | Minimum no. of reads to accept as a cluster. [Default = 30] |
|-ct CONSENSUS_THRESHOLD| Threshold to call a nucleotide to be consensus. [Default = 0.40] |
|-d, --dendrogram | Draw a dendrogram to help determine HD_Threshold (-ht) cutoff. |
|-hd, --keep_hdFile | Retain the Hamming Distances calculated. **Note:** Could take up a lot of space.|
|-kc, --keep_clusters | Retain the clustered reads. **Note:** Could take up alot of space.|

### Commandline
```sh
usage: indelRemover003G.py [-h] -b BAMFILE -c CODE_START -l READ_LENGTH -nr NUM_REF -q QUAL_THRESHOLD [-j JUMP] [-ht HD_THRESHOLD] [-mc MINREAD_CLUSTER] [-ct CONSENSUS_THRESHOLD] [-d] [-hd] [-kc]
```

#### Example
```sh
$  python indelRemover003G.py -b example.bam -c 1 -l 9000 -nr 1 -q 5 -j 10 
```


