# ONT-Tool
Python script for conservatively cleaning ONT reads from bam files and estimate variant frequencies.

### Brief Description
This python tool takes in a **bam file** performs a conservative cleaning procedure and then uses a hierachical clustering method based on hamming distance to identify potential variants. For each cluster formed using hierachical clustering, a consensus sequence is produced to represent that cluster and frequency of occurrence is measured based on reads per cluster over total reads extracted from the alignment.

This tool takes advantages of multiple cores if available, and in order to save RAM usage, temporary files will be created in the process. These temporary files will be deleted at the end of the run, however since storing pairwise hamming distance could take up a lot of space, the size and number of files could be large (100gb+).

##### Hardware used:
This algorithm was tested on a Linux machines: 

**Desktop**
*Operating System: Ubuntu 14.04
Memory: 32gb
Processor: Intel Core i7-4790 3.60Ghz*

**Internal Server**
*Operating System: RedHat Enterprise Linux Server Release 6.9 Santiago) 
Memory: 529 gb
Processors: 4x Intel Core Xeon E5-4650L 2.60GHz*

**Process of algorithm:**

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

[Minimap2](https://github.com/lh3/minimap2) and [GraphMap](https://github.com/isovic/graphmap) were used as aligners for sam/bam files. For Minimap2, *--MD* option needs to be turned on for the MD tag to appear in the alignment file, a requirement for ONT-Tool.

### Tool Arguments
 
| Arguments | Description |
| --------- | ----------- |
|-b BAMFILE | Filename of bam file. |
|-c CODE_START | Start codon position in the reference sequence |
| -l READ_LENGTH | Length cut off for read size |
|-nr NUM_REF | Number of references used in the alignment. |
|-q QUAL_THRESHOLD | Base quality score cut off. |
|-j JUMP | Increase this to make larger read intervals. Outputs less number of files but larger in size for the occasion when there's an upper limit to how many files are allowed to be opened for writing at the same time. |
| -ht HD_THRESHOLD | Hamming distance threshold used to call clusters. [Default = 234] |
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

### Algorithm Outputs
The algorithm will create a folder called **'Results'**, and within that folder there would be pairs of **.fa** files associated with the number of references (option -nr). For instance, if you've used 2 references in the bam file, then there should be 2 pairs of fasta files.

The first file will have suffix **'_ClusterConsensus.fa'**, which contains the estimated viral variants of the reference used (which will be the prefix of the fasta file). E.g. **Subject1_ClusterConsensus.fa** will contains viral variant estimation from reads aligned to Subject1 reference. The second file will simply be a collection of cleaned and trimmed reads used to identify clusters.

# Additional Help

##### How to set up Anaconda
1) [How](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/) to set up an environment in Anaconda with python.
2) How to install python packages after activatnig python environment in Anaconda
    a. [Video Tutorial](https://www.youtube.com/watch?v=Z_Kxg-EYvxM)
    b.  ```conda install [package_name]``` or ```pip install [package_name]```

##### Aligners
To align sequences and get a sam file, you can use either:

1) [Minimap2](https://github.com/lh3/minimap2)
2) [Graphmap](https://github.com/isovic/graphmap)

Here are some usage examples for:

1) [Minimap2 Examples](https://github.com/lh3/minimap2#getting-started)
2) [Graphmap Examples](https://github.com/isovic/graphmap#usage-examples)

Then you can transfer sam to bam files as well as any other alignment file manipulation through [samtools](http://www.htslib.org/download/).

##### Quick guide for samtools v1.9

Full guide [here](http://www.htslib.org/doc/samtools.html).


To convert sam file to bam file:

    samtools view -bt -o output_aln.bam input_aln.sam
    -b -> output bam format.
    -t -> tab delimited.

To sort a bam file:

    samtools sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam
    -T -> write temporary files to this location (if directory exists) using the given prefix.

To index a bam file:
    
    samtools index aln.sorted.bam 

This ReadMe file was written using [Dillinger](https://dillinger.io/).