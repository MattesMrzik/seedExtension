# seedExtension
This is the code of my bachelor's thesis, with the title "Extension of Seeds for Multiple Alignments in Genome Graphs".


Supervisor: [Prof. Dr. Mario Stanke](https://math-inf.uni-greifswald.de/en/department/about-us/employees/prof-dr-mario-stanke-english/), the inventor of the gene prediction tool [Augustus](http://bioinf.uni-greifswald.de/augustus/).

## Abstract
The number of sequenced genomes is rapidly increasing. However, the annotation process for these genomes is more elaborate and therefore slower. In order to reduce this discrepancy, new robust methods of genome annotation must be developed. The Geometric Hashing Algorithm ([GH](https://pubmed.ncbi.nlm.nih.gov/35689182/)), extended in the course of this work by the Seed Extension Algorithm (SE), addresses this task. The seeds created in GH are intended to be filtered, and this is implemented with a gapless seed extension. This idea can be easily realized in linearly stored sequences. In this work, however, a [colored de Bruijn graph](https://academic.oup.com/bioinformatics/article/33/20/3181/2995815) is utilized, which can have storage advantages. SE transfers the concept of gapless seed extension to the graph. The filtered seeds are then intended to reliably identify coding sequences and thus support the annotation process.

## Dependencies
[https://github.com/mabl3/metagraphInterface](https://github.com/mabl3/metagraphInterface)
