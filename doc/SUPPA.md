
## [SUPPA](https://github.com/comprna/SUPPA)

a flexible and powerful tool to study splicing at the transcript isoform or at the local alternative splicing event level


### install suppa
```
conda install -c bioconda suppa
```


### parameters
```
usage: suppa.py [-h]
                {generateEvents,psiPerEvent,psiPerIsoform,diffSplice,clusterEvents,joinFiles}
                ...

Description:

SUPPA allows you to generate all the possible Alternative Splicing events from an annotation file, 
calculate the PSI values per event, calculate differential splicing across multiple conditions 
with replicates, and cluster events across conditions 
For further information, see the help of each subcommand.

positional arguments:
  {generateEvents,psiPerEvent,psiPerIsoform,diffSplice,clusterEvents,joinFiles}
    generateEvents      Calculates the Alternative Splicing events from an annotation file.
    psiPerEvent         Calculates the PSI value for each event previously generated.
    psiPerIsoform       Calculates the PSI value for each isoform.
    diffSplice          Calculates differentially spliced events across multiple conditions.
    clusterEvents       Calculates clusters of events across conditions.
    joinFiles           Join multiple tab separated files into a single file.

optional arguments:
  -h, --help            show this help message and exit

```


