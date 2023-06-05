
### ultra parameters

```
$ ./ultra
No arguments recieved. 
usage](ultra [options] <input file>

available options:
-f  Output File Path](                  Path to output JSON data. [STDOUT]
-s  Score Threshold](                   Minimum score necessary for a repeat to be recorded. [-10000.000000]
-ml Region Length Trheshold](           Minimum total repeat length necessary for a repeat to be recorded. [10]
-mu Repeat Unit Threshold](             Minimum number of repeat units necessary for a repeat to be recorded. [3]
-at AT Richness](                       Frequency of As and Ts in the input sequence. [0.600000]
-atcg   ATCG Distribution](                 Frequency of As Ts Cs and Gs in the input sequence. [A=0.300000,T=0.300000,C=0.200000,G=0.200000]
-m  Match Probability](                 Probability of two nucleotides in consecutive repeat units being the same. [0.800000]

-n  Number of Threads](                 Number of threads used by ultra to analyze the input sequence. [1]
-p  Maximum Period](                    Largest repeat period detectable by Ultra. [15]
-mi Maximum Consecutive Insertions](    Maximum number of insertions that can occur in tandem. [8]
-md Maximum Consecutive Deletions](     Maximum number of deletions that can occur in tandem. [7]

-nr Repeat Probability](                Probability of transitioning from the nonrepetitive state to a repetitive state. [0.010000]
-rn Repeat Stop Probability](           Probability of transitioning from a repetitive state to the nonrepetitive state. [0.050000]
-pd Repeat probability decay](          **reminder to write a good description later. [0.850000]
-ri Insertion Probability](             Probability of transitioning from a repetitive state to an insertion state. [0.020000]
-rd Deletion Probability](              Probability of transitioning from a repetitive state to a deletion state. [0.020000]
-ii Consecutive Insertion Probability]( The probability of transitioning from an insertion state to another insertion state. [0.020000]
-dd Consecutive Deletion Probability](  The probability of transitioning from a deletion state to another deletion state. [0.020000]

-sr ***NOT FUNCTIONAL, DONT USE***](   Split repeats such as ATATATGCGCGCGC into subrepeats ATATAT and GCGCGCGC. [0]
-sd Split Depth](                       Number of repeat units to consider when splitting repeat. [5]
-sc Split Cutoff](                      Cutoff value used during splitting (smaller is more conservative). [2]
-msp    Split Max Period](                  Maximum repeat period that will be considered for splitting. [6]

-hs Hide Repeat Sequence](              Don t show the repetitive sequences in the results JSON. [Show Repeat Sequence]
-ss Show score Deltas](                 Output the score change per residue. [False]
-st Show traceback](                    Output the Viterbi traceback in the results JSON. [Hide traceback]
-wid    Show Window ID](                    Display the windowID corresponding to a repeat in the results JSON. [Hide Window ID]
-sl Show logo numbers](                 Output the corresponding logo annotation for a given repeat. [Hide logo numbers]
-json   Read JSON file](                    Process all passes in JSON file. [False]
-jpass  Process passes in JSON file](       Process selected passes in JSON file. [None]
-pid    Pass ID](                           Assigns a custom pass ID. [Smallest unused positive pass ID]
-R  Completely Read File](              Read the entire input file during initialization.. [Do not read whole file]
-ws Window Size](                       The number of nucleotides per sequence window. [8192]
-os Window Overlap](                    The number of nucleotides overlaped between two consecutive windows. [100]
-wn Number of Windows](                 Maximum number of windows stored in memory at once. [1024]
-v  Ultra Version](                     Shows Ultra s version. []

-doc    Debug overlap correction](          Report overlap correction in repeat information. [0]
```

```
$ /home/wuzhikun/software/ULTRA/ultra -mi 2 -md 2 -p 4001 -mu 2 -ws 100000 -os 10000  -n 20 -f /home/wuzhikun/Project/Vigna/scaffold/final3/centromere/ultra/cen_ultra_out.txt /home/wuzhikun/Project/Vigna/scaffold/final3/centromere/centromere_region.fasta
```


