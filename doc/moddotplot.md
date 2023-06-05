## [ModDotPlot](https://github.com/marbl/ModDotPlot)


### parameters

```
$ moddotplot -h

 _______  _______  ______          _______  _        _______ _________
(       )(  ___  )(  __  \        (  ____ )( \      (  ___  )\__   __/
| () () || (   ) || (  \  )       | (    )|| (      | (   ) |   ) (   
| || || || |   | || |   ) |       | (____)|| |      | |   | |   | |   
| |(_)| || |   | || |   | |       |  _____)| |      | |   | |   | |   
| |   | || |   | || |   ) |       | (      | |      | |   | |   | |   
| )   ( || (___) || (__/  )   _   | )      | (____/\| (___) |   | |   
|/     \|(_______)(______/   (_)  |/       (_______/(_______)   )_(   


usage: moddotplot [-h] -i INPUT [INPUT ...] [-k KMER] [-s SPARSITY] [-r RESOLUTION] [-id IDENTITY] [-o OUTPUT]
                  [-nc NON_CANONICAL] [--no-bed] [--no-plot] [--num-colors NUM_COLORS] [--bin-freq]
                  [--interactive] [--port PORT] [-q] [-v]

Mod.Plot, A Rapid and Interactive Visualization of Tandem Repeats.

optional arguments:
  -h, --help            show this help message and exit

Required input:
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Path to input fasta file(s)

Mod.Plot distance matrix commands:
  -k KMER, --kmer KMER  k-mer length. Must be < 32 (default: 21)
  -s SPARSITY, --sparsity SPARSITY
                        Modimizer sparsity value. Higher value will reduce the number of modimizers, but will
                        increase performance. Default set to 2 per Mbp of sequence, rounded up to nearest even
                        integer) (default: None)
  -r RESOLUTION, --resolution RESOLUTION
                        Dotplot resolution. (default: 1000)
  -id IDENTITY, --identity IDENTITY
                        Identity cutoff threshold. (default: 80)
  -o OUTPUT, --output OUTPUT
                        Name for bed file and plots. Will be set to input fasta file name if not provided.
                        (default: None)
  -nc NON_CANONICAL, --non-canonical NON_CANONICAL
                        Only consider forward strand when computing k-mers. (default: False)

Static plotting commands:
  --no-bed              Don't output bed file. (default: False)
  --no-plot             Don't output image plots. (default: False)
  --num-colors NUM_COLORS
                        Number of colors to map. Must be < 15. (default: 11)
  --bin-freq            By default, histograms are evenly spaced based on the number of colors and the
                        identity threshold. Select this argument to bin based on the frequency of observed
                        identity values. (default: False)

Interactive plotting commands:
  --interactive         Launch a interactive Dash application on localhost. (default: False)
  --port PORT           Port number for launching interactive mode on localhost. Only used in interactive
                        mode. (default: 8050)

Logging options:
  -q, --quiet           Supress help text when running. (default: False)
  -v, --verbose         Add additional logging info when running. (default: False)

```


### manual

```

$ source activate /home/zhangyun2/anaconda3/envs/python_latest/

$ source /home/zhangyun2/software/ModDotPlot-main/venv/bin/activate


$ moddotplot -i  /home/wuzhikun/Project/Centromere/Centromere/SRF/Models/G42_srf_aln_cluster_target.fasta  -o /home/wuzhikun/Project/Centromere/Centromere/SRF/Models/G42_srf_aln_cluster_target_moddotplot

```


out files:

```
-rw-rw-r-- 1 wuzhikun wuzhikun  44M May 27 13:23 /home/wuzhikun/Project/Centromere/Centromere/SRF/Models/G42_srf_aln_cluster_target_moddotplot.bed
-rw-rw-r-- 1 wuzhikun wuzhikun 9.4M May 27 13:23 /home/wuzhikun/Project/Centromere/Centromere/SRF/Models/G42_srf_aln_cluster_target_moddotplot.png

```




