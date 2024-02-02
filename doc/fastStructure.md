## [fastStructure](http://rajanil.github.io/fastStructure/)


### install

```
conda install -c bioconda faststructure
```

```
$ python /home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py -K 3 --input /home/wuzhikun/Project/PanSV2/NGS/Variant/Sample_variant_PruneLD --output test
Traceback (most recent call last):
  File "/home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py", line 4, in <module>
    import fastStructure 
  File "__init__.pxd", line 885, in init fastStructure
ValueError: numpy.ufunc has the wrong size, try recompiling. Expected 192, got 216

```

answer:
```
$ pip insatll -U numpy

```


### parameters

```
$ python /home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py --help
Incorrect options passed

Here is how you can use this script

Usage: python /home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py
     -K <int> (number of populations)
     --input=<file> (/path/to/input/file)
     --output=<file> (/path/to/output/file)
     --tol=<float> (convergence criterion; default: 10e-6)
     --prior={simple,logistic} (choice of prior; default: simple)
     --cv=<int> (number of test sets for cross-validation, 0 implies no CV step; default: 0)
     --format={bed,str} (format of input file; default: bed)
     --full (to output all variational parameters; optional)
     --seed=<int> (manually specify seed for random number generator; optional)

```


```
$ python /home/wuzhikun/anaconda3/envs/Assembly/bin/distruct.py --help
Incorrect options passed

Here is how you can use this script

Usage: python /home/wuzhikun/anaconda3/envs/Assembly/bin/distruct.py
     -K <int>  (number of populations)
     --input=<file>  (/path/to/input/file; same as output flag passed to structure.py)
     --output=<file> (/path/to/output/file)
     --popfile=<file> (file with known categorical labels; optional)
     --title=<figure title> (a title for the figure; optional)

```


```
$ python /home/wuzhikun/anaconda3/envs/Assembly/bin/chooseK.py --help
Incorrect options passed

Here is how you can use this script

Usage: python /home/wuzhikun/anaconda3/envs/Assembly/bin/chooseK.py
     --input=<file>

```



