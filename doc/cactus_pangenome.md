
## [cactus](https://github.com/ComparativeGenomicsToolkit/cactus)


### install cactus-pangenome


download pre-compiled binaries
```

wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.5.2/cactus-bin-v2.5.2.tar.gz
```

set the environment

```
$ cd /home/wuzhikun/software/cactus-bin-v2.5.2

$ virtualenv -p python3 cactus2.5
Running virtualenv with interpreter /home/wuzhikun/anaconda3/bin/python3
Using base prefix '/home/wuzhikun/anaconda3'
/usr/lib/python2.7/site-packages/virtualenv.py:1038: DeprecationWarning: the imp module is deprecated in favour of importlib; see the module s documentation for alternative uses
  import imp
New python executable in /home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/python3
Not overwriting existing python script /home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/python (you must use /home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/python3)
Installing setuptools, pip, wheel...done.


$ echo "export PATH=$(pwd)/bin:\$PATH" >> cactus2.5/bin/activate

$ echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus2.5/bin/activate

$ source cactus2.5/bin/activate
$ python3 -m pip install -U setuptools pip

$ python3 -m pip install -U .
$ python3 -m pip install -U -r ./toil-requirement.txt
```


### run  cactus

```
$ source /home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/activate

$ cactus-pangenome /home/wuzhikun/PublicData/Arabidopsis/genomes /home/wuzhikun/PublicData/Arabidopsis/genomes/genomeFile.txt  --outDir /home/wuzhikun/PublicData/Arabidopsis/genomes/cactus  --outName cactus_align --reference Col-0 --restart


[2023-06-07T17:21:10+0800] [MainThread] [I] [toil.statsAndLogging] Enabling realtime logging in Toil
[2023-06-07T17:21:10+0800] [MainThread] [I] [toil.statsAndLogging] Cactus Command: /home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/cactus-pangenome /home/wuzhikun/PublicData/Arabidopsis/genomes /home/wuzhikun/PublicData/Arabidopsis/genomes/genomeFile.txt --outDir /home/wuzhikun/PublicData/Arabidopsis/genomes/cactus --outName cactus_align --reference Col-0 --restart
[2023-06-07T17:21:10+0800] [MainThread] [I] [toil.statsAndLogging] Cactus Commit: a33d3eabb909873746ecd8e7e1528344e526d95b
Traceback (most recent call last):
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/toil/jobStores/fileJobStore.py", line 689, in read_shared_file_stream
    with open(self._get_shared_file_path(shared_file_name), 'rb' if encoding == None else 'rt', encoding=encoding, errors=errors) as f:
FileNotFoundError: [Errno 2] No such file or directory: '/home/wuzhikun/PublicData/Arabidopsis/genomes/files/shared/config.pickle'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/bin/cactus-pangenome", line 8, in <module>
    sys.exit(main())
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/cactus/refmap/cactus_pangenome.py", line 157, in main
    with Toil(options) as toil:
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/toil/common.py", line 939, in __enter__
    jobStore.resume()
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/toil/jobStores/fileJobStore.py", line 131, in resume
    super().resume()
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/toil/jobStores/abstractJobStore.py", line 212, in resume
    with self.read_shared_file_stream('config.pickle') as fileHandle:
  File "/home/wuzhikun/anaconda3/lib/python3.7/contextlib.py", line 112, in __enter__
    return next(self.gen)
  File "/home/wuzhikun/software/cactus-bin-v2.5.2/cactus2.5/lib/python3.7/site-packages/toil/jobStores/fileJobStore.py", line 694, in read_shared_file_stream
    raise NoSuchFileException(shared_file_name)
toil.jobStores.abstractJobStore.NoSuchFileException: File 'config.pickle' does not exist.

```
