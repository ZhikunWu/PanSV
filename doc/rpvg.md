
## [rpvg](https://github.com/jonassibbesen/)

### install
```
git clone --recursive https://github.com/jonassibbesen/rpvg.git
cd rpvg && mkdir build && cd build
cmake ..
make -j <threads>

```

Error
```
make[3]: Leaving directory '/home/wuzhikun/software/rpvg/deps/sdsl-lite/build'
SUCCESS: sdsl was installed successfully!
The sdsl include files are located in '/home/wuzhikun/software/rpvg/deps/sdsl-lite/include'.
The library files are located in '/home/wuzhikun/software/rpvg/deps/sdsl-lite/lib'.
 
Sample programs can be found in the examples-directory.
A program 'example.cpp' can be compiled with the command: 
g++ -std=c++14 -DNDEBUG -O3 [-msse4.2] \
   -I/home/wuzhikun/software/rpvg/deps/sdsl-lite/include -L/home/wuzhikun/software/rpvg/deps/sdsl-lite/lib \
   example.cpp -lsdsl -ldivsufsort -ldivsufsort64
 
Tests in the test-directory
A cheat sheet in the extras/cheatsheet-directory.
Have fun!
[  5%] Completed 'sdsl-lite-proj'
make[2]: Leaving directory '/home/wuzhikun/software/rpvg/build'
[ 15%] Built target sdsl-lite-proj
make[1]: Leaving directory '/home/wuzhikun/software/rpvg/build'
make: *** [Makefile:95: all] Error 2

```


```
docker run -d -ti -v /home/wuzhikun/Project/RNA/Test jsibbesen/rpvg > docker_id_rpvg.logl

```

