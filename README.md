# CR-WW
Study of color reconnection in WW->4q events at FCCee. Installing Pythia8 & Rivet is required. 

Clone repository from github

```
git clone git@github.com:izamiza10/FCCee_WW_CR.git
```
Install pythia & rivet. I have some example hepmc & yoda files generated for different Color reconnection (CR) schmes in ```HepMCExamples``` folder and programms used to generate the files can be found in ```PythiaPrograms```. Two simple rivet rountines can be found in ```rivet-ana```.


Installing Pythia8

```
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_106 x86_64-el9-gcc13-opt

tar -xzvf pythia8315.tar.gz

cd pythia8315

./configure \
    --with-lhapdf6=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-f90f6/x86_64-el9-gcc13-opt \
    --with-hepmc3=/cvmfs/sft.cern.ch/lcg/releases/hepmc3/3.2.7-77f68/x86_64-el9-gcc13-opt/ \
    --with-hepmc3-lib=/cvmfs/sft.cern.ch/lcg/releases/hepmc3/3.2.7-77f68/x86_64-el9-gcc13-opt/lib64 \
    --with-rivet=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/rivet/4.0.0-df6f7/x86_64-el9-gcc13-opt/ \
    --with-gzip
```

In the examples provided by Pyhtia main302.cc has some Color Reconnection (CR) checks. ```PythiaPrograms``` contains examples that can be used to study CR. To run the program update the Makefile if needed an simpy do:

```
make mainXXX
./mainXXX
```

Two examples of a simple rivet analysis are provided in ```rivet-ana```. To get a skelteon for a new rivet analysis one can simply run (after installing and setting up rivet)

```
rivet-mkanalysis MY_ANALYSIS
```
Then complie, run & plot:
```
rivet-build RivetMY_ANALYSIS.so MY_ANALYSIS.cc
rivet --analysis=MY_ANALYSIS examples/main305.hepmc -o QCDscheme.yoda
rivet-mkhtml QCDscheme.yoda
```
