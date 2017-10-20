# LiverICG-CompartmentModel
Liver ICG description with mathematical model based on ODE

### Requirements
librairies: lapack et libm installed
to plot the solution: python

Sundials

1/ install sundials libs:
```console
gunzip sundials-2.7.0.tar.gz
tar xvf  sundials-2.7.0.tar
mkdir  -p Librairies
cd Librairies
cmake -DCMAKE_INSTALL_PREFIX=. -DEXAMPLES_ENABLE=OFF ../sundials-2.7.0
make
make install
```
2/ make a root level to compile the code

### How to compile and run 

Once the requiered librairies the code can be compile and run. 

to compile :
```console
make clean
make
```

to run :
```console
./output
```

### Model parameters, sensitivity analsysis and observations
 * the result file `data.txt` is in the folder `./Res/`
to plot the output of the model a python script called `SolutionPlot.py`
to plot:
```console
python SolutionPlot.py
```

 * To change the parameters of the model the file `Model.c` needs to be modified. In the first lines of this file you can choose the parameters of the model and the result file name. Don't forger to re-compile the code before running it.

 * The files for sensitivity analysis are in `./SensitivityAnalysis` folder. You need to compile and run the sensitivity analysis code, then you can plot the resulting senstivity function with the python script `SensitivityPlot.py`. The parameters of the model can also be changed by modifying the file ```Model_Sensitivty.c```, don't forget to re-compile the code before running it. 

 * The observation files are from El-Desoky et. al. 1999, and are in `./Obs`



