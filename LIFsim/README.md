# Code for figures 4 and 4S

Supplied as a Jupyter Notebook: SpikingNeuronsSims.ipynb

The code is a python3 code using a custom C++ python extension.
The C++ extension is only used to efficiently generate Poisson spike-trains from
given input rates and to simulate the LIF dynamics.
The C++ extension relies on the Boost-Python library.
It would be possible to run the code without it by implementing
these two (relatively) simple functions.

## Building the C++ python extension (required for simulations):

```
cd source
make
mv _ST.so ../SpikingTempotron
make clean #(Optional)
```

The given makefile works on my setup (Ubuntu, with GCC and native Python3).
Compiler and linker parameters should probably be modified to fit different
systems.
