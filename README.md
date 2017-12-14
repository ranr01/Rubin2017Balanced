# Code repository for "Balanced excitation and inhibition are required for high-capacity, noise-robust neuronal selectivity"
## Ran Rubin, L. F. Abbott and Haim Sompolinsky

### Content

* QuadProg - Python code for finding maximal
$\kappa_\mathrm{in}$ and maximal $\kappa_\mathrm{out}$ weights for a sign and
norm constrained Perceptron.

* Theory - Python code for numerical solution of the saddle point equations.

* Recurrent - Python code for simulation of recurrent associative memory networks

* PerceptronLearning - Julia code for perceptron learning in the presence of
input and output noise.

* LIFsim - Simulations of LIF neurons performing classification with Balanced
and unbalances weights (Python w. custom C++ extention module).

Example code is given in Jupyter notebooks in QuadProg, Theory, Recurrent and
LIFsim. 

### Dependencies

* Python 3
* Jupyter
* Numpy
* CVXOPT
* Julia (0.6.1, For Perceptron Learning)
* Boost Python (For LIFsim)

Leagal notice:

THIS CODE MAY BE COPYRIGHTED BY THE AUTHOR (RAN RUBIN) AND/OR THE HEBREW UNIVERSITY IN JERUSALEM, AND/OR COLUMBIA UNIVERSITY. IT IS SUPPLIED, AS IS, FOR ACADEMIC PURPOSES WITHOUT ANY GUARANTY OF ANY KIND. THE AUTHOR OR ANY OTHER ORGANIZATION THE AUTHOR IS OR WAS AFFILIATED WITH, WILL NOT BARE ANY RESPONSIBILITY OR LIABILITY FOR ANY KIND OF DAMAGES RESULTING FROM THE USE OF THIS CODE.
