This repository contains code to perform Monte Carlo simulations of electrons
in a transition metal dichalcogenide heterobilayer (TMD) heterobilayer
Moire system in the strong coupling limit. For further details
see the associated manuscript: https://arxiv.org/abs/2112.08624

SYSTEM REQUIREMENTS:
  The software has the following dependencies, and has only been tested using the listed version numbers:
    Python >= 3.7.3,
    numpy >= 1.17.0,
    scipy >= 1.7.3,
    matplotlib >= 3.0.3.
  
  There is no required non-standard hardware.
  
INSTALLATION GUIDE:
  All one should have to do is clone the git repository.
  The dependencies can be installed by standard means e.g. pip or anaconda.
  Typical install time should not exceed a few minutes.

DEMO:
  A sample lattice for analysis is provided in test/test_lattices.py.
  Instantiating the NematicOP class in the Python REPL and calling the evaluate
  method on the sample lattices should return 0, and should take order milliseconds.
  The file test/test_mc.py calculates exact energies for a small lattice to compare to MC.
  Inputting the same l, n_particles, and t's into src/particle_mc.py should produce 
  almost identical sets of average energies, and take order 10's of minutes.
  
INSTRUCTIONS FOR USE: 

  RUNNING A SIMULATION:
    To do a simulation, run the file particle_mc.py It takes a system size,
    filling fraction numerator, and output filename prefix as command
    line arguments in that order. One may also wish to change the parameters
    detailed in the main function. The file run_sims.py auto-generates
    prefixes and runs simulations for multiple numerators simultaneously
    in the background.

  ANALYZING SIMULATION OUTPUT:
    The code outputs lattice configurations obtained at each MC step
    at which one collects data. Editing the 'prefix' variable in
    analysis.py will change which file is loaded, and MC averages
    of the nematic order parameter correlation function, orientational
    paameter, and structure factor are calculated. 
    One could add more operators as described below.

SYSTEM GEOMETRY:
  The simulations that this code performs take place on the triangular Moire
  lattice. The lattice is setup in a rhombus geometry. Due to the long-range
  interaction complicating simple periodic boundary conditions, I simulate
  a formally infinite system where particles interact within and between
  copies of the system out to some distance at which the interaction is small.
  This obeys the full point group symmetry of the lattice.

MONTE CARLO ALGORITHM:
  There are two types of Monte Carlo steps performed in this software.
  The first is a single-particle update that moves a particle to an
  unoccupied site according to standard Metropolis acceptance rules using
  the Boltzmann weight. The second is a cluster algorithm I made up based
  on the Geometric cluster algorithm and Wolff algorithm. For details about
  the algorithm and a proof of detailed balance, see appendix A of the 
  manuscript. The order two elements of the point group that I use
  are reflections.

  While the cluster algorithm is useful for escaping metastable
  states, it is much slower than single particle updates. As such I advise
  the user to perform mostly single-particle updates with some cluster updates
  mixed in.

OPERATOR CLASSES:
  The file operators.py contains the definitions of the operator classes
  that I use to compute the nematic order parameter and structure factor.
  A user could add additional operators by simply defining their own
  class with an evaluate method matching the call signature of those
  I've already defined.
