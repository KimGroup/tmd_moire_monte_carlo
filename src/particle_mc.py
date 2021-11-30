import sys, numpy as np
from fractions import Fraction
from operators import *
from lookup_tables import *
from string import ascii_lowercase
import mc_steps
import time

# function to initialize a lattice simulation cell in rhombus geometry
# ARGS:
#   lx <- [int] horizontal width of lattice in units of sites
#   ly <- [int] vertical width of lattice in units of sites
#   ff <- [Fraction instance] filling fraction of lattice
#   reuse <- [bool] optional, attempts to load lattice from file if
#            True, makes a fresh lattice if False, default False
#   faname <- [string,None] optiona, file name to load from, default None
# RETURNS:
#   2D Numpy int array of 0/1 indicating un/occupied lattice sites
def initialize_lattice(lx,ly,ff,reuse=False,fname = None):
  if reuse:
    try:
      lat = np.loadtxt(fname)
      return lat
    except:
      print("Requested file not found. Initializing fresh lattice")
  lat = np.zeros((ly,lx),dtype='int')
  n_sites = lx*ly
  # make sure filling fractor corresponds to integer number of electrons
  num = ff.numerator*n_sites
  assert((num % ff.denominator) == 0)
  n_filled = num // ff.denominator
  indices = [ind for ind in np.ndindex(lat.shape)]
  # occupy random sites according to ff
  for i in np.random.choice(len(indices),n_filled,replace=False):
    lat[indices[i]] = 1
  return lat

# function to perform Monte Carlo simulation
# ARGS:
#   lat <- [2D Numpy int array] array of 0/1 indicating un/occupied sites
#   temperature <- [float] simulation temperature in units of K_B
#   sweeps_per_site <- [int] number of MC sweeps per site to perform
#   collection_frequency <- [int] collect data every collection_frequency sweep
#   ops <- [list of instances of classes that inherit from Operator]
#          operators to evaluate when you collect data
#   field <- [float] strain field value, isn't working
# RETURNS:
#   Lattice after simulation and ops
def monte_carlo(lat, temperature, sweeps_per_site, collection_frequency, ops=[], field=None):
  beta = 1.0 / temperature
  l = lat.shape[0]
  n_sweeps = l*l*sweeps_per_site
  # Build lookup tables to avoid repeated computations
  energy_table = make_energy_table(lat.shape,field=field)
  pair_energy_table = make_pair_energy_table(l,field=field)
  transformation_table = make_transformation_table(l)
  neighbor_table = make_nn_table((l,l))
  lat = np.ravel(lat)
  
  sweep = 0
  collected = False # did you collect data this sweep?
  mode = "single_particle" # whether to do single_particle or cluster udpate
  cluster_frequency = 1001 # how often to do a cluster update

  while sweep < n_sweeps:
    # figure out what kind of sweep to do and do it. 
    # this is kind of gross right now, user probably shouldn't have to 
    # interact with this
    mode = "cluster" if (((sweep//(l*l)) % cluster_frequency) == 0) else "single_particle"
    #mode = "single_particle"
    if mode == "single_particle":
      #lat, log = mc_steps.single_particle_step(lat, l, beta, energy_table)
      lat, log = mc_steps.short_range_step(lat,l,beta,energy_table,neighbor_table)
    elif mode == "cluster":
      lat, log = mc_steps.cluster_step(lat, l, beta, pair_energy_table,\
                                       transformation_table)
    else: assert(False)

    if log["non_trivial"]: 
      sweep += 1
      collected = False
    if (not collected) and ((sweep % collection_frequency) == 0):
      for op in ops: op.evaluate(lat)
      collected = True

    #if sweep // (n_sweeps // 4) % 2 == 0: mode = "single_particle"
    #else: mode = "cluster"

  return lat.reshape((l,l)), ops

# main function to setup and run the simulation
# COMMAND LINE ARGS (IN ORDER):
#   l < - [int] linear dimension of lattice simulation cell
#   num <- [int] filling fraction numerator
#   prefix <- [string] output file name prefix
def main():
  l = int(sys.argv[1])
  num, denom = int(sys.argv[2]), l*l
  # user might want to change these
  collection_frequency = int(5e2)
  n_equil, n_sample = int(1e4), int(1e4)
  # MAKE SURE TEMPERATURES ARE DECREASING
  temperatures = np.linspace(0.08,0.001,20)

  # Set up the lattice
  filling_fraction = Fraction(num,denom)
  lat = initialize_lattice(l,l,filling_fraction)

  # Set up file
  prefix = sys.argv[3]
  operator_fname="saved_configs/"+prefix+"_l{l}_num{n}_denom{d}_temp{t}"
  
  # Do Monte Carlo
  for temp in temperatures:
    fname = operator_fname.format(l=l,n=num,d=denom,t=temp)
    # Equilibrate
    lat,_ = monte_carlo(lat,temp,n_equil,collection_frequency)
    # Sample
    config = LatticeConfig(l,sum_evals=False)
    ops = [config]
    lat,ops = monte_carlo(lat,temp,n_sample,collection_frequency,ops=ops)
    np.save(fname,ops[0].state)

if __name__ == '__main__': main()
