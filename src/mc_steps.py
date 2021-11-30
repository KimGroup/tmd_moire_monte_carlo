import numpy as np
import time
from collections import deque

# compute the energy contributed to the Hamiltonian by a single site
# ARGS:
#   lat <- Numpy array of 0/1, recording un/occupied lattice sites
#   site <- int index corresponding to lattice site
#   l <- linear dimension of lattice simulation cell
#   energy_table <- Lookup table constructed by make_energy_table
# RETURNS:
#   energy contributed to Hamiltonian by input site
def local_energy(lat,site,l,energy_table):
  sr,sc = site//l,site%l # site row and column
  site2 = np.arange(l*l)
  site2 = site2[lat[site2] != 0]
  site2 = ((sc - site2 % l)%l + l)%l + (((sr - site2 // l)%l + l)%l)*l
  energy = np.sum(energy_table[site2])
  return energy

# perform single cluster update
# ARGS:
#   lat <- Numpy array of 0/1, recording un/occupied lattice sites
#   l <- linear dimension of lattice simulation cell
#   beta <- inverse temperature
#   pair_energy_table <- Lookup table constructed by make_pair_energy_table
#   transformation_table <- Lookup table from make_transformation_table
# RETURNS:
#   updated lattice and dictionary recording whether step was trivial
def cluster_step(lat, l, beta, pair_energy_table, transformation_table):
  # pick a random transformation and cluster seed site
  transformation_ind = np.random.choice(transformation_table.shape[0])
  new_inds = transformation_table[transformation_ind,:]
  occ_diff = lat - lat[new_inds]
  seed_site = np.random.choice(l*l)
  if occ_diff[seed_site] == 0: 
    # if seed site is unoccupied update is trivial
    log_dict = {"non_trivial":False,"cluster_size":None}
    return lat,log_dict
  
  # Compute acceptance probabilities as laid out in Appendix A of MS
  avail = np.where(occ_diff != 0)[0]
  delta = pair_energy_table[new_inds,:] - pair_energy_table
  delta *= np.outer(occ_diff,occ_diff)
  delta /= 2
  delta = np.exp(-beta*delta)
  
  # Build cluster (see appendix A of MS)
  cluster = np.zeros(l*l,dtype='bool')
  cluster[seed_site] = True
  cluster[new_inds[seed_site]] = True
  stack = deque([seed_site]) 
  lat[seed_site],lat[new_inds[seed_site]] = lat[new_inds[seed_site]],lat[seed_site]
  while len(stack) > 0:
    i = stack.pop()
    # possible optimizations:
    # vectorize this loop. compute all random #'s at once
    # replace cluster with boolean array, do np.where to pick out relevant
    # sites
    for k in avail:
      if (not cluster[k]) and (np.random.random() < (1-delta[i,k])):
        cluster[k] = True
        cluster[new_inds[k]] = True
        lat[k],lat[new_inds[k]] = lat[new_inds[k]],lat[k]
        stack.append(k)
  log_dict = {"non_trivial":True,"cluster_size":len(cluster)} # TODO: Fix
  return lat,log_dict

# perform one single particle update
# ARGS:
#   lat <- Numpy array of 0/1, recording un/occupied lattice sites
#   l <- linear dimension of lattice simulation cell
#   beta <- inverse temperature
#   energy_table <- Lookup table constructed by make_energy_table
# RETURNS:
#   updated lattice and dictionary recordin whether step was trivial
def single_particle_step(lat,l,beta,energy_table):
  # pick random site to move from and move to
  site = np.random.choice(l*l)
  proposed_site = np.random.choice(l*l)
  if lat[proposed_site] == lat[site]: 
    # if both occupancies are same, move is trivial
    log_dict = {"non_trivial":False,"acceptance":None}
    return lat, log_dict
  
  # determine acceptance using standard Metropolis rules
  if lat[site] == 0:
    pre_move_local_energy = local_energy(lat,proposed_site,l,energy_table)
  else:
    pre_move_local_energy = local_energy(lat,site,l,energy_table)
  
  lat[proposed_site], lat[site] = lat[site], lat[proposed_site]
  
  if lat[site] == 1:
    post_move_local_energy = local_energy(lat,site,l,energy_table)
  else:
    post_move_local_energy = local_energy(lat,proposed_site,l,energy_table)

  delta_e = post_move_local_energy - pre_move_local_energy
  if delta_e < 0: accepted = True
  else: accepted = np.random.random() < np.exp(-beta*delta_e)
  if not accepted: 
    lat[proposed_site], lat[site] = lat[site], lat[proposed_site]
  log_dict = {"non_trivial":True,"acceptance":accepted}
  return lat, log_dict

# perform one single particle update, exchanging only within a few nn's
# ARGS:
#   lat <- Numpy array of 0/1, recording un/occupied lattice sties
#   l <- linear dimension of lattice simulation cell
#   beta <- inverse temperature
#   neighbor_table <- Lookup table constructed by make_nn_table
def short_range_step(lat,l,beta,energy_table,neighbor_table):
  # pick random ste to move from and move to
  site = np.random.choice(np.where(lat==1)[0])
  neighbors = neighbor_table[site]
  available_neighbors = neighbors[lat[neighbors] == 0]
  if len(available_neighbors) == 0:
    # If occupancies are same the move is trivial
    log_dict = {"non_trivial":False,"acceptance":None}
    return lat, log_dict
  
  # Determine acceptance based on standard Metropolis rules
  proposed_site = np.random.choice(available_neighbors)
  pre_move_local_energy = local_energy(lat,site,l,energy_table)
  lat[proposed_site],lat[site]=1,0
  post_move_local_energy = local_energy(lat,proposed_site,l,energy_table)
  
  delta_e = post_move_local_energy - pre_move_local_energy
  if delta_e < 0: accepted = True
  else: accepted = np.random.random() < np.exp(-beta*delta_e)
  if not accepted:
    lat[proposed_site],lat[site] = 0,1
  log_dict = {"non_trivial":True,"acceptance":accepted}
  return lat, log_dict
