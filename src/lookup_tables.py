import numpy as np
from scipy.special import kn

# builds a lookup table mapping a site index to an np array of its first
# few nearest neighbors for the simulated triangular latice
# ARGS:
#   shape <- length 2 tuple of ints corresponding to lattice shape
#   max_neighbor <- either 1,2,or 3, indicating to consider neighbors up to
#                   first, second, or third nearest neighbors
# RETURNS:
#   2D np array where row index is site and columns are neighbors
def make_nn_table(shape, max_neighbor = 2):
  assert(1 <= max_neighbor <= 3)
  nn_table = np.zeros((shape[0]*shape[1],6*max_neighbor),dtype='int')
  mod_index=lambda index: np.array([(i%s + s)%s for i,s in zip(index, shape)])
  # locations of first and  second nns relative to a site
  first_nns = [[0,-1],[1,0],[1,1],[0,1],[-1,0],[-1,-1]]
  second_nns = [[1,-1],[2,1],[1,2],[-1,1],[-2,-1],[-1,-2]]
  third_nns = [[0,-2],[2,0],[2,2],[0,2],[-2,0],[-2,-2]]
  all_nns = first_nns
  if max_neighbor > 1: all_nns += second_nns
  if max_neighbor > 2: all_nns += third_nns
  all_nns = [nn for nn in map(lambda x:np.array(x),all_nns)]

  for site in np.ndindex(shape):
    site_neighbors = [mod_index(site+nn) for nn in all_nns]
    site_neighbors = [idx[1]+idx[0]*shape[1] for idx in site_neighbors]
    nn_table[site[1]+site[0]*shape[1],:] = np.array(site_neighbors)
  return nn_table

# builds a lookup table corresponding to a vectorized matrix of the same
# shape as the lattice where each entry corresponds to the interaction
# energy between an electron at that site and an electron at the origin
# ARGS:
#   shape <- length 2 tuple of ints corresponding to lattice shape
#   d <- dielectric gate distance in units of Moire lattice constant
#   field <- strain field strength. Currently broken
# RETURNS:
#   np array of energies relative to origin
def make_energy_table(shape, d=10, field=None):
  if (shape[0] != shape[1]): assert(False)
  l = shape[0] # I'm just going to only do this for lx = ly right now
  # lattice vectors. I checked distances up to 3nn. Seems fine.
  a1, a2 = np.array([0,1]), np.array([np.sqrt(3)/2.0,-0.5])
  dist = lambda site: np.sqrt(np.sum(np.power(site[1]*a1+site[0]*a2,2)))
  # potential as in the MS
  v = lambda r: (4/d)*np.sum(kn(0,np.pi*r/d*np.arange(1,101,2)))
  #    if np.allclose(r,1) else 0
  #v = lambda r: 1 if np.allclose(r,1) else 0
  energies = np.zeros(shape)
  for site in np.ndindex(shape):
    if site == (0,0): continue
    energy = v(dist(site))

    # ADDED FOR STRAIN
    if not (field is None):
      if (site == (0,1)) or (site == (0,shape[1]-1)):
        energy += 0
      elif (site == (1,0)) or (site == (shape[0]-1,0)): 
        energy += field
      elif (site == (shape[0]-1,shape[1]-1)) or (site == (1,1)): 
        energy += field

    # this seems to be enough, I can put in a more formal tolerance eventually
    
    for n in range(-10,11):
      for m in range(-10,11):
        if(n == 0) and (m == 0): continue
        energy += v(dist((site[0]+n*l,site[1]+m*l)))
    
    energies[site] = energy
  return np.ravel(energies)

# calculate new indices after a reflection
# ARGS:
#   l <- linear dimension of lattice
#   direction <- either 1, 2, or 3, each corresponding to one of the three
#                possible reflection axes
#   shift <- int related to reflection center
# RETURNS:
#   vectorized 2d np array where the entry at each index corresponds to the 
#   index it gets mapped to by the reflection
def transformation(l, direction, shift):
  shifts = np.arange(l,dtype='int')
  shifts = np.roll(shifts,shift)
  inds = np.arange(l*l,dtype='int').reshape((l,l))
  if direction == 1:
    inds = inds[:,::-1]
    for i in range(inds.shape[0]):
      inds[i,:] = np.roll(inds[i,:],shifts[i])
  elif direction == 2:
    inds = inds[::-1,:]
    for i in range(inds.shape[1]):
      inds [:,i] = np.roll(inds[:,i],shifts[i])
  elif direction == 3:
    inds = inds[::-1,::-1]
    inds = np.transpose(inds)
    inds = np.roll(inds,shift,axis=1)
    inds = np.roll(inds,shift,axis=0)
  else:
    raise NotImplementedError()
  return np.ravel(inds)

# builds a lookup table corresponding to index maps for possible reflections
# ARGS:
#   l <- linear dimension of lattice simulation cell
# RETURNS:
#   2D np array where each row corresponds to a reflection and each column
#   corresponds to the index that that column is mapped to under the reflection
def make_transformation_table(l):
  directions = [1,2,3] # possible reflection directions
  shifts = [i for i in range(l)] # possible shift amounts
  transformations = np.zeros((len(directions)*len(shifts), l*l), dtype='int')
  row_counter = 0
  for direction in directions:
    for shift in shifts:
      transformations[row_counter,:] = transformation(l, direction, shift)
      row_counter += 1
  return transformations

# builds a lookup table corresponding to interaction energy between two sites
# ARGS:
#   l <- linear dimension of lattice simulation cell
#   energy_table <- optional starting energy table, default None
#   d <- distance between dielectric gates in units of Moire lattice constant
#   field <- strain field strength, doesn't work right now
# RETURNS:
#   2D np array where row and column indices correspond to lattice sites
#   and the entries correspond to the interaction energies between them
def make_pair_energy_table(l, energy_table=None, d=10, field=None):
  if energy_table is None: 
    energy_table = make_energy_table((l,l), d=d, field=field)

  def pair_energy(s1,s2):
    sr1, sc1 = s1 // l, s1 % l # site row and column
    sr2, sc2 = s2 // l, s2 % l
    relative_site = ((sc1 - sc2)%l + l)%l + (((sr1 - sr2)%l + l)%l)*l
    energy = energy_table[relative_site]
    return energy

  pair_energy_table = np.zeros((l*l,l*l))
  for site1 in range(l*l):
    for site2 in range(l*l):
      pair_energy_table[site1,site2] = pair_energy(site1,site2)
  return pair_energy_table

