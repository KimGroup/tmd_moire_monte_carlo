import numpy as np
from lookup_tables import make_energy_table
import structure_factor_paths as sf_pth

# Base class for operators
class Operator(object):
  # generic init method called upon instantiation of class Object
  # ARGS:
  #   l <- linear dimension of lattice simulation cell
  #   initial_value <- optional initial value for operator, default None
  #   sum_evals <- optional, whether to sum over Markov chain, default True
  # RETURNS:
  #   None
  def __init__(self,l,initial_value = None,sum_evals = True):
    self.sum_evals = sum_evals
    if not sum_evals: self.state = np.array([])
    else: self.state = 0
    # number of times evaluate has been called
    self.n_evals = 0
    return
  
  # delta argument corresponds to a known change since last evaluation
  # generic evaluate method
  # ARGS:
  #   lat <- Numpy array of 0/1, recording un/occupied lattice sites
  #   delta <- optional change in Operator since previous evaluate call
  # RETURNS:
  #   value of Operator evaluatd on lat
  def evaluate(self,lat,delta = None):
    self.n_evals += 1
    if self.sum_evals: self.state += 0
    else: np.append(self.state,0)
    return 0

# Operator class to record lattice configurations
class LatticeConfig(Operator):
  def __init__(self, l, initial_value=None, sum_evals=False):
    self.l = l
    self.sum_evals = sum_evals
    self.n_evals = 0
    if not (initial_value is None):
      if len(initial_value.shape) == 1:
        initial_value = initial_value[np.newaxis,:]
      self.state = initial_value
    else:
      if self.sum_evals: self.state = np.zeros(self.l*self.l, dtype='int')
      else: self.state = np.zeros((0,self.l*self.l), dtype='int')

  def evaluate(self, lat, delta=None):
    if self.sum_evals: 
      self.state += lat
    else: 
      self.state=np.append(self.state,lat[np.newaxis,:],axis=0)
    self.n_evals += 1
    return

# Operator class to evaluate the Hamiltonian from the MS on the lattice
class Energy(Operator):
  # Allowed kwargs:
  #   d <- dielectric gate spacing in units of Moire lattice constants
  #   field <- strain field, currently not working
  def __init__(self,l, initial_value=None, sum_evals=True, **kwargs):
    self.l = l
    d = 10 # gate spacing
    field = None # strain field, not working
    for key,value in kwargs.items():
      if key == "d": d = value
      if key == "field": field = value
    self.energy_table = make_energy_table((l,l),d=d,field=field)
    super(Energy,self).__init__(l,initial_value,sum_evals)
    return

  def evaluate(self,lat,delta = None):
    if not (delta is None): energy = self.state[-1]+delta
    else:
      l = self.l
      energy = 0
      site2 = np.arange(l*l)
      site2 = site2[lat[site2] != 0]
      for i in range(l*l):
        # calculate relative site position and lookup energy
        sr,sc = i//l,i%l
        relative_site = ((sc - site2 % l)%l + l)%l + (((sr-site2//l)%l+l)%l)*l
        energy += np.sum(self.energy_table[relative_site])*lat[i]
      # each interaction was double counted 
      energy /= 2
    self.state = np.append(self.state,energy)
    self.n_evals += 1
    return energy

# Operator class to evaluate the  nematic order parameter on the lattice
# See appendix B of MS for details
class NematicOP(Operator):
  # Allowed kwargs:
  #   spatial_average <- boolean whether to average OP over lattice
  def __init__(self, l, initial_value=None, sum_evals=True, **kwargs):
    self.l = l
    self.sum_evals = sum_evals
    self.spatial_average = True
    for key,value in kwargs.items():
      if key == "spatial_average":
        assert(isinstance(value,bool))
        self.spatial_average = value
    # Set up nn stuff for OP calculation. Note dx/dy/nn have consistent order
    delta_xs = np.zeros((self.l*self.l,6))
    delta_ys = np.zeros((self.l*self.l,6))
    delta_xs[:,:] = np.array([-0.5,-1,-0.5,0.5,1,0.5])
    delta_ys[:,:] = np.array([np.sqrt(3)/2,0,-np.sqrt(3)/2,-np.sqrt(3)/2,0,np.sqrt(3)/2])
    # Precompute possible OP component values
    self.n1_terms = delta_xs*delta_xs-delta_ys*delta_ys
    self.n2_terms = 2*delta_xs*delta_ys
    
    nns = np.vstack([[1,0],[0,-1],[-1,-1],[-1,0],[0,1],[1,1]])
    self.nn_table = np.zeros((self.l*self.l,len(nns)),dtype='int')
    mod_index = lambda index: np.array([(i%s + s)%s for i,s in zip(index,(self.l,self.l))])
    # precompute nearest neighbor indices for each site
    for site in np.ndindex((self.l,self.l)):
      site_neighbors = [mod_index(site+nn) for nn in nns]
      site_neighbors = [idx[1]+idx[0]*self.l for idx in site_neighbors]
      self.nn_table[site[1]+site[0]*self.l, :] = np.array(site_neighbors)

    # finish setting up
    self.n_evals = 0
    if not (initial_value is None):
      self.state = initial_value
    else:
      if self.spatial_average:
        if self.sum_evals: self.state = np.zeros(2)
        else: self.state = np.zeros((2,0))
      else:
        if self.sum_evals: self.state = np.zeros((2,self.l*self.l))
        else: self.state = np.zeros((2,self.l*self.l,0))
    return

  def evaluate(self,lat,delta = None):
    if not(delta is None):
      new_nematic_op = self.state[-1]+delta
      self.state = np.vstack([self.state,new_nematic_op])
    else:
      # occupied nearest neighbors
      occ = (lat[self.nn_table] == 1)
      if self.spatial_average:
        # calculate the OP contributions from occupied sites/nns
        n1 = np.sum(np.sum(self.n1_terms*occ,axis=1)*lat)
        n2 = np.sum(np.sum(self.n2_terms*occ,axis=1)*lat)
        if self.sum_evals: 
          self.state += np.array([n1,n2])
        else: 
          self.state = np.append(self.state,np.array([n1,n2])[:,np.newaxis],axis=1)
      else:
        n1 = np.sum(self.n1_terms*occ,axis=1)*lat
        n2 = np.sum(self.n2_terms*occ,axis=1)*lat
        if self.sum_evals:
          self.state += np.vstack([n1,n2])
        else:
          self.state = np.append(self.state,np.vstack([n1,n2])[:,:,np.newaxis],axis=2)
    self.n_evals += 1
    return

# Operator class to calculate the structure factor (see MS for definition)
class StructureFactor(Operator):
  # Allowed kwargs:
  #   mode <- "BZ" or "PATH" indicating to calculate structure factor on a
  #           mesh within the first Brillouin zone or along the paths in
  #           k-space defined in structure_factor_paths.py
  def __init__(self, l, initial_value = None, sum_evals = True, **kwargs):
    # lattice vectors
    a1, a2 = np.array([0,1]),np.array([np.sqrt(3)/2.0,-0.5])
    # reciprocal lattice vectors
    q1, q2 = np.array([2*np.pi/np.sqrt(3),2*np.pi]),np.array([4*np.pi/np.sqrt(3),0])
    self.mode = "PATH" #"BZ"
    for key,value in kwargs.items():
      if key == "mode":
        self.mode = value
    if self.mode == "BZ":
      # use 100 by 100 k-point mesh, kind of arbitrary
      qxs = np.linspace(-2*np.pi,2*np.pi,101)
      qys = np.linspace(-4*np.pi/np.sqrt(3),4*np.pi/np.sqrt(3),101)
      Qx,Qy = np.meshgrid(qxs,qys)
      Qx,Qy = Qx.ravel(),Qy.ravel()
    elif self.mode == "PATH":
      ts = np.linspace(0,1,200)
      qs = np.vstack([pth.coord(t) for pth in sf_pth.PATHS for t in ts])
      Qx,Qy = qs[:,0],qs[:,1]
    # precompute contributions to sum for occupied sites
    self.sf_terms = np.zeros((len(Qx),l*l),dtype='complex')
    for ind1 in range(len(Qx)):
      vec_q = np.array([Qy[ind1],Qx[ind1]])
      for ind2,(i,j) in enumerate(np.ndindex((l,l))):
        r = i*a2+j*a1
        self.sf_terms[ind1,ind2] = np.exp(-1j*np.dot(vec_q,r))/l
    # finish setting up    
    self.sum_evals = sum_evals
    if not(initial_value is None): self.state = initial_value
    else:
      if self.sum_evals: self.state = np.zeros(len(Qx))
      else: self.state = np.zeros((len(Qx),0))
    self.n_evals = 0
    return

  def evaluate(self, lat, delta = None):
    # compute which terms contribute to sum and evaluate
    sf = np.square(np.abs(np.matmul(self.sf_terms,lat)))
    if self.sum_evals:
      self.state += sf
    else:
      self.state = np.append(self.state,sf[:,np.newaxis],axis=1)
    self.n_evals += 1
    return
  
  # reset operator without redoing computations in init method
  def reset(self):
    if self.sum_evals: self.state = np.zeros(self.state.shape)
    else: self.state = np.zeros((self.state.shape[0],0))
    self.n_evals = 0
    return
    
