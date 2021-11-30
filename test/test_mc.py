from operators import *
from lookup_tables import make_energy_table
import numpy as np
from scipy.special import binom
import time

def ctz(m,n_bits):
  mask = sum([1 << i for i in range(n_bits)])
  tz_count = 0
  for i in range(n_bits):
    if((m & (mask ^ (1 << i))) != m): return tz_count
    tz_count += 1
  return 0

def state_perms(n_e,n_sites):
  fmt_str = '{0:0'+str(n_sites)+'b}'
  fmt = lambda n: [int(i) for i in fmt_str.format(n)]
  n_states = int(round(binom(n_sites,n_e)))
  s = sum([1 << i for i in range(n_e)])
  res = [fmt(s)]
  for j in range(1,n_states):
    t = s | (s-1)
    s = (t+1)|(((~t & -~t)-1)>>(ctz(s,n_sites)+1))
    res += [fmt(s)]
  return sorted(res)

def exact_average_energy(l,n_particles,temp,energy_vals):
  z = np.sum(np.exp(-energy_vals/temp))
  res = np.sum(energy_vals*np.exp(-energy_vals/temp))/z
  return res


def config_energy(l,lat,energy_table):
  energy = 0
  site2 = np.arange(l*l)
  site2 = site2[lat[site2] != 0]
  for i in range(l*l):
    sr,sc = i//l,i%l
    relative_site = ((sc - site2 % l)%l + l)%l + (((sr-site2//l)%l+l)%l)*l
    energy += np.sum(energy_table[relative_site])*lat[i]
  energy /= 2
  return energy

def find_ground_state(l,n_particles):
  energy_table = make_energy_table((l,l))
  n_sites = l*l
  n_e = n_particles
  fmt_str = '{0:0'+str(n_sites)+'b}'
  fmt = lambda n: np.array([int(i) for i in fmt_str.format(n)])
  n_states = int(round(binom(n_sites,n_e)))
  s = sum([1 << i for i in range(n_e)])

  config = fmt(s)
  lowest_configs = [config]
  t0 = time.time()
  lowest_energy = config_energy(l,config,energy_table)
  tf = time.time()
  print(tf-t0)
  return
  for j in range(1,n_states):
    t = s | (s-1)
    s = (t+1)|(((~t & -~t)-1)>>(ctz(s,n_sites)+1))
    config = fmt(s)
    nrg = config_energy(l,config,energy_table)
    if np.allclose(nrg,lowest_energy):
      lowest_configs += [config]
    elif (nrg < lowest_energy):
      lowest_energy = nrg
      lowest_configs = [config]
  return lowest_energy,lowest_configs


def main():
  l = 4
  n_particles = [6,7,8]
  ts = np.linspace(1,10,10)
  for n in n_particles:
    nrg = Energy(l,sum_evals=False)
    configs = np.array(state_perms(n,l*l))
    for config in configs: nrg.evaluate(config)
    vals = nrg.state
    res = []
    for t in ts:
      res += [exact_average_energy(l,n,t,vals)]
    print(res)

find_ground_state(5,12)
#main()
