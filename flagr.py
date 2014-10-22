#! /usr/bin/env python
# -*- coding: utf_8 -*-
###__________________________________________________________________________###
#####                              ::flagr::                              ######
#######                     a flake growth simulator                    ########
##########                                                           ###########
#######    Guilherme Stein (2014)                                       ########
#####     <guilherme.stein@physik.uni-wuerzburg.de>                       ######
###--------------------------------------------------------------------------###

from __future__ import print_function, division, generators
import time
import itertools as it
#import  math
#import  sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import flt
#import iters


#==============================================================================#
#                               class FlakeBase                                #
#==============================================================================#
class FlakeBase(object):

  """
  Base class defining common functions.
  """

  def attribute_setter(self, **kwargs):
    """
    Sets kwargs as instance attributes.
    """
    for key, value in kwargs.iteritems():
      setattr(self, key, value)

  def return_dict(self, return_list):
    """
    Extracts instance attributes into dictionary from return_list.
    """
    return {k: v for k, v in self.__dict__.iteritems() if k in return_list}

#==============================================================================#
#                                  class Grid                                  #
#==============================================================================#
class Grid(FlakeBase):

  """Here I write my docstring"""

  def __init__(self, **kwargs):
    """TODO: to be defined1.

    :**kwargs: TODO

    """
    attribute_setter(**kwargs)

#==============================================================================#
#                                 class Flake                                  #
#==============================================================================#


class Flake(FlakeBase):
  def __init__(self, **kwargs):
    """
    The "Flake" creates a lattice, more precise a mapping of a cubic- to a
    fcc-lattice. The lattice is then populated with a seed flake, containing a
    twin plane. As arguments a dictionary is expected, with following entries:

      --key--      --type--          --description--
    iteration    : int      | number of atoms initally deposited
    size         : int      | edge length of lattice cube
    twin_size    : int      | edge length of seed
    seed         : int      | height of twin planes
    e_base       : int      | base for binding energy =: BASE**BINDINGS
    nn_distance  : int      | radius including atoms as next neighbours
    """
    self.attribute_setter(**kwargs)

    # lattice specifications:
    self.a = np.zeros(4 * self.size**3).reshape(self.size,  # x'
                                                self.size,  # y'
                                                self.size,  # z'
                                                5)          # (x,y,z,t,u)
    self.seed_size = np.array(((self.size - self.seed) // 2,
                               (self.size + self.seed) // 2))
    self.nn = flt.next_neighbour('all')

    # attribute initialization:
    self.atom       = -1
    self.site       = (0, 0, 0)
    self.runtime    = 0
    self.inittime   = 0
    self.idx        = range(self.size)
    self.prob_sum   = 0
    self.prob_list  = []
    self.bindings   = []
    self.candidates = []

    # ok, what is this exactly ???
    for nn in range(len(self.nn)+1):
      self.candidates.append([])
      self.bindings.append(0)
      self.prob_list.append(0)

    # offset for plane borders, WIP!
    self.shift = self.seed_size % 2

  def lattice(self):
    """
    Creates a FCC-lattice based on flt.vector() and the size of the
    twin-planes.
    """
    # twin planes boundaries...
    twin_min = (self.size - self.twin_size) // 2
    twin_max = (self.size + self.twin_size) // 2
    # ... and index ranges
    twin_l_idx = range(twin_min)
    twin_m_idx = range(twin_min, twin_max)
    twin_u_idx = range(twin_max, self.size)

    for (i, j, k) in it.product(self.idx, self.idx, twin_l_idx):
      self.a[i, j, k] = flt.vector(i, j, k, regular=True)
    for (i, j, k) in it.product(self.idx, self.idx, twin_m_idx):
      self.a[i, j, k] = flt.vector(i, j, k, regular=False, shift=self.shift[0])
    for (i, j, k) in it.product(self.idx, self.idx, twin_u_idx):
      self.a[i, j, k] = flt.vector(i, j, k, regular=True, shift=self.shift[1])

  def make_seed(self):
    """
    Places a seed in the middle of the lattice.
    """
    seed_min = self.seed_size[0]
    seed_max = self.seed_size[1]
    idx = slice(seed_min, seed_max)
    self.seed_idx = list(range(seed_min, seed_max))
    self.a[idx, idx, idx, 3] = self.atom

  def candidate_check(self, site=None, check='single'):
    """
    Performs an iteration over the check-object and evaluates if this site
    is empty and has neighbouring atoms, it adds 2 for each atom in the
    neighbourhood. check-object can be 'single', 'seed' and 'all'.
  """
    if check == 'single':
      if site is None:
        site = self.site
      self.nn_check(site)
    else:
      if check == 'all':
        idx = self.idx
      elif check == 'seed':
        idx = self.seed_idx
      for (i, j, k) in it.product(idx, idx, idx):
        site = (i, j, k)     # index of site, contains (x,y,z,t)
        if self.a[site][3] == self.atom:
          self.nn_check(site)

  def neighbour(self, site):
    # zip together coordinate and next neighbour tuples
    pairs = [zip(site, nn) for nn in self.nn]
    # return list of tuples, containing added coordinates
    return [tuple(sum(y) for y in x) for x in pairs]

  def site_status(self, site_list=None, idx=3):
    if site_list == None:
      site_list = self.site
    return [self.a[site][idx] for site in site_list]


##  make get, set status, make class,  for lattice
  def nn_check(self, site=None):
    if site == None:
      site = self.site
    bind_idx = 0
    indices = []
    try:
      nb = self.neighbour(site)
    except IndexError:
      print("Lattice border reached")
    for space in nb:
      if self.site_status(space) == self.atom:
        bind_idx += 1
      #else: and self.site_status(4) == 1
        #self.candidates[bind_idx].append(space)
        #self.a[space][3] = bind_idx
    self.bindings[bind_idx] += 1
    if site_status(site) is not self.atom:
      site_status(site) = bind_idx
    indices.append((site, bind_idx))
    return indices

#         indices.append((nb, bind_idx))
#         self.candidates[bind_idx].append(nb)   # was list(nb)
#         self.bindings[bind_idx] += 1
#         self.a[nb][3] = bind_idx


#   def nn_check(self, site_idx, checkonly=False):
#     """
#     First loop iterates over next neighbours and collects them in nbhood if
#     they are empty. Second loop iterates over nbhood and adds a binding for
#     each of THOSE neighbours that are atoms.
#     """
#     nbhood = []        # neighbourhood
#     indices = []        # nb vectors with binding numbers
#     for step in self.nn:
#       neighbour_idx = tuple(np.array(site_idx) + np.array(step))
#       try:
#         neighbour = self.a[neighbour_idx][3]   # does it even exist?
#         dist = self.distance(site_idx, neighbour_idx)
#         if neighbour != self.atom and dist < self.nn_distance:
#           nbhood.append(neighbour_idx)
#       except Exception:
#         print("Lattice Border reached!")
#
#     for nb in nbhood:
#       # increase bindings counter, each entry counts the actual number of
#       # bindings represented by the list index
#       bind_idx = 0
#       for step in self.nn:
#         # NEW neighbour (of nbhood-item)
#         n_neighbour_idx = tuple(np.array(nb) + np.array(step))
#         # n_neighbour_idx = tuple(nb + step)    #neighbour's neighbour
#         try:
#           n_neighbour = self.a[n_neighbour_idx][3]
#           if n_neighbour == self.atom:
#             # now it counts sourrounding atoms
#             bind_idx += 1
#         except IndexError:
#           print("Neighbour is outside lattice border!")
#       # append lattice object to candidates list, according to binding number
#       # nb should always be free space here! index type
#       if checkonly:
#         indices.append((nb, bind_idx))
#         print("checkonly")
#       else:
#         indices.append((nb, bind_idx))
#         self.candidates[bind_idx].append(nb)   # was list(nb)
#         self.bindings[bind_idx] += 1
#         self.a[nb][3] = bind_idx
#     return indices

  def distance(self, atom1, atom2):
    """
    Returns the distance between two atoms.
    """
    (x1, y1, z1) = self.a[atom1][:3]
    (x2, y2, z2) = self.a[atom2][:3]
    dist_sq = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    return np.sqrt(dist_sq)

  def probability(self, index=None):
    # verify that this is correct !!
    def exp(x):
      return self.e_base**x
    if index is None:
      for idx, elem in enumerate(self.bindings):
        self.prob_list[idx] = elem * exp(idx)
        self.prob_sum = sum(self.prob_list)
    else:
      for elem in index:
        bindings = elem[1]
        new_p = exp(bindings)
        old_p = exp(bindings-1)
        self.prob_list[bindings] += new_p
        self.prob_list[bindings-1] -= old_p
        self.prob_sum += (new_p - old_p)

  def strike(self):
    """
    Create weighted probability distribution based on prob_list. After that draw
    a list and subsequently draw an element inside the list.
    """
    l_luck = np.random.randint(self.prob_sum)
    cand_list = []
    while not cand_list:
      adder = 0
      l_luck = np.random.randint(self.prob_sum)
      for idx, prob in enumerate(self.prob_list):
        adder += prob
        if adder >= l_luck:
          lucky_list = idx
          break
      cand_list = self.candidates[lucky_list]
    a_luck = np.random.randint(len(cand_list))
    try:
      lucky_atom = self.candidates[lucky_list].pop(a_luck)
    except Exception:
      print("candidates list empty, something wrong...")
    return lucky_atom

  def grow(self, rounds=1):
    starttime = time.time()
    for round in range(rounds):
      atom = self.strike()
      self.site = atom
      self.a[atom][3] = -1
      indices = self.candidate_check(site=atom)
      self.probability(indices)
    endtime = time.time()
    self.inittime = "%.3f" % (endtime - starttime)
    print("Grew ", rounds, " atoms in ", self.inittime, "seconds.")
    print("-------------------------------")

  def plot(self, save=False):
    # only add sites with status 1 or 2 to scatter plot:
    x, _x = [], self.a[:, ..., 0].flatten()
    y, _y = [], self.a[:, ..., 1].flatten()
    z, _z = [], self.a[:, ..., 2].flatten()
    t, _t = [], self.a[:, ..., 3].flatten()

    for i, elem in enumerate(_t):
      if elem:
        x.append(_x[i])
        y.append(_y[i])
        z.append(_z[i])
        t.append(_t[i])

    # options:
    font = {
        'family': 'serif',
        'color': 'darkred',
        'weight': 'normal',
        'size': 18}
    plt.xlabel('x', fontdict=font)
    plt.ylabel('y', fontdict=font)
    crange = (self.atom, len(self.nn) // 2)
    spectrum = 'spectral'                       # 'gist_stern'
    cbar = plt.cm.ScalarMappable(
        cmap=spectrum,
        norm=plt.normalize(
            vmin=crange[0],
            vmax=crange[1]))
    cbar._A = []

    # plotting:
    plt.clf()
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c=t, cmap=spectrum, vmin=crange[0], vmax=crange[1],
               s=self.sphere_size)
    if self.fixed_axes:
      ax.set_xlim(0, 2 * self.size)
      ax.set_ylim(0, 2 * self.size)
      ax.set_zlim(-self.size / 2 + 2 * self.seed, self.size / 2 + 2 * self.seed)
    plt.colorbar(cbar)

    # output
    if save:
      plt.savefig(
          'output/goldflake_base-{}_{}_{}_{}'.format(self.e_base, RSEED,
                                                     counter, thing) + '.png')
    else:
      plt.show()

  def info(self, file_out=False, grid_info=False):
    """
    Print out status information about the Flake.
    """
    if file_out:
      with open('output/timing.txt', mode='a') as logfile,\
      flt.RedirectStdoutTo(logfile):
        print(self.size, '\t\t', self.iteration, '\t\t', 'init:', self.inittime,
              'run:', self.runtime, '\t', self.comment, '\n')
    elif not file_out:
      print("-------------------------------")
      print("Total time: ", self.runtime, "seconds")
    if grid_info:
      print("Nearest neighbours:\n", self.nn)
      print()
      print("4th grid entry defines the status of site: \n",
            "0:\t empty site\n",
            "-1:\t gold atom\n",
            ">0:\t empty site next to atom => log(P)\n")

  def initialize(self):
    """
    Initialisation of array, seed and probabilities at beginning
    """
    print("---- Basis: ", self.e_base, "--------- Seed: ", RSEED)
    starttime = time.time()
    # self.array(self.size)
    self.lattice()
    self.make_seed()
    self.candidate_check(check='seed')
    self.probability()
    endtime = time.time()
    self.inittime = "%.3f" % (endtime - starttime)
    print("-------------------------------")
    print("Lattice Dimensions:", self.size, "x", self.size, "x", self.size)
    print("Initialized in ", self.inittime, "seconds.")
    print("-------------------------------")

  def main(self):
    """ main logic loop """
    starttime = time.time()
    self.initialize()
    self.grow(rounds=self.iteration)
    endtime = time.time()
    self.runtime = "%.3f" % (endtime - starttime)
    self.info()
    self.info(file_out=True)

# ::::: : :::::

flake_params  = {
              # lattice parameters
                'iteration': 100,
                'size': 20,
                'twin_size': 2,
                'seed': 5,
                'e_base': 10,
                'nn_distance': 2,
              # output parameters
                'sphere_size': 50,
                'fixed_axes': False,
                'comment': ''
                }

RSEED = np.random.randint(99999)
np.random.seed(RSEED)

# mass-production
if False:
  counter = 0
  GROW_ITER = 0
  gold = Flake()
  gold.initialize()
  MAX_GROW = gold.size**3//5000      # growth in thousands
  st_lst = [0, 100, 300, 600, 1000, 3000, 5000]
  st_lst.extend((MAX_GROW // 10 - 1) * [10**4])
  for thing in st_lst:
    GROW_ITER = thing
    gold.grow(rounds=thing)
    gold.plot(save=True)
    counter += 1
else:
  if __name__ == "__main__":
    gold = Flake(**flake_params)
    gold.main()
    #gold.plot()

go = Flake(**flake_params)
