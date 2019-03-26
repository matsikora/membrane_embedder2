import click
import sys,os
import pickle

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mdadist

# We should be able to have a way of creating a 
# connectivity graph, so that connected proteins are 
# put into the same domain for topology movements.
# we could have it by adding a "linked_to" list as a topology
# attribute. This would be sorted when assembling
# It could also be a way of getting any extra beads - 
# one could define them at the end of input file and assign connections.

class Builder:
   def __init__(self):
      self.topologies = []
      
   def add_topology(self,fname):
      """
      append new top to a list
      """
      self.topologies.append(self._readTOP(fname))
      
      
   def _readTOP(self,fname):
      """
      Read a pickle obj conaining a prepared tpoology
      """
      pickle_in = open(fname,'rb')
      topology = pickle.load(pickle_in)
      pickle_in.close
      return topology
   
   