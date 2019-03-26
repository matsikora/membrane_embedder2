import click
import sys,os
import pickle

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mdadist

class Topology:
   def __init__(self):
      """
      Keep the connectivity info for any single molecule type.
      This will correspond to a single "domain" in complexes
      """
      # At the least, add residue span per domain
      self.termini=[] 
      #~ self.box=box
      #~ self.pdb=pdb
      #~ self.topout='test.top'
      
   def add_domain_span(self,span):
      """
      Add residue range belonging to a single domain
      """
      self.termini.append(span)
   
   
   #~ def write_topology(self):
      #~ with open(self.topout,'w') as f:
         #~ f.write("tralala\n")
      
class Assembler:
   def __init__(self,box):
      """
      A general class that will handle molecule placement in the
      simulation box
      """
      self.box=box
      # Initialise empty array
      self.LEN = 10000 # To speed things up we do it statically and append once we exceed N
      self.global_coordinates = np.empty([self.LEN,3],dtype=np.float32)
      self.N = 0 # Keep current number of atoms
      self.segids=[] # Keep chain names
      self.residues=[] # Keep residue names
      self.topologies=[] #Keep topologies
      
      
   def add_molecule(self,options):
      """
      Process a single molecule
      """
      
      # Get relevant info
      molname=options["MolName"]
      PDB=options['PDB']
      Ncopies=options['Ncopies']
      Placement=options['Placement']
      Cutoff=options['Cutoff']
      coordinates,residues,segids=self._loadPDB(PDB)
      
      
      # There are 3 options: random, grid, coordinates 
      if Placement == "coordinates":
         # Directly take the coordinates as they are in the input file
         # Assume Ncopies=1
         if np.any(coordinates.max()>self.box) or np.any(coordinates.min()<0): # Check it, I think complexes can wrap a mol even if it is outside box
            raise BaseException("The molecule does not fit in the box!")
         
         # Now check if it does not collide with the rest of the box
         if self._has_clash(coordinates,Cutoff):
            raise BaseException("A collision detected in the system {} for PDB file {}".format(molname,PDB))
         self._append_mol(coordinates,residues,segids)
         
         
      elif Placement == "random":
         # place molecules randomly in the box
         
         # Move test molecule to the begining of the coord system
         coordinates-=np.mean(coordinates,axis=0) 

         # Iterate until a clash-free configuration is found
         imol=0
         while imol < Ncopies:
            pos = self._get_random_position()
            
            if not self._has_clash(coordinates+pos,Cutoff):
               self._append_mol(coordinates+pos,residues,segids)
               imol+=1
               print "accepted {} out of {}".format(imol,Ncopies)
            else:
               print "reroll"
               
               
               
      elif Placement == "grid":
         # Create a cubic grid in the box unless this is a membrane system
         # Move test molecule to the begining of the coord system
         coordinates-=np.mean(coordinates,axis=0) 
         positions = self._guess_grid(Ncopies,2)


         for imol in range(Ncopies):
            pos=positions[imol]
            if not self._has_clash(coordinates+pos,Cutoff):
               self._update_top(segids)
               self._append_mol(coordinates+pos,residues,segids)
               imol+=1
               print "accepted {} out of {}".format(imol,Ncopies)
            else:
               raise BaseException("Clash found in during grid formation in the system {} for PDB file {}  for molecule No {}".format(molname,PDB,imol))
      
   def _update_top(self,segids):
      """
      Update all topology features according to actual numbering
      """
      topology=Topology()
      topology.add_domain_span([self.N+1,self.N+len(segids)])
      self.topologies.append(topology)
      
   def _append_mol(self,coordinates,residues,segids):
      """
      Append mol at the end
      """
      self._paste_coordinates(coordinates)
      self.residues=self.residues+residues
      self.segids=self.segids+segids
      
   def _guess_grid(self,Ncopies,dimensions=3):
      """
      Guess spacing and distribute molecules on a grid
      Will guess the minimum number of planes and rows and fill with Ncopies
      
      This is a not well defined problem, so we make an approximation.
      
      If dimensions < 3 we will still return a 3D array but with zeros in the last 1 or 2 coords
      """
      
      points_per_dim=int(np.ceil(Ncopies**(1./dimensions)))
      
      # fill first dimension, then second, then third
      if dimensions==3:
         new_positions=list(np.ndindex(points_per_dim,points_per_dim,points_per_dim))
      elif dimensions==2:
         new_positions=list(np.ndindex(points_per_dim,points_per_dim))
      else:
         new_positions=list(np.ndindex(points_per_dim))
         
      # Get box fractions
      new_positions = np.array(new_positions[:Ncopies]).astype(np.float32)/(dimensions+1) # divide by N+1 to avoid touching via PBC
      
      # supplement zeros to have 3d array always
      new_positions=np.hstack([new_positions]+[np.zeros((Ncopies,1)) for i in xrange(3-dimensions)])

      return new_positions.astype(np.float32)*self.box

      
   def _get_random_position(self):
      """
      get a random position in the box
      """
      
      return (np.random.random(3)*self.box).astype(np.float32)
      
      
      
   def _has_clash(self,coordinates,cutoff):
      """
      Check for clashes in the system
      """
      if cutoff > 0:
         d=mdadist.distance_array(self.global_coordinates[:self.N,:],coordinates,box=self.box)
         if np.any(d<cutoff):
            return True
         else:
            return False
      else:
         return False

         
   def _paste_coordinates(self,coordinates):
      """
      Append coordinates to the global universe
      """
      # Determine how many atoms we paste
      M=coordinates.shape[0] 
      
      # Check if the host array is not too small
      if self.N+M>=self.global_coordinates.shape[0]:
         self.global_coordinates = np.append(self.global_coordinates,np.empty([self.LEN,3],dtype=np.float32),axis=0)
      
      if self.N == 0:
         # initialise
         self.global_coordinates[:M,:]=coordinates
         
      else:
         self.global_coordinates[self.N:(self.N+M),:] = coordinates
      self.N = self.N + M
      

      
   def _loadPDB(self,PDB):
      """
      Is it better to store it as MDA universe or rather as coordinates only?
      """
      u=mda.Universe(PDB)
      # copy relevant data and delete Universe
      coordinates=np.copy(u.atoms.positions)
      residues=list(u.atoms.resnames)
      segids=list(u.atoms.segids)
      del u
      return coordinates,residues,segids
   
   
   def _writePDB(self,pdbout):
      """
      Writes a current state of the assembled system
      
      It could be better to have n_segments corresponding to the real number of segments and then a correct
      index that would assign atoms to segments.
      """
      # create a universe
      atom_resindex=np.arange(0,self.N,1)
      residue_segindex=np.arange(0,self.N,1)
      u = mda.Universe.empty(n_atoms=self.N, trajectory=True, n_residues=self.N, n_segments=self.N, atom_resindex=atom_resindex,residue_segindex=residue_segindex)

      # add topology attributes to store the proteins properly in a pdb file
      u.add_TopologyAttr('resnames')
      u.add_TopologyAttr('segids')
      u.add_TopologyAttr('names')
      
      # fill in attributes
      u.atoms.residues.resnames = self.residues[:self.N]
      u.atoms.segments.segids = self.segids[:self.N]
      u.atoms.names = ["CA"]*self.N
      u.atoms.positions = self.global_coordinates[:self.N]
      # Write the pdb file
      u.atoms.write(pdbout)
      
   def _writeTOP(self,topout):
      """
      Write a list of topologies for further processing
      """
      
      #~ for top in self.topologies:
         #~ print top
      
      pickle_out=open(topout,'wb')
      pickle.dump(list(self.topologies),pickle_out)
      pickle_out.close()