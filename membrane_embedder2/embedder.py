import click
import sys,os

import numpy as np
import MDAnalysis as mda

# The idea is:
# 1. make it general and fast
# 2. Avoid mdanalysis untill it's neccessary i.e. at the end
# The problem with previous version was that it was creating
# and destroying universes whenever there was a protein added.
# We should be able to store coordinates as points in the 
# Universe class untill there is a need to write a pdb file 
# and make indices etc.

# Read in the structure for every protein. Syntax can stay as it was


# Each protein could also be a class, with all bonds, connections
# And additional beads assigned to this class. The class would first 
# create the full linkage, extra stuff before attempting to 
# place it in the box. Then roll the coordinates using just numpy array 
# and check with the big coordinate array for clashes. Reroll etc.
# When happy, append coordinats to the big array and append instance of 
# protein class to a list. The actual offsets for the bead 
# numbers could be added already when creating the class, 
# it will be known at the time. 

# Dump coordinates and atom names with correct indices to a pdb

# go through the list of Protein instances and assemble 
# a topology.

# Limitation! in this way we cannot have bonds between proteins... 
# We should therefore build the universe incrementally. But 
# bonded proteins should anyhow behave as a single domain. 
# Could it append to an existing domain when inter-protein 
# bond is detected? Or build Universe incrementally and then
# assembly into domains based on connectivity? But we need to make 
# the bonds such, that single proteins are placed in the correct 
# distance. 
# 
import src.parser as PRS
import src.assembler as ASM
import src.topology_builder as TPB

class Universe:
   def __init__(self):
      pass

@click.group()
@click.option('--debug/--no-debug', default=False)
@click.pass_context
def cli(ctx, debug):
   """Generate Complexes++ input"""
   ctx.obj['DEBUG'] = debug

@cli.command()
@click.option('--input','input',default=None,help="Input file for molecules")
@click.option('--box','box',nargs=3, default=(100,100,100), type=click.FLOAT, help='unit cell vectors a,b,c\n default is 100.0, 100.0,100.0')
@click.option('--pdbout','pdbout',nargs=1, default=None, type=click.STRING, help='name of an output pdb file',required=True)
@click.option('--topout','topout',nargs=1, default=None, type=click.STRING, help='name of an output topology file',required=True)
@click.pass_context
def assemble(ctx,input,box,pdbout,topout):
   if ctx.obj['DEBUG']:
      print('Debug is %s' % (ctx.obj['DEBUG'] and 'on' or 'off'))
   box=np.array(box,dtype=np.float32)
   P=PRS.Parse(input,ctx.obj['DEBUG'])
   A=ASM.Assembler(box)
   for mol in P.keys():
      #~ print "=",P[mol]
      A.add_molecule(P[mol])
   A._writePDB(pdbout)
   A._writeTOP(topout)
   
@cli.command()
@click.option('--topin','topin',nargs=1, default=None, type=click.STRING, help='name of an input topology file',required=True)
@click.pass_context
def build_topology(ctx,topin):
   print('Debug is %s' % (ctx.obj['DEBUG'] and 'on' or 'off'))
   B=TPB.Builder()
   B.add_topology(topin)

if __name__ == '__main__':
   cli(obj={})
    