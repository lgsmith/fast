from open_mm.app import *
from openmm import *
from openmm.unit import *
import json
from sys import argv
import mdtraj as md

with open(argv[1]) as f:
    conf = json.load(f)

sim = simulation.loadState(conf['sim'])
pos = md.load(conf['crds'])
sim.context.setPositions(pos.coords)
sim.minimizeEnergy()
sim.step(conf['n_steps'])
