from openmm.app import Simulation, PDBxFile,PDBFile, StateDataReporter, \
    GromacsTopFile, AmberPrmtopFile, CharmmPsfFile
from openmm import XmlSerializer, Platform
from mdtraj.reporters import XTCReporter
from pathlib import Path
import json
from sys import argv

top_loaders = {
    '.top': GromacsTopFile,
    '.pdb': PDBFile,
    '.cif': PDBxFile,
    '.prmtop': AmberPrmtopFile,
    '.psf': CharmmPsfFile
}

with open(argv[1]) as f:
    conf = json.load(f)

# use extension produced by 'suffix' to load topology
try:
    top_path = Path(conf['topology'])
    topology = top_loaders[top_path.suffix](str(top_path)).topology
except KeyError:
    print(top_path.suffix, 'is not a recognized topology extension. Choices are:',
          *top_loaders.keys())
    raise

# load the system components
integrator = XmlSerializer.deserialize(Path(conf['integrator']).read_text())
system = XmlSerializer.deserialize(Path(conf['system']).read_text())
simulation = Simulation(topology, system, integrator, Platform.getPlatformByName(conf['platform name']),
                        conf['platform properties'])
# Get positions from config

# Get new velocities at specified temp
simulation.context.setVelocitiesToTemperature(conf['temperature'])
# set up reporters
steps = conf['steps']
write_freq = conf['write_freq']
out_prefix = conf['output_prefix']
simulation.reporters.append(
    XTCReporter(
        out_prefix+'.xtc',
        write_freq,
        append=conf['append']
    )
)
simulation.reporters.append(
    StateDataReporter(
        out_prefix+'.out',
        write_freq,
        **(conf['out reporters']),
    )
)


simulation.step(steps)

