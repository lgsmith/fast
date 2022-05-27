from openmm.app import Simulation, PDBFile, DCDReporter, StateDataReporter
from openmm import XmlSerializer, Platform
from pathlib import Path
import json
from sys import argv

with open(argv[1]) as f:
    conf = json.load(f)

# load the system components
topology = PDBFile(conf['topology']).getTopology()
integrator = XmlSerializer.deserialize(Path(conf['integrator']).read_text())
system = XmlSerializer.deserialize(Path(conf['system']).read_text())
simulation = Simulation(topology, system, integrator, Platform.getPlatformByName(conf['platform name']),
                        conf['platform properties'])
simulation.loadState(conf['state'])
# set up reporters
steps = conf['steps']
write_freq = conf['write_freq']
out_prefix = conf['output_prefix']
simulation.reporters.append(
    DCDReporter(
        out_prefix+'.dcd',
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

# Get new velocities at specified temp
simulation.context.setVelocitiesToTemperature(conf['temperature'])

simulation.step(steps)

