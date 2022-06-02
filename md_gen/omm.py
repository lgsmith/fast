import mdtraj as md
import numpy as np
from .. import tools
from pathlib import Path
import os
import json
import submissions
from ..base import base


# read in default job & processing scripts, for use in Processing and Run objects below.
default_omm_script = Path('omm_jobber.py').absolute()


class OpenMMProcessing(base):
    """Writes a script for postprocessing simulations into directory,
    then passes out string to call script. Applies kwargs as key-value pairs to
    call to format script once it's read in as string, allowing customization."""

    def __init__(self, processing_template=None, **kwargs):
        if processing_template:
            self.processing_template = processing_template.read_text()
            self.processing_script = processing_template.format(**kwargs)
            self.processing_script_path = Path('post-process.sh')


    @property
    def class_name(self):
        return "OpenMMProcessing"

    @property
    def config(self):
        return {
            'processing_script': self.processing_script,
            'processing_script_path': self.processing_script_path
        }

    def update_processing_script(self, **kwargs):
        self.processing_script = self.processing_template.format(**kwargs)

    def run(self):
        if self.processing_script:
            with open('post-process.sh', 'w') as f:
                f.write(self.processing_script)
            return './post-process.sh'
        else:
            return ''


class OpenMM(base):
    """OpenMM wrapper for running md simulations or minimizing
    structures

    Parameters
    ----------
    topology : str,
        OpenMM compatible topology file.
    python_omm_script : str,
        python_OMM_script specifying simulation to be run.
    n_cpus : int, default = 1,
        The number of cpus to use with the simulation.
    n_gpus : int, default = None,
        The number of gpus to use with the simulation. If None, will
        only use cpus.
    processing_obj : object, default = None,
        Object that when run will write a provided script for processing a
        trajectory, then return a command to run it.
    submission_obj : object,
        Submission object used for running the simulation. Look into
        SlurmSub or OSSub.
    min_run : bool, default = False,
        Is this a minimization run? Helps with output naming.
    env_exports : str, default=None,
        A list of commands to submit before running a job.
    """
    def __init__(self, topology_fn, integrator_fn, system_fn, state_fn, steps,
                 temperature, submission_obj,
                 python_omm_script=default_omm_script,
                 write_freq=10000,
                 output_prefix='frame0',
                 plaform_name=None,
                 platform_properties=None,
                 output_reporters=None,
                 processing_obj=None,
                 min_run=False,
                 env_exports=None):
        # bind openmm script path
        self.output_dir = None
        self.start_name = None
        self.omm_script_p = python_omm_script

        # set up a bunch of reporters by default. These'll be written to a space delimited outfile.
        if not output_reporters:
            self.output_reporters = dict(
                step=True,
                time=True,
                potentialEnergy=True,
                totalEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                elapsed_time=True,
                totalSteps=steps,
                separator=' '
            )
        # set the default platform properties to be mixed prec.
        if not platform_properties:
            self.platform_properties = {'Precision': 'mixed'}


        # read in the simulation components for the study.
        self.topology_p = Path(topology_fn).absolute()
        self.integrator_p = Path(integrator_fn).absolute()
        self.system_p = Path(system_fn).absolute()
        self.state_p = Path(state_fn).absolute()
        self.temperature = temperature
        self.write_freq = write_freq
        self.output_prefix = output_prefix
        self.platform_name = plaform_name

        if env_exports is None:
            self.env_exports = ''
        else:
            self.env_exports = env_exports
        # If processing_obj is not specified use openMMProcessing as default.
        if processing_obj is None:
            self.processing_obj = OpenMMProcessing()
        else:
            self.processing_obj = processing_obj
        # bind submission obj
        self.submission_obj = submission_obj
        self.min_run = min_run

    @property
    def class_name(self):
        return 'OpenMM'

    @property
    def config(self):
        return {
            'topology': str(self.topology_p),
            'integrator': str(self.integrator_p),
            'system': str(self.system_p),
            'state': str(self.state_p),
            'out reporters': self.output_reporters,
            'platform name': self.platform_name,
            'platform properties': self.platform_properties,
            'python_OMM_script': self.omm_script_p,
            'write_freq': self.write_freq,
            'output_prefix': self.output_prefix,
            'processing_obj': self.processing_obj,
            'submission_obj': self.submission_obj,
            'min_run': self.min_run,
            'env_exports': self.env_exports,
            'start_name': self.start_name,
            'output_dir': self.output_dir
        }

    def setup_run(self, struct, output_dir=None):
        # set output directory
        self.output_dir = output_dir
        if self.output_dir is None:
            self.output_dir = "./"
        self.output_dir = os.path.abspath(self.output_dir)
        # generate directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            tools.run_commands('mkdir ' + self.output_dir)
        # determine starting structure filename
        if type(struct) is md.Trajectory:
            struct.save_gro(self.output_dir + '/start.gro')
            self.start_name = self.output_dir + '/start.gro'
        else:
            self.start_name = os.path.abspath(struct)
        # openmm jobber reads the config file, so dump it into the output dir.
        config_path = Path(self.output_dir)/'config.json'
        with config_path.open('w') as f:
            json.dump(self.config, f, indent=4)
        return config_path


    def run(self, struct, output_dir=None, check_continue=True):
        # setup_run
        config_path = self.setup_run(struct=struct, output_dir=output_dir)
        if self.min_run:
            base_output_name = 'em'
        else:
            base_output_name = 'md'
        # python run command to put in the jobscript
        cmd = 'python {jobber} {conf}'.format(jobber=self.omm_script_p,
                                              conf=config_path)
        job_id = self.submission_obj.run(cmd, output_dir=output_dir)
        return job_id

