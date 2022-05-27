import mdtraj as md
import numpy as np
from .. import tools
from pathlib import Path
import os
import submissions
from ..base import base


# read in default job & processing scripts, for use in Processing and Run objects below.
default_omm_script = Path('omm_jobber.py').absolute()
default_processing_script = Path('default-aligner.sh', 'r').absolute()


class OpenMMProcessing(base):
    """Writes a script for postprocessing simulations into directory, then passes out string that can call script."""

    def __init__(self, processing_script=default_processing_script, **kwargs):
        self.processing_script = processing_script.read_text()
        if processing_script:
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
    def __init__(self, topology_fn, integrator_fn, system_fn, state_fn, steps, temperature, submission_obj,
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
            self.processing_obj = OpenMMProcessing(self.topology_fn)
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
            'env_exports': self.env_exports
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
        # move over additional topology files
        if self.itp_files is not None:
            cmds = []
            for filename in self.itp_files:
                cmds.append('cp ' + filename + ' ' + self.output_dir + ' -r')
            tools.run_commands(cmds)
        return

    def run(self, struct, output_dir=None, check_continue=True):
        # setup_run
        self.setup_run(struct=struct, output_dir=output_dir)
        if self.min_run:
            base_output_name = 'em'
        else:
            base_output_name = 'md'
        # source command
        if self.source_file is None:
            source_cmd = ''
        else:
            source_cmd = 'source ' + self.source_file + '\n\n'
        # generate grompp command
        grompp_cmd = 'gmx grompp -f ' + self.mdp_file + ' -c ' + \
                     self.start_name + ' -p ' + self.topology_fn + ' -o ' + \
                     base_output_name + ' -maxwarn ' + self.max_warn
        # optionally add an index file
        if self.index_file is not None:
            grompp_cmd +=  ' -n ' + self.index_file
        grompp_cmd += '\n'
        # generate mdrun command
        # JRP added '-cpi state -g md' on 07-01-2019
        mdrun_cmd = 'gmx mdrun -cpi state -g md -s ' + base_output_name + ' -o ' + \
            base_output_name + ' -c after_' + base_output_name + ' -v -nt ' + \
            str(self.n_cpus)
        # if an MD run, make default name for trajectory
        if not self.min_run:
            mdrun_cmd += ' -x frame0'
        # add gpus to mdrun command
        if self.n_gpus is not None:
            if self.n_cpus%self.n_gpus != 0:
                raise
            mdrun_cmd += ' -ntmpi ' + str(self.n_gpus) + ' -ntomp ' + \
                str(int(self.n_cpus/self.n_gpus))
        # adds additional keywords that are not specified
        keys = list(self.kwargs.keys())
        values = list(self.kwargs.values())
        additions = " ".join(
            ['-' + i[0] + ' ' + i[1] for i in np.transpose([keys, values])])
        mdrun_cmd += ' ' + additions + '\n'
        # check for previous tpr for continuation
        if check_continue:
            tpr_filename = self.output_dir + "/md.tpr"
            bash_check_cmd = 'if [ ! -f "%s" ]; then\n' % tpr_filename
            bash_check_cmd += '    echo "Didn\'t find md.tpr, running grompp..."\n'
            bash_check_cmd += '    ls\n    pwd\n    %s' % grompp_cmd
            bash_check_cmd += 'else\n    echo "Found md.tpr, not running grompp"\nfi\n\n'
            grompp_cmd = bash_check_cmd
        # combine commands and submit to submission object
        cmds = [self.env_exports, source_cmd, grompp_cmd, mdrun_cmd]
        cmds.append(self.processing_obj.run())
        job_id = self.submission_obj.run(cmds, output_dir=output_dir)
        return job_id

