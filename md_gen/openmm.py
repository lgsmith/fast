import mdtraj as md
import numpy as np
import os
from .. import tools
from pathlib import Path
import submissions
from ..base import base

with open('openmm_jobber.py', 'r') as f:
    default_omm_script = f.read()

class OpenMMProcessing(base):
    """Writes a script for postprocessing simulations into directory, then passes out string that can call script."""

    def __init__(self, processing_script=None):
        self.processing_script = processing_script
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
    top_file : str,
        OpenMM compatible topology file.
    python_OMM_script : str,
        python_OMM_script specifying simulation to be run.
    n_cpus : int, default = 1,
        The number of cpus to use with the simulation.
    n_gpus : int, default = None,
        The number of gpus to use with the simulation. If None, will
        only use cpus.
    processing_obj : object, default = None,
        Object that when run will write a provided script for projcessing a
        trajectory, then return a command to run it.
    submission_obj : object,
        Submission object used for running the simulation. Look into
        SlurmSub or OSSub.
    min_run : bool, default = False,
        Is this a minimization run? Helps with output naming.
    env_exports : str, default=None,
        A list of commands to submit before running a job.
    """
    def __init__(self, top_file, submission_obj, python_OMM_script=default_omm_script, n_cpus=1, n_gpus=None,
                 processing_obj=None, index_file=None, itp_files=None,
                 min_run=False,
                 env_exports=None, **kwargs):
        self.top_file = os.path.abspath(top_file)
        self.mdp_file = os.path.abspath(mdp_file)
        self.n_cpus = n_cpus
        self.n_gpus = n_gpus
        if env_exports is None:
            self.env_exports = ''
        else:
            self.env_exports = env_exports
        # If processing_obj is not specified use GromaxProcessing as
        # default.
        if processing_obj is None:
            self.processing_obj = OpenMMProcessing()
        else:
            self.processing_obj = processing_obj
        # bind submission obj
        self.submission_obj = submission_obj
        # determine index file
        if type(index_file) is str:
            self.index_file = os.path.abspath(index_file)
        else:
            self.index_file = index_file
        # determine additional topology files
        if itp_files is None:
            self.itp_files = itp_files
        elif type(itp_files) is str:
            self.itp_files = np.array([os.path.abspath(itp_files)])
        else:
            self.itp_files = np.array(
                [os.path.abspath(top_file) for top_file in itp_files])
        # stringify max_warn
        self.max_warn = str(max_warn)
        self.min_run = min_run
        # get full path of source_file
        if source_file is None:
            self.source_file = source_file
        else:
            self.source_file = os.path.abspath(source_file)
        self.kwargs = kwargs

    @property
    def class_name(self):
        return 'OpenMM'

    @property
    def config(self):
        return {
            'top_file': self.top_file,
            'mdp_file': self.mdp_file,
            'n_cpus': self.n_cpus,
            'n_gpus': self.n_gpus,
            'processing_obj': self.processing_obj,
            'index_file': self.index_file,
            'itp_files': self.itp_files,
            'submission_obj': self.submission_obj,
            'max_warn': self.max_warn,
            'min_run': self.min_run,
            'source_file': self.source_file,
            'env_exports' : self.env_exports
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
            self.start_name + ' -p ' + self.top_file + ' -o ' + \
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

