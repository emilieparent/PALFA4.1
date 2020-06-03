import os.path
import queue_managers.slurm_beluga
import config.basic
################################################################
# JobPooler Configurations
################################################################
max_jobs_running = 650 # Maximum number of running jobs: 1000 on Beluga
max_jobs_queued = 850  # Can be kept small so that you don't hog the queue (>=1)
max_jobs_cluster = 1000
max_attempts = 15 #maximum number of times a job is attempted due to errors
submit_sleep = 2 #time for jobpooler to sleep after submitting a job in seconds (prevents overloading queue manager)
obstime_limit = 60 # lower limit in seconds of observation time for jobs to be submitted

# Arguments to moab.MoabManager are (node name, msub flags, walltime)

queue_manager = queue_managers.slurm_beluga.SLURMManager("%s_batchjob" % config.basic.survey, \
                                    walltime_per_gb=35.15, rapID="rrg-vkaspi-ad")


# Use the following to use an alternative script that gets submitted
# to the worker nodes. e.g. A script that sets things up before calling
# search.py.
#alternative_submit_script = '/sb/project/bgf-180-aa/pipeline/pipeline2.0/bin/search_onnode_python.sh'
alternative_submit_script = None

import jobpooler_check
jobpooler_check.jobpooler.populate_configs(locals())
jobpooler_check.jobpooler.check_sanity()
