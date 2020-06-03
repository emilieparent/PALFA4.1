import os.path

# All of the following 3 are simple strings with no
# pre-defined options
institution = 'McGill-beluga'
pipeline = "PRESTO4"
survey = "PALFA4"

# This is the root directory of the source for the pipeline as pulled
# from github
pipelinedir = "/project/rrg-vkaspi-ad/PALFA4/software/pipeline2.0" 

# This is the root directory of the source for the psrfits_utils code
# as pulled from github
psrfits_utilsdir = "/project/rrg-vkaspi-ad/software/src/psrfits_utils"

# A boolean value that determines if raw data is deleted when results
# for a job are successfully uploaded to the common DB, or if the
# maximum number of attempts for a job is reached.
delete_rawdata = True

# Should not need to change this unless you rearrange the pipeline filesystem
wapp_coords_table = os.path.join(pipelinedir, "lib", "PALFA_wapp_coords_table.txt")
mock_coords_table = os.path.join(pipelinedir, "lib", "PALFA_mock_coords_table.txt")
log_dir = "/project/rrg-vkaspi-ad/PALFA4/logs/"  
qsublog_dir = os.path.join(log_dir, "qsublog")

import basic_check
basic_check.basic.populate_configs(locals())
basic_check.basic.check_sanity()

