import os.path
import config.basic
################################################################
# Downloader Configuration
################################################################
# The url is from Adam and is unlikely to need changing
api_service_url = "http://arecibo.tc.cornell.edu/palfadataapi/dataflow.asmx"
# These are from Adam and are different from the database ones...
api_username = "mcgill"
api_password = "5u43RN0v@"
# You should not need to change the port or host for FTP
ftp_host = "arecibo.tc.cornell.edu"
ftp_port = 31001
# These are from Adam and are different from the database ones...
ftp_username = "palfadata"
ftp_password = "NAIC305m"

# This is the directory where the raw data is downloaded to
# it needs read/write perms for whoever is running the downloaded
datadir = "/scratch/eparent/eparent/PALFA4/rawdata/"
space_to_use = 18000 * 2**30    # max size of downloaded data in bytes
min_free_space = 50 * 2**30  # Minimum amount of disk space on the file system
                             # that must be kept free (bytes)

numdownloads = 10  # max number of files to download at once
numrestored = 2000   # max number of files waiting to be restored + 
                  # waiting to be downloaded at any given time
numretries = 3    # max number of times to try to download before failing

request_timeout = 72 # Integer number of hours before an un-restored request 
                    # is marked as failed

# The following control the data type and format requested from Cornell
request_numbits = 4 # or 16 (which we will likely never serve)
request_datatype = 'mock' # or 'wapp'

use_lftp = True # use lftp if True, if False use M2Crypto to retrieve data

import download_check
download_check.download.populate_configs(locals())
download_check.download.check_sanity()
