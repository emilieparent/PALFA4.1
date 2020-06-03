################################################################
# Common Database Configuration
################################################################
# These are the username / password and host at Cornell to be used to
# download data.  Adam will supply these.  Note that you need to have
# your IP address approved in order to pull data
username = 'mcgill'
password = 'M0nt43@l!'
host = 'narwell.tc.cornell.edu'

import commondb_check
commondb_check.commondb.populate_configs(locals())
commondb_check.commondb.check_sanity()
