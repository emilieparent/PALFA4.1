################################################################
# Email Notification Configuration
################################################################
"""
Note: admins on Beluga will not allow users to use SMTP for security reason.
Can't use the mailer in the pipeline at all: must use Slurm (e.g. #SBATCH --mail-type=FAIL) 
Emilie Parent, July 2019
"""

enabled = False   # whether error email is sent or not
#smtp_host = 'smtp.gmail.com' # None - For use of the local smtp server
#smtp_host = '10.241.128.1' # None - For use of the local smtp server

smtp_host = '10.74.73.3'  # SET FOR BELUGA
smtp_port = 25 # Port to use for connecting to SMTP mail server (should be 25 or 587)

smtp_username = 'mcgill.pipeline'
smtp_password = 'uni3600rphys227'
smtp_login = False # Whether username/password are used to log into SMTP server
smtp_usetls = False # Whether Transport Layer Security (TLS) is used
#recipient = 'pscholz0520+pipeline@gmail.com' # The address to send emails to
#sender = 'pscholz' # From address to show in email
recipient = 'emilieparent010@gmail.com' # The address to send emails to
sender = 'eparent' # From address to show in email
# Every "error" gives a failure...
send_on_failures = True
# After so many errors (determined in job pooler) you get a terminal failure
send_on_terminal_failures = True
# Crash is when one of the background scripts crash
send_on_crash = True
smtp_usessl = False


import email_check
email_check.email.populate_configs(locals())
email_check.email.check_sanity()
