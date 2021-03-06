import os.path
import sys
import os
import shutil
import time
import re
import threading
import traceback

import M2Crypto

import debug
import mailer
import OutStream
import datafile
import jobtracker
import CornellFTP
import CornellWebservice
import pipeline_utils
import config.background
import config.download
import config.email
import config.basic
dlm_cout = OutStream.OutStream("Download Module", \
                        os.path.join(config.basic.log_dir, "downloader.log"), \
                        config.background.screen_output)


def check_download_attempts():
    """For each download attempt with status 'downloading' check
        to see that its thread is still active. If not, mark it
        as 'unknown', and mark the file as 'unverified'.
    """
    attempts = jobtracker.query("SELECT * FROM download_attempts " \
                                "WHERE status='downloading'")

    active_ids = [int(t.getName()) for t in threading.enumerate() \
                            if isinstance(t, DownloadThread)]

    k=1
    for attempt in attempts:
        if attempt['id'] not in active_ids:
            dlm_cout.outs("Download attempt (ID: %d) is no longer running." % \
                            attempt['id'])
            queries = []
            queries.append("UPDATE files " \
                           "SET status='unverified', " \
                                "updated_at='%s', " \
                                "details='Download thread is no longer running' "
                           "WHERE id=%d" % (jobtracker.nowstr(), attempt['file_id']))
            queries.append("UPDATE download_attempts " \
                           "SET status='unknown', " \
                                "updated_at='%s', " \
                                "details='Download thread is no longer running' "
                           "WHERE id=%d" % (jobtracker.nowstr(), attempt['id']))
            jobtracker.query(queries)
	if k%10==0:
		print k,'/',len(attempts)
	k+=1


def can_request_more():
    """Returns whether Downloader can request more restores.
        This is based on took disk space allowed for downloaded
        file, disk space available on the file system, and maximum
        number of active requests allowed.

    Inputs:
        None
    Output:
        can_request: A boolean value. True if Downloader can make a request.
                        False otherwise.
    """
    # Note: Files are restored in pairs (so we multiply by 2)
    active_requests = jobtracker.query("SELECT IFNULL(SUM(numrequested)*2, 0) " \
                                       "FROM requests " \
                                       "WHERE status='waiting'", fetchone=True)[0]
    to_download = jobtracker.query("SELECT * FROM files " \
                                   "WHERE status NOT IN ('downloaded', " \
                                                        "'added', " \
                                                        "'deleted', " \
                                                        "'terminal_failure')")
    num_to_restore = active_requests
    num_to_download = len(to_download)
    used = get_space_used()
    avail = get_space_available()
    reserved = get_space_committed()

    can_request = ((num_to_restore+num_to_download) < config.download.numrestored) and \
            (avail-reserved > config.download.min_free_space) and \
            (used+reserved < config.download.space_to_use)
    return can_request


def get_space_used():
    """Return space used by the download directory (config.download.datadir)

    Inputs:
        None
    Output:
        used: Size of download directory (in bytes)
    """
    files = jobtracker.query("SELECT * FROM files " \
                             "WHERE status IN ('added', 'downloaded', 'unverified')")

    total_size = 0
    for file in files:
        if os.path.exists(file['filename']):
            total_size += file['size']
    return total_size


def get_space_available():
    """Return space available on the file system where files
        are being downloaded.
    
        Inputs:
            None
        Output:
            avail: Number of bytes available on the file system.
    """
    s = os.statvfs(os.path.abspath(config.download.datadir))
    total = s.f_bavail*s.f_bsize
    return total


def get_space_committed():
    """Return space reserved to files to be downloaded.

        Inputs:
            None
        Outputs:
            reserved: Number of bytes reserved by files to be downloaded.
    """
    reserved = jobtracker.query("SELECT SUM(size) FROM files " \
                                "WHERE status IN ('downloading', 'new', " \
                                                 "'retrying', 'failed')", \
                                fetchone=True)[0]
    if reserved is None:
        reserved = 0
    return reserved


def run():
    """Perform a single iteration of the downloader's loop.

        Inputs:
            None
        Outputs:
            numsuccess: The number of successfully downloaded files 
                        this iteration.
    """
    try:
        pipeline_utils.get_zaplist_tarball(verbose=True)
	print "Got zaplists"
    #except CornellFTP.M2Crypto.ftpslib.error_perm:
    except :
        exctype, excvalue, exctb = sys.exc_info()
        dlm_cout.outs("FTP error getting zaplist tarball.\n" \
                        "\tError: %s" % \
                        ("".join(traceback.format_exception_only(exctype, excvalue)).strip()))
        return 0
    check_active_requests2()
    start_downloads2()
    check_download_attempts()
    numsuccess = verify_files()
    recover_failed_downloads()
    check_downloading_requests()
    delete_downloaded_files()
    #if can_request_more():
    make_request()
    print "Number of sucessfull runs: ",numsuccess
    return numsuccess


def make_request(num_beams=None):
    """Make a request for data to be restored by connecting to the
        web services at Cornell.
    """
    if not num_beams:
        num_beams = get_num_to_request()
        if not num_beams:
            # Request size is 0
            return
    dlm_cout.outs("Requesting data\nIssuing a request of size %d" % num_beams)

    web_service = CornellWebservice.Client()
    guid = web_service.Restore(username=config.download.api_username, \
                               pw=config.download.api_password, \
                               pipeline=config.basic.pipeline.lower(), \
                               number=num_beams, \
                               bits=config.download.request_numbits, \
                               fileType=config.download.request_datatype)
    print "guid:",guid
    if guid == "fail":
        raise pipeline_utils.PipelineError("Request for restore returned 'fail'.")

    requests = jobtracker.query("SELECT * FROM requests " \
                             "WHERE guid='%s'" % guid)
    if requests:
        # Entries in the requests table exist with this GUID!?
        raise pipeline_utils.PipelineError("There are %d requests in the " \
                                           "job-tracker DB with this GUID %s" % \
                                           (len(requests), guid))

    jobtracker.query("INSERT INTO requests ( " \
                        "numbits, " \
                        "numrequested, " \
                        "file_type, " \
                        "guid, " \
                        "created_at, " \
                        "updated_at, " \
                        "status, " \
                        "details) " \
                     "VALUES (%d, %d, '%s', '%s', '%s', '%s', '%s', '%s')" % \
                     (config.download.request_numbits, num_beams, \
                        config.download.request_datatype, guid, \
                        jobtracker.nowstr(), jobtracker.nowstr(), 'waiting', \
                        'Newly created request'))
  

def check_active_requests():
    """Check for any requests with status='waiting'. If there are
        some, check if the files are ready for download.
    """
    active_requests = jobtracker.query("SELECT * FROM requests " \
                                       "WHERE status='waiting'")
    if not active_requests:
	print "Checking active requests: no active requests"
    
    web_service = CornellWebservice.Client()
    for request in active_requests:
        location = web_service.Location(guid=request['guid'], \
                                        username=config.download.api_username, \
                                        pw=config.download.api_password)
	print "location:",location,"\t\t"," guid:",request['guid']
        if location == "done":
            dlm_cout.outs("Restore (%s) is done. Will create file entries." % \
                            request['guid'])
            create_file_entries(request)
        else:

            query = "SELECT (julianday('%s')-julianday(created_at))*24 " \
                        "AS deltaT_hours " \
                    "FROM requests " \
                    "WHERE guid='%s'" % \
                        (jobtracker.nowstr(), request['guid'])

            row = jobtracker.query(query, fetchone=True)
            if row['deltaT_hours'] > config.download.request_timeout:
                dlm_cout.outs("Restore (%s) is over %d hr old " \
                                "and still not ready. Marking " \
                                "it as failed." % \
                        (request['guid'], config.download.request_timeout))
                jobtracker.query("UPDATE requests " \
                                 "SET status='failed', " \
                                    "details='Request took too long (> %d hr)', " \
                                    "updated_at='%s' " \
                                 "WHERE guid='%s'" % \
                    (config.download.request_timeout, jobtracker.nowstr(), \
                            request['guid']))

def check_active_requests2():
    """Check for any requests with status='waiting'. If there are
        some, check if the files are ready for download.
    """
    active_requests = jobtracker.query("SELECT * FROM requests " \
                                       "WHERE status='waiting'")
    if not active_requests:
	print "Checking active requests: no active requests"
    
    for request in active_requests:
	try:
		cftp = CornellFTP.CornellFTP()
		files = cftp.get_files(str(request['guid']))
	except CornellFTP.M2Crypto.ftpslib.error_perm:
	#except :
        	exctype, excvalue, exctb = sys.exc_info()
		dlm_cout.outs("FTP error getting file information.\n" \
        	                "\tGUID: %s\n\tError: %s" % (request['guid'], \
	                        "".join(traceback.format_exception_only(exctype, excvalue)).strip()))
        	files = []
	    # Couldn't get files on server 
            	print "Errors with restoring .."
            	query = "SELECT (julianday('%s')-julianday(created_at))*24 " \
                        "AS deltaT_hours " \
                    "FROM requests " \
                    "WHERE guid='%s'" % \
                        (jobtracker.nowstr(), request['guid'])

            	row = jobtracker.query(query, fetchone=True)
            	if row['deltaT_hours'] > config.download.request_timeout:
                	dlm_cout.outs("Restore (%s) is over %d hr old " \
                                "and still not ready. Marking " \
                                "it as failed." % \
                        (request['guid'], config.download.request_timeout))
                	jobtracker.query("UPDATE requests " \
                                 "SET status='failed', " \
                                    "details='Request took too long (> %d hr)', " \
                                    "updated_at='%s' " \
                                 "WHERE guid='%s'" % \
                    	(config.download.request_timeout, jobtracker.nowstr(), \
                            request['guid']))
		#cftp.close()


	else:
		cftp.close()
        	dlm_cout.outs("Restore (%s) is done. Will create file entries." % \
                            request['guid'])
		create_file_entries2(request,files)



def create_file_entries(request):
    """Given a row from the requests table in the job-tracker DB
        check the FTP server for its files and create entries in
        the files table.

        Input:
            request: A row from the requests table.
        Outputs:
            None
    """
    cftp = CornellFTP.CornellFTP()
    try:
        files = cftp.get_files(request['guid'])
    except CornellFTP.M2Crypto.ftpslib.error_perm:
        exctype, excvalue, exctb = sys.exc_info()
        dlm_cout.outs("FTP error getting file information.\n" \
                        "\tGUID: %s\n\tError: %s" % \
                        (request['guid'], \
                        "".join(traceback.format_exception_only(exctype, excvalue)).strip()))
        files = []
    print "Create_file_entries : %s new files "%str(len(files))
    total_size = 0
    num_files = 0
    queries = []
    kkk = 1 
    for fn, size in files:
	if kkk%10==0:
		print '\n',int(kkk),'/',len(files),'\n'
	kkk+=1
        if not pipeline_utils.can_add_file(fn,verbose=True):
            dlm_cout.outs("Skipping %s" % fn)
            continue

        # Insert entry into DB's files table
        queries.append("INSERT INTO files ( " \
                            "request_id, " \
                            "remote_filename, " \
                            "filename, " \
                            "status, " \
                            "created_at, " \
                            "updated_at, " \
                            "size) " \
                       "VALUES ('%s', '%s', '%s', '%s', '%s', '%s', %d)" % \
                       (request['id'], fn, os.path.join(config.download.datadir, fn), \
                        'new', jobtracker.nowstr(), jobtracker.nowstr(), size))
        total_size += size
        num_files += 1

    if num_files:
        dlm_cout.outs("Request (GUID: %s) has succeeded.\n" \
                        "\tNumber of files to be downloaded: %d" % \
                        (request['guid'], num_files))
        queries.append("UPDATE requests " \
                       "SET size=%d, " \
                            "updated_at='%s', " \
                            "status='downloading', " \
                            "details='Request has been filled' " \
                       "WHERE id=%d" % \
                       (total_size, jobtracker.nowstr(), request['id']))
    else:
        dlm_cout.outs("Request (GUID: %s) has failed.\n" \
                        "\tThere are no files to be downloaded." % \
                        request['guid'])

        # delete restore since there may be skipped files
        web_service = CornellWebservice.Client()
        delete_status = web_service.Deleter(guid=request['guid'], \
                                            username=config.download.api_username, \
                                            pw=config.download.api_password)
        if delete_status == "deletion successful":
            dlm_cout.outs("Deletion (%s) succeeded." % request['guid'])
	elif delete_status == "invalid user":
	    dlm_cout.outs("Deletion (%s) failed due to invalid user." % \
			  request['guid'])
	elif delete_status == "deletion failed":
	    dlm_cout.outs("Deletion (%s) failed for unknown reasons." % \
			  request['guid'])

	# redefine 'queries' because there are no files to update
	queries = ["UPDATE requests " \
		   "SET updated_at='%s', " \
			"status='failed', " \
			"details='No files to download.' " \
		   "WHERE id=%d" % \
		   (jobtracker.nowstr(), request['id'])]

    jobtracker.query(queries)


def create_file_entries2(request,files):
    """Given a row from the requests table in the job-tracker DB
        check the FTP server for its files and create entries in
        the files table.

        Input:
            request: A row from the requests table.
        Outputs:
            None
    """
    cftp = CornellFTP.CornellFTP()
    try:
        files = cftp.get_files(request['guid'])
    except CornellFTP.M2Crypto.ftpslib.error_perm:
        exctype, excvalue, exctb = sys.exc_info()
        dlm_cout.outs("FTP error getting file information.\n" \
                        "\tGUID: %s\n\tError: %s" % \
                        (request['guid'], \
                        "".join(traceback.format_exception_only(exctype, excvalue)).strip()))
        files = []
    print "Create_file_entries : %s new files "%str(len(files))
    total_size = 0
    num_files = 0
    queries = []
    k = 1 
    for fn, size in files:
	if k%10==0:
		print k,'/',len(files)
	k+=1
        if not pipeline_utils.can_add_file(fn,verbose=True):
            dlm_cout.outs("Skipping %s" % fn)
            continue

        # Insert entry into DB's files table
        queries.append("INSERT INTO files ( " \
                            "request_id, " \
                            "remote_filename, " \
                            "filename, " \
                            "status, " \
                            "created_at, " \
                            "updated_at, " \
                            "size) " \
                       "VALUES ('%s', '%s', '%s', '%s', '%s', '%s', %d)" % \
                       (request['id'], fn, os.path.join(config.download.datadir, fn), \
                        'new', jobtracker.nowstr(), jobtracker.nowstr(), size))
        total_size += size
        num_files += 1

    if num_files:
        dlm_cout.outs("Request (GUID: %s) has succeeded.\n" \
                        "\tNumber of files to be downloaded: %d" % \
                        (request['guid'], num_files))
        queries.append("UPDATE requests " \
                       "SET size=%d, " \
                            "updated_at='%s', " \
                            "status='downloading', " \
                            "details='Request has been filled' " \
                       "WHERE id=%d" % \
                       (total_size, jobtracker.nowstr(), request['id']))
    else:
        dlm_cout.outs("Request (GUID: %s) has failed.\n" \
                        "\tThere are no files to be downloaded." % \
                        request['guid'])

        # delete restore since there may be skipped files
	"""
        web_service = CornellWebservice.Client()
        delete_status = web_service.Deleter(guid=request['guid'], \
                                            username=config.download.api_username, \
                                            pw=config.download.api_password)
        if delete_status == "deletion successful":
            dlm_cout.outs("Deletion (%s) succeeded." % request['guid'])
	elif delete_status == "invalid user":
	    dlm_cout.outs("Deletion (%s) failed due to invalid user." % \
			  request['guid'])
	elif delete_status == "deletion failed":
	    dlm_cout.outs("Deletion (%s) failed for unknown reasons." % \
			  request['guid'])
	"""
	# redefine 'queries' because there are no files to update
	queries = ["UPDATE requests " \
		   "SET updated_at='%s', " \
			"status='failed', " \
			"details='No files to download.' " \
		   "WHERE id=%d" % \
		   (jobtracker.nowstr(), request['id'])]

    jobtracker.query(queries)

def start_downloads2():
    """Check for entries in the files table with status 'retrying'
        or 'new' and start the downloads.
    """
    todownload  = jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='retrying' " \
                                   "ORDER BY created_at ASC")
    todownload += jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='new' " \
                                   "ORDER BY created_at ASC")
    kkk=1
    for file in todownload:
	if kkk%10==0:
		print kkk,'/',len(todownload)
	kkk+=1
	dlm_cout.outs("Initiating download of %s" % \
                            os.path.split(file['filename'])[-1])

            # Update file status and insert entry into download_attempts
 	queries = []
        queries.append("UPDATE files " \
                           "SET status='downloading', " \
                                "details='Initiated download', " \
                                "updated_at='%s' " \
                            "WHERE id=%d" % \
                            (jobtracker.nowstr(), file['id']))
        queries.append("INSERT INTO download_attempts (" \
                                "status, " \
                                "details, " \
                                "updated_at, " \
                                "created_at, " \
                                "file_id) " \
                           "VALUES ('%s', '%s', '%s', '%s', %d)" % \
                           ('downloading', 'Initiated download', jobtracker.nowstr(), \
                                jobtracker.nowstr(), file['id']))
	insert_id = jobtracker.query(queries)
	attempt = jobtracker.query("SELECT * FROM download_attempts " \
                                       "WHERE id=%d" % insert_id, \
                                       fetchone=True)


def start_downloads():
    """Check for entries in the files table with status 'retrying'
        or 'new' and start the downloads.
    """
    todownload  = jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='retrying' " \
                                   "ORDER BY created_at ASC")
    todownload += jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='new' " \
                                   "ORDER BY created_at ASC")

    for file in todownload:
        if can_download():
            dlm_cout.outs("Initiating download of %s" % \
                            os.path.split(file['filename'])[-1])

            # Update file status and insert entry into download_attempts
            queries = []
            queries.append("UPDATE files " \
                           "SET status='downloading', " \
                                "details='Initiated download', " \
                                "updated_at='%s' " \
                            "WHERE id=%d" % \
                            (jobtracker.nowstr(), file['id']))
            queries.append("INSERT INTO download_attempts (" \
                                "status, " \
                                "details, " \
                                "updated_at, " \
                                "created_at, " \
                                "file_id) " \
                           "VALUES ('%s', '%s', '%s', '%s', %d)" % \
                           ('downloading', 'Initiated download', jobtracker.nowstr(), \
                                jobtracker.nowstr(), file['id']))
            insert_id = jobtracker.query(queries)
            attempt = jobtracker.query("SELECT * FROM download_attempts " \
                                       "WHERE id=%d" % insert_id, \
                                       fetchone=True)
    
            # download(attempt)
            DownloadThread(attempt).start()
        else:
            break


def get_num_to_request():
    """Return the number of files to request given the average
        time to download a file (including waiting time) and
        the amount of space available.

        Inputs:
            None

        Outputs:
            num_to_request: The size of the request.
    """
    ALLOWABLE_REQUEST_SIZES = [5,10,20,50,100,200]
    avgrate = jobtracker.query("SELECT AVG(files.size/" \
                                "(JULIANDAY(download_attempts.updated_at) - " \
                                "JULIANDAY(download_attempts.created_at))) " \
                               "FROM files, download_attempts " \
                               "WHERE files.id=download_attempts.file_id " \
                                    "AND download_attempts.status='downloaded'", \
                               fetchone=True)[0]
    avgsize = jobtracker.query("SELECT AVG(size/numrequested) FROM requests " \
                               "WHERE numbits=%d AND " \
                                    "file_type='%s'" % \
                                (config.download.request_numbits, \
                                    config.download.request_datatype.lower()), \
                                fetchone=True)[0]
    if avgrate is None or avgsize is None:
        return min(ALLOWABLE_REQUEST_SIZES)

    # Total number requested that can be downloaded per day (on average).
    max_to_request_per_day = avgrate/avgsize
    
    used = get_space_used()
    avail = get_space_available()
    reserved = get_space_committed()
    
    # Maximum number of bytes that we should request
    max_bytes = min([avail-reserved-config.download.min_free_space, \
                        config.download.space_to_use-reserved-used])
    # Maximum number to request
    max_num = max_bytes/avgsize

    ideal_num_to_request = min([max_num, max_to_request_per_day])

    if debug.DOWNLOAD:
        print "Average dl rate: %.2f bytes/day" % avgrate
        print "Average size per request unit: %d bytes" % avgsize
        print "Max can dl per day: %d" % max_to_request_per_day
        print "Max num to request: %d" % max_num
        print "Ideal to request: %d" % ideal_num_to_request

    # Return the closest allowable request size without exceeding
    # 'ideal_num_to_request'
    num_to_request = max([0]+[N for N in ALLOWABLE_REQUEST_SIZES \
                            if N <= ideal_num_to_request])
    return num_to_request


def can_download():
    """Return true if another download can be initiated.
        False otherwise.

        Inputs:
            None
        Output:
            can_dl: A boolean value. True if another download can be
                    initiated. False otherwise.
    """
    downloading = jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='downloading'")
    numdownload = len(downloading)
    used = get_space_used()
    avail = get_space_available()
    
    can_dl = (numdownload < config.download.numdownloads) and \
            (avail > config.download.min_free_space) and \
            (used < config.download.space_to_use)
    return can_dl 


def download(attempt):
    """Given a row from the job-tracker's download_attempts table,
        actually attempt the download.
    """
    file = jobtracker.query("SELECT * FROM files " \
                            "WHERE id=%d" % attempt['file_id'], \
                            fetchone=True)
    request = jobtracker.query("SELECT * FROM requests " \
                               "WHERE id=%d" % file['request_id'], \
                               fetchone=True)

    queries = []
    try:
        cftp = CornellFTP.CornellFTP()
        cftp.download(os.path.join(request['guid'], file['remote_filename']), guid=request['guid'], preserve_modtime=False)
    except Exception, e:
        queries.append("UPDATE files " \
                       "SET status='failed', " \
                            "updated_at='%s', " \
                            "details='Download failed - %s' " \
                       "WHERE id=%d" % \
                       (jobtracker.nowstr(), str(e), file['id']))
        queries.append("UPDATE download_attempts " \
                       "SET status='download_failed', " \
                            "details='Download failed - %s', " \
                            "updated_at='%s' " \
                       "WHERE id=%d" % \
                       (str(e), jobtracker.nowstr(), attempt['id']))
    else:
        queries.append("UPDATE files " \
                       "SET status='unverified', " \
                            "updated_at='%s', " \
                            "details='Download is complete - File is unverified' " \
                       "WHERE id=%d" % \
                       (jobtracker.nowstr(), file['id']))
        queries.append("UPDATE download_attempts " \
                       "SET status='complete', " \
                            "details='Download is complete', " \
                            "updated_at='%s' " \
                       "WHERE id=%d" % \
                       (jobtracker.nowstr(), attempt['id']))
    if queries:
	jobtracker.query(queries)


def verify_files():
    """For all downloaded files with status 'unverify' verify the files.
        
        Inputs:
            None
        Output:
            numverified: The number of files successfully verified.
    """
    toverify = jobtracker.query("SELECT * FROM files " \
                                "WHERE status='unverified'")
    k = 1
    numverified = 0
    for file in toverify:
	if k%10==0:
		print k,'/',len(toverify)
	k+=1
        if os.path.exists(file['filename']):
            actualsize = os.path.getsize(file['filename'])
        else:
            # Check if download.datadir has changed since file entry was created
            # and if that is why file missing.
            alt_path = os.path.join(config.download.datadir,os.path.split(file['filename'])[-1])
            if os.path.exists(alt_path):
                actualsize = os.path.getsize(alt_path)
                jobtracker.query("UPDATE files SET filename='%s' WHERE id=%d" % (alt_path, file['id']))
            else:
                actualsize = -1

        expectedsize = file['size']

        last_attempt_id = jobtracker.query("SELECT id " \
                                           "FROM download_attempts " \
                                           "WHERE file_id=%s " \
                                           "ORDER BY id DESC " %file['id'],fetchone=True)[0]
                                                
        queries = []
        if actualsize == expectedsize:
            dlm_cout.outs("Download of %s is complete and verified." % \
                            os.path.split(file['filename'])[-1])
            # Everything checks out!
            queries.append("UPDATE files " \
                           "SET status='downloaded', " \
                                "details='Download is complete and verified', " \
                                "updated_at='%s'" \
                           "WHERE id=%d" % \
                           (jobtracker.nowstr(), file['id']))
            queries.append("UPDATE download_attempts " \
                           "SET status='downloaded', " \
                                "details='Download is complete and verified', " \
                                "updated_at='%s'" \
                           "WHERE id=%d" % \
                           (jobtracker.nowstr(), last_attempt_id))
            numverified += 1
        else:
            dlm_cout.outs("Verification of %s failed. \n" \
                            "\tActual size (%d bytes) != Expected size (%d bytes)" % \
                            (os.path.split(file['filename'])[-1], actualsize, expectedsize))
            
            # Boo... verification failed.
            queries.append("UPDATE files " \
                           "SET status='failed', " \
                                "details='Downloaded file failed verification', " \
                                "updated_at='%s'" \
                           "WHERE id=%d" % \
                           (jobtracker.nowstr(), file['id']))
            queries.append("UPDATE download_attempts " \
                           "SET status='verification_failed', " \
                                "details='Downloaded file failed verification', " \
                                "updated_at='%s'" \
                           "WHERE id=%d" % \
                           (jobtracker.nowstr(), last_attempt_id))
        jobtracker.query(queries)
    return numverified


def recover_failed_downloads():
    """For each entry in the job-tracker DB's files table
        check if the download can be retried or not.
        Update status and clean up, as necessary.
    """
    failed_files = jobtracker.query("SELECT * FROM files " \
                                   "WHERE status='failed'")
    if not failed_files:
	print "Recovering failed files: no failed files"

    for file in failed_files:
        attempts = jobtracker.query("SELECT * FROM download_attempts " \
                                    "WHERE file_id=%d" % file['id'])
        if len(attempts) < config.download.numretries:
            # download can be retried
            jobtracker.query("UPDATE files " \
                             "SET status='retrying', " \
                                  "updated_at='%s', " \
                                  "details='Download will be attempted again' " \
                             "WHERE id=%s" % \
                             (jobtracker.nowstr(), file['id']))
        else:
            # Abandon this file
            if os.path.exists(file['filename']):
                os.remove(file['filename'])
            jobtracker.query("UPDATE files " \
                             "SET status='terminal_failure', " \
                                  "updated_at='%s', " \
                                  "details='This file has been abandoned' " \
                             "WHERE id=%s" % \
                             (jobtracker.nowstr(), file['id']))

    
def status():
    """Print downloader's status to screen.
    """
    used = get_space_used()
    avail = get_space_available()
    allowed = config.download.space_to_use
    print "Space used by downloaded files: %.2f GB of %.2f GB (%.2f%%)" % \
            (used/1024.0**3, allowed/1024.0**3, 100.0*used/allowed)
    print "Space available on file system: %.2f GB" % (avail/1024.0**3)

    numwait = jobtracker.query("SELECT COUNT(*) FROM requests " \
                               "WHERE status='waiting'", \
                               fetchone=True)[0]
    numfail = jobtracker.query("SELECT COUNT(*) FROM requests " \
                               "WHERE status='failed'", \
                               fetchone=True)[0]
    print "Number of requests waiting: %d" % numwait
    print "Number of failed requests: %d" % numfail

    numdlactive = jobtracker.query("SELECT COUNT(*) FROM files " \
                                   "WHERE status='downloading'", \
                                   fetchone=True)[0]
    numdlfail = jobtracker.query("SELECT COUNT(*) FROM files " \
                                 "WHERE status='failed'", \
                                 fetchone=True)[0]
    print "Number of active downloads: %d" % numdlactive
    print "Number of failed downloads: %d" % numdlfail


def check_downloading_requests():
    requests = jobtracker.query("SELECT * FROM requests "\
                                "WHERE status='downloading'")
    if len(requests) > 0:
        queries = []
	print "Number of downloading requests: %d \n"%len(requests)
	all_downloaded = 0
	all_files_requested = 0
        for request in requests:
            files_in_request = jobtracker.query("SELECT * FROM files "\
                                                "WHERE request_id=%d" % \
                                                request['id'])
            downloaded_files = 0
            all_files_requested += len(files_in_request)
            for f in files_in_request:
                if f['status'] in ('downloaded','deleted'): downloaded_files += 1
            if downloaded_files == len(files_in_request):
                queries.append("UPDATE requests " \
                               "SET status='finished', " \
                               "details='All files downloaded', " \
                               "updated_at='%s' " \
                               "WHERE id=%d" % \
                               (jobtracker.nowstr(), request['id']))

            print "\t Request guid %s:	downloaded %d out of %d"%(request['guid'],downloaded_files,len(files_in_request)) 
            all_downloaded+=downloaded_files
        jobtracker.query(queries)
        print "\nDownloaded %d files out of %d"%(all_downloaded,all_files_requested)
    else:
	print "Checking downloading requests: no downloading requests"
        pass




def delete_downloaded_files():
    requests_to_delete = jobtracker.query("SELECT * FROM requests " \
                                          "WHERE status='finished'")
    if len(requests_to_delete) > 0:
        web_service = CornellWebservice.Client()
        queries = []
        for request_to_delete in requests_to_delete:
            delete_status = web_service.Deleter(guid=request_to_delete['guid'], \
                                                username=config.download.api_username, \
                                                pw=config.download.api_password)
            if delete_status == "deletion successful":
                dlm_cout.outs("Deletion (%s) succeeded." % request_to_delete['guid'])
                queries.append("UPDATE requests " \
                               "SET status='cleaned_up', " \
                               "details='Files deleted from FTP server', " \
                               "updated_at='%s' " \
                               "WHERE id=%d" % \
                               (jobtracker.nowstr(), request_to_delete['id']))

            elif delete_status == "invalid user":
                dlm_cout.outs("Deletion (%s) failed due to invalid user." % \
                              request_to_delete['guid'])
            
            elif delete_status == "deletion failed":
                dlm_cout.outs("Deletion (%s) failed for unknown reasons." % \
                              request_to_delete['guid'])
        jobtracker.query(queries)
    else: pass
            


class DownloadThread(threading.Thread):
    """A sub-class of threading.Thread to download restored
        file from Cornell.
    """
    def __init__(self, attempt):
        """DownloadThread constructor.
            
            Input:
                attempt: A row from the job-tracker's download_attempts table.

            Output:
                self: The DownloadThread object constructed.
        """
        super(DownloadThread, self).__init__(name=attempt['id'])
        self.attempt = attempt

    def run(self):
        """Download data as a separate thread.
        """
        download(self.attempt)
