import os
import warnings
import traceback
import glob
import sys
import time
import shutil
import subprocess

import debug
import datafile
import header
import candidates
import sp_candidates
import diagnostics
import jobtracker
import database
import upload
import pipeline_utils
import CornellFTP
import config.upload
import config.basic
import ratings2.utils
import ratings2.database
import M2Crypto

# Suppress warnings produced by uploaders
# (typically because data, weights, scales, offsets are missing
#       from PSRFITS files)
warnings.filterwarnings("ignore", message="Can't find the .* column")
warnings.filterwarnings("ignore", message=".*NSUBOFFS reports 0 previous rows.*")
warnings.filterwarnings("ignore", message="Channel spacing changes in file 0!")

def run():
    """
    Drives the process of uploading results of the completed jobs.

    """
    query = "SELECT * FROM jobs " \
            "WHERE status='processed' AND details in ('Ready for upload','Processed without errors') order by updated_at desc"
    query = "SELECT * FROM jobs " \
            "WHERE status='processed' AND details in ('Ready for upload','Processed without errors','Processed with warnings') order by id desc"
    processed_jobs = jobtracker.query(query)
    print "Found %d processed jobs waiting for upload" % len(processed_jobs)
    for ii, job in enumerate(processed_jobs):
       	starttime0 = time.time()
    	if debug.UPLOAD:
        	upload.upload_timing_summary = {}
        # Get the job's most recent submit
        submit = jobtracker.query("SELECT * FROM job_submits WHERE job_id=%d "%job['id']+"AND status='processed'" \
        	" AND output_dir like '/project/rrg-vkaspi-ad/%' ORDER BY id DESC", fetchone=True)
        if submit is None:
            continue
        print "Upload %d of %d" % (ii+1, len(processed_jobs))
	make_fail = False
	try:
	        upload_results(submit)
		print "Total time for upload = %.2f s)"%(time.time()-starttime0)
	except:
		if len(glob.glob(submit['output_dir']+'/*'))==0:
			print "Empty results directory! Upload failed -> Updating jobtracking database"
			#cmd = "python ~/projects/rrg-vkaspi-ad/PALFA4/software/pipeline2.0/bin/stop_processing_jobs.py "
			#cmd+= " -s %s --upload-fail"%str(submit['id'])
			#subprocess.call(cmd,shell=True)
			pipeline_utils.change_job_status_jtdb(submit,'processing_failed')
		else:
			outdir = submit['output_dir']
                        outdir_part = outdir.split('/results/')[-1]
                        print "!! Upload_failed, skipping submit_id %s  (dir: %s)"%(str(submit['id']),outdir_part)
                        #cmd = "python /home/eparent/projects/rrg-vkaspi-ad/PALFA4/software/pipeline2.0/bin/stop_processing_jobs.py "
			#cmd+= " -s %s --upload-fail"%str(submit['id'])
                        #subprocess.call(cmd,shell=True)

                        f = open('/scratch/eparent/eparent/PALFA4/failed_uploads.txt','a')
                        line = str(submit['job_id'])+'\t'+str(submit['id'])+'\t'+outdir_part+'\n'
                        f.write(line)
                        f.close()

			fitsfile = get_fitsfiles(submit)[0]
			base = fitsfile.replace('.fits','')
			tar = base+'.tgz'
			report = base+'.report'
			to_keep = glob.glob(fitsfile)
			to_keep += glob.glob(tar)
			to_keep += glob.glob(report)
			all_files = glob.glob(outdir+'/*')

	    		if config.upload.upload_zerodm_periodicity or config.upload.upload_zerodm_singlepulse:
				outdir_zerodm = os.path.join(outdir, 'zerodm')
				to_keep += glob.glob(outdir_zerodm+'/*zerodm.tgz')
				to_keep += glob.glob(outdir_zerodm+'/*zerodm.report')
				all_files += glob.glob(outdir_zerodm+'/*')

				
			to_del = list(set(all_files)-set(to_keep))
	    		for f in to_del:
				if not f.endswith('/zerodm'):
					os.remove(f)
			if make_fail:
				pipeline_utils.change_job_status_jtdb(submit,'processing_failed')
				faildir_base = "/scratch/eparent/eparent/PALFA4/failed_uploads"
				faildir = os.path.join(faildir_base,outdir_part)
				if not os.path.exists(faildir):
					os.makedirs(faildir)
				for f in to_keep: 
					try:
						shutil.move(f,faildir)
					except:
						continue
				jobtracker.query("Update job_submits set output_dir='%s' where id=%d"%(str(faildir),submit['id']))
				
			#sys.exit()

def get_version_number(dir):
    """Given a directory containing results check to see if there is a file 
        containing the version number. If there is read the version number
        and return it. Otherwise get the current versions, create the file
        and return the version number.

        Input:
            dir: A directory containing results

        Output:
            version_number: The version number for the results contained in 'dir'.
    """
    vernum_fn = os.path.join(dir, "version_number.txt")
    if os.path.exists(vernum_fn):
        f = open(vernum_fn, 'r')
        version_number = f.readline()
        f.close()
    else:
        version_number = config.upload.version_num()
        f = open(vernum_fn, 'w')
        f.write(version_number+'\n')
        f.close()
    return version_number.strip()


def upload_results(job_submit):
    """
    Uploads Results for a given submit.

        Input:
            job_submit: A row from the job_submits table.
                Results from this job submission will be
                uploaded.

        Output:
            None
    """
    print "Attempting to upload results"

    print "\tJob ID: %d, Job submission ID: %d\n\tOutput Dir: %s" % \
            (job_submit['job_id'], job_submit['id'], job_submit['output_dir'])


    if debug.UPLOAD:
        upload.upload_timing_summary = {}
        starttime = time.time()
    try:
        # Connect to the DB
        db = database.Database('default', autocommit=False)

        # Prepare for upload
        dir = job_submit['output_dir']

	# NEW Beluga - Untar the tarball 
	import tarfile 
	to_keep = os.listdir(job_submit['output_dir'])
	tarball = glob.glob(job_submit['output_dir']+'/*00.tgz')[0]
	tar = tarfile.open(tarball,'r:gz')
	tar.extractall(path=job_submit['output_dir'])
	tar.close()

	all_files = os.listdir(job_submit['output_dir'])
	to_del = set(all_files)-set(to_keep)

	if config.upload.upload_zerodm_periodicity or config.upload.upload_zerodm_singlepulse:
		to_keep_zerodm = os.listdir(job_submit['output_dir']+'/zerodm')
		tarball = glob.glob(job_submit['output_dir']+'/zerodm/*zerodm.tgz')[0]
		tar = tarfile.open(tarball,'r:gz')
		tar.extractall(path=job_submit['output_dir']+'/zerodm')
		tar.close()
		all_files_zerodm = os.listdir(job_submit['output_dir']+'/zerodm')
		to_del_zerodm = set(all_files_zerodm)-set(to_keep_zerodm)

        pdm_dir = os.path.join(dir,"zerodm") if config.upload.upload_zerodm_periodicity else dir
        sp_dir = os.path.join(dir,"zerodm") if config.upload.upload_zerodm_singlepulse else dir

        if not os.path.exists(dir) or not os.listdir(dir):
            errormsg = 'ERROR: Results directory, %s, does not exist or is empty for job_id=%d' %\
                       (dir, job_submit['job_id'])
            raise upload.UploadNonFatalError(errormsg)
        elif len(os.listdir(dir)) == 1 and os.listdir(dir)[0] == 'zerodm' \
                                       and not os.listdir(os.path.join(dir,os.listdir(dir)[0])):
            errormsg = 'ERROR: Results directory, %s, does not exist or is empty for job_id=%d' %\
                       (dir, job_submit['job_id'])
            raise upload.UploadNonFatalError(errormsg)

        fitsfiles = get_fitsfiles(job_submit)
        try:
            data = datafile.autogen_dataobj(fitsfiles)
        except ValueError:
            raise upload.UploadNonFatalError
        version_number = get_version_number(dir)

        if debug.UPLOAD: 
            parsetime = time.time()
        # Upload results
        hdr = header.get_header(fitsfiles)
        
        print "\tHeader parsed."

        rat_inst_id_cache = ratings2.utils.RatingInstanceIDCache(dbname='common3')

        cands, tempdir = candidates.get_candidates(version_number, pdm_dir, \
                                                   timestamp_mjd=data.timestamp_mjd, \
                                                   inst_cache=rat_inst_id_cache)
        print "\tPeriodicity candidates parsed. (%d cands)"%len(cands) 
        sp_cands, tempdir_sp = sp_candidates.get_spcandidates(version_number, sp_dir, \
                                                              timestamp_mjd=data.timestamp_mjd, \
                                                              inst_cache=rat_inst_id_cache)
        print "\tSingle pulse candidates parsed. (%d cands)"%len(sp_cands)

        diags = diagnostics.get_diagnostics(data.obs_name, 
                                             data.beam_id, \
                                             data.obstype, \
                                             version_number, \
                                             pdm_dir, sp_dir)
        print "\tDiagnostics parsed."

        for c in (cands + sp_cands):
            hdr.add_dependent(c)
        
        if debug.UPLOAD: 
            upload.upload_timing_summary['Parsing'] = \
                upload.upload_timing_summary.setdefault('Parsing', 0) + \
                (time.time()-parsetime)

        # Perform the upload
        header_id = hdr.upload(db)
	print "Header ID: ",header_id
        for d in diags:
            d.upload(db)
        print "\tDB upload completed and checked successfully. header_id=%d" % \
                    header_id


    except (upload.UploadNonFatalError):
        # Parsing error caught. Job attempt has failed!
        exceptionmsgs = traceback.format_exception(*sys.exc_info())
        errormsg  = "Error while checking results!\n"
        errormsg += "\tJob ID: %d, Job submit ID: %d\n\n" % \
                        (job_submit['job_id'], job_submit['id'])
        errormsg += "".join(exceptionmsgs)
        
        sys.stderr.write("Error while checking results!\n")
        sys.stderr.write("Database transaction will not be committed.\n")
        sys.stderr.write("\t%s" % exceptionmsgs[-1])

        queries = []
        arglists = []
        queries.append("UPDATE job_submits " \
                       "SET status='upload_failed', " \
                            "details=?, " \
                            "updated_at=? " \
                       "WHERE id=?")
        arglists.append((errormsg, jobtracker.nowstr(), job_submit['id']))
        queries.append("UPDATE jobs " \
                       "SET status='failed', " \
                            "details='Error while uploading results', " \
                            "updated_at=? " \
                       "WHERE id=?")
        arglists.append((jobtracker.nowstr(), job_submit['job_id']))
        jobtracker.execute(queries, arglists)
        
        # Rolling back changes. 
        db.rollback()
    except (database.DatabaseConnectionError, ratings2.database.DatabaseConnectionError,\
               CornellFTP.CornellFTPTimeout, upload.UploadDeadlockError,\
               database.DatabaseDeadlockError), e:
        # Connection error while uploading. We will try again later.
        sys.stderr.write(str(e))
        sys.stderr.write("\tRolling back DB transaction and will re-try later.\n")
        
        # Rolling back changes. 
        db.rollback()
    except Exception, e:
        # Unexpected error!
        sys.stderr.write("Unexpected error!\n")
        sys.stderr.write("%s\n" % str(e))
        sys.stderr.write("\tRolling back DB transaction and re-raising.\n")
        
        # Rolling back changes. 
        db.rollback()
        raise
    else:
        # No errors encountered. Commit changes to the DB.
        db.commit()
	print " ..Committed to DB"

        #FTP any FTPables
        attempts = 0
        while attempts < 5:
            try:
                cftp = CornellFTP.CornellFTP()
                hdr.upload_FTP(cftp,db)
                cftp.quit()

            except (CornellFTP.CornellFTPTimeout, CornellFTP.CornellFTPError):
                # Connection error during FTP upload. Reconnect and try again.
                print "FTP connection lost. Reconnecting..."
                exceptionmsgs = traceback.format_exception(*sys.exc_info())
                print"".join(exceptionmsgs)
                attempts += 1
                try:
                    cftp.quit()
                except (EOFError, M2Crypto.SSL.SSLError):
                    pass
            except:
                # Unexpected error
                sys.stderr.write("Unexpected error during FTP upload!\n")
                sys.stderr.write("\tRolling back DB transaction and re-raising.\n")
        
                # Rolling back changes (just last uncommited FTP). 
                db.rollback()
                raise 
            else:
                print "\tFTP upload completed successfully. header_id=%d" % \
                        header_id
                break

        # remove temporary dir for PFDs
        shutil.rmtree(tempdir)
        # remove temporary dir for SPDs
        shutil.rmtree(tempdir_sp)

        if attempts >= 5:
            errmsg = "FTP upload failed after %d connection failures!\n" % attempts
            sys.stderr.write(errmsg)
            sys.stderr.write("\tRolling back DB transaction and raising error.\n")
           
            # Rolling back changes (just last uncommited FTP). 
            db.rollback()
            raise pipeline_utils.PipelineError(errmsg)
                        
        else:
	    # Update database statuses
	    queries = []
	    queries.append("UPDATE job_submits " \
			   "SET status='uploaded', " \
				"details='Upload successful (header_id=%d)', " \
				"updated_at='%s' " \
			   "WHERE id=%d" % 
			   (header_id, jobtracker.nowstr(), job_submit['id']))
	    queries.append("UPDATE jobs " \
			   "SET status='uploaded', " \
				"details='Upload successful (header_id=%d)', " \
				"updated_at='%s' " \
			   "WHERE id=%d" % \
			   (header_id, jobtracker.nowstr(), job_submit['job_id']))
	    jobtracker.query(queries)

	    print "Results successfully uploaded"

	    if config.basic.delete_rawdata:
		pipeline_utils.clean_up(job_submit['job_id'])

	    if debug.UPLOAD: 
		upload.upload_timing_summary['End-to-end'] = \
		    upload.upload_timing_summary.setdefault('End-to-end', 0) + \
		    (time.time()-starttime)
		#print "Upload timing summary:  (total = %.2f s)"%(time.time()-starttime0)
		for k in sorted(upload.upload_timing_summary.keys()):
		    print "    %s: %.2f s" % (k, upload.upload_timing_summary[k])
	    print "" # Just a blank line

	   # NEW Beluga - re-delete files
	    for f in to_del:
		os.remove(job_submit['output_dir']+'/'+f)
	   
	    if config.upload.upload_zerodm_periodicity or config.upload.upload_zerodm_singlepulse:
		for f in to_del_zerodm:
		   	os.remove(job_submit['output_dir']+'/zerodm/'+f)

def get_fitsfiles(job_submit):
    """Find the fits files associated with this job.
        There should be a single file in the job's result
        directory.

        Input:
            job_submit: A row from the job_submits table.
                A list of fits files corresponding to the submit
                are returned.
        Output:
            fitsfiles: list of paths to *.fits files in results
                directory.
    """
    return glob.glob(os.path.join(job_submit['output_dir'], "*.fits"))

