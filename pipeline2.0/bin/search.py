#!/usr/bin/env python
"""
A batch script to search pulsar data.

Patrick Lazarus, May 20, 2010
"""

import sys
import os
import os.path
import socket
import tempfile
import shutil
import subprocess
import tarfile
import warnings

import datafile
import astro_utils.calendar

import config.processing
import pipeline_utils
import getpass

warnings.filterwarnings("ignore", message="Can't find the .* column")
warnings.filterwarnings("ignore", message=".*NSUBOFFS reports 0 previous rows.*")
warnings.filterwarnings("ignore", message="Channel spacing changes in file 0!")

#Cluster users that operates the pipeline
allowed_usr = config.processing.users
current_usr = getpass.getuser()

def warn_to_stdout(message, category, filename, lineno, file=None, line=None):
    """A function to replace warnings.showwarning so that warnings are
        printed to STDOUT instead of STDERR.
        Usage: warnings.showwarning = warn_to_stdout
    """

    sys.stdout.write(warnings.formatwarning(message,category,filename,lineno))

warnings.showwarning = warn_to_stdout

def get_datafns():
    """Get data filenames from command line or environment variable. 
        Environment variable option is only checked if no files
        are provided on command line. 
        
        (NOTE: PBS does not provide batch scripts wtih command line 
                arguments, so we check for if "DATAFILES" environment 
                variable is set.)
    """
    if sys.argv[2:]:
        # First argument is results directory
        # Files provided on command line
        fns = sys.argv[2:]
    else:
        # Files provided with environment variable
        fns = os.getenv("DATAFILES", "").split(';')

    # Ensure all files exist
    # for fn in fns:
    #    if not os.path.exists(fn):
    #        raise ValueError("Data file %s doesn't exist!" % fn)

    # Ensure there are files
    if not fns:
        raise ValueError("No data files provided!")
    return fns


def get_outdir():
    """Get output directory from command line or environment variable.
        Environment variable option is only checked if no files
        are provided on command line. 
        
        (NOTE: PBS does not provide batch scripts wtih command line 
                arguments, so we check for if "OUTDIR" environment 
                variable is set.)
    """
    if sys.argv[1:]:
        # Check command line
        outdir = sys.argv[1]
    else:
        # Use environment variable
        outdir = os.getenv("OUTDIR", "")

    # Ensure output directory is defined
    if not outdir:
        raise ValueError("Output directory is not defined!")
    return outdir
        

def init_workspace():
    """Initialize workspace. 
        - Create working directory.
        - Create results directory.
        - Return 2-tuple (working directory, results directory).
    """
    if config.processing.use_slurm_subdir:
        slurm_job_id = os.getenv("SLURM_JOBID")         
        base_working_dir = os.path.join(config.processing.base_working_directory, \
                                        slurm_job_id) 	
    else:
        base_working_dir = config.processing.base_working_directory

    # Generate temporary working directory
    if not os.path.isdir(base_working_dir):
        print "Creating base work directory..."
        os.makedirs(base_working_dir)
    workdir = tempfile.mkdtemp(suffix="_tmp", prefix="PALFA_processing_", \
                        dir=base_working_dir)
    resultsdir = tempfile.mkdtemp(suffix="_tmp", prefix="PALFA_results_", \
                        dir=base_working_dir)
    return (workdir, resultsdir)


def system_call(cmd):
    retcode = subprocess.call(cmd, shell=True)
    if retcode < 0:
        raise SystemCallError("System call (%s) terminated by signal (%s)!" % \
                                (cmd, -retcode))
    elif retcode > 0:
        raise SystemCallError("System call (%s) failed with status (%s)!" % \
                                (cmd, retcode))
    else:
        # Exit code is 0, which is "Success". Do nothing.
        pass


class SystemCallError(Exception):
    """An exception to throw when a system call returns 
        with a non-zero exit code.
    """
    pass


def set_up():
    print "Running on ", socket.gethostname()
    fns = get_datafns()
    print "Searching %d files:" % len(fns)
    outdir = get_outdir()
    workdir, resultsdir = init_workspace()
   
    print "Local working directory:", workdir
    print "Local results directory:", resultsdir
    print "When finished results will be copied to: %s" % outdir
    
    # Copy data files locally
    for fn in fns:
        system_call("rsync -auvl %s %s" % (fn, workdir))

    fns = sorted([os.path.split(fn)[-1] for fn in fns])

    return fns, workdir, resultsdir, outdir


def search(fns, workdir, resultsdir):
    # Search the data
    print "Go-Go-Gadget pulsar search..."
    import PALFA2_presto_search
    PALFA2_presto_search.main(fns, workdir, resultsdir)
    
    # Remove data, weights, scales and offsets from fits files
    # and stash them in the results directory.
    print "Removing data, weights, scales and offsets."
    for fn in fns:
        system_call("fitsdelcol %s[SUBINT] DATA DAT_WTS DAT_SCL DAT_OFFS" % fn)
        system_call("rsync -auvl %s %s" % (fn, resultsdir))


def copy_zaplist(fns, workdir):
    # Copy zaplist to working directory
    filetype = datafile.get_datafile_type(fns)
    parsed = filetype.fnmatch(fns[0]).groupdict()
    if 'date' not in parsed.keys():
        parsed['date'] = "%04d%02d%02d" % \
                            astro_utils.calendar.MJD_to_date(int(parsed['mjd']))

    customzapfns = []
    # First, try to find a custom zaplist for this specific data file
    customzapfns.append(fns[0].replace(".fits", ".zaplist"))
    # Next, try to find custom zaplist for this beam
    customzapfns.append("%s.%s.b%s.zaplist" % \
                        (parsed['projid'], parsed['date'], parsed['beam']))
    # Try to find custom zaplist for this MJD
    customzapfns.append("%s.%s.all.zaplist" % (parsed['projid'], parsed['date']))

    zaptar_fn = pipeline_utils.get_zaplist_tarball(no_check=True)
    zaptar = tarfile.open(zaptar_fn, mode='r')

    members = zaptar.getmembers()
    for customzapfn in customzapfns:
        radar_samples_fn = os.path.splitext(customzapfn)[0] + '_merged_radar_samples.txt'
        matches = [mem for mem in members \
                    if mem.name.endswith(customzapfn)]
        radar_matches = [mem for mem in members \
                          if mem.name.endswith(radar_samples_fn)]
        if matches:
            ti = matches[0] # The first TarInfo object found 
                            # that matches the file name

            # Write custom zaplist to workdir
            localfn = os.path.join(workdir, customzapfn)
            f = open(localfn, 'w')
            f.write(zaptar.extractfile(ti).read())
            f.close()
            print "Copied custom zaplist: %s" % customzapfn

            if radar_matches:
                radar_ti = radar_matches[0]

                # Write radar samples list to workdir
                radar_localfn = os.path.join(workdir, radar_samples_fn)
                f = open(radar_localfn, 'w')
                f.write(zaptar.extractfile(radar_ti).read())
                f.close()
                print "Copied radar samples list: %s" % radar_samples_fn

            break
        else:
            # The member we searched for doesn't exist, try next one
            pass
    else:
        # Copy default zaplist
        if filetype == datafile.WappPsrfitsData:
            zapfn = config.processing.default_wapp_zaplist
        elif ( filetype == datafile.MockPsrfitsData or
             filetype == datafile.MergedMockPsrfitsData ):
            zapfn = config.processing.default_mock_zaplist
        else:
            raise ValueError("No default zaplist for data files of type %s" % \
                                filetype.__name__)
        shutil.copy(zapfn, workdir)
        print "No custom zaplist found. Copied default zaplist: %s" % \
                zapfn
    
    zaptar.close()


def copy_results(resultsdir, outdir):
    # Copy search results to outdir (only if no errors occurred)
    print "Copying contents of local results directory to", outdir
    system_call("mkdir -m 770 -p %s" % outdir)
    system_call("rsync -auvl --chmod=Dg+rX,Fg+r %s/ %s" % (resultsdir, outdir))

    path1 = os.path.split(outdir)[0]
    path2 = os.path.split(path1)[0]
    path3 = os.path.split(path2)[0]

    owner_newdir = pipeline_utils.get_file_ownership(outdir)
    owner1 = pipeline_utils.get_file_ownership(path1)
    owner2 = pipeline_utils.get_file_ownership(path2)
    owner3 = pipeline_utils.get_file_ownership(path3)

    paths  = [outdir,     path1,  path2,  path3]
    owners = [owner_newdir,owner1, owner2, owner3]
    for i in range(len(paths)):
	if owners[i]==current_usr:
		for user in allowed_usr:
			cmd = "setfacl -m u:%s:rwx %s"%(user,paths[i])
			subprocess.call(cmd,shell=True)
	else: 
		print "Warning: current user is not the owner of %s, can't change permissions"%paths[i]
 


def clean_up(workdir, resultsdir):
    print "Cleaning up..."
    if workdir is not None and os.path.isdir(workdir):
        print "Removing working directory:", workdir
        shutil.rmtree(workdir)
    if resultsdir is not None and os.path.isdir(resultsdir):
        print "Removing local results directory:", resultsdir
        shutil.rmtree(resultsdir)
    

def main():
    workdir = None
    resultsdir = None
    try:
        fns, workdir, resultsdir, outdir = set_up()
        os.chdir(workdir)
        ppfns = [os.path.split(fn)[-1] for fn in datafile.preprocess(fns)]
        copy_zaplist(ppfns, workdir)
        search(ppfns, workdir, resultsdir)
        copy_results(resultsdir, outdir)
    except:
        # Some error was encountered
        sys.stderr.write("\nProcessing errors! Job ran on %s\n\n" % socket.gethostname())
        # Now, simply re-raise the error so it gets reported in the error logs
        raise
    finally:
        # Remove working directory and output directory
        # even if an error occurred
        clean_up(workdir, resultsdir)

if __name__=='__main__':
    main()
