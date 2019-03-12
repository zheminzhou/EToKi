#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""

Convenience methods to make it easier to run external programs
and other os-related tools

sulsj (ssul@lbl.gov)


"""

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## libraries to use

from subprocess import Popen, call, PIPE
import os, glob, sys
import shlex
import unittest
import time
import grp
import errno
from threading import Timer ## for timer

from common import get_run_path

g_scale_inv = ((1024.*1024., "MB"), (1024., "KB"))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function definitions

defineCalledProcessError = False
try:
    from subprocess import CalledProcessError
except ImportError:
    defineCalledProcessError = True

if defineCalledProcessError:
    class CalledProcessError(OSError):
        def __init__(self, returncode, cmd, *l, **kw):
            OSError.__init__(self, *l, **kw)
            self.cmd = cmd
            self.returncode = returncode


"""
Run user command using subprocess.call
@param: popenargs: command and options to run
@param: kwargs: additional parameters
"""
def run(*popenargs, **kwargs):
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop('dryRun', False)

    if dryRun:
        print popenargs
    else:
        ## convert something like run("ls -l") into run("ls -l", shell=True)
        if isinstance(popenargs[0], str) and len(shlex.split(popenargs[0])) > 1:
            kw.setdefault("shell", True)

        ## > /dev/null 2>&1
        if kw.pop("supressAllOutput", False):
            stdnull = open(os.devnull, "w") ## incompat with close_fds on Windows
            kw.setdefault("stdout", stdnull)
            kw.setdefault("stderr", stdnull)
        else:
            stdnull = None

        returncode = call(*popenargs, **kw)
        if stdnull:
            stdnull.close()
        if returncode != 0:
            raise CalledProcessError(returncode=returncode, cmd=str(popenargs))

"""
Similar to shell backticks, e.g. a = `ls -1` <=> a = backticks(['ls','-1']).
If 'dryRun=True' is given as keyword argument, then 'dryRet' keyword must
provide a value to return from this function.
@param: popenargs: command and options to run
@param: kwargs: additional parameters
@return: command result (stdout)
"""
def back_ticks(*popenargs, **kwargs):
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop('dryRun', False)
    dryRet = kw.pop('dryRet', None)

    if dryRun:
        print popenargs
        return dryRet
    else:
        kw['stdout'] = PIPE
        p = Popen(*popenargs, **kw)
        retOut = p.communicate()[0]
        if p.returncode != 0:
            raise CalledProcessError(returncode=p.returncode, cmd=str(popenargs))
        return retOut

"""
Run a command, catch stdout and stderr and exitCode
@param: cmd
@param: live (boolean, default false - don't run the command but pretend we did)

@return: stdout, stderr, exit code
"""
# def run_command(cmd, live=False):
#     return run_sh_command(cmd, live=False)

def run_sh_command(cmd, live=False, log=None, runTime=False, stdoutPrint=True, timeoutSec=0):
    stdOut = None
    stdErr = None
    exitCode = None
    start = 0
    end = 0
    elapsedSec = 0

    if cmd:
        if not live:
            stdOut = "Not live: cmd = '%s'" % (cmd)
            exitCode = 0
        else:
            if log and stdoutPrint:
                log.info("cmd: %s" % (cmd))

            ##---------
            ## OLD
            if runTime:
                start = time.time()

            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            ## ref) http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
            if timeoutSec > 0:
                kill_proc = lambda proc: proc.kill()
                timer = Timer(timeoutSec, kill_proc, [p])

            #p.wait()

            try:
                stdOut, stdErr = p.communicate()
                exitCode = p.returncode
            finally:
                if timeoutSec > 0:
                    timer.cancel()
                    exitCode = 143
                else:
                    pass

            if runTime:
                end = time.time()
                elapsedSec = end - start
                if log:
                    log.info("*************************************")
                    if cmd.split(" ")[0].split("/")[-1]:
                        log.info(cmd.split(" ")[0].split("/")[-1])
                    log.info("Command took " + str(elapsedSec) + " sec.")

                    log.info("*************************************")

            if log and stdoutPrint:
                log.info("Return values: exitCode=" + str(exitCode) + ", stdOut=" + str(stdOut) + ", stdErr=" + str(stdErr))

            if exitCode != 0:
                if log:
                    log.warn("- The exit code has non-zero value.")

    else:
        if log:
            log.error("- No command to run.")
            return None, None, -1


    return stdOut, stdErr, exitCode


"""
Create one dir with pathname path or do nothing if it already exists.
Same as Linux 'mkdir -p'.
@param: path: path
@param: dryRun: dryrun directive
"""
def make_dir(path, perm=None, dryRun=False):
    if not dryRun:
        if not os.path.exists(path):
            if not perm:
                os.makedirs(path)
            else:
                os.makedirs(path, perm)
    else:
        print "make_dir %s" % (path, )

"""
The method make_dir_p() is recursive directory creation function.
Like mkdir(), but makes all intermediate-level directories needed to contain the leaf directory.
"""
def make_dir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

"""
Create muiltiple dirs with the same semantics as make_dir
@param: path: path
@param: dryRun: dryrun directive
"""
def make_dirs(paths, dryRun=False):
    for path in paths:
        make_dir(path=path, dryRun=dryRun)


"""
Assume that the argument is a file name and make all directories that are
part of it
@param: fileName: create dir to the file
"""
def make_file_path(fileName):
    dirName = os.path.dirname(fileName)
    if dirName not in ("", "."):
        make_dir(dirName)

"""
Remove dir
@param: path: path to delete
@param: dryRun: dryrun directive
"""
## To do: perhaps use shutil.rmtree instead?
def rm_dir(path, dryRun=False):
    run(["rm", "-rf", path], dryRun=dryRun)

## make alias
rmrf = rm_dir

"""
Remove file.
@param: path: path to delete
@param: dryRun: dryrun directive
"""
def remove_file(path, dryRun=False):
    for f in glob.iglob(path):
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError:
            pass

"""
Remove multiple files.
@param: path: path to delete
@param: dryRun: dryrun directive
"""
def remove_files(paths, dryRun=False):
    for f in paths:
        try:
            os.remove(f)
        except OSError:
            pass

"""
Create an empty dir with a given path.
If path already exists,  it will be removed first.
@param: path: path to delete
@param: dryRun: dryrun directive
"""
def remake_dir(path, dryRun=False):
    rmrf(path, dryRun=dryRun)
    make_dir(path, dryRun=dryRun)



"""
Change mode.
@param: path: path to chmod
@param: mode: the form `[ugoa]*([-+=]([rwxXst]*|[ugo]))+' OR change_mod(path, "0755")
@param: opts: additional chmod options
@param: dryRun: dryrun directive
"""
def change_mod(path, mode, opts='', dryRun=False):
    if isinstance(path, basestring):
        path = [path]
    else:
        path = list(path)
    run(["chmod"]+opts.split()+[mode]+path, dryRun=dryRun)



"""
Change grp.
"""
def change_grp(filepath, grpName):
    uid = os.stat(filepath).st_uid
    #uid = pwd.getpwnam("qc_user").pw_uid
    gid = grp.getgrnam(grpName).gr_gid
    os.chown(filepath, uid, gid)


"""
find files
"""
def find_files(patt):
    #return [os.path.join(d, f) if f.find(patt) != -1 for f in os.listdir(d)]
    return [f for f in glob.glob(patt)]


"""
Move file
"""
def move_file(source, dest, dryRun=False):
    try:
        if os.path.exists(source) and not os.path.exists(dest):
            run(["mv", source, dest], dryRun=dryRun)
    except OSError:
        pass

"""
get various mem usage properties of process with id pid in MB

@param VmKey
@param pid
"""
#-------------------------------------------------------------------------------
def _VmB(VmKey, pid):
#-------------------------------------------------------------------------------
    procStatus = '/proc/%d/status' % pid
    unitScale = {'kB': 1.0/1024.0, 'mB': 1.0,
                 'KB': 1.0/1024.0, 'MB': 1.0}

    ## get pseudo file /proc/<pid>/status
    try:
        if os.path.exists(procStatus):
            t = open(procStatus)
            v = t.read()
            t.close()
        else:
            return 0.0
    except OSError:
        #logger.exception("Failed to open /proc files.")
        print "Failed to open /proc files."
        return 0.0 # non-Linux?

    ## get VmKey line e.g. 'VmRSS: 9999 kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3) # by whitespace
    if len(v) < 3:
        return 0.0 # invalid format?

    ## convert Vm value to bytes
    return float(v[1]) * unitScale[v[2]]

"""
convert scale
"""
#-------------------------------------------------------------------------------
def to_scale(x):
#-------------------------------------------------------------------------------
    for sc in g_scale_inv:
        y = x/sc[0]
        if y >= 1:
            return "%.3f%s" % (y, sc[1])

    return "%.3f%s" % (y, "B")


"""
Return memory usage in bytes or as formatted string.

@param pid
@param since
@param asStr
"""
#-------------------------------------------------------------------------------
def get_virtual_memory_usage(pid, since=0.0, asStr=True):
#-------------------------------------------------------------------------------
    b = _VmB('VmSize:', pid) - since
    if asStr:
        return "VirtMem: " + to_scale(b)
    else:
        return b


"""
Return resident memory usage in bytes.

@param pid
@param since
@param asStr
"""
#-------------------------------------------------------------------------------
def get_resident_memory_usage(pid, since=0.0, asStr=True):
#-------------------------------------------------------------------------------
    b = _VmB('VmRSS:', pid) - since
    if asStr:
        return "ResMem: " + to_scale(b)
    else:
        return b


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## unit test
class TestOsUtility(unittest.TestCase):
    def testRun(self):
        try:
            run(["rm", "-rf", "./unittest"], dryRun=False)
        except CalledProcessError, msg:
            self.assertNotEqual(msg.returncode, 0)
        try:
            make_dir("./unittest", dryRun=False)
        except CalledProcessError, msg:
            self.assertEqual(msg.returncode, 0)
        try:
            rm_dir("./unittest", dryRun=False)
        except CalledProcessError, msg:
            self.assertEqual(msg.returncode, 0)

    def testBackTicks(self):
        cmd = "free"
        try:
            freeOut = back_ticks(cmd, shell=True)
        except CalledProcessError, msg:
            print >>sys.stderr, "Failed to call %s. Exit code=%s" % (msg.cmd, msg.returncode)
            sys.exit(1)

        ret = -1
        ret = float(freeOut.split('\n')[1].split()[2]) / \
              float(freeOut.split('\n')[1].split()[1]) * 100.0
        assert ret > -1



def sh(cmd):
    """ simple function to call a sh command """
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc.communicate()
    if err: print >>sys.stderr, err
    return out



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## main program


## EOF
