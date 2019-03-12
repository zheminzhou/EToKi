#!/usr/bin/env python

import os
import sys
import time # sleep
import re # pod_path

srcDir = os.path.dirname(__file__)
sys.path.append(os.path.join(srcDir, '../assemblyqc/lib'))

# from assemblyqc_constants import *
from common import *
from rqc_constants import *
from os_utility import run_sh_command, make_dir_p
from html_utility import html_link


'''
Run Stephan Trong's blast plus wrapper instead of megablast.
It specifies default parameters for megablast and megablast is the default aligner.
The dbs for which run_blastplus is called do not have masking and are not really indexed with makembindex.
TODO: Run all databases in parallel.
      Use the same blastn parameters as used in run_blastn, which are set as follows
      defaults = " -num_threads " + str(num_threads) + " " + RQCBlast.BLASTN_DEFAULTS;
'''
# def run_blastplus(queryFastaFile, db, outputPath, log, num_threads=16):
# REMOVED! 09092016


def get_cat_cmd(seqUnitName, log):
    zcatCmd = ""

    if seqUnitName.endswith(".bz2"):
        zcatCmd = RQCCommands.BZCAT_CMD
    elif seqUnitName.endswith(".gz"):
        zcatCmd = RQCCommands.ZCAT_CMD
    elif seqUnitName.endswith(".fastq") or seqUnitName.endswith(".fa"):
        zcatCmd = RQCCommands.CAT_CMD
    else:
        log.error("source_fastq should be either bzipped or gzipped fastq. " + str(seqUnitName))
        return zcatCmd, RQCExitCodes.JGI_FAILURE


    return zcatCmd, RQCExitCodes.JGI_SUCCESS



"""
Localize specified file to /scratch/rqc
- Bryce added the file check and wait

@param: fileName: full path to file
@param: log

sulsj (ssul@lbl.gov)
"""
def localize_file(fileNameFullPath, log):
    done = 0
    loopCount = 0

    while done == 0:
        loopCount += 1

        if os.path.isfile(fileNameFullPath):
            done = 1 # congratulations!  The file system seems to be working for the moment
        else:
            sleepTime = 60 * loopCount
            log.error("Filename doesn't exist on the system: %s, sleeping for %s seconds.", fileNameFullPath, sleepTime)
            time.sleep(sleepTime)

        if loopCount > 3:
            done = 1 # probably failed

    if not os.path.isfile(fileNameFullPath):
        return None

    destinationPath = "/scratch/rqc/localized-file"
    if not os.path.isdir(destinationPath):
        # make_dir_p(destinationPath)
        _, _, exitCode = run_sh_command("mkdir -p /scratch/rqc/localized-file && chmod 777 /scratch/rqc/localized-file", True, log, True)
        assert exitCode == 0

    fileName = safe_basename(fileNameFullPath, log)[0]
    localizeCmd = "rsync -av --omit-dir-times --no-perms " + fileNameFullPath + " " + os.path.join(destinationPath, fileName)

    _, _, exitCode = run_sh_command(localizeCmd, True, log, True)

    localizedFileNameFullPath = None

    if exitCode == RQCExitCodes.JGI_SUCCESS:
        localizedFileNameFullPath = os.path.join(destinationPath, fileName)
        log.info("File localization completed for %s" % (fileName))
    else:
        localizedFileNameFullPath = None
        log.error("File localization failed for %s" % (fileName))


    return localizedFileNameFullPath


# def localize_dir(dbFileOrDirName, log):
# REMOVED! 09092016


"""
New localize_dir

Revisions
  08022017 Enabled burst buffer on cori
  08072017 Added support for persistent bbuffer on cori for localized NCBI databases

"""
def localize_dir2(db, log, bbuffer=None):
    safeBaseName, exitCode = safe_basename(db, log)
    safeDirName, exitCode = safe_dirname(db, log)
    targetScratchDir = None

    nerscDb = "/scratch/blastdb/global/dna/shared/rqc/ref_databases/ncbi/CURRENT"
    # nerscDb = "/" ## temporarily do not use NERSC db
    ## nerscDbPersistentBBuffer = "/var/opt/cray/dws/mounts/batch/NCBI_DB_striped_scratch/ncbi"

    ## check if persistent burst buffer is ready
    if bbuffer is None and "DW_PERSISTENT_STRIPED_NCBI_DB" in os.environ and os.environ['DW_PERSISTENT_STRIPED_NCBI_DB'] is not None:
        nerscDb = os.path.join(os.environ['DW_PERSISTENT_STRIPED_NCBI_DB'], "ncbi")

    else:
        targetScratchDir = "/scratch/rqc"

        if bbuffer is not None:
            if not os.path.isdir(bbuffer):
                log.error("Burst Buffer does not initiated: %s", bbuffer)
                return None, RQCExitCodes.JGI_FAILURE
            else:
                targetScratchDir = bbuffer
                log.info("Localization will use Burst Buffer location at %s", targetScratchDir)

        elif not os.path.isdir(targetScratchDir):
            _, _, exitCode = run_sh_command("mkdir %s && chmod 777 %s" % (targetScratchDir, targetScratchDir), True, log, True)
            assert exitCode == 0

    rsyncOption = ""
    src = db
    dest = ""

    ## For
    # GREEN_GENES = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/green_genes16s.insa_gg16S.fasta"
    # LSU_REF = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/LSURef_115_tax_silva.fasta"
    # SSU_REF = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/SSURef_NR99_115_tax_silva.fasta"
    # LSSU_REF = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/LSSURef_115_tax_silva.fasta"
    # CONTAMINANTS = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/JGIContaminants.fa"
    # COLLAB16S = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/collab16s.fa"
    if db.endswith(".fa") or db.endswith(".fasta"):
        rsyncOption = "--include '*.n??' --exclude '*.fa' --exclude '*.fasta' --exclude '%s' --exclude '*.log'" % (safeBaseName)
        src = db + ".n??"
        dest = os.path.join(targetScratchDir, safeBaseName)
        blastDb = dest
        if os.path.isfile(dest): ## just in case for the file name already exists
            rmCmd = "rm -rf %s" % (dest)
            run_sh_command(rmCmd, True, log, True)

    ## For
    # NT_maskedYindexedN_BB = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/nt/bbtools_dedupe_mask/nt_bbdedupe_bbmasked_formatted"
    elif db.endswith("nt_bbdedupe_bbmasked_formatted"):
        if os.path.isdir(os.path.join(nerscDb, "nt/bbtools_dedupe_mask")):
            blastDb = os.path.join(nerscDb, "nt/bbtools_dedupe_mask")
            log.info("NERSC NCBI Database found: %s" % (blastDb))
            return blastDb, RQCExitCodes.JGI_SUCCESS

        rsyncOption = "--include '*.n??' --exclude '*.fna' --exclude 'nt_bbdedupe_bbmasked_formatted' --exclude '*.log'"
        src = safeDirName + '/'
        dest = os.path.join(targetScratchDir, safeBaseName)
        blastDb = dest

    ## For
    # NR                   = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/nr/nr"
    # REFSEQ_ARCHAEA       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.archaea/refseq.archaea"
    # REFSEQ_BACTERIA      = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.bacteria/refseq.bacteria"
    # REFSEQ_FUNGI         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.fungi/refseq.fungi"
    # REFSEQ_MITOCHONDRION = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.mitochondrion/refseq.mitochondrion"
    # REFSEQ_PLANT         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plant/refseq.plant"
    # REFSEQ_PLASMID       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plasmid/refseq.plasmid"
    # REFSEQ_PLASTID       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plastid/refseq.plastid"
    # REFSEQ_VIRAL         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.viral/refseq.viral"
    else:
        if os.path.isdir(os.path.join(nerscDb, safeBaseName)):
            blastDb = os.path.join(nerscDb, safeBaseName)
            log.info("NERSC NCBI Database found: %s" % (blastDb))
            return blastDb, RQCExitCodes.JGI_SUCCESS

        rsyncOption = "--include '*.n??' --exclude '*.fna' --exclude '%s' --exclude '*.log'" % (safeBaseName)
        src = safeDirName + '/'
        dest = os.path.join(targetScratchDir, safeBaseName)
        blastDb = os.path.join(dest, dest)


    rsyncCmd = "rsync -av --omit-dir-times --no-perms %s %s %s" % (rsyncOption, src, dest)
    dbDirNameLocalized, _, exitCode = run_sh_command(rsyncCmd, True, log, True)

    if exitCode != 0:
        log.error("rsync failed. Cannot localize " + str(db))
        return dbDirNameLocalized, RQCExitCodes.JGI_FAILURE

    else:
        cmd = "chmod -f -R 777 %s" % (dest)
        _, _, exitCode = run_sh_command(cmd, True, log)


    return blastDb, RQCExitCodes.JGI_SUCCESS


def safe_dirname(pathName, log):
    dirName = ""

    if pathName is None:
        return dirName, RQCExitCodes.JGI_FAILURE

    dirName = os.path.dirname(pathName)

    if dirName is None:
        log.error("Could not get basename for " + str(pathName))
        return dirName, RQCExitCodes.JGI_FAILURE

    return dirName, RQCExitCodes.JGI_SUCCESS


def safe_basename(pathName, log):
    baseName = ""

    if pathName is None:
        return baseName, RQCExitCodes.JGI_FAILURE

    baseName = os.path.basename(pathName)

    if baseName is None:
        log.error("Could not get basename for " + str(pathName))
        return baseName, RQCExitCodes.JGI_FAILURE

    return baseName, RQCExitCodes.JGI_SUCCESS


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## qsub functions

## module load uge
def submit_qsub(qsub_opt, qsub_cmd, qsub_name, output_path, log):
    job_id = 0
    job_share = 150

    cluster_project = "gentech-rqc.p"
    qsub_path = os.path.join(output_path, "qsub")

    if not os.path.isdir(qsub_path):
        os.makedirs(qsub_path)

    output_log = os.path.join(qsub_path, "qsub-%s.log" % qsub_name)

    my_qsub = "qsub -b y -j y -m n -w e -terse -N %s -P %s -o %s -js %s %s '%s'" % (qsub_name, cluster_project, output_log, job_share, qsub_opt, qsub_cmd)
    # append to qsub.txt

    cmd = "module load uge;%s" % my_qsub
    stdOut, stdErr, exitCode = run_sh_command(cmd, True, log)
    post_mortem_cmd(cmd, exitCode, stdOut, stdErr, log)

    if exitCode == 0:
        job_id = int(stdOut.strip())

    log.info("- cluster job id: %s", job_id)

    if job_id > 0:
        qsub_log = os.path.join(output_path, "qsub_list.txt")
        fh = open(qsub_log, "a")
        fh.write("%s,%s\n" % (job_id, "submitted"))
        fh.close()

    return job_id


## watch jobs running on the cluster
def watch_cluster(output_path, log, sleep_time=300):
    log.info("watch_cluster")

    is_job_complete = "/usr/common/usg/bin/isjobcomplete.new" # + job_id = nice way to ask if job is done
    qsub_log = os.path.join(output_path, "qsub_list.txt")

    #sleep_time = 300 # check isjobrunning every xx seconds

    hb_max = (180 * 3600) / sleep_time # total number of heartbeats before we give up, 180 hours worth of run time
    #hb_max = 5 # temp

    done = 0
    hb_cnt = 0 # heartbeat count

    if not os.path.isfile(qsub_log):
        done = 1

    while done == 0:

        hb_cnt += 1
        log.info("- heartbeat: %s", hb_cnt)

        qsub_list = []


        qsub_cnt = 0
        qsub_complete = 0
        qsub_err = 0
        fh = open(qsub_log, "r")
        for line in fh:

            qsub_list.append(line.strip())
            qsub_cnt += 1
        fh.close()

        fh = open(qsub_log, "w")


        for qsub in qsub_list:

            job_id, status = qsub.split(",")
            new_status = status

            if status in ("complete", "fail"):
                continue
            else:
                cmd = "%s %s" % (is_job_complete, job_id)
                stdOut, stdErr, exitCode = run_sh_command(cmd, True, log)
                #post_mortem_cmd(cmd, exitCode, stdOut, stdErr, log)

                running = "%s queued/running" % job_id
                not_running = "%s not queued/running" % job_id
                error_qw = "%s queued/running/error" % job_id

                if stdOut.strip == running:
                    new_status = "running"

                elif stdOut.strip() == not_running:
                    new_status = "complete" # might have failed
                    qsub_complete += 1

                elif stdOut.strip() == error_qw:
                    new_status = "error"
                    qsub_err += 1

                if stdOut.strip() == "not queued/running":
                    new_status = "complete"

            fh.write("%s,%s\n" % (job_id, new_status))
            log.info("- job_id: %s, status: %s", job_id, new_status)

        fh.close()

        qsub_running = qsub_cnt - (qsub_complete + qsub_err)
        log.info("- job count: %s, running: %s, err: %s", qsub_cnt, qsub_running, qsub_err)

        if qsub_cnt == (qsub_complete + qsub_err):
            done = 1

        if hb_cnt >= hb_max:
            done = 1

        if done == 0:
            time.sleep(sleep_time)


'''
Open the file, count number of lines/hits
- skip header (#)

Used in sag.py and sag_decontam.py, ce.py and iso.py too
* out of memory when reading big files ...
'''
def get_blast_hit_count(blast_file):
    hit_count = 0

    if os.path.isfile(blast_file):
        # if more than a gig we get a memory error doing this, even though its not supposed to cause a problem ...
        # - maybe the hit_count var grows too big? no, not counting the hits still causes issues
        if os.path.getsize(blast_file) < 100000000:
            fh_blast = open(blast_file, "r")
            for line in fh_blast:
                if line.startswith("#"):
                    continue
                hit_count += 1
            fh_blast.close()
        else:
            hit_count = 9999
    return hit_count


'''
Pads a string with 0's and splits into groups of 2
e.g. padStringPath(44) returns "00/00/00/44"

@param myString: string to pad
@return: padded string
'''

def pad_string_path(myString, padLength=8, depth=None):
    myString = str(myString)
    padLength = int(padLength)

    if padLength > 8 or padLength <= 0:
        padLength = 8

    # left-pad with 0's
    myString = myString.zfill(padLength)

    # use re.findall function to split into pairs of strings
    stringList = re.findall("..", myString)

    # create ss/ss/ss/ss
    if not depth:
        padString = "/".join(stringList)
    else:
        padString = "/".join(stringList[:depth])

    padString = padString + "/"

    return padString

def get_dict_obj(sfile):
    kvmap = {}
    if os.path.isfile(sfile):
        with open(sfile, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                tks = line.strip().split('=')
                if len(tks) == 2:
                    kvmap[tks[0].strip().lower()] = tks[1].strip()
    return kvmap

# helper function
def pipeline_val(keyname, kwargs, stats, files=None, file_path_prefix=None):
    keyname = keyname.lower()
    form = kwargs.get('type')
    if form == 'file' and files:
        val = files.get(keyname)
        if val:
            if kwargs.get('filter') == 'base':
                val = os.path.basename(val)
            elif kwargs.get('filter') != 'full' and file_path_prefix and val.startswith(file_path_prefix):    # make the file path relative to run dir
                val = val[len(file_path_prefix)+1:]

            if kwargs.get('filter') == 'link':
                val = html_link(val, kwargs.get('label'))
    elif stats:
        val = stats.get(keyname, 'n.a.')
        if val != 'n.a.':
            if val.endswith('%'):   # remove last "%"
                val = val[:-1]
            if form == 'bigint' or form == 'int':
                try:
                    tmp = int(val)
                    if form == 'bigint':
                        val = format(tmp, ',')  # e.g. 123,456,789
                    else:
                        val = tmp # int
                except:
                    pass
            elif form == 'pct':
                try:
                    tmp = 100.0 * float(val)
                    if 'filter' in kwargs and type(kwargs['filter']) == int:
                        fformat = '%.' + str(kwargs['filter']) + 'f'
                    else:
                        fformat = '%.2f'
                    
                    val = fformat % tmp
                except Exception as e:
                    print('ERROR pipeline_val - form=pct: %s' % e)
                    
            elif form == 'floatstr' or form == 'float':
                try:
                    tmp = float(val)
                    if form == 'floatstr':
                        val = '%.4f' % tmp
                    elif form == 'float':
                        val = tmp # float
                        if 'filter' in kwargs and type(kwargs['filter']) == int:
                            fformat = '%.' + str(kwargs['filter']) + 'f'
                            val = fformat % val
                except Exception as e:
                    print('ERROR pipeline_val - form=floatstr|float: %s' % e)
                    
            
        if val == 'n.a.' and kwargs.get('filter') == 0:
            if form == 'int':
                val = 0
            else:
                val = '0'
    else:
        val = None
    
    if 'vtype' not in kwargs or kwargs['vtype'] != 'numeric':
        val = str(val)
        
    return val # this is for string replacement, so return string

## EOF

if __name__ == '__main__':
    print('unit test ...')
    tok_map = {
            'overall bases Q score mean' : {'token' : '[_BASE-QUALITY-SCORE_]', 'type': 'float', 'filter': 1},
            'overall bases Q score std' : {'token' : '[_BASE-QUALITY-SCORE-STD_]', 'type': 'float'},
            'Q30 bases Q score mean' : {'token' : '[_Q30-BASE-QUALITY-SCORE_]', 'type': 'float'},
            'Q30 bases Q score std' : {'token' : '[_Q30-BASE-QUALITY-SCORE-STD_]', 'type': 'float'},
    }
    
    odir = '/global/projectb/scratch/syao/standalone/readqc_bb'
    stats = get_dict_obj(os.path.join(odir, 'readqc_stats.txt'))
    # rawCnt = pipeline_val('inputReads', {'type': 'int'}, stats)
    # print(rawCnt)
    # print(type(rawCnt))
    # from pprint import pprint
    # pprint(stats)
 
    for key in tok_map:
        dat = tok_map[key]
        key = key.lower()
        print(key)
        print(pipeline_val(key, dat, stats))
        print(stats[key])
        print('\n')
        
        
        