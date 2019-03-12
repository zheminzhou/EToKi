#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Function definitions common to all programs.


"""

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## libraries to use

#import re
import os
import time
import sys
#import getpass
import logging
#from colorlog import ColoredFormatter
# import EnvironmentModules # get_read_count_fastq
from subprocess import Popen, PIPE
from email.mime.text import MIMEText


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function definitions


'''
creates a logging instance
https://docs.python.org/2/howto/logging.html
https://pypi.python.org/pypi/colorlog
'''
def get_logger(log_name, log_file, log_level = "INFO", stdout = False, color = False):
    log = logging.getLogger(log_name)
    handler = None
    if stdout:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(log_file)

    formatter = logging.Formatter('%(filename)-15s:%(process)d %(asctime)s %(levelname)s: %(message)s')

    if color and 1==2:
        """
        formatter = ColoredFormatter("%(filename)-15s:%(process)d %(asctime)s %(log_color)s%(levelname)s: %(message)s", datefmt=None, reset=True,
                                log_colors={
                                    'DEBUG': 'blue',
                                    'INFO': 'green',
                                    'WARNING': 'yellow',
                                    'ERROR': 'red',
                                    'CRITICAL': 'red, bg_white',
                                },
                                secondary_log_colors={},
                                style='%')
        Not working in conda - 2017-04-29
        """
    handler.setFormatter(formatter)

    log.addHandler(handler)
    log.setLevel(log_level)

    return log


'''
Checkpoint the status plus a timestamp
- appends the status

@param status_log: /path/to/status.log (or whatever you name it)
@param status: status to append to status.log

'''
def checkpoint_step(status_log, status):
    status_line = "%s,%s\n" % (status, time.strftime("%Y-%m-%d %H:%M:%S"))

    with open(status_log, "a") as myfile:
        myfile.write(status_line)

'''
returns the last step (status) from the pipeline
@param status_log: /path/to/status.log (or whatever you name it)
@param log: logger object

@return last status in the status log, "start" if nothing there
'''
def get_status(status_log, log = None):
    #status_log = "%s/%s" % (output_path, "test_status.log")

    status = "start"
    timestamp = str(time.strftime("%Y-%m-%d %H:%M:%S"))

    if os.path.isfile(status_log):
        fh = open(status_log, 'r')
        lines = fh.readlines()
        fh.close()

        for line in lines:
            if line.startswith('#'): continue
            line_list = line.split(",")
            assert len(line_list) == 2
            status = str(line_list[0]).strip()
            timestamp = str(line_list[1]).strip()

        if not status:
            status = "start"

        if log:
            log.info("Last checkpointed step: %s (%s)", status, timestamp)

    else:
        if log:
            log.info("Cannot find status.log (%s), assuming new run", status_log)

    status = status.strip().lower()

    return status



'''
run a command from python

@param cmd: command to run
@param live: False = run in dry mode (print command), True = run normally
@param log: logger object

@return std_out, std_err, exit_code
'''
def run_command(cmd, live=False, log=None):
    stdOut = None
    stdErr = None
    exitCode = None
    #end = 0
    #elapsedSec = 0

    if cmd:
        if not live:
            stdOut = "Not live: cmd = '%s'" % (cmd)
            exitCode = 0
        else:
            if log: log.info("cmd: %s" % (cmd))

            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            stdOut, stdErr = p.communicate()
            exitCode = p.returncode

        if log:
            log.info("Return values: exitCode=" + str(exitCode) + ", stdOut=" + str(stdOut) + ", stdErr=" + str(stdErr))
            if exitCode != 0:
                log.warn("- The exit code has non-zero value.")

    else:
        if log:
            log.error("- No command to run.")
            return None, None, -1

    return stdOut, stdErr, exitCode


'''
replacement for run_command
- includes logging, convert_cmd & post_mortem
'''

def run_cmd(cmd, log=None):

    std_out = None
    std_err = None
    exit_code = 0

    if cmd:
        # convert to work on genepool/denovo
        cmd = convert_cmd(cmd)

        if log:
            log.info("- cmd: %s", cmd)

        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        std_out, std_err = p.communicate()
        exit_code = p.returncode

        post_mortem_cmd(cmd, exit_code, std_out, std_err, log)

    return std_out, std_err, exit_code


'''
Simple function to output to the log what happened only if exit code > 0

Typical usage:
    std_out, std_err, exit_code = run_command(cmd, True)
    post_mortem_cmd(cmd, exit_code, std_out, std_err)

'''
def post_mortem_cmd(cmd, exit_code, std_out, std_err, log = None):
    if exit_code > 0:
        if log:
            log.error("- cmd failed: %s", cmd)
            log.error("- exit code: %s", exit_code)

        else:
            print "- cmd failed: %s" % (cmd)
            print "- exit code: %s" % (exit_code)


        if std_out:
            if log:
                log.error("- std_out: %s", std_out)
            else:
                print "- std_out: %s" % (std_out)

        if std_err:
            if log:
                log.error("- std_err: %s", std_err)
            else:
                print "- std_err: %s" % (std_err)


'''
Convert command to use genepool or denovo (shifter) to run
replace #placeholder; with shifter or module load command
#placeholder.v; should specify the version to use

This should be the only place in the pipelines that specifies the images/modules translation
'''
def convert_cmd(cmd):
    new_cmd = cmd

    shifter_img = {
        "#bbtools" : "shifter --image=bryce911/bbtools ",
        "#pigz" : "module load pigz;",
        "#jamo" : "shifter --image=registry.services.nersc.gov/htandra/jamo_dev:1.0 ", # works, but would like simple module to use - have one on Denovo but not Cori
        "#gnuplot" : "shifter --image=bryce911/bbtools ", # (1)
        "#spades/3.9.0" : "shifter --image=bryce911/spades3.9.0 ",
        "#spades/3.10.1" : "shifter --image=bryce911/spades3.10.1 ",
        "#spades/3.11.0" : "shifter --image=bryce911/spades-3.11.0 ", # GAA-3383
        "#spades/3.11.1-check" : "shifter --image=bryce911/spades3.11.1-check ", # development
        "#prodigal/2.6.3" : "shifter --image=registry.services.nersc.gov/jgi/prodigal ", # RQCSUPPORT-1318
        "#prodigal/2.5.0" : "shifter --image=registry.services.nersc.gov/jgi/prodigal ",
        "#prodigal/2.50" : "shifter --image=registry.services.nersc.gov/jgi/prodigal ", 
        "#lastal/869" : "shifter --image=bryce911/lastal:869 ",
        "#lastal/828" : "shifter --image=bryce911/lastal:869 ",
        #"#lastal" : "shifter --image=bryce911/lastal:869 ",
        "#R/3.3.2" : "module load R/3.3.2;",
        "#texlive" : "shifter --image=bryce911/bbtools ", # (1)
        "#java" : "shifter --image=bryce911/bbtools ", # (1)
        "#blast+/2.6.0" : "shifter --image=sulsj/ncbi-blastplus:2.6.0 ",
        "#blast" : "shifter --image=sulsj/ncbi-blastplus:2.7.0 ",
        "#megahit-1.1.1" : "shifter --image=foster505/megahit:v1.1.1-2-g02102e1 ",
        "#smrtanalysis/2.3.0_p5" : "shifter --image=registry.services.nersc.gov/jgi/smrtanalysis:2.3.0_p5 ", # meth - need more memory
        "#mummer/3.23" : "shifter --image=bryce911/mummer3.23 ", # 3.23
        "#hmmer" : "shifter --image=registry.services.nersc.gov/jgi/hmmer:latest ", # 3.1b2
        "#samtools/1.4" : "shifter --image=rmonti/samtools ",
        "#mothur/1.39.5" : "shifter --image=bryce911/mothur1.39.5 ",
        "#vsearch/2.4.3" : "shifter --image=bryce911/vsearch2.4.3 ",
        "#graphviz" : "shifter --image=bryce911/bbtools ",
        "#ssu-align/0.1.1" : "shifter --image=bryce911/ssu-align0.1.1 ", # openmpi/1.10 included in docker container
        "#smrtlink/4.0.0.190159" : "shifter --image=registry.services.nersc.gov/jgi/smrtlink:4.0.0.190159 /smrtlink/smrtcmds/bin/", # progs not in path
        "#smrtlink/5.0.1.9585" : "shifter --image=registry.services.nersc.gov/jgi/smrtlink:5.0.1.9585 /smrtlink/smrtcmds/bin/", # progs not in path, Tony created 2017-10-16
        "#smrtlink" : "shifter --image=registry.services.nersc.gov/jgi/smrtlink:5.0.1.9585 /smrtlink/smrtcmds/bin/", # progs not in path
        "#prodege" : "shifter --image=bryce911/prodege ", # 2.2.1
        #"#hmmer" : "shifter --image=registry.services.nersc.gov/jgi/hmmer ", # 3.1b2 - Feb 2015, latest as of Oct 2017
        "#checkm" : "shifter --image=registry.services.nersc.gov/jgi/checkm ",
    }
    

    # (1) -  added as part of the bryce911 bbtools package

    #cmd = "#bbtools-shijie;bbmap...."
    # this dict will be deprecated as of March 2018 when genepool passes into legend
    genepool_mod = {
        "#bbtools" : "module load bbtools",
        "#pigz" : "module load pigz",
        "#jamo" : "module load jamo",
        "#gnuplot" : "module load gnuplot/4.6.2", # sag,iso,sps,ce:gc_cov, gc_histogram, contig_gc
        "#spades/3.9.0" : "module load spades/3.9.0",
        "#spades/3.10.1" : "module load spades/3.10.1",
        "#spades/3.11.1" : "module load spades/3.11.1-check",
        "#prodigal/2.6.3" : "module load prodigal/2.50", # aka 2.50, also 2.60 is available
        "#prodigal/2.5.0" : "module load prodigal/2.50",
        "#prodigal/2.50" : "module load prodigal/2.50",
        #"#lastal" : "module load last/828",
        "#lastal/828" : "module load last/828",
        "#R/3.3.2" : "module unload R;module load R/3.3.1", # 3.3.2 not on genepool - RQCSUPPORT-1516 unload R for Kecia
        "#texlive" : "module load texlive",
        "#blast+/2.6.0" : "module load blast+/2.6.0",
        #"#blast+/2.7.0" : "module load blast+/2.7.0", # not created
        "#blast" : "module load blast+/2.6.0",
        "#java" : "", # java runs natively on genepool
        "#megahit-1.1.1" : "module load megahit/1.1.1",
        "#smrtanalysis/2.3.0_p5" : "module load smrtanalysis/2.3.0_p5",
        "#smrtanalysis/2.3.0_p5_xmx32g" : "module load smrtanalysis/2.3.0_p5;export _JAVA_OPTIONS='-Xmx32g'",
        "#mummer/3.23" : "module load mummer/3.23",
        "#hmmer" : "module load hmmer/3.1b2",
        "#samtools/1.4" : "module load samtools/1.4",
        "#mothur/1.39.5" : "module load mothur/1.32.1", # 1.26.0 default, 1.32.1
        "#vsearch/2.4.3" : "module load vsearch/2.3.0", # 2.3.0
        "#graphviz" : "module load graphviz",
        "#ssu-align/0.1.1" : "module load ssu-align",
        "#smrtlink/4.0.0.190159" : "module load smrtlink/4.0.0.190159",
        "#smrtlink" : "module load smrtlink/5.0.1.9585",
        "#smrtlink/5.0.1.9585" : "module load smrtlink/5.0.1.9585",
        "#prodege" : "module load R;/projectb/sandbox/rqc/prod/pipelines/external_tools/sag_decontam/prodege-2.2/bin/",
        "#checkm" : "module load hmmer prodigal pplacer", # checkm installed in python by default on genepool
    }

    #bbtools;stats.sh
    
    if cmd.startswith("#"):

        cluster = "genepool"
        # any other env ids to use?

        # cori, denovo, genepool
        cluster = os.environ.get('NERSC_HOST', 'unknown')


        f = cmd.find(";")
        mod = "" # command to replace
        if f > -1:
            mod = cmd[0:f]

        if mod:
            
            
            # use module load jamo on denovo
            if mod == "#jamo" and cluster == "denovo":
                shifter_img[mod] = "module load jamo;"
                
            if cluster in ("denovo", "cori"):
                
                if mod in shifter_img:
                    new_cmd = new_cmd.replace(mod + ";", shifter_img[mod])

            else:
                if mod in genepool_mod:
                    if genepool_mod[mod] == "":
                        new_cmd = new_cmd.replace(mod + ";", "")
                    else:
                        new_cmd = new_cmd.replace(mod, genepool_mod[mod])

    if new_cmd.startswith("#"):
        print "Command not found!  %s" % new_cmd
        sys.exit(18)

    #print new_cmd
    return new_cmd


'''
returns human readable file size
@param num = file size (e.g. 1000)

@return: readable float e.g. 1.5 KB
'''
def human_size(num):
    if not num:
        num = 0.0

    for x in ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'XB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

    return "%3.1f %s" % (num, 'ZB')



'''
send out email
@param emailTo: email receipient (e.g. bryce@lbl.gov)
@param emailSubject: subject line for the email
@param emailBody: content of the email
@param emailFrom: optional email from

'''
def send_email(email_to, email_subject, email_body, email_from = 'rqc@jgi-psf.org', log = None):
    msg = ""
    err_flag = 0

    if not email_to:
        msg = "- send_email: email_to parameter missing!"

    if not email_subject:
        msg = "- send_email: email_subject parameter missing!"

    if not email_body:
        msg = "- send_email: email_body parameter missing!"


    if err_flag == 0:
        msg = "- sending email to: %s" % (email_to)

    if log:
        log.info(msg)
    else:
        print msg

    if err_flag == 1:
        return 0

    # assume html
    email_msg = MIMEText(email_body, "html") # vs "plain"
    email_msg['Subject'] = email_subject
    email_msg['From'] = email_from
    email_msg['To'] = email_to

    p = Popen(["/usr/sbin/sendmail", "-t"], stdin = PIPE)
    p.communicate(email_msg.as_string())

    return err_flag

'''
Write to rqc_file (e.g. rqc-files.tmp) the file_key and file_value

@param rqc_file_log: full path to file containing key=file
@param file_key: key for the entry
@param file_value: value for the entry

'''
def append_rqc_file(rqc_file_log, file_key, file_value, log=None):
    if file_key:

        buffer = "%s = %s\n" % (file_key, file_value)

        with open(rqc_file_log, "a") as myfile:
            myfile.write(buffer)

        if log: log.info("append_rqc_file: %s:%s" % (file_key, file_value))

    else:

        if log: log.warning("key or value error: %s:%s" % (file_key, file_value))

'''
Write to rqc_stats (e.g. rqc-stats.tmp) the stats_key and stats_value
@param rqc_file_log: full path to file containing key=file
@param file_key: key for the entry
@param file_value: value for the entry

'''
def append_rqc_stats(rqc_stats_log, stats_key, stats_value, log=None):
    if stats_key:

        buffer = "%s = %s\n" % (stats_key, stats_value)

        with open(rqc_stats_log, "a") as myfile:
            myfile.write(buffer)

        if log: log.info("append_rqc_stats: %s:%s" % (stats_key, stats_value))

    else:

        if log: log.warning("key or value error: %s:%s" % (stats_key, stats_value))

'''
Return the file system path to jgi-rqc-pipeline so we can use */tools and */lib

@return /path/to/jgi-rqc-pipelines
'''
def get_run_path():
    current_path = os.path.dirname(os.path.abspath(__file__))
    run_path = os.path.abspath(os.path.join(current_path, os.pardir))

    return run_path



'''
Simple read count using bbtools n_contigs field
- slightly different than in rqc_utility
n_scaffolds     n_contigs       scaf_bp contig_bp       gap_pct scaf_N50        scaf_L50        ctg_N50 ctg_L50 scaf_N90        scaf_L90        ctg_N90 ctg_L90 scaf_max        ctg_max scaf_n_gt50K     scaf_pct_gt50K  gc_avg  gc_std
1346616 1346616 405331416       405331415       0.000   1346616 301     1346615 301     1346616 301     1346615 301     301     301     0       0.000   0.44824 0.02675

'''
def get_read_count_fastq(fastq, log = None):
    read_cnt = 0

    if os.path.isfile(fastq):

        # EnvironmentModules.module(["load", "bbtools"])
        # bbtools faster than zcat | wc because bbtools uses pigz
        # cmd = "stats.sh format=3 in=%s" % fastq
        cmd = "#bbtools;stats.sh format=3 in=%s" % fastq
        cmd = convert_cmd(cmd)


        if log:
            log.info("- cmd: %s", cmd)

        std_out, std_err, exit_code = run_command(cmd, True)

        # EnvironmentModules.module(["unload", "bbtools"])

        if exit_code == 0 and std_out:

            line_list = std_out.split("\n")
            #print line_list
            val_list = str(line_list[1]).split() #.split('\t')
            #print "v = %s" % val_list

            read_cnt = int(val_list[1])

            if log:
                log.info("- read count: %s", read_cnt)

        else:
            if log:
                post_mortem_cmd(cmd, exit_code, std_out, std_err, log)


    else:
        log.error("- fastq: %s does not exist!", fastq)

    return read_cnt



'''
Subsampling calculation
0 .. 250k reads = 100%
250k .. 25m = 100% to 1%
25m .. 600m = 1%
600m+ .. oo < 1%

July 2014 - 15 runs > 600m (HiSeq-2500 Rapid) - 4 actual libraries / 85325 seq units
- returns new subsampling rate
'''
def get_subsample_rate(read_count):
    subsample = 0
    subsample_rate = 0.01
    max_subsample = 6000000 # 4 hours of blast time

    new_subsample_rate = 250000.0/read_count
    subsample_rate = max(new_subsample_rate, subsample_rate)
    subsample_rate = min(1, subsample_rate) # if subsample_rate > 1, then set to 1

    subsample = int(read_count * subsample_rate)

    if subsample > max_subsample:
        subsample = max_subsample

    subsample_rate = subsample / float(read_count)

    return subsample_rate


'''
Set color hash
- need to update to remove "c" parameter - used in too many places
'''
def set_colors(c, use_color = False):

    if use_color == False:

        color = {
            'black' : "",
            'red' : "",
            'green' : "",
            'yellow' : "",
            'blue' : "",
            'pink' : "",
            'cyan' : "",
            'white' : "",
            '' : ""
        }

    else:

        color = {
            'black' : "\033[1;30m",
            'red' : "\033[1;31m",
            'green' : "\033[1;32m",
            'yellow' : "\033[1;33m",
            'blue' : "\033[1;34m",
            'pink' : "\033[1;35m",
            'cyan' : "\033[1;36m",
            'white' : "\033[1;37m",
            '' : "\033[m"
        }


    return color

'''
New function that just returns colors
'''
def get_colors():
    
    color = {
        'black' : "\033[1;30m",
        'red' : "\033[1;31m",
        'green' : "\033[1;32m",
        'yellow' : "\033[1;33m",
        'blue' : "\033[1;34m",
        'pink' : "\033[1;35m",
        'cyan' : "\033[1;36m",
        'white' : "\033[1;37m",
        '' : "\033[m"
    }


    return color



'''
Returns msg_ok, msg_fail, msg_warn colored or not colored
'''
def get_msg_settings(color):

    msg_ok = "[  "+color['green']+"OK"+color['']+"  ]"
    msg_fail = "[ "+color['red']+"FAIL"+color['']+" ]"
    msg_warn = "[ "+color['yellow']+"WARN"+color['']+" ]"

    return msg_ok, msg_fail, msg_warn


'''
Use RQC's ap_tool to get the status
set mode = "-sa" to show all, even completed
'''
def get_analysis_project_id(seq_proj_id, target_analysis_project_id, target_analysis_task_id, output_path, log = None, mode = ""):

    if log:
        log.info("get_analysis_project_id: spid = %s, tapid = %s, tatid = %s", seq_proj_id, target_analysis_project_id, target_analysis_task_id)

    analysis_project_id = 0
    analysis_task_id = 0
    project_type = None
    task_type = None

    ap_list = os.path.join(output_path, "ap-info.txt")
    AP_TOOL = "/global/dna/projectdirs/PI/rqc/prod/jgi-rqc-pipeline/tools/ap_tool.py"
    #AP_TOOL = "/global/homes/b/brycef/git/jgi-rqc-pipeline/tools/ap_tool.py"
    cmd = "%s -spid %s -m psv -tapid %s -tatid %s %s > %s 2>&1" % (AP_TOOL, seq_proj_id, target_analysis_project_id, target_analysis_task_id, mode, ap_list)
    if log:
        log.info("- cmd: %s", cmd)
    else:
        print "- cmd: %s" % cmd
    std_out, std_err, exit_code = run_command(cmd, True)

    post_mortem_cmd(cmd, exit_code, std_out, std_err, log)

    if os.path.isfile(ap_list):

        ap_dict = {} # header = value
        cnt = 0
        fh = open(ap_list, "r")
        for line in fh:
            arr = line.strip().split("|")
            if cnt == 0:
                c2 = 0 # position of title in header
                for a in arr:
                    ap_dict[a.lower()] = c2
                    c2 += 1

            else:

                for a in ap_dict:

                    if ap_dict[a] + 1 > len(arr):
                        pass
                    else:
                        ap_dict[a] = arr[ap_dict[a]]


            cnt += 1

        fh.close()


        analysis_project_id = ap_dict.get("analysis project id")
        analysis_task_id = ap_dict.get("analysis task id")
        project_type = ap_dict.get("analysis product name")
        task_type = ap_dict.get("analysis task name")

    # nno such project
    if cnt == 1:
        analysis_project_id = 0
        analysis_task_id = 0

    if log:
        log.info("- project type: %s, task type: %s", project_type, task_type)
        log.info("- analysis_project_id: %s, analysis_task_id: %s", analysis_project_id, analysis_task_id)

    try:
        analysis_project_id = int(analysis_project_id)
        analysis_task_id = int(analysis_task_id)
    except:
        analysis_project_id = 0
        analysis_task_id = 0


    # ap = 4, at = 8 means its using the column names but didn't find anything
    if analysis_project_id < 100:
        analysis_project_id = 0
    if analysis_task_id < 100:
        analysis_task_id = 0

    return analysis_project_id, analysis_task_id


'''
For creating a dot file from the pipeline flow
'''
def append_flow(flow_file, orig_node, orig_label, next_node, next_label, link_label):

    fh = open(flow_file, "a")
    fh.write("%s|%s|%s|%s|%s\n" % (orig_node, orig_label, next_node, next_label, link_label))
    fh.close()

'''
Flow file format:
# comment
*label|PBMID Pipeline run for BTXXY<br><font point-size="10">Run Date: 2017-09-28 14:22:50</font>
# origin node, origin label, next node, next label, link label
input_h5|BTXXY H5<br><font point-size="10">3 smrtcells</font>|assembly|HGAP Assembly<FONT POINT-SIZE="10"><br>3 contigs, 13,283,382bp</FONT>|HGAP v4.0.1

nodes should be the output of the transformation between the nodes
e.g. input fastq (25m reads) --[ bbtools subsampling ]--> subsampled fastq (10m reads)

creates a dot file, to convert to png use:
$ module load graphviz
$ dot -T png (dot file) > (png file)

More info on formatting the labels
http://www.graphviz.org/content/node-shapes#html
'''
def dot_flow(flow_file, dot_file, log = None):

    if not os.path.isfile(flow_file):
        if log:
            log.info("- cannot find flow file:  %s", flow_file)
        else:
            print "Cannot find flow file: %s" % flow_file
        return


    fhw = open(dot_file, "w")

    fhw.write("// dot file\n")
    fhw.write("digraph rqc {\n") # directed graph
    fhw.write("    node [shape=box];\n")
    fhw.write("    rankdir=LR;\n")

    fh = open(flow_file, "r")

    for line in fh:
        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        # graph label
        if line.startswith("*label"):
            arr = line.split("|")
            label = flow_replace(str(arr[1]))

            fhw.write("    label=<%s>;\n" % label)
            fhw.write("    labelloc=top;\n")

        else:

            arr = line.split("|")
            #print arr

            if len(arr) == 5:

                org_node = arr[0]
                org_label = str(arr[1])
                next_node = arr[2]
                next_label = str(arr[3])
                link_label = str(arr[4])

                # must be <br/> in the dot file, I have a habit of using <br>
                org_label = flow_replace(org_label)
                next_label = flow_replace(next_label)
                link_label = flow_replace(link_label)

                # label are enclosed by < > instead of " " to handle html-ish markups

                if next_node:
                    link = "    %s -> %s;\n" % (org_node, next_node)
                    if link_label:
                        link = "    %s -> %s [label=<%s>];\n" % (org_node, next_node, link_label)
                    fhw.write(link)

                if org_label:
                    label = "    %s [label=<%s>];\n" % (org_node, org_label)
                    fhw.write(label)

                if next_label:
                    label = "    %s [label=<%s>];\n" % (next_node, next_label)
                    fhw.write(label)



    fh.close()

    fhw.write("}\n")

    fhw.close()
    if log:
        log.info("- created dot file: %s", dot_file)

    return dot_file

'''
simple replacements
'''
def flow_replace(my_str):

    new_str = my_str.replace("<br>", "<br/>").replace("<smf>", "<font point-size=\"10\">").replace("</f>", "</font>")

    return new_str


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## main program


if __name__ == "__main__":
    # unit tests

    print human_size(102192203)
    print human_size(250000000000)
    #print get_read_count_fastq("/global/projectb/scratch/brycef/sag/phix/11185.1.195330.UNKNOWN_matched.fastq.gz")

    cmd = "#bbtools;bbduk.sh in=/global/dna/dm_archive/sdm/illumina//01/14/88/11488.1.208132.UNKNOWN.fastq.gz ref=/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa outm=/global/projectb/scratch/brycef/phix/11488/11488.1.208132.UNKNOWN_matched.fastq.gz outu=/global/projectb/scratch/brycef/phix/11488/11488.1.208132.UNKNOWN_unmatched.fastq.gz"
    print convert_cmd(cmd)
    cmd = "#pigz;pigz /global/projectb/scratch/brycef/align/BTOYH/genome/11463.6.208000.CAAGGTC-AGACCTT.filter-RNA.fastq.gz-genome.sam"
    print convert_cmd(cmd)

    cmd = "#java;java -version"
    print convert_cmd(cmd)

    dot_flow("/global/projectb/scratch/brycef/pbmid/BWOAU/f2.flow", "/global/projectb/scratch/brycef/pbmid/BWOAU/BWOUAx.dot")

    sys.exit(0)


