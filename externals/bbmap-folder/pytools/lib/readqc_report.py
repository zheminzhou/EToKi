#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
 Readqc report: record stat key-value in readqc-stats.txt

 ### JGI_Analysis_Utility_Illumina::illumina_read_level_report

 Created: Jul 24 2013

 sulsj (ssul@lbl.gov)

"""


import os
import sys

## custom libs in "../lib/"
srcDir = os.path.dirname(__file__)
sys.path.append(os.path.join(srcDir, 'tools'))     ## ./tools
sys.path.append(os.path.join(srcDir, '../lib'))    ## rqc-pipeline/lib
sys.path.append(os.path.join(srcDir, '../tools'))  ## rqc-pipeline/tools

from readqc_constants import RQCReadQcConfig, ReadqcStats
from rqc_constants import RQCExitCodes
from os_utility import run_sh_command
from common import append_rqc_stats, append_rqc_file


statsFile = RQCReadQcConfig.CFG["stats_file"]
filesFile = RQCReadQcConfig.CFG["files_file"]


"""
 Title      : read_megablast_hits
 Function   : This function generates tophit list of megablast against different databases.
 Usage      : read_megablast_hits(db_name, log)
 Args       : blast db name or full path
 Returns    : SUCCESS
              FAILURE
 Comments   :

"""
def read_megablast_hits(db, log):
    currentDir = RQCReadQcConfig.CFG["output_path"]
    megablastDir = "megablast"
    megablastPath = os.path.join(currentDir, megablastDir)

    statsFile = RQCReadQcConfig.CFG["stats_file"]
    filesFile = RQCReadQcConfig.CFG["files_file"]


    ##
    ## Process blast output files
    ##
    matchings = 0
    hitCount = 0
    parsedFile = os.path.join(megablastPath, "megablast.*.%s*.parsed" % (db))
    matchings, _, exitCode = run_sh_command("grep -v '^#' %s 2>/dev/null | wc -l " % (parsedFile), True, log)

    if exitCode == 0: ## if parsed file found.
        t = matchings.split()

        if len(t) == 1 and t[0].isdigit():
            hitCount = int(t[0])

        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_MATCHING_HITS + " " + db, hitCount, log)

        ##
        ## add .parsed file
        ##
        parsedFileFound, _, exitCode = run_sh_command("ls %s" % (parsedFile), True, log)

        if parsedFileFound:
            parsedFileFound = parsedFileFound.strip()
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_PARSED_FILE + " " + db, os.path.join(megablastPath, parsedFileFound), log)
        else:
            log.error("- Failed to add megablast parsed file of %s." % (db))
            return RQCExitCodes.JGI_FAILURE

        ##
        ## wc the top hits
        ##
        topHit = 0
        tophitFile = os.path.join(megablastPath, "megablast.*.%s*.parsed.tophit" % (db))
        tophits, _, exitCode = run_sh_command("grep -v '^#' %s 2>/dev/null | wc -l " % (tophitFile), True, log)

        t = tophits.split()

        if len(t) == 1 and t[0].isdigit():
            topHit = int(t[0])

        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_TOP_HITS + " " + db, topHit, log)

        ##
        ## wc the taxonomic species
        ##
        spe = 0
        taxlistFile = os.path.join(megablastPath, "megablast.*.%s*.parsed.taxlist" % (db))
        species, _, exitCode = run_sh_command("grep -v '^#' %s 2>/dev/null | wc -l " % (taxlistFile), True, log)

        t = species.split()

        if len(t) == 1 and t[0].isdigit():
            spe = int(t[0])

        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_TAX_SPECIES + " " + db, spe, log)

        ##
        ## wc the top 100 hit
        ##
        top100hits = 0
        top100hitFile = os.path.join(megablastPath, "megablast.*.%s*.parsed.top100hit" % (db))
        species, _, exitCode = run_sh_command("grep -v '^#' %s 2>/dev/null | wc -l " % (top100hitFile), True, log)

        t = species.split()

        if len(t) == 1 and t[0].isdigit():
            top100hits = int(t[0])

        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_TOP_100HITS + " " + db, top100hits, log)

        ##
        ## Find and add taxlist file
        ##
        taxListFound, _, exitCode = run_sh_command("ls %s" % (taxlistFile), True, log)
        taxListFound = taxListFound.strip()

        if taxListFound:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_TAXLIST_FILE + " " + db, os.path.join(megablastPath, taxListFound), log)
        else:
            log.error("- Failed to add megablast taxlist file of %s." % (db))
            return RQCExitCodes.JGI_FAILURE

        ##
        ## Find and add tophit file
        ##
        tophitFound, _, exitCode = run_sh_command("ls %s" % (tophitFile), True, log)
        tophitFound = tophitFound.strip()

        if tophitFound:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_TOPHIT_FILE + " " + db, os.path.join(megablastPath, tophitFound), log)
        else:
            log.error("- Failed to add megablast tophit file of %s." % (db))
            return RQCExitCodes.JGI_FAILURE

        ##
        ## Find and add top100hit file
        ##
        top100hitFound, _, exitCode = run_sh_command("ls %s" % (top100hitFile), True, log)
        top100hitFound = top100hitFound.strip()

        if top100hitFound:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_TOP100HIT_FILE + " " + db, os.path.join(megablastPath, top100hitFound), log)
        else:
            log.error("- Failed to add megablast top100hit file of %s." % (db))
            return RQCExitCodes.JGI_FAILURE

    else:
        log.info("- No blast hits for %s." % (db))


    return RQCExitCodes.JGI_SUCCESS



"""
 Title      : read_level_qual_stats
 Function   : Generate qual scores and plots of read level 20mer sampling
 Usage      : read_level_mer_sampling($analysis, $summary_file_dir)
 Args       : 1) A reference to an JGI_Analysis object
              2) current working folder wkdir/uniqueness
 Returns    : JGI_SUCCESS: Illumina read level report could be successfully generated.
              JGI_FAILURE: Illumina read level report could not be generated.
 Comments   : This function is intended to be called at the very end of the illumina read level data processing script.

"""
def read_level_mer_sampling(dataToRecordDict, dataFile, log):
    retCode = RQCExitCodes.JGI_FAILURE

    ## Old data
    #nSeq   nStartUniMer    fracStartUniMer nRandUniMer fracRandUniMer
    ##  0         1                 2            3             4
    ##25000     2500               0.1         9704         0.3882

    ## New data
    #count  first   rand    first_cnt   rand_cnt
    #  0     1        2       3        4
    #25000  66.400  76.088  16600   19022
    #50000  52.148  59.480  13037   14870
    #75000  46.592  53.444  11648   13361
    #100000 43.072  49.184  10768   12296 ...

    if os.path.isfile(dataFile):
        with open(dataFile, "r") as merFH:
            lines = merFH.readlines()

            ## last line
            t = lines[-1].split('\t')
            # breaks 2016-09-07
            #assert len(t) == 5

            totalMers = int(t[0])

            ## new by bbcountunique
            uniqStartMerPer = float("%.2f" % (float(t[1])))
            uniqRandtMerPer = float("%.2f" % (float(t[2])))

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_20MER_SAMPLE_SIZE] = totalMers
            dataToRecordDict[ReadqcStats.ILLUMINA_READ_20MER_PERCENTAGE_STARTING_MERS] = uniqStartMerPer
            dataToRecordDict[ReadqcStats.ILLUMINA_READ_20MER_PERCENTAGE_RANDOM_MERS] = uniqRandtMerPer

            retCode = RQCExitCodes.JGI_SUCCESS

    else:
        log.error("- qhist file not found: %s" % (dataFile))


    return retCode




"""
 Title      : base_level_qual_stats
 Function   : Generate qual scores and plots of read level QC
 Usage      : base_level_qual_stats($analysis, $)
 Args       : 1) A reference to an JGI_Analysis object
              2) current working folder wkdir/qual
 Returns    : JGI_SUCCESS: Illumina read level report could be successfully generated.
              JGI_FAILURE: Illumina read level report could not be generated.
 Comments   : This function is intended to be called at the very end of the illumina base level data processing script.

"""
def base_level_qual_stats(dataToRecordDict, reformatObqhistFile, log):
    cummlatPer = 0
    cummlatBase = 0

    statsPerc = {30:0, 25:0, 20:0, 15:0, 10:0, 5:0}
    statsBase = {30:0, 25:0, 20:0, 15:0, 10:0, 5:0}

    Q30_seen = 0
    Q25_seen = 0
    Q20_seen = 0
    Q15_seen = 0
    Q10_seen = 0
    Q5_seen = 0

    ## New format
    ##Median    38
    ##Mean  37.061
    ##STDev 4.631
    ##Mean_30   37.823
    ##STDev_30  1.699
    ##Quality   bases   fraction
    #0  159 0.00008
    #1  0   0.00000
    #2  12175   0.00593
    #3  0   0.00000
    #4  0   0.00000
    #5  0   0.00000
    #6  0   0.00000

    allLines = open(reformatObqhistFile).readlines()
    for l in allLines[::-1]:
        l = l.strip()

        ##
        ## obqhist file format example
        ##
        # #Median   36
        # #Mean 33.298
        # #STDev    5.890
        # #Mean_30  35.303
        # #STDev_30 1.517
        # #Quality  bases   fraction
        # 0 77098   0.00043
        # 1 0   0.00000
        # 2 0   0.00000
        # 3 0   0.00000
        # 4 0   0.00000
        # 5 0   0.00000
        # 6 0   0.00000

        if len(l) > 0:
            if l.startswith("#"):
                if l.startswith("#Mean_30"):
                    dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q30_SCORE_MEAN] = l.split('\t')[1]

                elif l.startswith("#STDev_30"):
                    dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q30_SCORE_STD] = l.split('\t')[1]

                elif l.startswith("#Mean"):
                    dataToRecordDict[ReadqcStats.ILLUMINA_BASE_OVERALL_BASES_Q_SCORE_MEAN] = l.split('\t')[1]

                elif l.startswith("#STDev"):
                    dataToRecordDict[ReadqcStats.ILLUMINA_BASE_OVERALL_BASES_Q_SCORE_STD] = l.split('\t')[1]

                continue


            qavg = None
            nbase = None
            percent = None

            t = l.split()

            try:
                qavg = int(t[0])
                nbase = int(t[1])
                percent = float(t[2])
            except IndexError:
                log.warn("parse error in base_level_qual_stats: %s %s %s %s" % (l, qavg, nbase, percent))
                continue

            log.debug("base_level_qual_stats(): qavg and nbase and percent: %s %s %s" % (qavg, nbase, percent))

            cummlatPer += percent * 100.0
            cummlatPer = float("%.f" % (cummlatPer))

            if cummlatPer > 100:
                cummlatPer = 100.0 ## RQC-621

            cummlatBase += nbase

            if qavg == 30:
                Q30_seen = 1
                statsPerc[30] = cummlatPer
                statsBase[30] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q30] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C30] = cummlatBase

            elif qavg == 25:
                Q25_seen = 1
                statsPerc[25] = cummlatPer
                statsBase[25] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q25] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C25] = cummlatBase

            elif qavg == 20:
                Q20_seen = 1
                statsPerc[20] = cummlatPer
                statsBase[20] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q20] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C20] = cummlatBase

            elif qavg == 15:
                Q15_seen = 1
                statsPerc[15] = cummlatPer
                statsBase[15] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q15] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C15] = cummlatBase

            elif qavg == 10:
                Q10_seen = 1
                statsPerc[10] = cummlatPer
                statsBase[10] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q10] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C10] = cummlatBase

            elif qavg == 5:
                Q5_seen = 1
                statsPerc[5] = cummlatPer
                statsBase[5] = cummlatBase

                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q5] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C5] = cummlatBase

    ## Double check that no value is missing.
    if Q25_seen == 0 and Q30_seen != 0:
        Q25_seen = 1
        statsPerc[25] = statsPerc[30]
        statsBase[25] = statsBase[30]
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q25] = cummlatPer
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C25] = cummlatBase

    if Q20_seen == 0 and Q25_seen != 0:
        Q20_seen = 1
        statsPerc[20] = statsPerc[25]
        statsBase[20] = statsBase[25]
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q20] = cummlatPer
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C20] = cummlatBase

    if Q15_seen == 0 and Q20_seen != 0:
        Q15_seen = 1
        statsPerc[15] = statsPerc[20]
        statsBase[15] = statsBase[20]
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q15] = cummlatPer
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C15] = cummlatBase

    if Q10_seen == 0 and Q15_seen != 0:
        Q10_seen = 1
        statsPerc[10] = statsPerc[15]
        statsBase[10] = statsBase[15]
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q10] = cummlatPer
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C10] = cummlatBase

    if Q5_seen == 0 and Q10_seen != 0:
        Q5_seen = 1
        statsPerc[5] = statsPerc[10]
        statsBase[5] = statsBase[10]
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_Q5] = cummlatPer
        dataToRecordDict[ReadqcStats.ILLUMINA_BASE_C5] = cummlatBase

    if Q30_seen == 0:
        log.error("Q30 is 0. Base quality values are ZERO.")

    log.debug("Q and C values: %s" % (dataToRecordDict))


    return RQCExitCodes.JGI_SUCCESS


"""

 Title      : q20_score
 Function   : this method returns Q20 using a qrpt file as input
 Usage      : JGI_QC_Utility::qc20_score($qrpt)
 Args       : $_[0] : qrpt file.
 Returns    : a number of Q20 score
 Comments   :

"""
# def q20_score(qrpt, log):
#     log.debug("qrpt file %s" % (qrpt))
# 
#     q20 = None
#     num = 0
# 
#     if os.path.isfile(qrpt):
#         with open(qrpt, "r") as qrptFH:
#             for l in qrptFH:
#                 num += 1
# 
#                 if num == 1:
#                     continue
# 
#                 ##############
#                 ## Old format
#                 ## READ1.qrpt
#                 ## column   count   min max sum       mean  Q1  med Q3  IQR lW  rW  A_Count C_Count G_Count T_Count N_Count Max_count
#                 ## 1        378701  2   34  12447306  32.87 31  34  34  3   27  34  108573  83917   81999   104127  85      378701
#                 ## 2        378701  2   34  12515957  33.05 33  34  34  1   32  34  112178  83555   84449   98519   0       378701
#                 ## 3        378701  2   34  12519460  33.06 33  34  34  1   32  34  104668  72341   80992   120700  0       378701
#                 ## 4        378701  2   37  13807944  36.46 37  37  37  0   37  37  96935   95322   83958   102440  46      378701
#                 ## 5        378701  2   37  13790443  36.42 37  37  37  0   37  37  114586  68297   78020   117740  58      378701
#                 ##
#                 ## or
#                 ##
#                 ## READ2.qrpt
#                 ## column   count   min max sum       mean  Q1  med Q3  IQR lW  rW  A_Count C_Count G_Count T_Count N_Count Max_count
#                 ## 1        378701  2   34  8875097   23.44 25  26  28  3   21  32  106904  84046   81795   105956  0       378701
#                 ## 2        378701  2   34  6543224   17.28 15  16  26  11  2   34  107573  77148   97953   88998   7029    378701
#                 ## 3        378701  2   34  7131741   18.83 16  16  26  10  2   34  96452   83003   107891  91355   0       378701
#                 ## 4        378701  2   37  9686653   25.58 19  32  33  14  2   37  97835   78304   87944   114618  0       378701
#                 ## 5        378701  2   37  10208226  26.96 25  33  35  10  10  37  98021   90611   89040   101029  0       378701
# 
#                 pos = None
#                 mean = None
#                 t = l.split("\t")
#                 assert len(t) > 6
#                 pos = int(t[0])
#                 mean = float(t[5])
# 
#                 if mean and pos:
#                     if mean < 20:
#                         return pos - 1
#                     else:
#                         q20 = pos
# 
#     else:
#         log.error("- qhist file not found: %s" % (qrpt))
#         return None
# 
# 
#     return q20


def q20_score_new(bqHist, readNum, log):
    log.debug("q20_score_new(): bqHist file = %s" % (bqHist))

    q20 = None

    if os.path.isfile(bqHist):
        with open(bqHist, "r") as qrptFH:
            for l in qrptFH:
                if l.startswith('#'):
                    continue

                ## New data
                #  0        1       2        3       4      5        6      7       8       9       10      11      12     13       14     15       16      17      18
                ##BaseNum   count_1 min_1   max_1   mean_1  Q1_1    med_1   Q3_1    LW_1    RW_1    count_2 min_2   max_2   mean_2  Q1_2    med_2   Q3_2    LW_2    RW_2
                #  0        6900    0       36      33.48   33      34      34      29      36      6900    0       36      33.48   33      34      34      29       36

                pos = None
                mean = None
                t = l.split("\t")
                pos = int(t[0]) + 1

                if readNum == 1:
                    mean = float(t[4])
                else:
                    mean = float(t[13])

                if mean and pos:
                    if mean < 20:
                        return pos - 1
                    else:
                        q20 = pos
                        
    else:
        log.error("- bqHist file not found: %s" % (bqHist))
        return None    

    return q20


"""
 Title      : read_level_qual_stats
 Function   : Generate qual scores and plots of read level QC
 Usage      : read_level_qual_stats($analysis, $)
 Args       : 1) A reference to an JGI_Analysis object
              2) current working folder wkdir/qual
 Returns    : JGI_SUCCESS: Illumina read level report could be successfully generated.
              JGI_FAILURE: Illumina read level report could not be generated.
 Comments   : This function is intended to be called at the very end of the illumina read level data processing script.


"""
def read_level_qual_stats(dataToRecordDict, qhistTxtFullPath, log):
    retCode = RQCExitCodes.JGI_FAILURE

    cummlatPer = 0.0
    Q30_seen = 0
    Q25_seen = 0
    Q20_seen = 0
    Q15_seen = 0
    Q10_seen = 0
    Q5_seen = 0

    if os.path.isfile(qhistTxtFullPath):
        stats = {30:0, 25:0, 20:0, 15:0, 10:0, 5:0}

        allLines = open(qhistTxtFullPath).readlines()

        for l in allLines[::-1]:
            if not l:
                break
            if l.startswith('#'):
                continue

            t = l.split()
            assert len(t) == 3
            qavg = int(t[0])
            percent = float(t[2]) * 100.0 ## 20140826 Changed for bbtools

            cummlatPer = cummlatPer + percent
            cummlatPer = float("%.2f" % cummlatPer)

            if qavg <= 30 and qavg > 25 and Q30_seen == 0:
                Q30_seen = 1
                stats[30] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q30] = cummlatPer

            elif qavg <= 25 and qavg > 20 and Q25_seen == 0:
                Q25_seen = 1
                stats[25] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q25] = cummlatPer

            elif qavg <= 20 and qavg > 15 and Q20_seen == 0:
                Q20_seen = 1
                stats[20] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q20] = cummlatPer

            elif qavg <= 15 and qavg > 10 and Q15_seen == 0:
                Q15_seen = 1
                stats[15] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q15] = cummlatPer

            elif qavg <= 10 and qavg > 5 and Q10_seen == 0:
                Q10_seen = 1
                stats[10] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q10] = cummlatPer

            elif qavg <= 5 and Q5_seen == 0:
                Q5_seen = 1
                stats[5] = cummlatPer
                dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q5] = cummlatPer

        ### Double check that no value is missing.
        if Q25_seen == 0 and Q30_seen != 0:
            Q25_seen = 1
            stats[25] = stats[30]

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q25] = cummlatPer
        if Q20_seen == 0 and Q25_seen != 0:
            Q20_seen = 1
            stats[20] = stats[25]

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q20] = cummlatPer
        if Q15_seen == 0 and Q20_seen != 0:
            Q15_seen = 1
            stats[15] = stats[20]

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q15] = cummlatPer
        if Q10_seen == 0 and Q15_seen != 0:
            Q10_seen = 1
            stats[10] = stats[15]

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q10] = cummlatPer
        if Q5_seen == 0 and Q10_seen != 0:
            Q5_seen = 1
            stats[5] = stats[10]

            dataToRecordDict[ReadqcStats.ILLUMINA_READ_Q5] = cummlatPer
        if Q30_seen == 0:
            log.error("Q30 is 0 . Read quality values are ZERO.")

        log.debug("Q30 %s, Q25 %s, Q20 %s, Q15 %s, Q10 %s, Q5 %s" % \
                  (stats[30], stats[25], stats[20], stats[15], stats[10], stats[5]))


        retCode = RQCExitCodes.JGI_SUCCESS

    else:
        log.error("- qhist file not found: %s" % (qhistTxtFullPath))


    return retCode


"""

 Title      : read_gc_mean
 Function   : This function generates average GC content % and its standard deviation and put them into database.
 Usage      : read_gc_mean($analysis)
 Args       : 1) A reference to an JGI_Analysis object
 Returns    : JGI_SUCCESS:
              JGI_FAILURE:
 Comments   :

"""
def read_gc_mean(histFile, log):
    mean = 0.0
    stdev = 0.0
    retCode = RQCExitCodes.JGI_FAILURE

    if os.path.isfile(histFile):
        with open(histFile, "r") as histFH:
            line = histFH.readline() ## we only need the first line

            # Ex) #Found 1086 total values totalling 420.3971. <0.387106 +/- 0.112691>
            if len(line) == 0 or not line.startswith("#Found"):
                log.error("- GC content hist text file does not contains right results: %s, %s" % (histFile, line))
                retCode = RQCExitCodes.JGI_FAILURE

            else:
                toks = line.split()
                assert len(toks) == 9
                mean = float(toks[6][1:]) * 100.0
                stdev = float(toks[8][:-1])  * 100.0
                log.debug("mean, stdev = %.2f, %.2f" % (mean, stdev))

            retCode = RQCExitCodes.JGI_SUCCESS
    else:
        log.error("- gc hist file not found: %s" % (histFile))


    return retCode, mean, stdev


if __name__ == "__main__":

    exit(0)


## EOF
