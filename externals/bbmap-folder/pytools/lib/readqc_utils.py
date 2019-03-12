#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
 Read qc utilities

 Created: Jul 24 2013
 sulsj

"""

import os
import subprocess
import matplotlib
import numpy as np

from common import checkpoint_step
from os_utility import make_dir, change_mod, run_sh_command, rm_dir
from readqc_constants import RQCReadQcConfig, RQCContamDb, RQCReadQcReferenceDatabases, RQCReadQc
from rqc_utility import safe_basename, get_cat_cmd, localize_file, safe_dirname
from rqc_constants import RQCExitCodes

matplotlib.use("Agg")  ## This needs to skip the DISPLAY env var checking
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import mpld3
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


""" STEP1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fast_subsample_fastq_sequences

 Title    : fast_subsample_fastq_sequences

 Function : This function subsamples the data from a specified fastq file.

 Usage    : fast_subsample_fastq_sequences( $seq_unit, subsampledFile, $subsamplingRate, $max_count, \$totalBaseCount, $totalReadNum, log)

 Args     : 1) The source fastq file.
            2) The destination fastq file.
            3) The percentage of read subsampling.
            4) The maximum number of reads at which to stop subsampling.
            5) A reference to the variable that will store the
               basecount.
            6) A reference to the variable that will store the
               number of reads.
            7) A reference to a JGI_Log object.

 Returns  : JGI_SUCCESS: The fastq data was successfully sampled.
            JGI_FAILURE: The fastq data could not be sampled.

 Comments : Pass as parameters both the subsample_rate and the read_count
            in order to stop subsampling at the read_count.
            The read_count can also be null, in which case the number
            of reads corresponding to the percentage subsample_rate will be subsampled.

@param fastq:  source fastq file (full path)
@param outputFileName: subsampled output file name (basename)
@param subsamplingRate: sample rate < 1.0
@param isStatGenOnly: boolean -> generate stats output or not
@param log

@return retCode: success or failure
@return subsampledFile: subsampled output file name (full path)
@return totalBaseCount: total #bases (to be added to readqc_stats.txt)
@return totalReadNum: total #reads (to be added to readqc_stats.txt)
@return subsampledReadNum: total #reads sampled (to be added to readqc_stats.txt)

"""
def fast_subsample_fastq_sequences(sourceFastq, outputFileName, subsamplingRate, isStatGenOnly, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbtoolsReformatShCmd = os.path.join(cdir, '../../reformat.sh')  #RQCReadQcCommands.BBTOOLS_REFORMAT_CMD

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
    log.info("Sampling %s at %.2f rate", sourceFastq, subsamplingRate)

    retCode = None
    totalBaseCount = 0
    totalReadNum = 0
    subsampledReadNum = 0
    bIsPaired = False
    readLength = 0

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)
    make_dir(subsamplePath)
    change_mod(subsamplePath, "0755")

    subsampledFile = os.path.join(subsamplePath, outputFileName)

    fileSize = os.path.getsize(sourceFastq)
    log.info("Source fastq file size = %s", fileSize)

    ## Replaced subsampler with reformat.sh
    ## subsample with bbtoolsReformatShCmd:
    ## $ reformat.sh in=7348.8.68143.fastq out=subsample.fastq samplerate=0.01 qout=33
    ## - 21G == 180.399 seconds ~ 6x faster than subsample_fastq_pl
    ## new subampler from BBTOOLS
    ## without qin=33 then it uses auto detect, Illumina is phread64 but we need to convert to phred33
    ##reformat.sh in=7257.1.64419.CACATTGTGAG.s1.0.fastq out=temp.out samplerate=0.02 qin=33 qout=33 overwrite

    ## 20140820
    ## bhist=<file>     Write a base composition histogram to file.             ## Cycle Nucleotide Composition
    ## gchist=<file>    Write a gc content histogram to file.                   ## Read GC, mean, std
    ## qhist=<file>     Write a quality histogram to file.                      ## Average Base Position Quality
    ## bqhist=<file>    Write a quality histogram designed for box plots.       ## Average Base Position Quality Boxplot
    ## obqhist=<file>   Write a base quality histogram to file.                 ## Base quality histogram; *.base_qual.stats

    reformatPrefix = os.path.basename(subsampledFile).replace(".fastq", "")

    reformatLogFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.log")
    reformatGchistFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.gchist.txt")  ## Read GC
    reformatBhistFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.bhist.txt")  ## Cycle Nucleotide Composition
    reformatQhistFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.qhist.txt")  ## Average Base Position Quality
    reformatBqhistFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.bqhist.txt")  ## Average Base Position Quality Boxplot
    reformatObqhistFile = os.path.join(subsamplePath, reformatPrefix + ".reformat.obqhist.txt")  ## Base quality histogram

    if not isStatGenOnly:  ## if subsampling for blast, do not generate the stat files
        subsampleCmd = "%s in=%s out=%s samplerate=%s qin=33 qout=33 ow=t > %s 2>&1 " % \
                       (bbtoolsReformatShCmd, sourceFastq, subsampledFile, subsamplingRate, reformatLogFile)
    else:
        subsampleCmd = "%s in=%s out=%s samplerate=%s qin=33 qout=33 ow=t gcplot=t bhist=%s qhist=%s gchist=%s gcbins=auto bqhist=%s obqhist=%s > %s 2>&1 " % \
                       (bbtoolsReformatShCmd, sourceFastq, subsampledFile, subsamplingRate, reformatBhistFile,
                        reformatQhistFile, reformatGchistFile, reformatBqhistFile, reformatObqhistFile, reformatLogFile)

    _, _, exitCode = run_sh_command(subsampleCmd, True, log, True)

    if exitCode == 0:
        ##java -ea -Xmx200m -cp /usr/common/jgi/utilities/bbtools/prod-33.18/lib/BBTools.jar jgi.ReformatReads in=7257.1.64419.CACATTGTGAG.s1.0.fastq out=temp.out samplerate=0.02 qin=33 qout=33 overwrite
        ##Executing jgi.ReformatReads [in=7257.1.64419.CACATTGTGAG.s1.0.fastq, out=temp.out, samplerate=0.02, qin=33, qout=33, overwrite]
        ##
        ##Unspecified format for output temp.out; defaulting to fastq.
        ##Input is being processed as paired
        ##Writing interleaved.
        ##Input:                    6661 reads              1671911 bases
        ##Processed:                278 reads           69778 bases
        ##Output:                   278 reads (4.17%)   69778 bases (4.17%)
        ##
        ##Time:                             0.181 seconds.
        ##Reads Processed:         278  1.54k reads/sec
        ##Bases Processed:       69778  0.39m bases/sec

        ## NEW
        if os.path.isfile(reformatLogFile):
            with open(reformatLogFile) as STAT_FH:
                for l in STAT_FH.readlines():
                    if l.startswith("Input:"):
                        toks = l.split()
                        totalBaseCount = int(toks[3])
                        totalReadNum = int(toks[1])
                    # elif l.startswith("Processed:") or l.startswith("Output:"):
                    elif l.startswith("Output:"):
                        toks = l.split()
                        subsampledReadNum = int(toks[1])
                    elif l.startswith("Input is being processed as"):
                        toks = l.split()
                        if toks[-1].strip() == "paired":
                            bIsPaired = True

            log.info("Total base count of input fastq = %s", totalBaseCount)
            log.info("Total num reads of input fastq = %s", totalReadNum)
            log.info("Total num reads of sampled = %s", subsampledReadNum)

            readLength = int(totalBaseCount / totalReadNum)

            log.info("Read length = %d", readLength)
            log.info("Paired = %s", bIsPaired)

            if totalReadNum > 0 and subsampledReadNum > 0:
                ##
                ## TODO: deal with sequnits with small number of reads
                ## How to record the status in readqc.log
                ##
                retCode = RQCExitCodes.JGI_SUCCESS
                log.info("illumina_readqc_subsampling complete: output file = %s", subsampledFile)

            elif totalReadNum > 0 and subsampledReadNum <= 0:
                retCode = RQCExitCodes.JGI_FAILURE
                log.error("illumina_readqc_subsampling failure. subsampledReadNum <= 0.")

            else:
                retCode = RQCExitCodes.JGI_FAILURE
                log.error("illumina_readqc_subsampling failure. totalReadNum <= 0 and subsampledReadNum <= 0.")

        else:
            retCode = RQCExitCodes.JGI_FAILURE
            log.error("illumina_readqc_subsampling failure. Can't find stat file from subsampling.")

    else:
        retCode = RQCExitCodes.JGI_FAILURE
        log.error("illumina_readqc_subsampling failure. Failed to run bbtoolsReformatShCmd. Exit code != 0")
        with open(reformatLogFile, 'r') as f:
            log.error(f.read())
        retCode = -2

    return retCode, subsampledFile, totalBaseCount, totalReadNum, subsampledReadNum, bIsPaired, readLength


""" STEP2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_unique_20_mers

 Title:        write_unique_k_mers (k=20 or 25)
 Function:     Given a fastq file, finds unique 20/25 mers from the start
               of the read and along a random position of the read
 Usage:        write_unique_20_mers(\@seq_files, $log)
 Args:         1) ref to an array containing bz2 zipped fastq file path(s)
               2) output directory
               3) log file object
 Returns:      exitCode, merSamplerOutFile, pngPlotFile, htmlPlotFile
 Comments:     Using bbcountunique.sh's output file named merSampler.<fastq_name>.m20.e25000
               create a plot png file, merSampler.<fastq_name>.m20.e25000.png

               bbcountunique.sh: Generates a kmer uniqueness histogram, binned by file position.
               There are 3 columns for single reads, 6 columns for paired:
               count        number of reads or pairs processed
               r1_first     percent unique 1st kmer of read 1
               r1_rand      percent unique random kmer of read 1
               r2_first     percent unique 1st kmer of read 2
               r2_rand      percent unique random kmer of read 2
               pair         percent unique concatenated kmer from read 1 and 2

@param fastq:  source fastq file (full path)
@param log

@return retCode: success or failure
@return mersampler_out_file: plot data (to be added to readqc_files.txt)
@return pngPlotFile: output plot (to be added to readqc_files.txt)
@return htmlPlotFile: output d3 interactive plot (to be added to readqc_files.txt)

"""
def write_unique_20_mers(fastq, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbcountuniqueShCmd = os.path.join(cdir, '../../bbcountunique.sh')   #RQCReadQcCommands.BBCOUNTUNIQUE_SH_CMD

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
    uniqMerDir = "uniqueness"
    uniqMerPath = os.path.join(READ_OUTPUT_PATH, uniqMerDir)
    make_dir(uniqMerPath)
    change_mod(uniqMerPath, "0755")

    ## cmd line for merSampler
    uniqMerSize = RQCReadQc.ILLUMINA_MER_SAMPLE_MER_SIZE  ## 20 ==> 25 RQC-823 08102016
    reportFreq = RQCReadQc.ILLUMINA_MER_SAMPLE_REPORT_FRQ  ## 25000

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## bbcountunique.sh, a new version of sampler from bbtools
    ## bbcountunique.sh in=$FILENAME out=out.txt percent=t count=t cumulative=t
    ## ex) bbcountunique.sh in=$SEQDIR/dna/$SEQFILE.fastq.gz out=7378.1.69281.CGATG-2.txt percent=t count=t cumulative=f int=f
    ## ex2)
    ## cmd:  bbcountunique.sh k=20 interval=25000 in=7601.1.77813.CTTGTA.fastq.gz out=7601.1.77813.CTTGTA.merSampler.m20.e25000_2 percent=t count=t cumulative=f int=f

    log.info("bbcountunique.sh started.")

    merSamplerOutFile = os.path.join(uniqMerPath, sequnitFileNamePrefix + ".merSampler.m" + str(uniqMerSize) + ".e" + str(reportFreq) + "_2")

    ## RQC-823
    ## Adding shuffling before bbcountunique
    ## 08302016 Reverted to no shuffling
    ##
    ## shuffle.sh in=input.fastq.gz out=stdout.fq -Xmx40g | bbcountunique.sh in=stdin.fq -Xmx40g ==> not working
    # shuffledFastqFile = os.path.join(uniqMerPath, sequnitFileNamePrefix + ".shuffled.fq")
    # suffleCmd = "%s in=%s out=%s" % (shuffleShCmd, fastq, shuffledFastqFile)
    # stdOut, stdErr, exitCode = run_sh_command(suffleCmd, True, log, True)
    # if exitCode != 0:
    #     log.error("failed to suffle fastq for unique mer analysis.")
    #     return RQCExitCodes.JGI_FAILURE, None, None, None

    bbcountuniqCmd = "%s in=%s out=%s k=%s interval=%s percent=t count=t cumulative=f int=f ow=t" \
                     % (bbcountuniqueShCmd, fastq, merSamplerOutFile, uniqMerSize, reportFreq)

    _, _, exitCode = run_sh_command(bbcountuniqCmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to sample unique %s mers by bbcountunique.sh.", uniqMerSize)
        return RQCExitCodes.JGI_FAILURE, None, None, None

    log.info("bbcountunique.sh completed.")

    ## Old plotting data
    # nSeq   nStartUniMer    fracStartUniMer nRandUniMer fracRandUniMer
    ##  0         1                 2            3             4
    ##25000     2500               0.1         9704         0.3882

    ## New plotting data from bbcountunique
    # count  first   rand    first_cnt   rand_cnt
    #  0     1        2       3        4
    # 25000  66.400  76.088  16600   19022
    # 50000  52.148  59.480  13037   14870
    # 75000  46.592  53.444  11648   13361
    # 100000 43.072  49.184  10768   12296 ...

    pngPlotFile = None
    htmlPlotFile = None

    if os.path.isfile(merSamplerOutFile):
        ## sanity check

        ## OLD
        ## #nSeq    nStartUniMer    fracStartUniMer nRandUniMer fracRandUniMer
        ## ex) 25000    16594   0.6638  18986   0.7594
        ##     50000    29622   0.5211  33822   0.5934
        ##     75000    41263   0.4656  47228   0.5362
        ##    100000    52026   0.4305  59545   0.4927 ...
        """
        2016-09-07
        #count  first   rand    first_cnt       rand_cnt        avg_quality     perfect_prob
        25000   96.480  98.636  24120   24659   30.36   80.94
        50000   96.204  97.996  24051   24499   30.41   81.17
        75000   95.512  97.568  23878   24392   29.99   80.06
        100000  95.408  97.588  23852   24397   30.24   80.78
        125000  95.176  97.240  23794   24310   30.23   80.86
        """

        line = None
        numData = 0

        with open(merSamplerOutFile, "r") as FH:
            lines = FH.readlines()
            line = lines[-1]  ## get the last line
            numData = sum(1 for l in lines)

        toks = line.split()
        assert len(toks) == 7, "ERROR: merSamplerOutFile format error: %s " % (merSamplerOutFile)

        if numData < 3:
            log.error("Not enough data in merSamplerOutFile: %s", merSamplerOutFile)
            return RQCExitCodes.JGI_SUCCESS, None, None, None


        ## Generating plots
        rawDataMatrix = np.loadtxt(merSamplerOutFile, delimiter='\t', comments='#')
        # Bryce: 2016-09-07, its 7 now.  Its failed 622 pipelines ...
        # assert len(rawDataMatrix[1][:]) == 5

        fig, ax = plt.subplots()

        markerSize = 5.0
        lineWidth = 1.5

        ## Note: no need to show all the data points
        ## If the number of data points > 5k, get only 5k data.
        jump = 1
        if len(rawDataMatrix[:, 0]) > 10000:
            jump = int(len(rawDataMatrix[:, 0]) / 5000)

        xData = rawDataMatrix[:, 0][0::jump]
        yData = rawDataMatrix[:, 1][0::jump]
        yData2 = rawDataMatrix[:, 2][0::jump]

        totalReadNum = rawDataMatrix[-1, 0]  ## sampled read num from the last line of the data file
        assert int(totalReadNum) > 0
        maxX = int(totalReadNum) * 3

        p1 = ax.plot(xData, yData, 'g', marker='x', markersize=markerSize, linewidth=lineWidth, label="Starting %s Mer Uniqueness" % (str(uniqMerSize)), alpha=0.5)
        p2 = ax.plot(xData, yData2, 'b', marker='x', markersize=markerSize, linewidth=lineWidth, label="Random %s Mer Uniqueness" % (str(uniqMerSize)), alpha=0.5)

        ## Curve-fitting
        from scipy.optimize import curve_fit

        ## fit function: f(x)=a*log(x)+b
        def fit_func(x, a, b):
            return a * np.log(x) + b

        fitpars, _ = curve_fit(fit_func, rawDataMatrix[:, 0], rawDataMatrix[:, 1])

        fix_x = [i for i in range(25000, maxX, 25000)]
        ax.plot(fix_x, fit_func(fix_x, *fitpars), 'r', linewidth=lineWidth, label="fit", alpha=0.5)

        ax.set_xlabel("Read Sampled", fontsize=12, alpha=0.5)
        ax.set_ylabel("Percentage Unique", fontsize=12, alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        fontProp = FontProperties()
        fontProp.set_size("small")
        fontProp.set_family("Bitstream Vera Sans")
        ax.legend(loc=1, prop=fontProp)
        ax.set_xlim([0, maxX])
        ax.set_ylim([0, 100])
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        labels = ["%.2f" % i for i in rawDataMatrix[:, 1]]
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=labels))

        labels = ["%.2f" % i for i in rawDataMatrix[:, 2]]
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=labels))

        ## Create both dynamic and static plots
        pngPlotFile = merSamplerOutFile + "_mer_sampler_plot.png"
        plt.savefig(pngPlotFile, dpi=fig.dpi)

        htmlPlotFile = merSamplerOutFile + "_mer_sampler_plot_d3.html"
        mpld3.save_html(fig, htmlPlotFile)

        log.info("New data file from bbcountunique: %s", merSamplerOutFile)
        log.info("New png plot: %s", pngPlotFile)
        log.info("New D3 plot: %s", htmlPlotFile)

    else:
        log.error("Cannot find merSamplerOutFile by bbcountunique.sh, %s", merSamplerOutFile)
        return RQCExitCodes.JGI_FAILURE, None, None, None

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # if os.path.isfile(merSamplerOutFile):
    #     line = None
    #     numData = 0
    #
    #     with open(merSamplerOutFile, "r") as FH:
    #         line = FH.readlines()[-1]  ## get the last line
    #         numData = len(line)
    #
    #     ## #nSeq    nStartUniMer    fracStartUniMer nRandUniMer fracRandUniMer
    #     ## ex) 25000    16594   0.6638  18986   0.7594
    #     ##     50000    29622   0.5211  33822   0.5934
    #     ##     75000    41263   0.4656  47228   0.5362
    #     ##    100000    52026   0.4305  59545   0.4927 ...
    #     """
    #     2016-09-07
    #     #count  first   rand    first_cnt       rand_cnt        avg_quality     perfect_prob
    #     25000   96.480  98.636  24120   24659   30.36   80.94
    #     50000   96.204  97.996  24051   24499   30.41   81.17
    #     75000   95.512  97.568  23878   24392   29.99   80.06
    #     100000  95.408  97.588  23852   24397   30.24   80.78
    #     125000  95.176  97.240  23794   24310   30.23   80.86
    #     """
    #
    #     toks = line.split()
    #     assert len(toks)==7, "ERROR: merSamplerOutFile format error: %s " % (merSamplerOutFile)
    #
    #     # totalReadNum = int(toks[0])
    #     # log.info("Total number of reads = %s." % (totalReadNum))
    #     # assert totalReadNum > 0
    #
    # else:
    #     log.error("cannot find mersampler_out_file, %s" % (merSamplerOutFile))
    #     return RQCExitCodes.JGI_FAILURE, None, None, None
    #
    # if numData < 3:
    #     log.error("not enough data in %s" % (merSamplerOutFile))
    #     return RQCExitCodes.JGI_FAILURE, None, None, None

    ## verify that the mer sampler output file was created
    # if not os.path.isfile(merSamplerOutFile):
    #     log.error("failed to find output file for %s Mer Uniqueness: %s" % (str(uniqMerSize), merSamplerOutFile))
    # else:
    #     log.info("MerSampler output file successfully generated (%s)." % (merSamplerOutFile))

    ## verify that the mer sampler plot png file was created
    if not os.path.isfile(pngPlotFile):
        log.warning("Failed to find output plot png file for %s Mer Uniqueness", str(uniqMerSize))
    else:
        log.info("MerSampler output png file successfully generated (%s)", pngPlotFile)

    if not os.path.isfile(htmlPlotFile):
        log.warning("Failed to find output d3 plot html file for %s Mer Uniqueness", str(uniqMerSize))
    else:
        log.info("MerSampler output png file successfully generated (%s)", htmlPlotFile)


    return RQCExitCodes.JGI_SUCCESS, merSamplerOutFile, pngPlotFile, htmlPlotFile


""" STEP3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_read_gc

 Title:        illumina_read_gc
 Function:     Takes path to fastq file and generates
               read gc histograms (txt and png)
 Usage:        illumina_read_gc($fastq_path, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.

@param fastq:  source fastq file (full path)
@param log

@return reformatGchistFile: hist text data (to be added to readqc_files.txt)
@return pngFile: output plot (to be added to readqc_files.txt)
@return htmlFile: output d3 interactive plot (to be added to readqc_files.txt)
@return meanVal: gc mean (to be added to readqc_stats.txt)
@return stdevVal: gc stdev (to be added to readqc_stats.txt)

"""
def illumina_read_gc(fastq, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatGchistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.gchist.txt")  ## gc hist

    log.debug("gchist file: %s", reformatGchistFile)

    ## Gen Average Base Position Quality plot
    if not os.path.isfile(reformatGchistFile):
        log.error("Gchist file not found: %s", reformatGchistFile)
        return None, None, None, None, None, None, None

    ## File format
    ## #Mean    41.647
    ## #Median  42.000
    ## #Mode    42.000
    ## #STDev   4.541
    ## #GC  Count
    ## 0.0  0
    ## 1.0  0
    ## 2.0  0
    ## 3.0  0
    ## 4.0  0
    meanVal = None
    medVal = None
    modeVal = None
    stdevVal = None

    with open(reformatGchistFile, "r") as STAT_FH:
        for l in STAT_FH.readlines():
            if l.startswith("#Mean"):
                meanVal = l.strip().split('\t')[1]
            elif l.startswith("#Median"):
                medVal = l.strip().split('\t')[1]
            elif l.startswith("#Mode"):
                modeVal = l.strip().split('\t')[1]
            elif l.startswith("#STDev"):
                stdevVal = l.strip().split('\t')[1]

    rawDataMatrix = np.loadtxt(reformatGchistFile, comments='#', usecols=(0, 1, 2))  ## only use 3 colums: GC, Count, Cumulative

    ## In addition to the %GC and # reads, the cumulative read % is added.
    assert len(rawDataMatrix[1][:]) == 3

    fig, ax = plt.subplots()

    markerSize = 5.0
    lineWidth = 1.5

    p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, alpha=0.5)

    ax.set_xlabel("%GC", fontsize=12, alpha=0.5)
    ax.set_ylabel("Read count", fontsize=12, alpha=0.5)
    ax.grid(color="gray", linestyle=':')

    ## Add tooltip
    toolTipStrReadCnt = ["Read count=%d" % i for i in rawDataMatrix[:, 1]]
    toolTipStrGcPerc = ["GC percent=%.1f" % i for i in rawDataMatrix[:, 0]]
    toolTipStrReadPerc = ["Read percent=%.1f" % (i * 100.0) for i in rawDataMatrix[:, 2]]
    toolTipStr = ["%s, %s, %s" % (i, j, k) for (i, j, k) in
                  zip(toolTipStrGcPerc, toolTipStrReadCnt, toolTipStrReadPerc)]
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=toolTipStr))

    pngFile = os.path.join(qualPath, sequnitFileNamePrefix + ".gchist.png")
    htmlFile = os.path.join(qualPath, sequnitFileNamePrefix + ".gchist.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFile)

    ## Save Matplotlib plot in png format
    plt.savefig(pngFile, dpi=fig.dpi)


    return reformatGchistFile, pngFile, htmlFile, meanVal, stdevVal, medVal, modeVal


""" STEP5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_avg_base_quality_stats

 Title:        write_base_quality_stats
 Function:     Takes path to fastq file and generates
               quality plots for each read
 Usage:        write_base_quality_stats($fastq_path, $analysis, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      None.
 Comments:     None.


@param fastq:  source fastq file (full path)
@param log

@return reformatObqhistFile: output data file (to be added to readqc_files.txt)

"""
def write_avg_base_quality_stats(fastq, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatObqhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.obqhist.txt")  ## base composition histogram

    log.debug("obqhist file: %s", reformatObqhistFile)

    ## Gen base composition histogram
    if not os.path.isfile(reformatObqhistFile):
        log.error("Obqhist file not found: %s", reformatObqhistFile)
        return None
    else:
        return reformatObqhistFile


""" STEP6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_count_q_score

 Title:        count_q_score
 Function:     Given a fastq (bz2 zipped or unzipped)
               file path, creates a histogram of
               the quality scores in the file
 Usage:        count_q_score($fastq, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     Generates a file named <fastq_path>.qhist that has
               the quality score histogram.

@param fastq: source fastq file (full path)
@param log

@return reformatObqhistFile: output data file (to be added to readqc_files.txt)
@return pngFile: output plot file (to be added to readqc_files.txt)
@return htmlFile: output d3 interactive plot file (to be added to readqc_files.txt)

"""
def illumina_count_q_score(fastq, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatObqhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.obqhist.txt")  ## base composition histogram

    log.debug("obqhist file: %s", reformatObqhistFile)

    rawDataMatrix = np.loadtxt(reformatObqhistFile, delimiter='\t', comments='#')
    assert len(rawDataMatrix[1][:]) == 3

    ## Qavg nrd percent
    ##    0    1     2
    fig, ax = plt.subplots()

    markerSize = 5.0
    lineWidth = 1.5

    p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 2], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, alpha=0.5)
    ax.set_xlabel("Average Read Quality", fontsize=12, alpha=0.5)
    ax.set_ylabel("Fraction of Reads", fontsize=12, alpha=0.5)
    ax.grid(color="gray", linestyle=':')

    ## Add tooltip
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=list(rawDataMatrix[:, 2])))

    pngFile = os.path.join(qualPath, sequnitFileNamePrefix + ".avg_read_quality_histogram.png")
    htmlFile = os.path.join(qualPath, sequnitFileNamePrefix + ".avg_read_quality_histogram.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFile)

    ## Save Matplotlib plot in png format
    plt.savefig(pngFile, dpi=fig.dpi)

    return reformatObqhistFile, pngFile, htmlFile


""" STEP7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_calculate_average_quality

 Title:        illumina_calculate_average_quality
 Function:     Given a fastq (subsampled) file, calculates average quality in 21 mer windows.
 Usage:        illumina_calculate_average_quality($fastq_path, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     Several output files are generated in the directory that
               the script is run.
               1) Text output file with 21mer start position, number of mers read,
               total mers, and average accuracy of the bin
               2) A gnuplot png file named <fastq_name>.21mer.qual.png

               The function assumes that the fastq file exists.
               The 21mer qual script was writtten by mli


@return retCode: success or failure
@return stat_file: plot data (to be added to readqc_files.txt)
@return pngFile: output plot (to be added to readqc_files.txt)

"""
## Removed!
##def illumina_calculate_average_quality(fastq, log):



""" STEP8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_find_common_motifs

 Title:        illumina_find_common_motifs
 Function:     Given a fastq (subsampled) file, finds most common N-string motifs.
 Usage:        illumina_find_common_motifs($fastq_path, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     An output file is generated in the directory
               1) Text output summary file most common motifs and perecent total motifs.

               The function assumes that the fastq file exists.
               The nstutter script was writtten by jel


@param fastq: source fastq file (full path)
@param log

@return retCode: success or failure
@return nstutterStatFile: output stutter data file (to be added to readqc_files.txt)


"""
def illumina_find_common_motifs(fastq, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    patterNFastqPlCmd = os.path.join(cdir, '../tools/patterN_fastq.pl') #RQCReadQcCommands.PATTERN_FASTQ_PL

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
    stutterDir = "stutter"
    stutterPath = os.path.join(READ_OUTPUT_PATH, stutterDir)
    make_dir(stutterPath)
    change_mod(stutterPath, "0755")

    nstutterStatFile = os.path.join(stutterPath, sequnitFileNamePrefix + ".nstutter.stat")

    ## ex) patterN_fastq.pl -analog -PCT 0.1 -in 7601.1.77813.CTTGTA.s0.01.fastq > 7601.1.77813.CTTGTA.s0.01.nstutter.stat ; wait;
    makeStatFileCmd = "%s -analog -PCT 0.1 -in %s > %s " % (patterNFastqPlCmd, fastq, nstutterStatFile)
    combinedCmd = "%s; wait; " % (makeStatFileCmd)

    _, _, exitCode = run_sh_command(combinedCmd, True, log, True)

    if exitCode != 0:
        log.error("failed to run patterNFastqPlCmd. Exit code != 0.")
        return RQCExitCodes.JGI_FAILURE, None

    if os.path.isfile(nstutterStatFile):
        log.info("N stutter stat file successfully created (%s)", nstutterStatFile)
    else:
        log.warning("Could not locate N stutter stat file %s", nstutterStatFile)
        nstutterStatFile = "failed to generate"

    return RQCExitCodes.JGI_SUCCESS, nstutterStatFile


""" STEP9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_run_bwa

 Title:        illumina_run_bwa
 Function:     Given a fastq (subsampled) file path, runs bwa aligner
               Reads are aligned to each other.
 Usage:        illumina_run_bwa($fastq_file_path, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.


@param fastq: source fastq file (full path)
@param log

@return retCode: success or failure
@return summary_file: output bwa summary file  (to be added to readqc_files.txt)

"""

## REMOVED!
##def illumina_run_bwa(fastq, log):

def illumina_run_dedupe(fastq, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbdedupeShCmd = os.path.join(cdir, '../../bbdedupe.sh') #RQCReadQcCommands.BBDEDUPE_SH

    sequnitFileName, exitCode = safe_basename(fastq, log)
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    dupDir = "dupes"
    dupPath = os.path.join(READ_OUTPUT_PATH, dupDir)
    make_dir(dupPath)
    change_mod(dupPath, "0755")

    dedupeSummaryFile = os.path.join(dupPath, sequnitFileName + ".bwa.summary")

    ## dedupe.sh in=reads.fq s=0 ftr=49 ac=f int=f
    xmx = "-Xmx23G"
    bbdedupeShCmd = "%s %s in=%s out=null qin=33 ow=t s=0 ftr=49 ac=f int=f> %s 2>&1 " % (bbdedupeShCmd, xmx, fastq, dedupeSummaryFile)

    _, _, exitCode = run_sh_command(bbdedupeShCmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run bbdedupeShCmd.sh")
        return RQCExitCodes.JGI_FAILURE, None

    return RQCExitCodes.JGI_SUCCESS, dedupeSummaryFile


""" STEP10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_run_tagdust

 Title:        illumina_run_tagdust
 Function:     Given a fastq (subsampled) file path, runs tag dust to
               find common illumina artifacts.
 Usage:        illumina_run_tagdust($fastq_file_path, $log)
 Args:         1) path to subsampled fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.


@param fastq: source fastq file (full path)
@param log

@return retCode: success or failure
@return tagdust_out: output tagdust file (to be added to readqc_files.txt)

"""
## No longer needed!!
##def illumina_run_tagdust(fastq, log):



""" STEP11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_detect_read_contam

@param fastq: source fastq file (full path)
@param firstBp: first bp length to cut for read contam detection
@param log

@return retCode: success or failure
@return outFileList: output duk stat file list (to be added to readqc_files.txt)
@return ratioResultDict: output stat value dict (to be added to readqc_stats.txt)

"""
##def illumina_detect_read_contam(fastq, log):
## REMOVED!

# # def illumina_detect_read_contam2(fastq, firstBp, log):
## REMOVED! 08302016


"""
Contam removal by seal.sh

"""
def illumina_detect_read_contam3(fastq, firstBp, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    sealShCmd = os.path.join(cdir, '../../seal.sh')     #RQCReadQcCommands.SEAL_SH_CMD

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    numBadFiles = 0
    ratio = 0
    outFileDict = {}
    ratioResultDict = {}
    contamStatDict = {}

    catCmd, exitCode = get_cat_cmd(fastq, log)

    ## TODO: remove JGI Contaminants from CONTAM_DBS
    # CONTAM_DBS['artifact'] = ARTIFACT_FILE_NO_SPIKEIN
    # CONTAM_DBS['artifact_50bp'] = ARTIFACT_FILE_NO_SPIKEIN ## 20131203 Added for 50bp contam
    # CONTAM_DBS['DNA_spikein'] = ARTIFACT_FILE_DNA_SPIKEIN
    # CONTAM_DBS['RNA_spikein'] = ARTIFACT_FILE_RNA_SPIKEIN
    # CONTAM_DBS['contaminants'] = CONTAMINANTS
    # CONTAM_DBS['fosmid'] = FOSMID_VECTOR
    # CONTAM_DBS['mitochondrion'] = MITOCHONDRION_NCBI_REFSEQ
    # CONTAM_DBS['phix'] = PHIX
    # CONTAM_DBS['plastid'] = CHLOROPLAST_NCBI_REFSEQ
    # CONTAM_DBS['rrna'] = GENERAL_RRNA_FILE
    # CONTAM_DBS['microbes'] = MICROBES ## non-synthetic
    # CONTAM_DBS['synthetic'] = SYNTHETIC
    # CONTAM_DBS['adapters'] = ADAPTERS
    for db in RQCContamDb.CONTAM_DBS.iterkeys():
        if db == "artifact_50bp" and int(firstBp) == 20:
            sealStatsFile = os.path.join(qualPath, sequnitFileNamePrefix + ".artifact_20bp.seal.stats")
        else:
            sealStatsFile = os.path.join(qualPath, sequnitFileNamePrefix + "." + db + ".seal.stats")

        ## Localization file to /scratch/rqc
        # log.info("Contam DB localization started for %s", db)
        # localizedDb = localize_file(RQCContamDb.CONTAM_DBS[db], log)

        ## 04262017 Skip localization temporarily until /scratch can be mounted
        ## in shifter container
        # if os.environ['NERSC_HOST'] == "genepool":
        #     localizedDb = localize_file(RQCContamDb.CONTAM_DBS[db], log)
        # else:
        #     localizedDb = None
        #
        # if localizedDb is None:
        #     localizedDb = RQCContamDb.CONTAM_DBS[db]  ## use the orig location
        # else:
        #     log.info("Use the localized file, %s", localizedDb)
        localizedDb = RQCContamDb.CONTAM_DBS[db]

        ## 09112017 Manually add -Xmx23G
        xmx = "-Xmx23G"
        if db == "artifact_50bp":
            cutCmd = "cut -c 1-%s | " % (firstBp)
            cmd = "set -e; %s %s | %s %s in=stdin.fq out=null ref=%s k=22 hdist=0 stats=%s ow=t statscolumns=3 %s " % \
                  (catCmd, fastq, cutCmd, sealShCmd, localizedDb, sealStatsFile, xmx)
        elif db == "microbes":
            cmd = "%s in=%s out=null ref=%s hdist=0 mm=f mkf=0.5 ambig=random minlen=120 qtrim=rl trimq=10 stats=%s ow=t statscolumns=3 %s " % \
                  (sealShCmd, fastq, localizedDb, sealStatsFile, xmx)
        else:
            cmd = "%s in=%s out=null ref=%s k=22 hdist=0 stats=%s ow=t statscolumns=3 %s " % \
                  (sealShCmd, fastq, localizedDb, sealStatsFile, xmx)

        _, _, exitCode = run_sh_command(cmd, True, log, True)

        if exitCode != 0:
            log.error("Failed to run seal.sh cmd")
            return RQCExitCodes.JGI_FAILURE, None, None, None

        ## Parsing seal output
        if not os.path.isfile(sealStatsFile):
            log.warning("Cannot open contam output file %s", sealStatsFile)
            numBadFiles += 1
            continue  ## continue to next contam db

        maxStatCount = 0
        with open(sealStatsFile, "r") as sealFH:
            for line in sealFH:
                line.strip()
                if line.find("#Matched") != -1:
                    ## ex) #Matched  1123123 77.31231%
                    toks = line.split()
                    assert len(toks) == 3
                    ratio = toks[-1].replace('%', '')

                ## contamintaion stat
                if not line.startswith('#'):
                    t = line.rstrip().split('\t')
                    # contamStatDict["%s:%s" % (db, t[0])] = t[2].replace('%', '')
                    if maxStatCount < 10 and t[0].startswith("gi|"):
                        # contamStatDict["contam:%s:%s" % (db, t[0])] = t[2].replace('%', '')
                        contamStatDict["contam:%s:%s" % (db, "|".join(t[0].split('|')[:2]))] = t[2].replace('%', '') ## save only gi part (RQC-906)
                        maxStatCount += 1

        ## RQC-743
        if db == "artifact_50bp" and int(firstBp) == 20:
            db = "artifact_20bp"

        outFileDict[db + ".seal.stats"] = sealStatsFile
        log.debug("Contam db and matched ratio: %s = %f", RQCContamDb.CONTAM_KEYS[db], float(ratio))
        ratioResultDict[RQCContamDb.CONTAM_KEYS[db]] = float(ratio)

    if numBadFiles:
        log.info("Number of bad I/O cases = %s", numBadFiles)
        return RQCExitCodes.JGI_FAILURE, None, None, None


    return RQCExitCodes.JGI_SUCCESS, outFileDict, ratioResultDict, contamStatDict


""" STEP12 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_sciclone_analysis

 Title:        illumina_sciclone_analysis
 Function:     Takes path to fastq file and determines
               if it is from a multiplexed run or not
 Usage:        illumina_sciclone_analysis($subfastq, $log)
 Args:         1) fastq file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.


@param origFastq: source fastq file (full path)
@param isPairedEnd: pair- or single-ended
@param log

@return retCode: success or failure
@return ratioResultDict: output stat value dict (to be added to readqc_stats.txt)
@return dnaCountFile: output sam stat file (to be added to readqc_files.txt)
@return rnaCountFile: output sam stat file (to be added to readqc_files.txt)


"""

## Removed!
##def illumina_sciclone_analysis(origFastq, isPairedEnd, log, libName=None, isRna=None):

def illumina_sciclone_analysis2(origFastq, isPairedEnd, log, libName=None, isRna=None):
    ## detect lib is rna or not
    sequnitFileName, exitCode = safe_basename(origFastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    # if libName is None and isRna is None:
    #     _, _, libName, isRna = get_lib_info(sequnitFileNamePrefix, log)  ## seqProjectId, ncbiOrganismName not used
    #
    # if isRna == "N/A":
    #     log.error("Failed to get lib info for %s", sequnitFileNamePrefix)
    #     return RQCExitCodes.JGI_FAILURE, -1, -1

    if isRna == '1':
        isRna = True

    elif isRna == '0':
        isRna = False

    if isRna:
        log.debug("The lib is RNA (%s)", libName)
    else:
        log.debug("The lib is DNA (%s)", libName)

    if isPairedEnd:
        log.debug("It's pair-ended.")
    else:
        log.debug("It's single-ended.")

    ## output dir
    sciclone_dir = "sciclone_analysis"
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
    sciclone_path = os.path.join(READ_OUTPUT_PATH, sciclone_dir)

    ## Define the subdirectory.  If it exists already, remove it
    if os.path.isdir(sciclone_path):
        rm_dir(sciclone_path)

    make_dir(sciclone_path)
    change_mod(sciclone_path, "0755")

    ## NOTE: Save count file in analysis object
    dnaCountFile = None
    rnaCountFile = None

    if isRna:
        rnaCountFile = os.path.join(sciclone_path, sequnitFileNamePrefix + "_bbduk_sciclone_rna_count.txt")
    else:
        dnaCountFile = os.path.join(sciclone_path, sequnitFileNamePrefix + "_bbduk_sciclone_dna_count.txt")

    cdir = os.path.dirname(__file__)
    bbdukShCmd = os.path.join(cdir, '../../bbduk.sh')       #RQCReadQcCommands.BBDUK_SH_CMD
    bbdukRnaDb = RQCReadQcReferenceDatabases.SCICLONE_RNA2
    # bbdukDnaDb = RQCReadQcReferenceDatabases.SCICLONE_DNA2

    cmd = None

    if isRna:
        ## Localization file to /scratch/rqc
        # log.info("Sciclone RNA ref DB localization started for %s", bbdukRnaDb)
        # localizedDb = localize_file(bbdukRnaDb, log)

        ## 04262017 Skip localization temporarily until /scratch can be mounted
        ## in shifter container
        # if os.environ['NERSC_HOST'] == "genepool":
        #     localizedDb = localize_file(bbdukRnaDb, log)
        # else:
        #     localizedDb = None
        #
        # if localizedDb is None:
        #     localizedDb = bbdukRnaDb  ## use the orig location
        # else:
        #     log.info("Use the localized file, %s", localizedDb)

        localizedDb = bbdukRnaDb  ## use the orig location

        ## bbduk.sh in=7365.2.69553.AGTTCC.fastq.gz ref=/global/projectb/sandbox/gaag/bbtools/data/sciclone_rna.fa out=null fbm=t k=31 mbk=0 stats=sciclone2.txt
        cmd = "%s in=%s ref=%s out=null fbm=t k=31 mbk=0 stats=%s statscolumns=3 " % (bbdukShCmd, origFastq, localizedDb, rnaCountFile)

    else:
        ## Localization file to /scratch/rqc
        # log.info("Sciclone DNA ref DB localization started for %s", bbdukDnaDb)
        # localizedDb = localize_file(bbdukDnaDb, log)

        ## 04262017 Skip localization temporarily until /scratch can be mounted
        ## in shifter container
        # if os.environ['NERSC_HOST'] == "genepool":
        #     localizedDb = localize_file(bbdukRnaDb, log)
        # else:
        #     localizedDb = None
        #
        # if localizedDb is None:
        #     localizedDb = bbdukRnaDb  ## use the orig location
        # else:
        #     log.info("Use the localized file, %s", localizedDb)

        localizedDb = bbdukRnaDb  ## use the orig location

        ## bbduk.sh in=7257.1.64419.CACATTGTGAG.fastq.gz ref=/global/projectb/sandbox/gaag/bbtools/data/sciclone_dna.fa out=null fbm=t k=31 mbk=0 stats=sciclone1.txt
        cmd = "%s in=%s ref=%s out=null fbm=t k=31 mbk=0 stats=%s statscolumns=3 " % (bbdukShCmd, origFastq, localizedDb, dnaCountFile)

    _, _, exitCode = run_sh_command(cmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run bbduk.sh cmd")
        return RQCExitCodes.JGI_FAILURE, None, None

    log.debug("rnaCountFile = %s", rnaCountFile)
    log.debug("dnaCountFile = %s", dnaCountFile)


    return RQCExitCodes.JGI_SUCCESS, dnaCountFile, rnaCountFile


""" STEP13 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

illumina_read_megablast

 Title:        illumina_read_megablast
 Function:     Takes path(s) to bz2 zipped or gzipped fastq file
               and runs megablast against the reads.
 Usage:        illumina_read_megablast(\@seq_files, $subsampledFile, $read_length, $log)
 Args:         1) reference to an array containing bz2 zipped or gzipped fastq file path(s);
                  the files should all be compressed the same way
               2) subsampled fastq file
               3) read length
               4) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.


@param subsampledFile: source fastq file (full path)
@param log

@return retCode: success or failure

"""
## No longer needed. Removed.
##def illumina_read_megablast(subsampledFile, read_num_to_pass, log, blastDbPath=None):



""" STEP14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
##def illumina_read_blastn_refseq_microbial(subsampledFile, log, blastDbPath=None):
## 12212015 sulsj REMOVED!



""" STEP14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def illumina_read_blastn_refseq_archaea(subsampledFile, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbtoolsReformatShCmd = os.path.join(cdir, '../../reformat.sh')      #RQCReadQcCommands.BBTOOLS_REFORMAT_CMD

    ## verify the fastq file
    if not os.path.isfile(subsampledFile):
        log.error("Failed to find fastq file")
        return RQCExitCodes.JGI_FAILURE

    log.info("Read level contamination analysis using blastn")

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    ## output dir
    megablastDir = "megablast"
    megablastPath = os.path.join(READ_OUTPUT_PATH, megablastDir)
    make_dir(megablastPath)
    change_mod(megablastPath, "0755")

    ## 20140929 Replaced with reformat.sh
    queryFastaFileName = "reads.fa"

    ## reformat.sh for converting fastq to fasta
    cmd = "%s in=%s out=%s qin=33 qout=33 ow=t " % (bbtoolsReformatShCmd, subsampledFile, os.path.join(megablastPath, queryFastaFileName))

    _, _, exitCode = run_sh_command(cmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run reformat.sh to convert fastq to fasta: %s", cmd)
        return RQCExitCodes.JGI_FAILURE

    megablastOutputFile = None
    db = "refseq.archaea"

    log.info("---------------------------------------------")
    log.info("Start blastn search against %s", db)
    log.info("---------------------------------------------")

    ## final output ==> READ_OUTPUT_PATH/megablast
    retCode, megablastOutputFile = run_blastplus_py(os.path.join(megablastPath, queryFastaFileName), db, log)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        if megablastOutputFile is None:
            log.error("Failed to run blastn against %s. Ret = %s", db, retCode)
            retCode = RQCExitCodes.JGI_FAILURE
        elif megablastOutputFile == -143:
            log.warning("Blast overtime. Skip the search against %s.", db)
            retCode = -143  ## blast overtime
    else:
        log.info("Successfully ran blastn of reads against %s", db)
        retCode = RQCExitCodes.JGI_SUCCESS


    return retCode


""" STEP15 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def illumina_read_blastn_refseq_bacteria(subsampledFile, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbtoolsReformatShCmd = os.path.join(cdir, '../../reformat.sh')      #RQCReadQcCommands.BBTOOLS_REFORMAT_CMD

    ## verify the fastq file
    if not os.path.isfile(subsampledFile):
        log.error("Failed to find fastq file for blastn")
        return RQCExitCodes.JGI_FAILURE

    log.info("Read level contamination analysis using blastn")

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    ## output dir
    megablastDir = "megablast"
    megablastPath = os.path.join(READ_OUTPUT_PATH, megablastDir)
    make_dir(megablastPath)
    change_mod(megablastPath, "0755")

    ## 20140929 Replaced with reformat.sh
    queryFastaFileName = "reads.fa"

    ## reformat.sh for converting fastq to fasta
    cmd = "%s in=%s out=%s qin=33 qout=33 ow=t " % (bbtoolsReformatShCmd, subsampledFile, os.path.join(megablastPath, queryFastaFileName))

    _, _, exitCode = run_sh_command(cmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run reformat.sh to convert fastq to fasta: %s", cmd)
        return RQCExitCodes.JGI_FAILURE

    megablastOutputFile = None
    db = "refseq.bacteria"

    log.info("---------------------------------------------")
    log.info("Start blastn search against %s", db)
    log.info("---------------------------------------------")

    ## final output ==> READ_OUTPUT_PATH/megablast
    retCode, megablastOutputFile = run_blastplus_py(os.path.join(megablastPath, queryFastaFileName), db, log)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        if megablastOutputFile is None:
            log.error("Failed to run blastn against %s. Ret = %s", db, retCode)
            retCode = RQCExitCodes.JGI_FAILURE
        elif megablastOutputFile == -143:
            log.warning("Blast overtime. Skip the search against %s.", db)
            retCode = -143  ## blast overtime
    else:
        log.info("Successfully ran blastn of reads against %s", db)
        retCode = RQCExitCodes.JGI_SUCCESS

    return retCode


""" STEP16 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def illumina_read_blastn_nt(subsampledFile, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbtoolsReformatShCmd = os.path.join(cdir, '../../reformat.sh')      #RQCReadQcCommands.BBTOOLS_REFORMAT_CMD

    ## verify the fastq file
    if not os.path.isfile(subsampledFile):
        log.error("Failed to find fastq file for blastn")
        return RQCExitCodes.JGI_FAILURE

    log.info("Read level contamination analysis using blastn")

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    ## output dir
    megablastDir = "megablast"
    megablastPath = os.path.join(READ_OUTPUT_PATH, megablastDir)
    make_dir(megablastPath)
    change_mod(megablastPath, "0755")

    ## 20140929 Replaced with reformat.sh
    queryFastaFileName = "reads.fa"

    ## reformat.sh
    cmd = "%s in=%s out=%s qin=33 qout=33 ow=t " % (bbtoolsReformatShCmd, subsampledFile, os.path.join(megablastPath, queryFastaFileName))

    _, _, exitCode = run_sh_command(cmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run reformat.sh to convert fastq to fasta: %s", cmd)
        return RQCExitCodes.JGI_FAILURE

    megablastOutputFile = None
    db = "nt"

    log.info("----------------------------------")
    log.info("Start blastn search against %s", db)
    log.info("----------------------------------")

    ## final output ==> READ_OUTPUT_PATH/megablast
    retCode, megablastOutputFile = run_blastplus_py(os.path.join(megablastPath, queryFastaFileName), db, log)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        if megablastOutputFile is None:
            log.error("Failed to run blastn against %s. Ret = %s", db, retCode)
            retCode = RQCExitCodes.JGI_FAILURE
        elif megablastOutputFile == -143:
            log.warning("Blast overtime. Skip the search against %s.", db)
            retCode = -143  ## blast overtime
    else:
        log.info("Successfully ran blastn of reads against %s", db)
        retCode = RQCExitCodes.JGI_SUCCESS

    return retCode


""" STEP17 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
illumina_generate_index_sequence_detection_plot

"""
def illumina_generate_index_sequence_detection_plot(fastq, log, isMultiplexed=None):    ## TO BE REMOVED!
    isMultiplexed = 0
    if not os.path.isfile(fastq):
        log.error("Failed to find the input fastq file, %s", fastq)
        return RQCExitCodes.JGI_FAILURE, None, None, None
    else:
        log.info("fastq file for index sequence analysis: %s", fastq)

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    # if isMultiplexed is None:
    #     isMultiplexed = get_multiplex_info(sequnitFileNamePrefix, log)

    #retCode = None

    demultiplexStatsFile = None
    demultiplexPlotDataFile = None
    detectionPlotPngFile = None
    storedDemultiplexStatsFile = None

    if int(isMultiplexed) == 1:
        log.info("Multiplexed - start analyzing...")

        READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]
        demul_dir = "demul"
        demulPath = os.path.join(READ_OUTPUT_PATH, demul_dir)
        make_dir(demulPath)
        change_mod(demulPath, "0755")

        ## This version is sorted by percent for readability of stats
        demultiplexStatsFile = os.path.join(demulPath, sequnitFileNamePrefix + ".demultiplex_stats")

        ## This version has index column and sort by index for plot
        demultiplexPlotDataFile = os.path.join(demulPath, sequnitFileNamePrefix + ".demultiplex_stats.tmp")

        ## This path is relative to final qual location to be stored in analysis obj.
        detectionPlotPngFile = os.path.join(demulPath, sequnitFileNamePrefix + ".index_sequence_detection.png")

        storedDemultiplexStatsFile = os.path.join(demulPath, sequnitFileNamePrefix + ".demultiplex_stats")

        if not os.path.isfile(demultiplexStatsFile):
            indexSeq = None
            line = None
            header = None
            indexSeqCounter = {}

            catCmd, exitCode = get_cat_cmd(fastq, log)

            if fastq.endswith(".gz"):
                catCmd = "zcat"  ## pigz does not work with subprocess

            seqCount = 0

            try:
                proc = subprocess.Popen([catCmd, fastq], bufsize=2 ** 16, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

                while 1:
                    line = proc.stdout.readline()
                    if not line:
                        break

                    ## First line is header
                    line.strip()
                    header = line

                    ## Get second line of record - sequence
                    line = proc.stdout.readline()

                    ## Get third line - junk
                    line = proc.stdout.readline()

                    ## Get the final line (4th line) - quality
                    line = proc.stdout.readline()

                    ## Parse the header
                    headerFields = header.split(":")

                    ## The last index is the index
                    indexSeq = headerFields[-1].strip()
                    assert indexSeq

                    ## Increment the counts and store the index
                    if indexSeq in indexSeqCounter:
                        indexSeqCounter[indexSeq] += 1
                    else:
                        indexSeqCounter[indexSeq] = 1

                    seqCount += 1

            except Exception as e:
                if log:
                    log.error("Exception in file reading: %s", e)
                    log.error("Failed to read the given fastq file [%s]", fastq)
                    log.error("Fastq header doesn't have the index sequence: %s", header)
                    log.error("Index sequence analysis is skipped!")

                return RQCExitCodes.JGI_SUCCESS, None, None, None

            ## Open the output file handles for writing
            log.info("demultiplexPlotDataFile = %s", demultiplexPlotDataFile)
            log.info("detectionPlotPngFile = %s", detectionPlotPngFile)
            log.info("demultiplexStatsFile = %s", demultiplexStatsFile)
            log.info("storedDemultiplexStatsFile = %s", storedDemultiplexStatsFile)

            plotDataFH = open(demultiplexPlotDataFile, "w")
            statsFH = open(demultiplexStatsFile, "w")

            ## Count the total number of indexes found
            numIndexesFound = len(indexSeqCounter)

            ## Store the data header information for printing
            reportHeader = """# Demultiplexing Summary
#
# Seq unit name: %s
# Total sequences: %s
# Total indexes found: %s
# 1=indexSeq 2=index_sequence_count 3=percent_of_total
#
""" % (sequnitFileName, seqCount, numIndexesFound)

            statsFH.write(reportHeader)

            ## Sort by value, descending
            log.debug("Sort by value of indexSeqCounter")
            for indexSeq in sorted(indexSeqCounter, key=indexSeqCounter.get, reverse=True):
                perc = float(indexSeqCounter[indexSeq]) / float(seqCount) * 100
                l = "%s\t%s\t%.6f\n" % (indexSeq, indexSeqCounter[indexSeq], perc)
                statsFH.write(l)

            ## Sort by index and add id column for plotting
            log.debug("Sort by index of indexSeqCounter")
            i = 1
            for indexSeq in sorted(indexSeqCounter.iterkeys()):
                perc = float(indexSeqCounter[indexSeq]) / float(seqCount) * 100
                l = "%s\t%s\t%s\t%.6f\n" % (i, indexSeq, indexSeqCounter[indexSeq], perc)
                plotDataFH.write(l)
                i += 1

            plotDataFH.close()
            statsFH.close()

        log.debug("demultiplex plotting...")
        ## matplotlib plotting
        # data
        # Index_seq_id Index_seq indexSeqCounter percent
        # 1 AAAAAAAAAAAA    320 0.000549
        # 2 AAAAAAAAAAAC    16  0.000027
        # 3 AAAAAAAAAAAG    8   0.000014
        # 4 AAAAAAAAAACA    4   0.000007
        # 5 AAAAAAAAAACG    2   0.000003
        # 6 AAAAAAAAAAGA    6   0.000010
        # 7 AAAAAAAAAATA    6   0.000010

        # rawDataMatrix = np.loadtxt(demultiplexPlotDataFile, delimiter='\t', comments='#')
        # assert len(rawDataMatrix[1][:]) == 4

        ## For a textfile with 4000x4000 words this is about 10 times faster than loadtxt.
        ## http://stackoverflow.com/questions/14985233/load-text-file-as-strings-using-numpy-loadtxt
        def load_data_file(fname):
            data = []

            with open(fname, 'r') as FH:
                lineCnt = 0
                for line in FH:
                    if lineCnt > 10000: ## experienced out of mem with a stat file with 9860976 index sequences
                        log.warning("Too many index sequences. Only 10000 index sequences will be used for plotting.")
                        break
                    data.append(line.replace('\n', '').split('\t'))
                    lineCnt += 1

            return data

        def column(matrix, i, opt):
            if opt == "int":
                return [int(row[i]) for row in matrix]
            elif opt == "float":
                return [float(row[i]) for row in matrix]
            else:
                return [row[i] for row in matrix]

        rawDataMatrix = load_data_file(demultiplexPlotDataFile)

        fig, ax = plt.subplots()

        markerSize = 6.0
        lineWidth = 1.5

        p1 = ax.plot(column(rawDataMatrix, 0, "int"), column(rawDataMatrix, 3, "float"), 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Index Sequence ID", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        labels = ["%s" % i for i in column(rawDataMatrix, 1, "str")]

        ## Show index_seq in the plot
        for i in rawDataMatrix:
            ax.text(i[0], float(i[3]) + .2, "%s" % i[1])

        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=labels))

        detectionPlotHtml = os.path.join(demulPath, sequnitFileNamePrefix + ".index_sequence_detection.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, detectionPlotHtml)

        ## Save Matplotlib plot in png format
        plt.savefig(detectionPlotPngFile, dpi=fig.dpi)

        if exitCode != 0:
            log.error("Failed to create demulplex plot")
            return RQCExitCodes.JGI_FAILURE, None, None, None

        log.info("demulplex stats and plot generation completed!")


    else:
        log.info("Not multiplexed - skip this analysis.")
        return RQCExitCodes.JGI_SUCCESS, None, None, None

    if detectionPlotPngFile is not None and storedDemultiplexStatsFile is not None and os.path.isfile(
            detectionPlotPngFile) and os.path.isfile(storedDemultiplexStatsFile):
        return RQCExitCodes.JGI_SUCCESS, demultiplexStatsFile, detectionPlotPngFile, detectionPlotHtml

    else:
        return RQCExitCodes.JGI_FAILURE, None, None, None


""" STEP18 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end_of_read_illumina_adapter_check

    usage: kmercount_pos.py [-h] [-k <int>] [-c <int>] [-t <int>] [-p <file>]
                        fastaFile fastqFile [fastqFile ...]

    Count occurance of database kmers in reads

    positional arguments:
      fastaFile Input FASTA file(s). Text or gzip
      fastqFile Input FASTQ file(s). Text or gzip

    optional arguments:
      -h, --help show this help message and exit
      -k <int> kmer length (default: 16)
      -c <int> minimum allowed coverage (default: 2)
      -t <int> target coverage (default: 30)
      -p <file>, --plot <file> plot data and save as png to <file> (default: None)

    * NOTE
      - RQC-383 04082014: Updated to the newest version of kmercount_pos.py (readqc ver 5.0.4)

"""
def end_of_read_illumina_adapter_check(firstSubsampledFastqFile, log):
    cdir = os.path.dirname(__file__)
    kmercountPosCmd = os.path.join(cdir, 'kmercount_pos.py')    #RQCReadQcCommands.KMERCOUNT_POS_CMD
    adapterDbName = RQCReadQcReferenceDatabases.END_OF_READ_ILLUMINA_ADAPTER_CHECK_DB

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, exitCode = safe_basename(firstSubsampledFastqFile, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    plotFile = ""
    dataFile = ""

    adapterCheckDir = "adapter"
    adapterCheckPath = os.path.join(READ_OUTPUT_PATH, adapterCheckDir)
    make_dir(adapterCheckPath)
    change_mod(adapterCheckPath, "0755")

    plotFile = os.path.join(adapterCheckPath, sequnitFileNamePrefix + ".end_of_read_adapter_check.png")  ## ignored
    dataFile = os.path.join(adapterCheckPath, sequnitFileNamePrefix + ".end_of_read_adapter_check.txt")

    ## Localization file to /scratch/rqc
    log.info("illumina_adapter_check DB localization started for %s", adapterDbName)

    # if os.environ['NERSC_HOST'] == "genepool":
    #     localizedDb = localize_file(adapterDbName, log)
    # else:
    #     localizedDb = None
    #
    # if localizedDb is None:
    #     localizedDb = adapterDbName  ## use the orig location
    # else:
    #     log.info("Use the localized file, %s", localizedDb)

    localizedDb = adapterDbName  ## use the orig location

    ## ex) kmercount_pos.py --plot plot.png Artifacts.adapters_primers_only.fa subsample.fastq > counts_by_pos.txt
    cmd = "%s --plot %s %s %s > %s " % (kmercountPosCmd, plotFile, localizedDb, firstSubsampledFastqFile, dataFile)

    ## Run cmd
    _, _, exitCode = run_sh_command(cmd, True, log, True)
    assert exitCode == 0

    ## mpld3 plots gen
    rawDataMatrix = np.loadtxt(dataFile, delimiter='\t', comments='#', skiprows=0)
    assert len(rawDataMatrix[1][:]) == 3 or len(rawDataMatrix[1][:]) == 5

    ##pos    read1   read2
    ## 0       1       2

    ## This output file format is changed on 2013.06.26 (RQC-442)
    ##
    ## pos    read1_count     read1_perc      read2_count     read2_perc
    ##
    fig, ax = plt.subplots()

    markerSize = 3.5
    lineWidth = 1.0

    if len(rawDataMatrix[1][:]) != 5:  ## support for old file
        p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label="read1", alpha=0.5)
        p2 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 2], 'g', marker='d', markersize=markerSize, linewidth=lineWidth, label="read2", alpha=0.5)
        ax.set_ylabel("Read Count with Database K-mer", fontsize=12, alpha=0.5)

    else:
        p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 2], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label="read1", alpha=0.5)
        p2 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 4], 'g', marker='d', markersize=markerSize, linewidth=lineWidth, label="read2", alpha=0.5)
        ax.set_ylim([0, 100])
        ax.set_ylabel("Percent Reads with Database K-mer", fontsize=12, alpha=0.5)

    ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
    ax.yaxis.set_label_coords(-0.095, 0.75)
    fontProp = FontProperties()
    fontProp.set_size("small")
    fontProp.set_family("Bitstream Vera Sans")
    ax.legend(loc=1, prop=fontProp)
    ax.grid(color="gray", linestyle=':')

    ## Add tooltip
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=list(rawDataMatrix[:, 1])))
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=list(rawDataMatrix[:, 2])))

    pngFile = os.path.join(adapterCheckPath, sequnitFileNamePrefix + ".end_of_read_adapter_check.png")
    htmlFile = os.path.join(adapterCheckPath, sequnitFileNamePrefix + ".end_of_read_adapter_check.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFile)

    ## Save Matplotlib plot in png format
    plt.savefig(pngFile, dpi=fig.dpi)

    if exitCode != 0:
        log.error("Failed to run kmercountPosCmd")
        return RQCExitCodes.JGI_FAILURE, None, None, None

    if os.path.isfile(plotFile) and os.path.isfile(dataFile):
        log.info("kmercount_pos completed.")
        return RQCExitCodes.JGI_SUCCESS, dataFile, pngFile, htmlFile

    else:
        log.error("cannot find the output files from kmercount_pos. kmercount_pos failed.")
        return RQCExitCodes.JGI_FAILURE, None, None, None


""" STEP19 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

insert_size_analysis

Using bbmerge.sh from bbtools, create insert size histogram (static/interactive) plots using D3 and data file

e.q.
    java -ea -Xmx200m -cp /usr/common/jgi/utilities/bbtools/prod-v32.28/lib/BBTools.jar jgi.BBMerge
    in=/global/projectb/scratch/brycef/rqc-dev/staging/00/00/66/26/6626.2.48981.TTCTCC.fastq ihist=ihist.txt
    Executing jgi.BBMerge [in=/global/projectb/scratch/brycef/rqc-dev/staging/00/00/66/26/6626.2.48981.TTCTCC.fastq, ihist=ihist.txt]

    e.g.
    bbmerge.sh in=[path-to-fastq] hist=hist.txt
    - you should use the whole fastq, not just the subsampled fastq

"""
def insert_size_analysis(fastq, log):
    ## Tools
    cdir = os.path.dirname(__file__)
    bbmergeShCmd = os.path.join(cdir, '../../bbmerge.sh')       #RQCReadQcCommands.BBMERGE_SH_CMD

    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    plotFile = ""
    dataFile = ""

    insertSizeOutDir = "insert_size_analysis"
    insertSizeOutPath = os.path.join(READ_OUTPUT_PATH, insertSizeOutDir)
    make_dir(insertSizeOutPath)
    change_mod(insertSizeOutPath, "0755")

    plotFile = os.path.join(insertSizeOutPath, sequnitFileNamePrefix + ".insert_size_histo.png")
    htmlFile = os.path.join(insertSizeOutPath, sequnitFileNamePrefix + ".insert_size_histo.html")
    dataFile = os.path.join(insertSizeOutPath, sequnitFileNamePrefix + ".insert_size_histo.txt")

    ## TODO
    ## if it's single ended
    ## 1. rqcfilter.sh for adapter trim
    ## 2. reformat.sh for getting lhist
    ## 3. analyze lhist.txt

    ## ex) bbmerge.sh in=7601.1.77813.CTTGTA.fastq.gz hist=.../insert_size_analysis/7601.1.77813.CTTGTA.insert_size_histo.txt
    ## reads=1000000 --> 1M reads are enough for insert size analysis
    cmd = "%s in=%s hist=%s reads=1000000 " % (bbmergeShCmd, fastq, dataFile)

    ## Run cmd
    _, stdErr, exitCode = run_sh_command(cmd, True, log, True)

    if exitCode != 0:
        log.error("Failed to run bbmerge_sh_cmd.")
        return RQCExitCodes.JGI_FAILURE, None, None, None, None

    retCode = {}

    ## File format
    ## BBMerge version 5.0
    ## Finished reading
    ## Total time: 8.410 seconds.
    ##
    ## Pairs:           1000000
    ## Joined:          556805      55.681%
    ## Ambiguous:       9665        0.967%
    ## No Solution:     433474      43.347%
    ## Too Short:       56          0.006%
    ## Avg Insert:                  234.6
    ## Standard Deviation:          33.9
    ## Mode:                    250
    ##
    ## Insert range:            26 - 290
    ## 90th percentile:         277
    ## 75th percentile:         262
    ## 50th percentile:         238
    ## 25th percentile:         211
    ## 10th percentile:         188

    for l in stdErr.split('\n'):
        toks = l.split()
        if l.startswith("Total time"):
            retCode["total_time"] = toks[2]
        elif l.startswith("Reads"):
            retCode["num_reads"] = toks[1]
        elif l.startswith("Pairs"):
            retCode["num_reads"] = toks[1]
        elif l.startswith("Joined"):
            retCode["joined_num"] = toks[1]
            retCode["joined_perc"] = toks[2]
        elif l.startswith("Ambiguous"):
            retCode["ambiguous_num"] = toks[1]
            retCode["ambiguous_perc"] = toks[2]
        elif l.startswith("No Solution"):
            retCode["no_solution_num"] = toks[2]
            retCode["no_solution_perc"] = toks[3]
        elif l.startswith("Too Short"):
            retCode["too_short_num"] = toks[2]
            retCode["too_short_perc"] = toks[3]
        elif l.startswith("Avg Insert"):
            retCode["avg_insert"] = toks[2]
        elif l.startswith("Standard Deviation"):
            retCode["std_insert"] = toks[2]
        elif l.startswith("Mode"):
            retCode["mode_insert"] = toks[1]
        elif l.startswith("Insert range"):
            retCode["insert_range_start"] = toks[2]
            retCode["insert_range_end"] = toks[4]
        elif l.startswith("90th"):
            retCode["perc_90th"] = toks[2]
        elif l.startswith("50th"):
            retCode["perc_50th"] = toks[2]
        elif l.startswith("10th"):
            retCode["perc_10th"] = toks[2]
        elif l.startswith("75th"):
            retCode["perc_75th"] = toks[2]
        elif l.startswith("25th"):
            retCode["perc_25th"] = toks[2]

    log.debug("Insert size stats: %s", str(retCode))

    ## -------------------------------------------------------------------------
    ## plotting
    rawDataMatrix = np.loadtxt(dataFile, delimiter='\t', comments='#')

    ## File format
    ## loc val
    ## 1 11
    try:
        rawDataX = rawDataMatrix[:, 0]
        rawDataY = rawDataMatrix[:, 1]

    except IndexError:
        log.info("No solution from bbmerge.")
        return RQCExitCodes.JGI_SUCCESS, dataFile, None, None, retCode

    fig, ax = plt.subplots()

    markerSize = 5.0
    lineWidth = 1.5

    p1 = ax.plot(rawDataX, rawDataY, 'r', marker='o', markersize=markerSize, linewidth=lineWidth, alpha=0.5)

    ax.set_xlabel("Insert Size", fontsize=12, alpha=0.5)
    ax.set_ylabel("Count", fontsize=12, alpha=0.5)

    ax.grid(color="gray", linestyle=':')

    ## Add tooltip
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=list(rawDataY)))

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFile)

    ## Save Matplotlib plot in png format
    plt.savefig(plotFile, dpi=fig.dpi)

    ## Checking outputs
    if os.path.isfile(plotFile) and os.path.isfile(dataFile) and os.path.isfile(htmlFile):
        log.info("insert_size_analysis completed.")
        return RQCExitCodes.JGI_SUCCESS, dataFile, plotFile, htmlFile, retCode

    else:
        log.error("cannot find the output files. insert_size_analysis failed.")
        return RQCExitCodes.JGI_FAILURE, None, None, None, None


""" STEP20 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gc_divergence_analysis


"""
def gc_divergence_analysis(fastq, bIsPaired, srcDir, log):
    ## Tools
    # rscriptCmd = RQCReadQcCommands.RSCRIPT_CMD

    READ_FILES_FILE = RQCReadQcConfig.CFG["files_file"]
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, exitCode = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    ## get the bhist.txt produced by reformat.sh
    reformatBhistFile = None
    with open(READ_FILES_FILE, "r") as FFH:
        for l in FFH.readlines():
            if l.startswith("read base count text 1"):
                reformatBhistFile = l.split("=")[1].strip()
                assert os.path.isfile(reformatBhistFile), "reformatBhistFile does not exist: %s" % (l)
                break

    assert reformatBhistFile, "ERROR: reformatBhistFile cannot be found in %s." % (READ_FILES_FILE)
    log.debug("gc_divergence_analysis(): bhist file = %s", reformatBhistFile)

    gcDivergenceTransformedFile = os.path.join(qualPath, sequnitFileNamePrefix + ".gc.divergence.transformed.csv")  ## gc divergence transformed csv
    gcDivergenceTransformedPlot = os.path.join(qualPath, sequnitFileNamePrefix + ".gc.divergence.transformed.png")  ## gc divergence transformed plot
    gcDivergenceCoefficientsFile = os.path.join(qualPath, sequnitFileNamePrefix + ".gc.divergence.coefficients.csv")  ## gc divergence coefficients csv

    ## Check base composition histogram
    if not os.path.isfile(reformatBhistFile):
        log.error("Bhist file not found: %s", reformatBhistFile)
        return RQCExitCodes.JGI_FAILURE, None, None, None, None

    ##
    ## Need R/3.1.2 b/c of library(dplyr) and library(tidyr) are only installed for the R/3.1.2 version
    ##

    ## Transform bhist output from reformat.sh to csv
    if bIsPaired:
        # transformCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; format_signal_data --input %s --output %s --read both --type composition " %\
        transformCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; format_signal_data --input %s --output %s --read both --type composition " %\
                       (reformatBhistFile, gcDivergenceTransformedFile)

        # cmd = " ".join(["module load R; Rscript", "--vanilla", os.path.join(SRCDIR, "tools", "jgi-fastq-signal-processing", "format_signal_data"), ])
        # transformCmd = "module unload R; module load R; module load R/3.3.2; Rscript --vanilla %s --input %s --output %s --read both --type composition" %\
        if os.environ['NERSC_HOST'] in ("denovo", "cori"):
            transformCmd = "Rscript --vanilla %s --input %s --output %s --read both --type composition" %\
                           (os.path.join(srcDir, "tools/jgi-fastq-signal-processing/bin", "format_signal_data"), reformatBhistFile, gcDivergenceTransformedFile)
    else:
        # transformCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; format_signal_data --input %s --output %s --read 1 --type composition " %\
        transformCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; format_signal_data --input %s --output %s --read 1 --type composition " %\
                       (reformatBhistFile, gcDivergenceTransformedFile)

        # transformCmd = "module unload R; module load R; module load R/3.3.2; Rscript --vanilla %s --input %s --output %s --read 1 --type composition" %\
        if os.environ['NERSC_HOST'] in ("denovo", "cori"):
            transformCmd = "Rscript --vanilla %s --input %s --output %s --read 1 --type composition" %\
                           (os.path.join(srcDir, "tools/jgi-fastq-signal-processing/bin", "format_signal_data"), reformatBhistFile, gcDivergenceTransformedFile)

    log.debug("Transform cmd = %s", transformCmd)

    _, _, exitCode = run_sh_command(transformCmd, True, log, True)

    if exitCode != 0:
        log.info("Failed to run GC_DIVERGENCE_TRANSFORM.")
        return -1, None, None, None, None

    ## Compute divergence value
    coeff = []
    # modelCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; model_read_signal --input %s --output %s " %\
    modelCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; model_read_signal --input %s --output %s " %\
               (gcDivergenceTransformedFile, gcDivergenceCoefficientsFile)

    # modelCmd = "module unload python; module unload R; module load R/3.3.2; Rscript --vanilla %s --input %s --output %s" %\
    if os.environ['NERSC_HOST'] in ("denovo", "cori"):
        modelCmd = "Rscript --vanilla %s --input %s --output %s" %\
                   (os.path.join(srcDir, "tools/jgi-fastq-signal-processing/bin", "model_read_signal"), gcDivergenceTransformedFile, gcDivergenceCoefficientsFile)

    log.debug("Model cmd = %s", modelCmd)

    _, _, exitCode = run_sh_command(modelCmd, True, log, True)

    if exitCode != 0:
        log.info("Failed to run GC_DIVERGENCE_MODEL.")
        return -1, None, None, None, None

    ## Parsing coefficients.csv
    ## ex)
    ##  "read","variable","coefficient"
    ##  "Read 1","AT",2
    ##  "Read 1","AT+CG",2.6
    ##  "Read 1","CG",0.6
    ##  "Read 2","AT",1.7
    ##  "Read 2","AT+CG",2.3
    ##  "Read 2","CG",0.7
    assert os.path.isfile(gcDivergenceCoefficientsFile), "GC divergence coefficient file not found."

    with open(gcDivergenceCoefficientsFile) as COEFF_FH:
        for l in COEFF_FH.readlines():
            if l.startswith("\"read\""):
                continue  ## skip header line
            toks = l.strip().split(',')
            assert len(toks) == 3, "Unexpected GC divergence coefficient file format."
            coeff.append({"read": toks[0], "variable": toks[1], "coefficient": toks[2]})

    ## Plotting
    # plotCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; plot_read_signal --input %s --output %s --type composition " %\
    plotCmd = "module unload R; module load R/3.2.4; module load jgi-fastq-signal-processing/2.x; plot_read_signal --input %s --output %s --type composition " %\
              (gcDivergenceTransformedFile, gcDivergenceTransformedPlot)

    # plotCmd = "module unload python; module unload R; module load R/3.3.2; Rscript --vanilla %s --input %s --output %s --type composition" %\
    if os.environ['NERSC_HOST'] in ("denovo", "cori"):
        plotCmd = "Rscript --vanilla %s --input %s --output %s --type composition" %\
                  (os.path.join(srcDir, "tools/jgi-fastq-signal-processing/bin", "plot_read_signal"), gcDivergenceTransformedFile, gcDivergenceTransformedPlot)

    log.debug("Plot cmd = %s", plotCmd)

    _, _, exitCode = run_sh_command(plotCmd, True, log, True)

    if exitCode != 0:
        log.info("Failed to run GC_DIVERGENCE_PLOT.")
        return -1, None, None, None, None


    return RQCExitCodes.JGI_SUCCESS, gcDivergenceTransformedFile, gcDivergenceTransformedPlot, gcDivergenceCoefficientsFile, coeff


""" STEP22 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cleanup_readqc

Cleaning up the ReadQC analysis directory with unwanted files.

@param log

@return retCode: always return success

"""
def cleanup_readqc(log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    ## Purge FASTQs
    cmd = "rm -f %s/%s/*.fastq " % (READ_OUTPUT_PATH, "subsample")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    cmd = "rm -f %s/%s/*.fastq " % (READ_OUTPUT_PATH, "qual")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    ## Purge FASTA qual and Fasta file
    cmd = "rm -f %s/%s/reads.fa " % (READ_OUTPUT_PATH, "megablast")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    cmd = "rm -f %s/%s/reads.qual " % (READ_OUTPUT_PATH, "megablast")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    ## Delete Sciclone files
    cmd = "rm -f %s/%s/*.fastq " % (READ_OUTPUT_PATH, "sciclone_analysis")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    cmd = "rm -f %s/%s/*.fq " % (READ_OUTPUT_PATH, "sciclone_analysis")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    cmd = "rm -f %s/%s/*.sam " % (READ_OUTPUT_PATH, "sciclone_analysis")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    cmd = "rm -f %s/%s/*.sai " % (READ_OUTPUT_PATH, "sciclone_analysis")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    ## Purge Blast files
    cmd = "rm -f %s/%s/megablast*v*JFfTIT " % (READ_OUTPUT_PATH, "megablast")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)

    ## purge megablast.reads.fa.v.nt.FmLD2a10p90E30JFfTITW45; megablast.reads.fa.v.refseq.microbial.FmLD2a10p90E30JFfTITW45
    cmd = "rm -f %s/%s/megablast*v*FfTITW45 " % (READ_OUTPUT_PATH, "megablast")
    _, _, exitCode = run_sh_command(cmd, True, log)

    if exitCode != 0:
        log.error("Failed to execute %s; may be already purged.", cmd)


    return RQCExitCodes.JGI_SUCCESS


## ===========================================================================================================================


""" For NEW STEP4

QC plot generation using the outputs from reformat.sh

"""
def gen_average_base_position_quality_plot(fastq, bIsPaired, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatQhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.qhist.txt")  ## Average Base Position Quality

    log.debug("qhist file: %s", reformatQhistFile)

    ## Gen Average Base Position Quality Plot
    if not os.path.isfile(reformatQhistFile):
        log.error("Qhist file not found: %s", reformatQhistFile)
        return None, None, None

    ## New data format
    ## Load data from txt
    rawDataMatrix = np.loadtxt(reformatQhistFile, delimiter='\t', comments='#', skiprows=0)
    assert len(rawDataMatrix[1][:]) == 5 or len(rawDataMatrix[1][:]) == 3

    ## New data (paired)
    # BaseNum    Read1_linear    Read1_log   Read2_linear    Read2_log
    # 1  33.469  30.347  32.459  29.127
    # 2  33.600  32.236  32.663  29.532
    # 3  33.377  30.759  32.768  29.719

    fig, ax = plt.subplots()

    markerSize = 3.5
    lineWidth = 1.0
    p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='read1')

    if bIsPaired:
        p2 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 3], 'g', marker='d', markersize=markerSize, linewidth=lineWidth, label='read2')

    ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
    ax.set_ylabel("Average Quality Score", fontsize=12, alpha=0.5)
    ax.set_ylim([0, 45])
    fontProp = FontProperties()
    fontProp.set_size("small")
    fontProp.set_family("Bitstream Vera Sans")
    ax.legend(loc=1, prop=fontProp)
    ax.grid(color="gray", linestyle=':')

    ## Add tooltip
    mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=list(rawDataMatrix[:, 1])))

    if bIsPaired:
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=list(rawDataMatrix[:, 3])))

    pngFile = os.path.join(qualPath, sequnitFileNamePrefix + ".r1_r2_baseposqual.png")
    htmlFile = os.path.join(qualPath, sequnitFileNamePrefix + ".r1_r2_baseposqual.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFile)

    ## Save Matplotlib plot in png format
    plt.savefig(pngFile, dpi=fig.dpi)

    return reformatQhistFile, pngFile, htmlFile


def gen_average_base_position_quality_boxplot(fastq, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatBqhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.bqhist.txt")  ## Average Base Position Quality Boxplot

    log.debug("qhist file: %s", reformatBqhistFile)

    ## Gen Average Base Position Quality Boxplot
    if not os.path.isfile(reformatBqhistFile):
        log.error("Bqhist file not found: %s", reformatBqhistFile)
        return None, None, None, None, None

    ## New data format (paired)
    #  0        1       2        3       4      5        6      7       8       9       10      11      12     13       14     15       16      17      18
    ##BaseNum    count_1 min_1   max_1   mean_1  Q1_1    med_1   Q3_1    LW_1    RW_1    count_2 min_2   max_2   mean_2  Q1_2    med_2   Q3_2    LW_2    RW_2
    #  0        6900    0       36      33.48   33      34      34      29      36      6900    0       36      33.48   33      34      34      29       36

    rawDataMatrix = np.loadtxt(reformatBqhistFile, delimiter='\t', comments='#', skiprows=0)
    assert len(rawDataMatrix[1][:]) == 19 or len(rawDataMatrix[1][:]) == 10

    bIsPaired = True

    if len(rawDataMatrix[1][:]) == 10:
        bIsPaired = False

    ## create data for boxplot
    boxplot_data_r1 = []
    boxplot_data_r2 = []

    for i in rawDataMatrix:
        idx = int(i[0]) - 1  ## read base loc

        spread = [rawDataMatrix[idx, 5], rawDataMatrix[idx, 7]]  # Q1 ~ Q3
        center = [rawDataMatrix[idx, 6]]  # median
        flier_high = [rawDataMatrix[idx, 8]]  # whisker lW 2%
        flier_low = [rawDataMatrix[idx, 9]]  # whisker rW 98%
        boxplot_data_r1.append(np.concatenate((spread, center, flier_high, flier_low), 0))

        if bIsPaired:
            spread = [rawDataMatrix[idx, 14], rawDataMatrix[idx, 16]]  # Q1 ~ Q3
            center = [rawDataMatrix[idx, 15]]  # median
            flier_high = [rawDataMatrix[idx, 17]]  # whisker lW 2%
            flier_low = [rawDataMatrix[idx, 18]]  # whisker rW 98%
            boxplot_data_r2.append(np.concatenate((spread, center, flier_high, flier_low), 0))

    fig, ax = plt.subplots()
    ax.boxplot(boxplot_data_r1)
    plt.subplots_adjust(left=0.06, right=0.9, top=0.9, bottom=0.1)

    F = plt.gcf()

    ## How to get the current size?
    ## DPI = F.get_dpi() ## = 80
    ## DefaultSize = F.get_size_inches() ## = (8, 6)
    F.set_size_inches(18, 6)  ## plot size (w, h)

    ax.set_xlabel("Read Position", fontsize=11, alpha=0.5)
    ax.set_ylabel("Quality Score (Solexa Scale: 40=Highest, -15=Lowest)", fontsize=11, alpha=0.5)

    ax.set_ylim([-20, 45])
    ax.yaxis.grid(True, linestyle=':', which="major")

    majorLocator_x = MultipleLocator(5)
    majorFormatter = FormatStrFormatter("%d")
    minorLocator = MultipleLocator(1)

    ax.xaxis.set_major_locator(majorLocator_x)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)

    majorLocator_y = MultipleLocator(5)
    minorLocator = MultipleLocator(1)
    ax.yaxis.set_major_locator(majorLocator_y)
    ax.yaxis.set_minor_locator(minorLocator)

    pngFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".r1_average_base_position_quality_boxplot.png")
    htmlFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".r1_average_base_position_quality_boxplot.html")

    ## Save D3 interactive plot in html format
    mpld3.save_html(fig, htmlFileR1)

    ## Save Matplotlib plot in png format
    plt.savefig(pngFileR1, dpi=fig.dpi)

    plotPngPlotFileR2 = None
    plotHtmlPlotFileR2 = None

    if bIsPaired:
        fig, ax = plt.subplots()
        ax.boxplot(boxplot_data_r2)
        plt.subplots_adjust(left=0.06, right=0.9, top=0.9, bottom=0.1)

        F = plt.gcf()

        ## How to get the current size?
        ## DPI = F.get_dpi() ## = 80
        ## DefaultSize = F.get_size_inches() ## = (8, 6)
        F.set_size_inches(18, 6)  ## plot size (w, h)

        ax.set_xlabel("Read Position", fontsize=11, alpha=0.5)
        ax.set_ylabel("Quality Score (Solexa Scale: 40=Highest, -15=Lowest)", fontsize=11, alpha=0.5)

        ax.set_ylim([-20, 45])
        ax.yaxis.grid(True, linestyle=':', which='major')

        majorLocator_x = MultipleLocator(5)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1)

        ax.xaxis.set_major_locator(majorLocator_x)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        majorLocator_y = MultipleLocator(5)
        minorLocator = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator_y)
        ax.yaxis.set_minor_locator(minorLocator)

        plotPngPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".r2_average_base_position_quality_boxplot.png")
        plotHtmlPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".r2_average_base_position_quality_boxplot.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, plotHtmlPlotFileR2)

        ## Save Matplotlib plot in png format
        plt.savefig(plotPngPlotFileR2, dpi=fig.dpi)


    return reformatBqhistFile, pngFileR1, plotPngPlotFileR2, htmlFileR1, plotHtmlPlotFileR2


def gen_cycle_nucleotide_composition_plot(fastq, readLength, isPairedEnd, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatBhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.bhist.txt")  ## base composition histogram

    log.debug("gen_cycle_nucleotide_composition_plot(): bhist file = %s", reformatBhistFile)

    ## Genbase composition histogram
    if not os.path.isfile(reformatBhistFile):
        log.error("Bhist file not found: %s", reformatBhistFile)
        return None, None, None, None, None

    ## data
    #  0    1   2   3   4   5
    ##Pos   A   C   G   T   N
    # 0  0.15111 0.26714 0.51707 0.06412 0.00056
    # 1  0.20822 0.20773 0.25543 0.32795 0.00068

    rawDataMatrix = np.loadtxt(reformatBhistFile, delimiter='\t', comments='#')
    assert len(rawDataMatrix[1][:]) == 6

    if isPairedEnd:
        rawDataR1 = rawDataMatrix[:readLength - 1][:]
        rawDataR2 = rawDataMatrix[readLength:][:]

        ## r1 ------------------------------------------------------------------
        fig, ax = plt.subplots()

        markerSize = 2.5
        lineWidth = 1.0
        p2 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='A', alpha=0.5)
        p3 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 4], 'g', marker='s', markersize=markerSize, linewidth=lineWidth, label='T', alpha=0.5)
        p4 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 3], 'b', marker='*', markersize=markerSize, linewidth=lineWidth, label='G', alpha=0.5)
        p5 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 2], 'm', marker='d', markersize=markerSize, linewidth=lineWidth, label='C', alpha=0.5)
        p6 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 5], 'c', marker='v', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        fontProp = FontProperties()
        fontProp.set_size("small")
        fontProp.set_family("Bitstream Vera Sans")
        ax.legend(loc=1, prop=fontProp)
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=list(rawDataR1[:, 1])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p3[0], labels=list(rawDataR1[:, 4])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p4[0], labels=list(rawDataR1[:, 3])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p5[0], labels=list(rawDataR1[:, 2])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p6[0], labels=list(rawDataR1[:, 5])))

        pngFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition_r1.png")
        htmlFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition_r1.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, htmlFileR1)

        ## Save Matplotlib plot in png format
        plt.savefig(pngFileR1, dpi=fig.dpi)

        ## r2 ------------------------------------------------------------------
        fig, ax = plt.subplots()

        markerSize = 2.5
        lineWidth = 1.0
        p2 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='A', alpha=0.5)
        p3 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 4], 'g', marker='s', markersize=markerSize, linewidth=lineWidth, label='T', alpha=0.5)
        p4 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 3], 'b', marker='*', markersize=markerSize, linewidth=lineWidth, label='G', alpha=0.5)
        p5 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 2], 'm', marker='d', markersize=markerSize, linewidth=lineWidth, label='C', alpha=0.5)
        p6 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 5], 'c', marker='v', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        fontProp = FontProperties()
        fontProp.set_size("small")
        fontProp.set_family("Bitstream Vera Sans")
        ax.legend(loc=1, prop=fontProp)
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=list(rawDataR2[:, 1])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p3[0], labels=list(rawDataR2[:, 4])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p4[0], labels=list(rawDataR2[:, 3])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p5[0], labels=list(rawDataR2[:, 2])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p6[0], labels=list(rawDataR2[:, 5])))

        plotPngPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition_r2.png")
        plotHtmlPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition_r2.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, plotHtmlPlotFileR2)

        ## Save Matplotlib plot in png format
        plt.savefig(plotPngPlotFileR2, dpi=fig.dpi)

        return reformatBhistFile, pngFileR1, htmlFileR1, plotPngPlotFileR2, plotHtmlPlotFileR2

    else:
        fig, ax = plt.subplots()

        markerSize = 2.5
        lineWidth = 1.0
        p2 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 1], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='A', alpha=0.5)
        p3 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 4], 'g', marker='s', markersize=markerSize, linewidth=lineWidth, label='T', alpha=0.5)
        p4 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 3], 'b', marker='*', markersize=markerSize, linewidth=lineWidth, label='G', alpha=0.5)
        p5 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 2], 'm', marker='d', markersize=markerSize, linewidth=lineWidth, label='C', alpha=0.5)
        p6 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 5], 'c', marker='v', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        fontProp = FontProperties()
        fontProp.set_size("small")
        fontProp.set_family("Bitstream Vera Sans")
        ax.legend(loc=1, prop=fontProp)
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p2[0], labels=list(rawDataMatrix[:, 1])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p3[0], labels=list(rawDataMatrix[:, 4])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p4[0], labels=list(rawDataMatrix[:, 3])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p5[0], labels=list(rawDataMatrix[:, 2])))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p6[0], labels=list(rawDataMatrix[:, 5])))

        pngFile = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition.png")
        htmlFile = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_nucl_composition.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, htmlFile)

        ## Save Matplotlib plot in png format
        plt.savefig(pngFile, dpi=fig.dpi)

        return reformatBhistFile, pngFile, htmlFile, None, None


def gen_cycle_n_base_percent_plot(fastq, readLength, isPairedEnd, log):
    READ_OUTPUT_PATH = RQCReadQcConfig.CFG["output_path"]

    sequnitFileName, _ = safe_basename(fastq, log)
    sequnitFileNamePrefix = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    subsampleDir = "subsample"
    subsamplePath = os.path.join(READ_OUTPUT_PATH, subsampleDir)

    qualDir = "qual"
    qualPath = os.path.join(READ_OUTPUT_PATH, qualDir)
    make_dir(qualPath)
    change_mod(qualPath, "0755")

    reformatBhistFile = os.path.join(subsamplePath, sequnitFileNamePrefix + ".reformat.bhist.txt")  ## base composition histogram

    log.debug("gen_cycle_n_base_percent_plot(): bhist file = %s", reformatBhistFile)

    ## Genbase composition histogram
    if not os.path.isfile(reformatBhistFile):
        log.error("Bhist file not found: %s", reformatBhistFile)
        return None, None, None, None, None

    ## data
    #  0    1   2   3   4   5
    ##Pos   A   C   G   T   N
    # 0  0.15111 0.26714 0.51707 0.06412 0.00056
    # 1  0.20822 0.20773 0.25543 0.32795 0.00068

    rawDataMatrix = np.loadtxt(reformatBhistFile, delimiter='\t', comments='#')
    assert len(rawDataMatrix[1][:]) == 6

    if isPairedEnd:
        rawDataR1 = rawDataMatrix[:readLength - 1][:]
        rawDataR2 = rawDataMatrix[readLength:][:]

        ## r1 ------------------------------------------------------------------
        fig, ax = plt.subplots()

        markerSize = 5.0
        lineWidth = 1.5
        p1 = ax.plot(rawDataR1[:, 0], rawDataR1[:, 5], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        ax.grid(color="gray", linestyle=':')
        ax.set_xlim([0, readLength])

        ## Add tooltip
        labels = ["%.5f" % i for i in rawDataR1[:, 5]]
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=labels))

        pngFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent_r1.png")
        htmlFileR1 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent_r1.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, htmlFileR1)

        ## Save Matplotlib plot in png format
        plt.savefig(pngFileR1, dpi=fig.dpi)

        ## r2 ------------------------------------------------------------------
        fig, ax = plt.subplots()

        markerSize = 5.0
        lineWidth = 1.5
        p1 = ax.plot([(x - readLength) for x in rawDataR2[:, 0]], rawDataR2[:, 5], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        ax.grid(color="gray", linestyle=':')
        ax.set_xlim([0, readLength])

        ## Add tooltip
        labels = ["%.5f" % i for i in rawDataR2[:, 5]]
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=labels))

        plotPngPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent_r2.png")
        plotHtmlPlotFileR2 = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent_r2.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, plotHtmlPlotFileR2)

        ## Save Matplotlib plot in png format
        plt.savefig(plotPngPlotFileR2, dpi=fig.dpi)

        return reformatBhistFile, pngFileR1, htmlFileR1, plotPngPlotFileR2, plotHtmlPlotFileR2

    else:
        fig, ax = plt.subplots()

        markerSize = 5.0
        lineWidth = 1.5
        p1 = ax.plot(rawDataMatrix[:, 0], rawDataMatrix[:, 5], 'r', marker='o', markersize=markerSize, linewidth=lineWidth, label='N', alpha=0.5)

        ax.set_xlabel("Read Position", fontsize=12, alpha=0.5)
        ax.set_ylabel("Fraction", fontsize=12, alpha=0.5)
        ax.grid(color="gray", linestyle=':')

        ## Add tooltip
        labels = ["%.5f" % i for i in rawDataMatrix[:, 5]]
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(p1[0], labels=labels))

        pngFile = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent.png")
        htmlFile = os.path.join(qualPath, sequnitFileNamePrefix + ".cycle_n_base_percent.html")

        ## Save D3 interactive plot in html format
        mpld3.save_html(fig, htmlFile)

        ## Save Matplotlib plot in png format
        plt.savefig(pngFile, dpi=fig.dpi)

        return reformatBhistFile, pngFile, htmlFile, None, None


""" For STEP12
sciclone_sam2summary

@param sam_file: input sam file
@param count_file: stat file for writing

@return retCode: success or failure
@return count_file: output count file

"""
## Removed!
##def sciclone_sam2summary(sam_file, log):



""" For STEP12
run_rna_strandedness

 Title:        run_rna_strandedness
 Function:     Takes sam file generated from rna data set and determines
               mapping to sense and antisense strand
 Usage:        run_rna_strandedness($sam_file, $log)
 Args:         1) sam file
               2) log file object
 Returns:      JGI_SUCCESS
               JGI_FAILURE
 Comments:     None.

@param sam_file

@return retCode
@return outputLogFile: log file
@return outResultDict: stat value dict (to be added to readqc_stats.txt)

"""
## Removed!



""" For STEP12
separate_paired_end_fq

 Title      : separate_paired_end_fq
 Function   : Given a fastq file, this function splits the fastq into
              read1 and read2 fastq file.
 Usage      : sequence_lengths( $fastq, $read1_fq, $read2_fq, $log )
 Args       : 1) The name of a fastq sequence file.
              2) The name of the read1 fastq output file
              3) The name of the read2 fastq output file
              4) A JGI_Log object.
 Returns    : JGI_SUCCESS: The fastq file was successfully separated.
              JGI_FAILURE: The fastq file could not be separated.
 Comments   : For paired end fastq files only.


 sulsj

 - Added gzip'd fastq support
 ## TODO: need to write a fastq IO class

@param fastq
@param read1_outfile
@param read2_outfile
@param log

@return retCode

"""
## Removed!



""" For STEP12

NOTE: Borrowed from alignment.py and updated to get lib name and isRna

get_lib_info

Look up the seqProjectId, bio_name from the library_info table
- to look up the references

@param sequnitFileName: name of the seq unit in the seq_units.sequnitFileName field

@return seqProjectId
@return ncbiOrganismName
@return libraryName
@return isRna

"""

""" For STEP14 & STEP15
run_blastplus_py

Call jgi-rqc-pipeline/tools/run_blastplus.py

"""
def run_blastplus_py(queryFastaFile, db, log):
    timeoutCmd = 'timeout'  #RQCReadQcCommands.TIMEOUT_CMD
    blastOutFileNamePrefix = None
    # retCode = None

    outDir, exitCode = safe_dirname(queryFastaFile, log)
    queryFastaFileBaseName, exitCode = safe_basename(queryFastaFile, log)
    dbFileBaseName, exitCode = safe_basename(db, log)

    # runBlastnCmd = "/global/homes/s/sulsj/work/bitbucket-repo/jgi-rqc-pipeline/tools/run_blastplus.py" ## debug
    # runBlastnCmd = "/global/dna/projectdirs/PI/rqc/prod/jgi-rqc-pipeline/tools/run_blastplus.py"
    # runBlastnCmd = "/global/homes/s/sulsj/work/bitbucket-repo/jgi-rqc-pipeline/tools/run_blastplus_taxserver.py" ## debug
    runBlastnCmd = "/global/dna/projectdirs/PI/rqc/prod/jgi-rqc-pipeline/tools/run_blastplus_taxserver.py"

    blastOutFileNamePrefix = outDir + "/megablast." + queryFastaFileBaseName + ".vs." + dbFileBaseName

    ## Should use  jigsaw/2.6.0 for not checking database reference fasta file
    ## 07212016 Added -s to add lineage to the subject field
    cmd = "%s 21600s %s -d %s -o %s -q %s -s > %s.log 2>&1 " % (timeoutCmd, runBlastnCmd, db, outDir, queryFastaFile, blastOutFileNamePrefix)

    _, _, exitCode = run_sh_command(cmd, True, log, True)

    ## Added timeout to terminate blast run manually after 6hrs
    ## If exitCode == 124 or exitCode = 143, this means the process exits with timeout.
    ## Timeout exits with 128 plus the signal number. 143 = 128 + 15 (SGITERM)
    ## Ref) http://stackoverflow.com/questions/4189136/waiting-for-a-command-to-return-in-a-bash-script
    ##      timeout man page ==> If the command times out, and --preserve-status is not set, then exit with status 124.
    ##
    if exitCode in (124, 143):
        ## BLAST timeout
        ## Exit with success so that the blast step can be skipped.
        log.warning("##################################")
        log.warning("BLAST TIMEOUT. JUST SKIP THE STEP.")
        log.warning("##################################")
        return RQCExitCodes.JGI_FAILURE, -143

    elif exitCode != 0:
        log.error("Failed to run_blastplus_py. Exit code != 0")
        return RQCExitCodes.JGI_FAILURE, None

    else:
        log.info("run_blastplus_py complete.")

    return RQCExitCodes.JGI_SUCCESS, blastOutFileNamePrefix


"""===========================================================================
    checkpoint_step_wrapper

"""
def checkpoint_step_wrapper(status):
    assert RQCReadQcConfig.CFG["status_file"]
    checkpoint_step(RQCReadQcConfig.CFG["status_file"], status)


"""===========================================================================
    get the file content
"""
def get_analysis_file_contents(fullPath):
    retCode = ""

    if os.path.isfile(fullPath):
        with open(fullPath, "r") as FH:
            retCode = FH.readlines()
        return retCode
    else:
        return "file not found"


'''===========================================================================
    file_name_trim

'''
def file_name_trim(fname):
    return fname.replace(".gz", "").replace(".fastq", "").replace(".fasta", "")


## EOF
