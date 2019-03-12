package hiseq;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.Dedupe;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.ScheduleMaker;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Analyzes a flow cell for low-quality areas.
 * Removes reads in the low-quality areas.
 * 
 * @author Brian Bushnell
 * @date August 31, 2016
 *
 */
public class AnalyzeFlowCell {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		AnalyzeFlowCell x=new AnalyzeFlowCell(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeFlowCell(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		Parser parser=parse(args);
		
		if(gToN || discardG){MicroTile.TRACK_CYCLES=true;}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
			

			trimq=parser.trimq;
			trimE=parser.trimE();
			minlen=parser.minReadLength;
			trimLeft=parser.qtrimLeft;
			trimRight=parser.qtrimRight;
		}
		
		checkFiles();
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutbad=FileFormat.testOutput(outbad, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	private Parser parse(String[] args){
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimRight=trimRight;
		parser.trimq=trimq;
		parser.minReadLength=minlen;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("divisor") || a.equals("size")){
				Tile.xSize=Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("xdivisor") || a.equals("xsize")){
				Tile.xSize=Integer.parseInt(b);
			}else if(a.equals("ydivisor") || a.equals("ysize")){
				Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("target")){
				targetAverageReads=Integer.parseInt(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("indump") || a.equals("ind") || a.equals("dumpin")){
				dumpIn=b;
			}else if(a.equals("pound")){
				pound=Tools.parseBoolean(b);
			}else if(a.equals("loadkmers") || a.equals("usekmers")){
				loadKmers=Tools.parseBoolean(b);
			}else if(a.equals("lqo") || a.equals("lowqualityonly")){
				discardOnlyLowQuality=Tools.parseBoolean(b);
			}else if(a.equals("dl") || a.equals("discardlevel")){
				discardLevel=Integer.parseInt(b);
			}else if(a.equals("outbad") || a.equals("outb") || a.equals("outtoss") || a.equals("outt") || a.equals("outunwanted") || a.equals("outu")){
				outbad=b;
			}
			
			else if(a.equals("deviations") || a.equals("d")){
				qDeviations=uDeviations=eDeviations=gDeviations=Float.parseFloat(b);
			}else if(a.equals("qdeviations") || a.equals("qd") || a.equals("dq")){
				qDeviations=Float.parseFloat(b);
			}else if(a.equals("udeviations") || a.equals("ud") || a.equals("du")){
				uDeviations=Float.parseFloat(b);
			}else if(a.equals("edeviations") || a.equals("ed") || a.equals("de")){
				eDeviations=Float.parseFloat(b);
			}else if(a.equals("gdeviations") || a.equals("gd") || a.equals("dg")){
				gDeviations=Float.parseFloat(b);
			}
			
			else if(a.equals("qfraction") || a.equals("qf")){
				qualFraction=Float.parseFloat(b);
			}else if(a.equals("ufraction") || a.equals("uf")){
				uniqueFraction=Float.parseFloat(b);
			}else if(a.equals("efraction") || a.equals("ef")){
				errorFreeFraction=Float.parseFloat(b);
			}else if(a.equals("gfraction") || a.equals("gf")){
				gFraction=Float.parseFloat(b);
			}
			
			else if(a.equals("qabsolute") || a.equals("qa")){
				qualAbs=Float.parseFloat(b);
			}else if(a.equals("uabsolute") || a.equals("ua")){
				uniqueAbs=Float.parseFloat(b);
			}else if(a.equals("eabsolute") || a.equals("ea")){
				errorFreeAbs=Float.parseFloat(b);
			}else if(a.equals("gabsolute") || a.equals("ga")){
				gAbs=Float.parseFloat(b);
			}
			
			else if(a.equals("gton")){
				gToN=Tools.parseBoolean(b);
			}else if(a.equals("discardg")){
				discardG=Tools.parseBoolean(b);
			}
			
			else if(a.equals("minpolyg")){
				MicroTile.MIN_POLY_G=Integer.parseInt(b);
			}else if(a.equals("trackcycles")){
				MicroTile.TRACK_CYCLES=Tools.parseBoolean(b);
			}
			
			else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		return parser;
	}
	
	private void checkFiles(){
		doPoundReplacement();
		adjustInterleaving();
		checkFileExistence();
		checkStatics();
	}
	
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	private void checkFileExistence(){
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outbad, dump)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outbad);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outbad+", "+dump+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outbad, dump)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
	}
	
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		if(dumpIn==null){
			flowcell=new FlowCell();
			if(loadKmers){loadKmers();}
			fillTiles();
			keySets=null;
		}else{
			flowcell=new FlowCell(dumpIn);
			
			if(flowcell.avgReads<targetAverageReads){
				flowcell=flowcell.widen(targetAverageReads);
			}
			
			avgQuality=flowcell.avgQuality;
			avgUnique=flowcell.avgUnique;
			avgErrorFree=flowcell.avgErrorFree;
			avgG=flowcell.avgG;
			stdQuality=flowcell.stdQuality;
			stdUnique=flowcell.stdUnique;
			stdErrorFree=flowcell.stdErrorFree;
			stdG=flowcell.stdG;

			long readsToDiscard=markTiles(flowcell.toList(), flowcell.avgReads);
		}
		processReads(t);
	}

	/** Create read streams and process all data */
	void loadKmers(){
		Timer t2=new Timer();
		outstream.print("Loading kmers:  \t");
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		loadKmersInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		t2.stop();
		outstream.println(t2);
	}

	/** Create read streams and process all data */
	void fillTiles(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		fillTilesInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
	}

	/** Create read streams and process all data */
	void processReads(Timer t){
		
		if(ffout1!=null || ffoutbad!=null){
			Timer t2=new Timer();
			outstream.print("Filtering reads:\t");

			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			boolean paired=cris.paired();
			//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

			//Optionally read output streams
			final ConcurrentReadOutputStream ros, rosb;
			final int buff=4;
			if(ffout1!=null){
				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
				ros.start(); //Start the stream
			}else{ros=null;}
			
			if(ffoutbad!=null){
				rosb=ConcurrentReadOutputStream.getStream(ffoutbad, null, null, null, buff, null, false);
				rosb.start(); //Start the stream
			}else{rosb=null;}

			//Process the read stream
			processInner(cris, ros, rosb);

			if(verbose){outstream.println("Finished; closing streams.");}

			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris, ros, rosb);

			t2.stop();
			outstream.println(t2);
		}
		
		//Report timing and results
		{
			t.stop();
			lastReadsOut=readsProcessed-readsDiscarded;
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			
			if(ffout1!=null || ffoutbad!=null){

				String rpstring=Tools.padKM(readsDiscarded, 8);
				String bpstring=Tools.padKM(basesDiscarded, 8);
				String gpstring=Tools.padKM(gsTransformedToN, 8);
				outstream.println();
				outstream.println("Reads Discarded:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.3f%%", readsDiscarded*100.0/readsProcessed));
				outstream.println("Bases Discarded:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.3f%%", basesDiscarded*100.0/basesProcessed));
				if(gToN){outstream.println("Gs Masked By N:     "+gpstring+" \t"+String.format(Locale.ROOT, "%.3f%%", gsTransformedToN*100.0/basesProcessed));}
				outstream.println();
			}
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		readsProcessed=0;
		basesProcessed=0;
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				ArrayList<Read> keepList=new ArrayList<Read>(reads.size());
				ArrayList<Read> tossList=new ArrayList<Read>(4);
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					boolean keep=processReadPair(r1, r2);
					if(keep){
						keepList.add(r1);
					}else{
						tossList.add(r1);
						readsDiscarded+=r1.pairCount();
						basesDiscarded+=initialLength1+initialLength2;
					}
				}
				
				//Output reads to the output stream
				if(ros!=null){ros.add(keepList, ln.id);}
				if(rosb!=null){rosb.add(tossList, ln.id);}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/** Iterate through the reads */
	public void loadKmersInner(final ConcurrentReadInputStream cris){
		
		keySets=new AbstractKmerTable[WAYS];

		//Initialize tables
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 0.8);
		int[] schedule=scheduleMaker.makeSchedule();
		for(int j=0; j<WAYS; j++){
			keySets[j]=new HashArray1D(schedule, -1L);
		}
//		for(int j=0; j<WAYS; j++){
//			keySets[j]=new HashArray1D(512000, -1, -1L, true); //TODO: Set maxSize
//		}
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					if(initialLength1>=k && randy.nextBoolean()){
						final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer, 1);
						}
					}
					
					if(initialLength2>=k && randy.nextBoolean()){
						final long kmer=toKmer(r2.bases, r2.quality, randy.nextInt(initialLength2-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer, 1);
						}
					}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
	}
	
	/** Iterate through the reads */
	public void fillTilesInner(final ConcurrentReadInputStream cris){


		Timer t2=new Timer();
		outstream.print("Filling tiles:  \t");
		
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					final MicroTile mt=flowcell.getMicroTile(r1.id);
					
					if(loadKmers){
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
						
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, r1.quality, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
					}

					mt.add(r1);
					mt.add(r2);
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		t2.stop();
		outstream.println(t2);
		
		ArrayList<MicroTile> mtList=flowcell.calcStats();
		if(flowcell.avgReads<targetAverageReads){
			flowcell=flowcell.widen(targetAverageReads);
			mtList=flowcell.toList();
		}
		avgQuality=flowcell.avgQuality;
		avgUnique=flowcell.avgUnique;
		avgErrorFree=flowcell.avgErrorFree;
		avgG=flowcell.avgG;
		
		stdQuality=flowcell.stdQuality;
		stdUnique=flowcell.stdUnique;
		stdErrorFree=flowcell.stdErrorFree;
		stdG=flowcell.stdG;
		
		long readsToDiscard=markTiles(mtList, flowcell.avgReads);
		
		if(dump!=null){
			TextStreamWriter tsw=new TextStreamWriter(dump, overwrite, append, false);
			tsw.start();

			tsw.println("#xSize\t"+Tile.xSize);
			tsw.println("#ySize\t"+Tile.ySize);
			tsw.println("#reads\t"+String.format(Locale.ROOT, "%d", flowcell.readsProcessed));
			tsw.println("#avgReads\t"+String.format(Locale.ROOT, "%.1f", flowcell.avgReads));
			
			tsw.println("#avgQuality\t"+String.format(Locale.ROOT, "%.3f", avgQuality));
			tsw.println("#avgUnique\t"+String.format(Locale.ROOT, "%.3f", avgUnique));
			tsw.println("#avgErrorFree\t"+String.format(Locale.ROOT, "%.3f", avgErrorFree));
			tsw.println("#avgG\t"+String.format(Locale.ROOT, "%.3f", avgG));
			
			tsw.println("#stdQuality\t"+String.format(Locale.ROOT, "%.5f", stdQuality));
			tsw.println("#stdUnique\t"+String.format(Locale.ROOT, "%.5f", stdUnique));
			tsw.println("#stdErrorFree\t"+String.format(Locale.ROOT, "%.5f", stdErrorFree));
			tsw.println("#stdG\t"+String.format(Locale.ROOT, "%.5f", stdG));
			
			tsw.println((pound ? "#" : "") + "lane\ttile\tx1\tx2\ty1\ty2\treads\tunique\tquality\tprobErrorFree\tdiscard");
			
			for(Lane lane : flowcell.lanes){
				if(lane!=null){
					for(Tile tile : lane.tiles){
						if(tile!=null){
							tsw.print(tile.toString());
						}
					}
				}
			}
			tsw.poisonAndWait();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean processReadPair(final Read r1, final Read r2){
		boolean passes=processReadPair_inner(r1, r2);
		if(passes){return true;}
		if(trimq>0){
			TrimRead.trimFast(r1, trimLeft, trimRight, trimq, trimE, 0);
			if(r2!=null){TrimRead.trimFast(r2, trimLeft, trimRight, trimq, trimE, 0);}
			return r1.length()>=minlen && (r2==null || r2.length()>=minlen);
		}else{
			return false;
		}
	}
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	boolean processReadPair_inner(final Read r1, final Read r2){
		
		MicroTile mt=flowcell.getMicroTile(r1.id);
		if(mt==null){
			if(!warned){
				outstream.println("\nWarning - a read was found with no corresponding MicroTile:\n"+r1.id);
				warned=true;
			}
			return true;
		}
		if(mt.discard<discardLevel){return true;}
		if(!discardOnlyLowQuality){return false;}
		int len1=r1.length(), len2=r1.mateLength();
		if(len1>0){
			double qual=r1.avgQualityByProbabilityDouble(true, len1);
			double prob=100*r1.probabilityErrorFree(true, len1);
			if(qual<avgQuality-qDeviations*stdQuality){return false;}
			if(prob<avgErrorFree-eDeviations*stdErrorFree){return false;}
			if(discardG && shouldDiscardG(r1, mt)){return false;}
			if(gToN){gsTransformedToN+=doGToN(r1, mt);}
		}
		if(len2>0){
			double qual=r2.avgQualityByProbabilityDouble(true, len2);
			double prob=100*r2.probabilityErrorFree(true, len2);
			if(qual<avgQuality-qDeviations*stdQuality){return false;}
			if(prob<avgErrorFree-eDeviations*stdErrorFree){return false;}
			if(discardG && shouldDiscardG(r2, mt)){return false;}
			if(gToN){gsTransformedToN+=doGToN(r2, mt);}
		}
		return true;
	}
	
	private boolean shouldDiscardG(Read r, MicroTile mt){
		final byte[] bases=r.bases;
		final float[] gArray=mt.tracker.cycleAverages[2];
		
		final float thresh=(float)(avgG+Tools.max(gDeviations*stdG, avgG*gFraction, gAbs));
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='G' && gArray[i]>thresh){
				return true;
			}
		}
		return false;
	}
	
	private int doGToN(Read r, MicroTile mt){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final float[] gArray=mt.tracker.cycleAverages[2];
		
		final float thresh=(float)(avgG+Tools.max(gDeviations*stdG, avgG*gFraction, gAbs));
		int changes=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='G' && gArray[i]>thresh){
				bases[i]='N';
				changes++;
				if(quals!=null){quals[i]=0;}
			}
		}
		return changes;
	}
		
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private static final long toKmer(final byte[] bases, byte[] quals, final int start, final int klen){
//		if(minprob>0 && quals!=null){
//			float prob=toProb(quals, start, klen);
//			if(prob<minprob){return -1;}
//		}
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	private final long markTiles(ArrayList<MicroTile> mtList, double avgReads){
		for(MicroTile mt : mtList){
			mt.discard=0;
		}
		long readsToDiscard=0;
		
		cDiscards=qDiscards=eDiscards=uDiscards=mtDiscards=gDiscards=mtRetained=0;
		
		for(MicroTile mt : mtList){
			double q=mt.averageQuality();
			double e=mt.percentErrorFree();
			double u=mt.uniquePercent();
			double g=mt.maxG();

			double dq=avgQuality-q;
			double de=avgErrorFree-e;
			double du=u-avgUnique;
			double dg=g-avgG;
			
			if(mt.readCount<10 && mt.readCount<0.02f*avgReads){
				mt.discard++;
				cDiscards++;
			}
			
			if(dq>qDeviations*stdQuality && dq>avgQuality*qualFraction && dq>qualAbs){
				mt.discard++;
				qDiscards++;
			}
			if(de>eDeviations*stdErrorFree && de>avgErrorFree*errorFreeFraction && de>errorFreeAbs){
				mt.discard++;
				eDiscards++;
			}
			if(avgUnique>2 && avgUnique<98){
				if(du>uDeviations*stdUnique && du>avgUnique*uniqueFraction && du>uniqueAbs){
					mt.discard++;
					uDiscards++;
				}
			}
			
			//Or mode
//			if(dq>qDeviations*stdQuality || dq>avgQuality*qualFraction || dq>qualAbs){
//				mt.discard++;
//				qDiscards++;
//			}
//			if(de>eDeviations*stdErrorFree || de>avgErrorFree*errorFreeFraction || de>errorFreeAbs){
//				mt.discard++;
//				eDiscards++;
//			}
//			if(avgUnique>2 && avgUnique<98){
//				if(du>uDeviations*stdUnique || du>avgUnique*uniqueFraction || du>uniqueAbs){
//					mt.discard++;
//					uDiscards++;
//				}
//			}
			
			if((discardG || gToN) && (dg>gDeviations*stdG && dg>avgG*gFraction && dg>gAbs)){
				mt.discard++;
				gDiscards++;
			}
			if(mt.discard>0){
				mtDiscards++;
				readsToDiscard+=mt.readCount;
			}
			else{mtRetained++;}
		}
		
		outstream.println();
		outstream.println("Flagged "+mtDiscards+" of "+(mtDiscards+mtRetained)+" micro-tiles, containing "+readsToDiscard+" reads:");
		outstream.println(uDiscards+" exceeded uniqueness thresholds.");
		outstream.println(qDiscards+" exceeded quality thresholds.");
		outstream.println(eDiscards+" exceeded error-free probability thresholds.");
		outstream.println(gDiscards+" contained G spikes.");
		outstream.println(cDiscards+" had too few reads to calculate statistics.");
		outstream.println();
		
		return readsToDiscard;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Discard output file path */
	private String outbad=null;

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private boolean pound=true;
	private String dump=null;
	private String dumpIn=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;

	/** Number of reads discarded */
	public long readsDiscarded=0;
	/** Number of bases discarded */
	public long basesDiscarded=0;

	protected long cDiscards=0;
	protected long uDiscards=0;
	protected long qDiscards=0;
	protected long eDiscards=0;
	protected long gDiscards=0;
	protected long mtDiscards=0;
	protected long mtRetained=0;
	
	protected long gsTransformedToN=0; 
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private AbstractKmerTable[] keySets;
	
	private int targetAverageReads=800;

	private float minprob=0;
	private static final int WAYS=31;
	private static final int k=31, k2=30;
	
	private final Random randy=Shared.threadLocalRandom();
	private FlowCell flowcell;
	
//	private int[] avgQualityArray;
//	private int[] avgUniqueArray;
	
	private long minCountToUse=0;

	private float qDeviations=2f;
	private float uDeviations=1.5f;
	private float eDeviations=2f;
	private float gDeviations=2f;
	
	private float qualFraction=0.01f;
	private float uniqueFraction=0.01f;
	private float errorFreeFraction=0.01f;
	private float gFraction=0.1f;
	
	private float qualAbs=1f;
	private float uniqueAbs=1f;
	private float errorFreeAbs=1f;
	private float gAbs=0.05f;
	
	private double avgQuality;
	private double avgUnique;
	private double avgErrorFree;
	private double avgG;
	
	private double stdQuality;
	private double stdUnique;
	private double stdErrorFree;
	private double stdG;

	private boolean loadKmers=true;
	private boolean discardOnlyLowQuality=true;
	private int discardLevel=1;
	private boolean gToN=false;
	private boolean discardG=false;
	
	private int minlen=30;
	private float trimq=-1;
	private final float trimE;
	private boolean trimLeft=false;
	private boolean trimRight=true;
	
	private boolean warned=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	/** Output for discarded reads */
	private final FileFormat ffoutbad;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of reads output in the last run */
	public static long lastReadsOut;
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
