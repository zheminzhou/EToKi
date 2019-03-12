package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import stream.SamStreamer;
import structures.ListNum;

/**
 * Removes lines with unsupported variations from a sam file.
 * 
 * @author Brian Bushnell
 * @date January 25, 2018
 *
 */
public class FilterSam {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		FilterSam x=new FilterSam(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public FilterSam(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		SamLine.SET_FROM_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}else if(a.equals("ss") || a.equals("streamer")){
				useStreamer=Tools.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("outbad") || a.equals("outb")){
				outBad=b;
			}else if(a.equals("vars") || a.equals("variants") || a.equals("varfile") || a.equals("inv")){
				varFile=b;
			}else if(a.equals("vcf") || a.equals("vcffile")){
				vcfFile=b;
			}else if(a.equals("maxbadsubs") || a.equals("maxbadbars")){
				maxBadSubs=Integer.parseInt(b);
			}else if(a.equals("maxbadsubdepth") || a.equals("maxbadvardepth") || a.equals("maxbadsuballeledepth") || a.equals("maxbadvaralleledepth") || a.equals("mbsad")){
				maxBadSubAlleleDepth=Integer.parseInt(b);
			}else if(a.equals("minbadsubreaddepth") || a.equals("minbadvarreaddepth") || a.equals("mbsrd")){
				minBadSubReadDepth=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			outGood=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null || (vcfFile==null && varFile==null)){throw new RuntimeException("Error - an input file and a VCF file are required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGood, outBad)){
			outstream.println((outGood==null)+", "+(outBad==null));
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outGood+", "+outBad+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, vcfFile, varFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, vcfFile, varFile, outBad, outGood)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffoutGood=FileFormat.testOutput(outGood, FileFormat.SAM, null, true, overwrite, append, ordered);
		ffoutBad=FileFormat.testOutput(outBad, FileFormat.SAM, null, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		
		Timer t=new Timer(outstream, true);
		{
			t.start();
			outstream.print("Loading scaffolds:  ");
			if(ref==null){
				scafMap=ScafMap.loadSamHeader(in1);
			}else{
				scafMap=ScafMap.loadReference(ref, true);
			}
			assert(scafMap!=null && scafMap.size()>0) : "No scaffold names were loaded.";
			t.stop(pad(scafMap.size(), 12)+" \t");
		}
		t.start();
		if(varFile!=null){
			outstream.print("Loading vars:       ");
			varMap=VcfLoader.loadVars(varFile, scafMap);
			t.stop(pad(varMap.size(), 12)+" \t");
		}else if(vcfFile!=null){t.start();
			outstream.print("Loading vcf:        ");
			varMap=VcfLoader.loadVcf(vcfFile, scafMap, true, false);
			t.stop(pad(varMap.size(), 12)+" \t");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		
		final SamReadStreamer ss;
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		if(useStreamer && (maxReads<0 || maxReads==Long.MAX_VALUE)){
			cris=null;
			ss=new SamReadStreamer(ffin1, streamerThreads, true);
			ss.start();
			if(verbose){outstream.println("Started streamer");}
		}else{
			ss=null;
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros, rosb;
		if(ffoutGood!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			ros=ConcurrentReadOutputStream.getStream(ffoutGood, null, null, null, buff, null, false);
			ros.start(); //Start the stream
		}else{ros=null;}
		if(ffoutBad!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutBad, null, null, null, buff, null, false);
			rosb.start(); //Start the stream
		}else{rosb=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ss, ros, rosb);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			final long rg=readsOut, rb=readsProcessed-readsOut;
			final long bg=basesOut, bb=basesProcessed-basesOut;
			
			//Add "k" and "m" for large numbers
//			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
//			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
//			String rgstring=(rg<100000 ? ""+rg : rg<100000000 ? (rg/1000)+"k" : (rg/1000000)+"m");
//			String bgstring=(bg<100000 ? ""+bg : bg<100000000 ? (bg/1000)+"k" : (bg/1000000)+"m");
//			String rbstring=(rb<100000 ? ""+rb : rb<100000000 ? (rb/1000)+"k" : (rb/1000000)+"m");
//			String bbstring=(bb<100000 ? ""+bb : bb<100000000 ? (bb/1000)+"k" : (bb/1000000)+"m");
			
			final int len=12;
			
			final String rpstring=pad(readsProcessed, len);
			final String bpstring=pad(basesProcessed, len);
			final String rgstring=pad(rg, len);
			final String bgstring=pad(bg, len);
			final String rbstring=pad(rb, len);
			final String bbstring=pad(bb, len);
			
			final double mappedReadsDiscarded=mappedReadsProcessed-mappedReadsRetained;
			double avgQGood=(qSumGood/(double)mappedReadsRetained);
			double avgQBad=(qSumBad/(double)mappedReadsDiscarded);
			double avgMapQGood=(mapqSumGood/(double)mappedReadsRetained);
			double avgMapQBad=(mapqSumBad/(double)mappedReadsDiscarded);
			double avgSubsGood=(subSumGood/(double)mappedReadsRetained);
			double avgSubsBad=(subSumBad/(double)mappedReadsDiscarded);
			final String avgQGoodS=pad(String.format(Locale.ROOT, "%.2f", avgQGood), len);
			final String avgQBadS=pad(String.format(Locale.ROOT, "%.2f", avgQBad), len);
			final String avgMapQGoodS=pad(String.format(Locale.ROOT, "%.2f", avgMapQGood), len);
			final String avgMapQBadS=pad(String.format(Locale.ROOT, "%.2f", avgMapQBad), len);
			final String avgSubsGoodS=pad(String.format(Locale.ROOT, "%.2f", avgSubsGood), len);
			final String avgSubsBadS=pad(String.format(Locale.ROOT, "%.2f", avgSubsBad), len);
			
			outstream.println("Time:                         \t\t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
			outstream.println();
			outstream.println("Reads Retained:     "+rgstring+" \t"+String.format(Locale.ROOT, "%.2f%%", (rg*100.0/readsProcessed)));
			outstream.println("Bases Retained:     "+bgstring+" \t"+String.format(Locale.ROOT, "%.2f%%", (bg*100.0/basesProcessed)));
			outstream.println("Avg. Qual Retained: "+avgQGoodS);
			outstream.println("Avg. MapQ Retained: "+avgMapQGoodS);
			outstream.println("Avg. Subs Retained: "+avgSubsGoodS);
			outstream.println();
			outstream.println("Reads Discarded:    "+rbstring+" \t"+String.format(Locale.ROOT, "%.2f%%", (rb*100.0/readsProcessed)));
			outstream.println("Bases Discarded:    "+bbstring+" \t"+String.format(Locale.ROOT, "%.2f%%", (bb*100.0/basesProcessed)));
			outstream.println("Avg. Qual Discarded:"+avgQBadS);
			outstream.println("Avg. MapQ Discarded:"+avgMapQBadS);
			outstream.println("Avg. Subs Discarded:"+avgSubsBadS);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private static String pad(long s, int len){
		return pad(""+s, len);
	}
	
	private static String pad(String s, int len){
		while(s.length()<len){s=" "+s;}
		return s;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final SamStreamer ss, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ss, ros, rosb, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			mappedReadsProcessed+=pt.mappedReadsProcessedT;
			mappedBasesProcessed+=pt.mappedBasesProcessedT;
			mappedReadsRetained+=pt.mappedReadsRetainedT;
			mappedBasesRetained+=pt.mappedBasesRetainedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			qSumGood+=pt.qSumGoodT;
			qSumBad+=pt.qSumBadT;
			subSumGood+=pt.subSumGoodT;
			subSumBad+=pt.subSumBadT;
			mapqSumGood+=pt.mapqSumGoodT;
			mapqSumBad+=pt.mapqSumBadT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final SamStreamer ss_, final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosb_, final int tid_){
			cris=cris_;
			ss=ss_;
			ros=ros_;
			rosb=rosb_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			if(useStreamer){
				processSS();
			}else{
				processCris();
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processCris(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				ArrayList<Read> good=new ArrayList<Read>(reads.size());
				ArrayList<Read> bad=new ArrayList<Read>(Tools.max(1, reads.size()/4));
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength1;
					
					{
						//Reads are processed in this block.
						boolean keep=processRead(r1);
						if(keep){
							//Increment counters
							readsOutT++;
							basesOutT+=initialLength1;
							good.add(r1);
						}else{
							bad.add(r1);
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(good, ln.id);}
				if(rosb!=null){rosb.add(bad, ln.id);}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/** Iterate through the reads */
		void processSS(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=ss.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				ArrayList<Read> good=new ArrayList<Read>(reads.size());
				ArrayList<Read> bad=new ArrayList<Read>(Tools.max(1, reads.size()/4));
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength1;
					
					{
						//Reads are processed in this block.
						boolean keep=processRead(r1);
						if(keep){
							//Increment counters
							readsOutT++;
							basesOutT+=initialLength1;
							good.add(r1);
						}else{
							bad.add(r1);
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(good, ln.id);}
				if(rosb!=null){rosb.add(bad, ln.id);}

				//Fetch a new list
				ln=ss.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Process a read.
		 * @param r1 Read 1
		 * @return True if the read should be kept, false if it should be discarded.
		 */
		boolean processRead(final Read r1){
			return passesVariantFilter(r1);
		}
		
		private final boolean passesVariantFilter(Read r){
			if(!r.mapped() || r.bases==null || r.obj==null || r.match==null){return true;}
			final int subs=(Read.countSubs(r.match));
			final int len=r.length();
			final double q=r.avgQualityByProbabilityDouble(false, r.length());
			final SamLine sl=(SamLine) r.obj;
			mappedReadsProcessedT++;
			mappedBasesProcessedT+=len;
			
			if(subs<=maxBadSubs){
				subSumGoodT+=subs;
				qSumGoodT+=q;
				mapqSumGoodT+=sl.mapq;
				mappedReadsRetainedT++;
				mappedBasesRetainedT+=len;
				return true;
			}
			ArrayList<Var> list=CallVariants.findUniqueSubs(r, sl, varMap, scafMap, maxBadSubAlleleDepth, minBadSubReadDepth);
			if(list==null || list.size()<=maxBadSubs){
				subSumGoodT+=subs;
				qSumGoodT+=q;
				mapqSumGoodT+=sl.mapq;
				mappedReadsRetainedT++;
				mappedBasesRetainedT+=len;
				return true;
			}else{
				subSumBadT+=subs;
				qSumBadT+=q;
				mapqSumBadT+=sl.mapq;
				return false;
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of reads processed by this thread */
		protected long basesProcessedT=0;
		/** Number of mapped reads processed by this thread */
		protected long mappedReadsProcessedT=0;
		/** Number of mapped bases processed by this thread */
		protected long mappedBasesProcessedT=0;
		/** Number of mapped reads retained by this thread */
		protected long mappedReadsRetainedT=0;
		/** Number of mapped bases retained by this thread */
		protected long mappedBasesRetainedT=0;
		/** Number of good reads processed by this thread */
		protected long readsOutT=0;
		/** Number of good bases processed by this thread */
		protected long basesOutT=0;
		protected double qSumGoodT=0;
		protected double qSumBadT=0;
		protected long subSumGoodT=0;
		protected long subSumBadT=0;
		protected long mapqSumGoodT=0;
		protected long mapqSumBadT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared input stream */
		private final SamStreamer ss;
		/** Good output stream */
		private final ConcurrentReadOutputStream ros;
		/** Bad output stream */
		private final ConcurrentReadOutputStream rosb;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	
	/** Optional reference */
	private String ref=null;

	/** Good output file path */
	private String outGood=null;
	/** Bad output file path */
	private String outBad=null;

	/** VCF file path */
	private String varFile=null;
	/** Variant file path */
	private String vcfFile=null;
	private VarMap varMap=null;
	private ScafMap scafMap=null;
	
	/*--------------------------------------------------------------*/
	
	/** Maximum allowed unsupported substitutions in a read */
	private int maxBadSubs=1;
	/** Maximum variant depth for a variant to be considered unsupported */
	private int maxBadSubAlleleDepth=2;
	/** Minimum read depth for a variant to be considered unsupported */
	private int minBadSubReadDepth=2;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of mapped reads processed by this thread */
	protected long mappedReadsProcessed=0;
	/** Number of mapped bases processed by this thread */
	protected long mappedBasesProcessed=0;
	/** Number of mapped reads retained by this thread */
	protected long mappedReadsRetained=0;
	/** Number of mapped bases retained by this thread */
	protected long mappedBasesRetained=0;
	/** Number of good reads processed */
	protected long readsOut=0;
	/** Number of good bases processed */
	protected long basesOut=0;
	
	protected double qSumGood=0;
	protected double qSumBad=0;
	protected long subSumGood=0;
	protected long subSumBad=0;
	protected long mapqSumGood=0;
	protected long mapqSumBad=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	static boolean useStreamer=true;
	static int streamerThreads=3;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffoutGood;
	/** Secondary output file */
	private final FileFormat ffoutBad;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
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
	/** Reads are output in input order */
	private boolean ordered=true;
	
}
