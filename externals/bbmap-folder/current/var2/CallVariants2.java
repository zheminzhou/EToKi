package var2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import bloom.KCountArray7MTA;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import stream.SamStreamer;
import structures.ListNum;

/**
 * Calls variants from one or more sam or bam files.
 * 
 * @author Brian Bushnell
 * @date December 18, 2016
 *
 */
public class CallVariants2 {
	
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
		CallVariants2 x=new CallVariants2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CallVariants2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		SamLine.PARSE_0=false;
//		SamLine.PARSE_2=false;
//		SamLine.PARSE_5=false;
//		SamLine.PARSE_6=false;
//		SamLine.PARSE_7=false;
		SamLine.PARSE_8=false;
//		SamLine.PARSE_10=false;
//		SamLine.PARSE_OPTIONAL=false;
		SamLine.PARSE_OPTIONAL_MD_ONLY=true; //I only need the MD tag..
		
		SamLine.RNAME_AS_BYTES=false;
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimLeft=qtrimLeft;
		parser.qtrimRight=qtrimRight;
		parser.trimq=trimq;
		Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=true;

		samFilter.includeUnmapped=false;
		samFilter.includeSupplimentary=false;
		samFilter.includeDuplicate=false;
		samFilter.includeNonPrimary=false;
		samFilter.includeQfail=false;
		samFilter.minMapq=4;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("multi") || a.equals("multisample")){
				boolean multi=Tools.parseBoolean(b);
				assert(multi) : "\nThis program is for multisample variant calling.  Please use CallVariants for single-sample variant calling.\n";
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("ss") || a.equals("samstreamer")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					useStreamer=true;
					streamerThreads=Tools.max(1, Integer.parseInt(b));
				}else{
					useStreamer=Tools.parseBoolean(b);
				}
			}else if(a.equals("parsename")){
				SamLine.PARSE_0=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("noPassDotGenotype") || a.equalsIgnoreCase("noPassDot")){
				Var.noPassDotGenotype=Tools.parseBoolean(b);
			}else if(a.equals("extended")){
				Var.extendedText=Tools.parseBoolean(b);
			}else if(a.equals("useidentity")){
				Var.useIdentity=Tools.parseBoolean(b);
			}else if(a.equals("usehomopolymer") || a.equals("homopolymer")){
				Var.useHomopolymer=Tools.parseBoolean(b);
			}else if(a.equals("usepairing")){
				Var.usePairing=Tools.parseBoolean(b);
			}else if(a.equals("usebias")){
				Var.useBias=Tools.parseBoolean(b);
			}else if(a.equals("nscan") || a.equals("donscan")){
				Var.doNscan=Tools.parseBoolean(b);
			}else if(a.equals("useedist")){
				Var.useEdist=Tools.parseBoolean(b);
			}else if(a.equals("prefilter")){
				prefilter=Tools.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("vcf") || a.equals("vcfout") || a.equals("outvcf") || a.equals("out")){
				vcf=b;
			}else if(a.equals("vcf0") || a.equals("vcfout0") || a.equals("outvcf0")){
				vcf0=b;
			}else if(a.equals("scorehist") || a.equals("qualhist") || a.equals("qhist") || a.equals("shist")){
				scoreHistFile=b;
			}else if(a.equals("border")){
				border=Integer.parseInt(b);
			}else if(a.equals("sample") || a.equals("samplename")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(new File(b).exists()){sampleNames.add(b);}
				else{
					for(String s : b.split(",")){sampleNames.add(s);}
				}
			}
			
			else if(a.equals("ca3") || a.equals("32bit")){
				Scaffold.setCA3(Tools.parseBoolean(b));
			}else if(a.equals("strandedcov") || a.equals("trackstrand")){
				Scaffold.setTrackStrand(Tools.parseBoolean(b));
			}
			
			else if(a.equals("realign")){
				realign=Tools.parseBoolean(b);
			}else if(a.equals("unclip")){
				unclip=Tools.parseBoolean(b);
			}else if(a.equals("realignrows") || a.equals("rerows")){
				Realigner.defaultMaxrows=Integer.parseInt(b);
			}else if(a.equals("realigncols") || a.equals("recols")){
				Realigner.defaultColumns=Integer.parseInt(b);
			}else if(a.equals("realignpadding") || a.equals("repadding") || a.equals("padding")){
				Realigner.defaultPadding=Integer.parseInt(b);
			}else if(a.equals("msa")){
				Realigner.defaultMsaType=b;
			}
			
			else if(samFilter.parse(arg, a, b)){
				//do nothing
			}
			
			else if(a.equals("in") || a.equals("in1") || a.equals("in2")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(new File(b).exists()){in.add(b);}
				else{
					for(String s : b.split(",")){in.add(s);}
				}
			}else if(a.equals("list")){
				for(String line : TextFile.toStringLines(b)){
					in.add(line);
				}
			}else if(a.equals("clearfilters")){
				if(Tools.parseBoolean(b)){
					varFilter.clear();
					samFilter.clear();
				}
			}else if(varFilter.parse(a, b, arg)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(arg.indexOf('=')<0 && (new File(arg).exists() || arg.indexOf(',')>0)){
				if(new File(arg).exists()){in.add(arg);}
				else{
					for(String s : arg.split(",")){in.add(s);}
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		if(ploidy<1){System.err.println("WARNING: ploidy not set; assuming ploidy=1."); ploidy=1;}
		samFilter.setSamtoolsFilter();
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=parser.overwrite;
			
			extin=parser.extin;
			extout=parser.extout;

			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();
			
			trimWhitespace=Shared.TRIM_READ_COMMENTS;
		}
		if(vcf==null){Scaffold.setTrackStrand(false);}
		
		assert(FastaReadInputStream.settingsOK());
		
		ploidyArray=new long[ploidy+1];
		
		//Ensure there is an input file
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, vcf)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+vcf+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Create input FileFormat objects
		for(String s : in){
			FileFormat ff=FileFormat.testInput(s, FileFormat.SAM, extin, true, false);
			ffin.add(ff);
		}
		
		fixSampleNames();
		assert(sampleNames.size()==in.size()) : "Number of sample names and file names must match.";
		
		assert(ref!=null) : "Please specify a reference fasta.";
	}
	
	private void fixSampleNames(){
		if(sampleNames.size()!=0){assert(sampleNames.size()==in.size()) : "Different number of input files ("+in.size()+") and sample names ("+sampleNames.size()+")";}
		if(sampleNames.size()==0){
			HashMap<String, Integer> map=new HashMap<String, Integer>();
			for(String s : in){
				String core=ReadWrite.stripToCore(s);
				if(map.containsKey(core)){
					int x=map.get(core)+1;
					map.put(core, x);
					sampleNames.add(core+"_copy_"+x);
				}else{
					map.put(core, 1);
					sampleNames.add(core);
				}
			}
		}
		
//		assert(false) : sampleNames;
		assert(sampleNames.size()==in.size()) : "Different number of input files ("+in.size()+") and sample names ("+sampleNames.size()+")";
		
		HashSet<String> set=new HashSet<String>();
		for(String s : sampleNames){
			assert(!set.contains(s)) : "Duplicate sample name "+s;
			set.add(s);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void loadReference(){
		if(loadedRef){return;}
		assert(ref!=null);
		ScafMap.loadReference(ref, scafMap, samFilter, true);
		if(realign){Realigner.map=scafMap;}
		loadedRef=true;
	}
	
	/** Create read streams and process all data */
	public void process(Timer t){

		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		Timer t2=new Timer();
		
		if(ref!=null){
			loadReference();
		}
		
		VarMap sharedVarMap=new VarMap(scafMap);
		
		ArrayList<Sample> samples=new ArrayList<Sample>(ffin.size());
		
		for(int i=0; i<ffin.size(); i++){
			FileFormat ff=ffin.get(i);
			String sname=sampleNames.get(i);
			Sample sample=new Sample(ff, sname);
			samples.add(sample);
		}
		
		t2.start("Calculating which variants pass filters.");
		
		long loadedVars=0;
		long varsProcessed0=0;
		for(Sample sample : samples){
			loadedVars+=sample.process1(sharedVarMap);
			varsProcessed0+=sample.varsProcessed;
			sample.clear();
			scafMap.clearCoverage();
		}
		
		t2.stop(loadedVars+" variants passed filters.");
		
		t2.start("Processing second pass.");
		
		long readsProcessed=0;
		long basesProcessed=0;
		long pairedInSequencingReadsProcessed=0;
		long properlyPairedReadsProcessed=0;
		long trimmedBasesProcessed=0;
		long realignmentsAttempted=0;
		long realignmentsSucceeded=0;
		long realignmentsImproved=0;
		long realignmentsRetained=0;
		long varsPrefiltered=0;
		long varsProcessed=0;
		
		for(Sample sample : samples){
			sample.process2(sharedVarMap);
			
			if(sample.vcfName!=null){
				VcfWriter vw=new VcfWriter(sample.varMap, varFilter, sample.readsProcessed-sample.readsDiscarded, 
						sample.pairedInSequencingReadsProcessed, sample.properlyPairedReadsProcessed,
						sample.trimmedBasesProcessed, ref, trimWhitespace, sample.name);
				vw.writeVcfFile(sample.vcfName);
			}
			
			readsProcessed+=sample.readsProcessed;
			basesProcessed+=sample.basesProcessed;
			pairedInSequencingReadsProcessed+=sample.pairedInSequencingReadsProcessed;
			properlyPairedReadsProcessed+=sample.properlyPairedReadsProcessed;
			trimmedBasesProcessed+=sample.trimmedBasesProcessed;
			realignmentsAttempted+=sample.realignmentsAttempted;
			realignmentsSucceeded+=sample.realignmentsSucceeded;
			realignmentsImproved+=sample.realignmentsImproved;
			realignmentsRetained+=sample.realignmentsRetained;
			varsPrefiltered+=sample.varsPrefiltered;
			varsProcessed+=sample.varsProcessed;
			
			sample.clear();
			scafMap.clearCoverage();
		}
		

		t2.start("Finished second pass.");
		
		long[] types=sharedVarMap.countTypes();
		
		if(vcf!=null){
			t2.start("Writing output.");
			MergeSamples merger=new MergeSamples();
			merger.filter=varFilter;
			merger.mergeSamples(samples, scafMap, vcf, scoreHistFile);
			t2.stop("Time: ");
		}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		{
			t.stop();
			
			long size=scafMap.lengthSum();
			long a=varsProcessed0, b=loadedVars, c=varsPrefiltered, d=varsProcessed;
			double amult=100.0/a;
			double bmult=100.0/b;
			outstream.println();
			if(prefilter){
				outstream.println(c+" of "+d+" events were screened by the prefilter ("+String.format(Locale.ROOT, "%.4f%%", c*100.0/d)+").");
			}
			outstream.println(b+" of "+a+" variants passed filters ("+String.format(Locale.ROOT, "%.4f%%", b*amult)+").");
			outstream.println();
			outstream.println("Substitutions: \t"+types[Var.SUB]+String.format(Locale.ROOT, "\t%.1f%%", types[Var.SUB]*bmult));
			outstream.println("Deletions:     \t"+types[Var.DEL]+String.format(Locale.ROOT, "\t%.1f%%", types[Var.DEL]*bmult));
			outstream.println("Insertions:    \t"+types[Var.INS]+String.format(Locale.ROOT, "\t%.1f%%", types[Var.INS]*bmult));
			outstream.println("Variation Rate:\t"+(b==0 ? 0 : 1)+"/"+(size/Tools.max(1,b))+"\n");
			
			if(realign){
				outstream.println("Realignments:  \t"+realignmentsAttempted);
				outstream.println("Successes:     \t"+realignmentsSucceeded);
				outstream.println("Improvements:  \t"+realignmentsImproved);
				outstream.println("Retained:      \t"+realignmentsRetained);
				outstream.println();
			}
			
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorStateOverall){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	class Sample{

		public Sample(FileFormat ff_, String sname_){
			ff=ff_;
			name=sname_;
			vcfName=vcf0==null ? null : vcf0.replaceFirst("%", name);
		}

		/** Find all variants that pass filters for this sample and add them to the shared map */
		public long process1(VarMap sharedVarMap) {
			
			Timer t2=new Timer();
			outstream.println("Processing sample "+name+".");
			
			final KCountArray7MTA kca;
			if(prefilter){
				t2.start("Loading the prefilter.");
				kca=prefilter(varFilter.minAlleleDepth);
				double used=(100.0*kca.cellsUsed())/kca.cells;
				outstream.println("Added "+varsProcessed+" events to prefilter; approximately "+(long)(kca.estimateUniqueKmers(2))+" were unique.");
				outstream.println(String.format(Locale.ROOT, "The prefilter is %.2f%% full.", used));
				varsProcessed=0;
				t2.stop("Time: ");
				outstream.println();
			}else{
				kca=null;
			}
			
			t2.start("Loading variants.");
			
			assert(varMap==null);
			varMap=new VarMap(scafMap);
			processInput(ff, kca, null, true);
			final double properPairRate=properlyPairedReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double pairedInSequencingRate=pairedInSequencingReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double totalQualityAvg=totalQualitySum/(double)Tools.max(1, trimmedBasesProcessed);
			final double totalMapqAvg=totalMapqSum/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double readLengthAvg=trimmedBasesProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
//			System.err.println(properlyPairedReadsProcessed+", "+readsProcessed+", "+readsDiscarded);
			varMap.ploidy=ploidy;
			varMap.properPairRate=properPairRate;
			varMap.pairedInSequencingRate=pairedInSequencingRate;
			varMap.totalQualityAvg=totalQualityAvg;
			varMap.totalMapqAvg=totalMapqAvg;
			varMap.readLengthAvg=readLengthAvg;
			t2.stop("Time: ");
			outstream.println();
			
//			assert(false) : filter.toString(properPairRate, ploidy);
			
			long initialCount=varMap.size();
			t2.start("Processing variants.");
			final long[] types=processVariants();
			t2.stop("Time: ");
			outstream.println();
			
			long initialSharedCount=0, finalSharedCount=0;
			for(int i=0; i<=VarMap.MASK; i++){
				synchronized(sharedVarMap.maps[i]){
					initialSharedCount+=sharedVarMap.maps[i].size();
					sharedVarMap.maps[i].putAll(varMap.maps[i]);
					finalSharedCount+=sharedVarMap.maps[i].size();
				}
			}
			return finalSharedCount-initialSharedCount;
		}

		public long process2(VarMap sharedVarMap) {
			Timer t2=new Timer();
			outstream.println("Processing sample "+name+".");
			
			t2.start("Loading variants.");
			
			assert(varMap==null);
			varMap=new VarMap(scafMap);
			processInput(ff, null, sharedVarMap, true);
			final double properPairRate=properlyPairedReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double pairedInSequencingRate=pairedInSequencingReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double totalQualityAvg=totalQualitySum/(double)Tools.max(1, trimmedBasesProcessed);
			final double totalMapqAvg=totalMapqSum/(double)Tools.max(1, readsProcessed-readsDiscarded);
			final double readLengthAvg=trimmedBasesProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
//			System.err.println(properlyPairedReadsProcessed+", "+readsProcessed+", "+readsDiscarded);
			varMap.ploidy=ploidy;
			varMap.properPairRate=properPairRate;
			varMap.pairedInSequencingRate=pairedInSequencingRate;
			varMap.totalQualityAvg=totalQualityAvg;
			varMap.totalMapqAvg=totalMapqAvg;
			varMap.readLengthAvg=readLengthAvg;
			t2.stop("Time: ");
			outstream.println();
			
//			assert(false) : filter.toString(properPairRate, ploidy);
			
			long initialCount=varMap.size();
			t2.start("Processing variants.");
			final long[] types=addSharedVariants(sharedVarMap);
			varMap.calcCoverage(scafMap);
			t2.stop("Time: ");
			outstream.println();
			
			long initialSharedCount=0, finalSharedCount=0;
			for(int i=0; i<=VarMap.MASK; i++){
				synchronized(sharedVarMap.maps[i]){
					initialSharedCount+=sharedVarMap.maps[i].size();
					sharedVarMap.maps[i].putAll(varMap.maps[i]);
					finalSharedCount+=sharedVarMap.maps[i].size();
				}
			}
			return finalSharedCount-initialSharedCount;
		}

		/*--------------------------------------------------------------*/
		
		private KCountArray7MTA prefilter(int minReads){
			int cbits=2;
			while((1L<<cbits)-1<minReads){
				cbits*=2;
			}
			
			long mem=Shared.memAvailable(4);
			long prebits=mem; //1 bit per byte; 1/8th of the memory

			long precells=prebits/cbits;
			if(precells<100000){ //Not enough memory - no point.
				return null;
			}

			KCountArray7MTA kca=new KCountArray7MTA(precells, cbits, 0, 2, null, minReads);

			if(ref==null){
				ScafMap.loadSamHeader(ff, scafMap);
			}

			final SamReadStreamer ss;
			//Create a read input stream
			final ConcurrentReadInputStream cris;
			if(useStreamer && (maxReads<0 || maxReads==Long.MAX_VALUE)){
				cris=null;
				ss=new SamReadStreamer(ff, streamerThreads, false);
				ss.start();
				if(verbose){outstream.println("Started streamer");}
			}else{
				ss=null;
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}

			final int threads=Shared.threads();

			//Fill a list with ProcessThreads
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){
				alpt.add(new ProcessThread(cris, ss, i, kca, null, true, false));
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
				varsProcessed+=pt.varsProcessedT;

				//Accumulate per-thread statistics
				success&=pt.success;
			}

			//Track whether any threads failed
			if(!success){errorState=true;}
			
			kca.shutdown();
			return kca;
		}

		/** Create read streams and process all data */
		void processInput(FileFormat ff, KCountArray7MTA kca, VarMap sharedVarMap, boolean calcCoverage){
			assert(ff.samOrBam()) : ff.name();
			
			if(ref==null){
				ScafMap.loadSamHeader(ff, scafMap);
			}
			
			final SamReadStreamer ss;
			//Create a read input stream
			final ConcurrentReadInputStream cris;
			if(useStreamer && (maxReads<0 || maxReads==Long.MAX_VALUE)){
				cris=null;
				ss=new SamReadStreamer(ff, streamerThreads, false);
				ss.start();
				if(verbose){outstream.println("Started streamer");}
			}else{
				ss=null;
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			
			//Process the reads in separate threads
			spawnThreads(cris, ss, kca, sharedVarMap, calcCoverage);
			
			if(verbose){outstream.println("Finished; closing streams.");}
			
			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris);
		}
		
		private long[] processVariants(){
			return varMap.processVariantsMT(varFilter, scoreArray, ploidyArray);
		}
		
		private long[] addSharedVariants(VarMap sharedVarMap){
			return varMap.addSharedVariantsST(varFilter, sharedVarMap);
		}
		
		/** Spawn process threads */
		private void spawnThreads(final ConcurrentReadInputStream cris, final SamReadStreamer ss,
				final KCountArray7MTA kca, final VarMap sharedVarMap, final boolean calcCoverage){
			
			//Do anything necessary prior to processing
			if(calcCoverage){
				scafMap.clearCoverage();
			}
			
			//Determine how many threads may be used
			final int threads=Shared.threads();
			
			//Fill a list with ProcessThreads
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){
				alpt.add(new ProcessThread(cris, ss, i, null, sharedVarMap, false, calcCoverage));
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
				trimmedBasesProcessed+=pt.trimmedBasesProcessedT;
				readsDiscarded+=pt.readsDiscardedT;
				pairedInSequencingReadsProcessed+=pt.pairedInSequencingReadsProcessedT;
				properlyPairedReadsProcessed+=pt.properlyPairedReadsProcessedT;
				varsPrefiltered+=pt.prefilteredT;
				varsProcessed+=pt.varsProcessedT;
				totalQualitySum+=pt.totalQualitySumT;
				totalMapqSum+=pt.totalMapqSumT;
				success&=pt.success;
				if(pt.realigner!=null){
					realignmentsAttempted+=pt.realigner.realignmentsAttempted;
					realignmentsImproved+=pt.realigner.realignmentsImproved;
					realignmentsSucceeded+=pt.realigner.realignmentsSucceeded;
					realignmentsRetained+=pt.realigner.realignmentsRetained;
				}
			}
			
			//Track whether any threads failed
			if(!success){errorState=true;}
		}
		
		private int dumpVars(HashMap<Var, Var> mapT){
			int added=varMap.dumpVars(mapT);
			assert(mapT.size()==0);
			return added;
		}
		
		public int fixVars(Read r, SamLine sl){
			assert(false);
			return -1;//fixVars(r, sl, varMap, scafMap);
		}
		
		private void clear(){
			readsProcessed=0;
			basesProcessed=0;
			trimmedBasesProcessed=0;
			readsDiscarded=0;
			pairedInSequencingReadsProcessed=0;
			properlyPairedReadsProcessed=0;
			varsPrefiltered=0;
			varsProcessed=0;
			totalQualitySum=0;
			totalMapqSum=0;

			realignmentsAttempted=0;
			realignmentsImproved=0;
			realignmentsSucceeded=0;
			realignmentsRetained=0;

			varMap=null;
//			varMap.clear();
			
			//errorState=false;
		}
		
		/*--------------------------------------------------------------*/
		
		final FileFormat ff;
		final String name;
		final String vcfName;
		
		/** Number of reads processed */
		protected long readsProcessed=0;
		/** Number of bases processed */
		protected long basesProcessed=0;
		/** Number of trimmed, mapped bases processed */
		protected long trimmedBasesProcessed=0;
		/** Number of reads discarded */
		protected long readsDiscarded=0;
		/** Number of paired reads processed by this thread, whether or not they mapped as pairs */
		protected long pairedInSequencingReadsProcessed=0;
		/** Number of properly paired reads processed */
		protected long properlyPairedReadsProcessed=0;
		/** Number of vars ignored via prefilter */
		protected long varsPrefiltered=0;
		/** Number of vars processed */
		protected long varsProcessed=0;

		/** Sum of trimmed, mapped base qualities */
		protected long totalQualitySum=0;
		/** Sum of mapqs */
		protected long totalMapqSum=0;

		protected long realignmentsAttempted=0;
		protected long realignmentsImproved=0;
		protected long realignmentsSucceeded=0;
		protected long realignmentsRetained=0;

		public VarMap varMap;
		
		boolean errorState=false;
		
		

		
		/*--------------------------------------------------------------*/
		
		/** This class is static to prevent accidental writing to shared variables.
		 * It is safe to remove the static modifier. */
		private class ProcessThread extends Thread {
			
			//Constructor
			ProcessThread(final ConcurrentReadInputStream cris_, final SamReadStreamer ss_, final int tid_,
					final KCountArray7MTA kca_, final VarMap sharedVarMap_, final boolean prefilterOnly_,
					final boolean calcCoverage_){
				cris=cris_;
				ss=ss_;
				tid=tid_;
				kca=kca_;
				prefilterOnly=prefilterOnly_;
				realigner=(realign ? new Realigner() : null);
				sharedVarMap=sharedVarMap_;
				calcCoverage=calcCoverage_;
			}
			
			//Called by start()
			@Override
			public void run(){
				//Do anything necessary prior to processing
				
				//Process the reads
				if(cris==null){
					processInner_ss();
				}else{
					processInner_cris();
				}
				
				//Do anything necessary after processing
				if(!varMapT.isEmpty()){
					dumpVars(varMapT);
				}
				assert(varMapT.isEmpty());
				
				//Indicate successful exit status
				success=true;
			}
			
			/** Iterate through the reads */
			void processInner_cris(){
				
				//Grab the first ListNum of reads
				ListNum<Read> ln=cris.nextList();
				//Grab the actual read list from the ListNum
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				//As long as there is a nonempty read list...
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//					if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

					//Loop through each read in the list
					for(int idx=0; idx<reads.size(); idx++){
						final Read r=reads.get(idx);
						assert(r.mate==null);
						
						if(!r.validated()){r.validate(true);}
						
						//Track the initial length for statistics
						final int initialLength=r.length();

						//Increment counters
						readsProcessedT++;
						basesProcessedT+=initialLength;
						
						boolean b=processRead(r);
						
						if(!b){
							readsDiscardedT++;
						}
					}

					//Notify the input stream that the list was used
					cris.returnList(ln);
//					if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

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
			void processInner_ss(){
				
				//Grab the actual read list from the ListNum
				ListNum<Read> ln=ss.nextList();

				//As long as there is a nonempty read list...
				while(ln!=null && ln.size()>0){
					ArrayList<Read> reads=ln.list;
//					if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

					//Loop through each read in the list
					for(int idx=0; idx<reads.size(); idx++){
						final Read r=reads.get(idx);
						assert(r.mate==null);
						
						if(!r.validated()){r.validate(true);}
						
						//Track the initial length for statistics
						final int initialLength=r.length();

						//Increment counters
						readsProcessedT++;
						basesProcessedT+=initialLength;
						
						boolean b=processRead(r);
						if(!b){
							readsDiscardedT++;
						}
					}
					
					ln=ss.nextList();
				}
			}
			
			/**
			 * Process a read.
			 * @param r Read 1
			 * @return True if the reads should be kept, false if they should be discarded.
			 */
			boolean processRead(final Read r){
				if(r.bases==null || r.length()<=1){return false;}
				SamLine sl=(SamLine) r.obj;
				
//				final SamLine oldSL=new SamLine(sl);
//				final Read oldRead=r.clone();
				
				if(samFilter!=null && !samFilter.passesFilter(sl)){return false;}
				if(sl.properPair()){properlyPairedReadsProcessedT++;}
				if(sl.hasMate()){pairedInSequencingReadsProcessedT++;}
				final Scaffold scaf=scafMap.getScaffold(sl);
				final int scafnum=scaf.number;
				
				r.toLongMatchString(false); //Not necessary if scoring can be done on short match string
				if(realign){
					realigner.realign(r, sl, scaf, unclip);
				}
//				System.err.println(sl);
//				System.err.println(new String(r.match));
				
				int leftTrimAmount=border, rightTrimAmount=border;
				if(qtrimLeft || qtrimRight){
					long packed=TrimRead.testOptimal(r.bases, r.quality, trimE);
					if(qtrimLeft){leftTrimAmount=Tools.max(leftTrimAmount, (int)((packed>>32)&0xFFFFFFFFL));}
					if(qtrimRight){rightTrimAmount=Tools.max(rightTrimAmount, (int)((packed)&0xFFFFFFFFL));}
				}
				
				int trimmed=TrimRead.trimReadWithMatch(r, sl, leftTrimAmount, rightTrimAmount, 0, scaf.length, false);
				if(trimmed<0){return false;}
				int extra=(qtrimLeft || qtrimRight) ? trimmed/2 : Tools.min(border, trimmed/2);
//				System.err.println(sl);
//				System.err.println(new String(r.match));
//				System.err.println();
				
				ArrayList<Var> vars=null;
				//				try {
				vars = Var.toVars(r, sl, callNs, scafnum);
				//				} catch (Throwable e) {
				//					// TODO Auto-generated catch block
				//					System.err.println("Bad line:");
				//					System.err.println(oldRead.toString());
				//					System.err.println(r.toString());
				//					System.err.println(oldSL.toString());
				//					System.err.println(sl.toString());
				//					System.err.println("\n");
				//				}

				if(prefilterOnly){
					if(vars==null){return true;}
					for(Var v : vars){
						long key=v.toKey();
						kca.increment(key);
					}
				}else{
					trimmedBasesProcessedT+=r.length();
					totalQualitySumT+=Tools.sum(r.quality);
					totalMapqSumT+=sl.mapq;
					if(calcCoverage){scaf.add(sl);}
					if(vars==null){return true;}

					for(Var v : vars){
						int depth=Integer.MAX_VALUE;
						if(kca!=null){
							depth=kca.read(v.toKey());
						}
						if(depth>=varFilter.minAlleleDepth){
							if(sharedVarMap==null || sharedVarMap.containsKey(v)){
								v.endDistMax+=extra;
								v.endDistSum+=extra;

								Var old=varMapT.get(v);
								if(old==null){varMapT.put(v, v);}
								else{old.add(v);}
							}
						}else{
							prefilteredT++;
						}
					}
					if(varMapT.size()>vmtSizeLimit){
						dumpVars(varMapT);
					}
				}
				varsProcessedT+=vars.size();
				return true;
			}

			private final KCountArray7MTA kca;
			private final boolean prefilterOnly;
			private final VarMap sharedVarMap;
			
			HashMap<Var, Var> varMapT=new HashMap<Var, Var>();
			
			/** Number of vars blocked by the prefilter */
			protected long prefilteredT=0;
			/** Number of vars processed */
			protected long varsProcessedT=0;
			
			/** Sum of trimmed, mapped base qualities */
			protected long totalQualitySumT=0;
			/** Sum of mapqs */
			protected long totalMapqSumT=0;
			
			/** Number of reads processed by this thread */
			protected long readsProcessedT=0;
			/** Number of bases processed by this thread */
			protected long basesProcessedT=0;
			/** Number of trimmed, mapped bases processed. */
			protected long trimmedBasesProcessedT=0;
			/** Number of reads discarded by this thread */
			protected long readsDiscardedT=0;
			/** Number of paired reads processed by this thread, whether or not they mapped as pairs */
			protected long pairedInSequencingReadsProcessedT=0;
			/** Number of properly paired reads processed by this thread */
			protected long properlyPairedReadsProcessedT=0;
			
			/** True only if this thread has completed successfully */
			boolean success=false;
			
			/** Shared input stream */
			private final ConcurrentReadInputStream cris;
			/** Optional SamReadStreamer for high throughput */
			private final SamReadStreamer ss;
			/** For realigning reads */
			final Realigner realigner;
			
			final boolean calcCoverage;
			
			/** Thread ID */
			final int tid;
		}
	}
	
	public static int fixVars(Read r, VarMap varMap, ScafMap scafMap){
		if(r==null || r.bases==null || r.match==null || r.obj==null){return 0;}
		SamLine sl=(SamLine) r.obj;
		if(!sl.mapped()){return 0;}
		return fixVars(r, sl, varMap, scafMap);
	}
	
	public static void unfixVars(Read r){
		if(r==null || r.bases==null || r.match==null || r.obj==null){return;}
		for(int i=0; i<r.match.length; i++){
			if(r.match[i]=='V'){r.match[i]='S';}
		}
	}
	
	public static int fixVars(Read r, SamLine sl, VarMap varMap, ScafMap scafMap){
		if(r==null || r.bases==null || r.match==null){return 0;}
		assert(r.mapped());
		
		if(r.match!=null && r.shortmatch()){
			r.toLongMatchString(false);
		}
		
		int varsFound=0;
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		
		if(r.strand()==Shared.MINUS){r.reverseComplement();}
		
		int rpos=sl.pos-1-SamLine.countLeadingClip(sl.cigar, true, true);
		final int scafnum=scafMap.getNumber(sl.rnameS());
		assert(scafnum>=0) : "Can't find scaffold "+sl.rnameS();

		for(int qpos=0, mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
			final byte b=bases[qpos];
			
			if(m=='S' && scafnum>=0){
				Var v=new Var(scafnum, rpos, rpos+1, b, Var.SUB);
				if(varMap.containsKey(v)){
					varsFound++;
					match[mpos]='V';
				}
			}
			
			if(m!='D'){qpos++;}
			if(m!='I'){rpos++;}
		}
		if(r.strand()==Shared.MINUS){r.reverseComplement();}
		return varsFound;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in=new ArrayList<String>();

	/** VCF output file path */
	private String vcf=null;
	
	/** Individual vcf files */
	private String vcf0="individual_%.vcf.gz";
	
	private String scoreHistFile=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	/** A fasta file. */
	private String ref=null;
	
	private boolean loadedRef=false;
	
	private boolean qtrimLeft=false;
	private boolean qtrimRight=true;
	private float trimq=10;
	private final float trimE;
	
	/*--------------------------------------------------------------*/
	
	public final ScafMap scafMap=new ScafMap();

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	public int ploidy=-1;
	
	public int border=5;

	public boolean realign=false;
	public boolean unclip=false;

	public boolean prefilter=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final ArrayList<FileFormat> ffin=new ArrayList<FileFormat>();
	
	/** Sample names */
	private final ArrayList<String> sampleNames=new ArrayList<String>();
	
	public final VarFilter varFilter=new VarFilter();
	public final SamFilter samFilter=new SamFilter();
	public final long[] scoreArray=new long[200]; //Not used.
	public final long[] ploidyArray; //Not used.
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int vmtSizeLimit=10000;
	
	static boolean callNs=false;
	static boolean trimWhitespace=true;
	
	static boolean useStreamer=true;
	static int streamerThreads=SamStreamer.DEFAULT_THREADS;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorStateOverall=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	
}
