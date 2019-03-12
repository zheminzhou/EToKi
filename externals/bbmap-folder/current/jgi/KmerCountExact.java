package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import assemble.Shaver;
import assemble.Tadpole;
import bloom.KmerCountAbstract;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import kmer.DumpThread;
import kmer.KmerTableSet;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.FastaReadInputStream;
import ukmer.KmerTableSetU;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class KmerCountExact {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		KmerCountExact x=new KmerCountExact(args);
		t2.stop();
//		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerCountExact(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		/* Initialize local variables with defaults */
		boolean useForest_=false, useTable_=false, useArray_=true;
		Parser parser=new Parser();
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
			}else if(a.equals("mincounttodump") || a.equals("mindump") || a.equals("mincount")){
				minToDump=Integer.parseInt(b);
			}else if(a.equals("maxcounttodump") || a.equals("maxdump") || a.equals("maxcount")){
				maxToDump=Integer.parseInt(b);
			}else if(a.equals("dumpthreads")){
				DumpThread.NUM_THREADS=Integer.parseInt(b);
			}else if(a.equals("hist") || a.equals("khist")){
				outHist=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("forest")){
				useForest_=Tools.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("table")){
				useTable_=Tools.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("array")){
				useArray_=Tools.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Tools.parseBoolean(b);
			}else if(a.equals("shave")){
				shave=Tools.parseBoolean(b);
			}else if(a.equals("rinse")){
				rinse=Tools.parseBoolean(b);
			}else if(a.equals("shavedepth")){
				shaveDepth=Integer.parseInt(b);
			}else if(a.equals("histcolumns")){
				histColumns=Integer.parseInt(b);
			}else if(a.equals("histmax") || a.equals("histlen") || a.equals("khistlen") || a.equals("histsize") || a.equals("khistsize")){
				histMax=Tools.parseIntKMG(b);
			}else if(a.equals("histheader")){
				histHeader=Tools.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				histZeros=!Tools.parseBoolean(b);
			}else if(a.equals("gchist")){
				gcHist=Tools.parseBoolean(b);
			}else if(a.equals("logscale")){
				doLogScale=Tools.parseBoolean(b);
			}else if(a.equals("logwidth")){
				logWidth=Double.parseDouble(b);
			}else if(a.equals("logpasses")){
				logPasses=Integer.parseInt(b);
			}
			
			else if(a.equals("minheight")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("peaks") || a.equals("peaksout")){
				outPeaks=b;
			}else if(a.equals("smooth") || a.equals("smoothe")){
				smoothKhist=smoothPeaks=Tools.parseBoolean(b);
			}else if(a.equals("smoothkhist") || a.equals("smoothhist")){
				smoothKhist=Tools.parseBoolean(b);
			}else if(a.equals("smoothpeaks")){
				smoothPeaks=Tools.parseBoolean(b);
			}else if(a.equals("smoothradius") || a.equals("smootheradius")){
				smoothRadius=Integer.parseInt(b);
			}else if(a.equals("maxradius")){
				CallPeaks.maxRadius=Integer.parseInt(b);
			}else if(a.equals("progressivemult")){
				CallPeaks.progressiveMult=Float.parseFloat(b);
			}
			
			else if(KmerTableSet.isValidArgument(a)){
				//Do nothing
			}else if(a.equals("decimals")){
				decimals=Integer.parseInt(b);
			}
			
			else if(a.equals("sketchmode")){
				KmerCountAbstract.SKETCH_MODE=Tools.parseBoolean(b);
			}else if(a.equals("sketch")){
				sketchPath=b;
			}else if(a.equals("sketchlen") || a.equals("sketchlength") || a.equals("sketchsize")){
				sketchLength=Tools.parseIntKMG(b);
			}else if(a.equals("sketchname")){
				sketchName=b;
			}else if(a.equals("sketchid")){
				sketchID=Integer.parseInt(b);
			}else if(SketchObject.parseSketchFlags(arg, a, b)){
				//Do nothing
			}
			
			else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(sketchPath!=null){
			SketchObject.postParse();
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		/* Adjust I/O settings and filenames */
		
		assert(FastaReadInputStream.settingsOK());

		if(outKmers!=null && !Tools.canWrite(outKmers, overwrite)){throw new RuntimeException("Output file "+outKmers+" already exists, and overwrite="+overwrite);}
		if(sketchPath!=null && !Tools.canWrite(sketchPath, overwrite)){throw new RuntimeException("Output file "+sketchPath+" already exists, and overwrite="+overwrite);}
		
		assert(THREADS>0);
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			Shared.printMemory();
			outstream.println();
		}
		
//		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		k=Tadpole.preparseK(args);
		
		if(k<=31){//TODO: 123 add "false" to the clause to force KmerTableSetU usage.
			tables=new KmerTableSet(args, 12);
		}else{
			tables=new KmerTableSetU(args, 0);
		}
		if(tables.prefilter){tables.minProbMain=false;}
		
		ffSketch=FileFormat.testOutput(sketchPath, FileFormat.TXT, null, true, overwrite, append, false);
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, outKmers, outHist, outPeaks, sketchPath);
		
		/* Count kmers */
		process2();
		
		makeKhistAndPeaks();
		
		if(ffSketch!=null){
			makeSketch();
		}
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(){
		
		/* Start phase timer */
		Timer t=new Timer();
		
		AbstractKmerTableSet.DISPLAY_STATS=false;
		
		/* Fill tables with kmers */
		tables.process(t);
		
		if(DISPLAY_PROGRESS){
			outstream.println("After loading:");
			Shared.printMemory();
			outstream.println();
		}
		
		errorState|=tables.errorState;
		
		t.stop();
		outstream.println("Input:                      \t"+tables.readsIn+" reads \t\t"+tables.basesIn+" bases.");
		
		if(tables.qtrimLeft() || tables.qtrimRight()){
			outstream.println("QTrimmed:               \t"+tables.readsTrimmed+" reads ("+String.format(Locale.ROOT, "%.2f",tables.readsTrimmed*100.0/tables.readsIn)+"%) \t"+
					tables.basesTrimmed+" bases ("+String.format(Locale.ROOT, "%.2f",tables.basesTrimmed*100.0/tables.basesIn)+"%)");
		}
		if(tables.minAvgQuality()>0){
			outstream.println("Low quality discards:   \t"+tables.lowqReads+" reads ("+String.format(Locale.ROOT, "%.2f",tables.lowqReads*100.0/tables.readsIn)+"%) \t"+
					tables.lowqBases+" bases ("+String.format(Locale.ROOT, "%.2f",tables.lowqBases*100.0/tables.basesIn)+"%)");
		}
		
		if(shave || rinse){
			kmersRemoved=shave(shave, rinse, shaveDepth);
		}
		
		outstream.println("\nFor K="+tables.kbig());
		outstream.println("Unique Kmers:               \t"+tables.kmersLoaded);
		if(shave || rinse){
			outstream.println("After Shaving:              \t"+(tables.kmersLoaded-kmersRemoved));
		}

		averageCount=tables.kmersIn*1.0/tables.kmersLoaded;
		double actualDepth=Tools.observedToActualCoverage(averageCount);
		double readDepth=(actualDepth*tables.basesIn)/(tables.kmersIn);
		
		outstream.println("Average Kmer Count:         \t"+String.format("%."+decimals+"f", averageCount));
		outstream.println("Estimated Kmer Depth:       \t"+String.format("%."+decimals+"f", actualDepth));
		outstream.println("Estimated Read Depth:       \t"+String.format("%."+decimals+"f", readDepth));
		outstream.println();
		
		outstream.println("Load Time:                  \t"+t);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	long shave(boolean shave, boolean rinse, int maxShaveDepth){
		long sum=0;

		for(int i=0; i<maxShaveDepth; i++){
			int a=i+1, b=maxShaveDepth, c=i+1;
			//				if(i>3){Shaver2.verbose2=true;}
			outstream.println("\nShave("+a+", "+b+", "+c+")");
			final Shaver shaver=Shaver.makeShaver(tables, THREADS, a, b, c, 1, 3, 100, 100, shave, rinse);
			long removed=shaver.shave(a, b);
			sum+=removed;
		}

		System.err.println();
		return sum;
	}
	
	private double makeKhist(String fname, String peaks, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime, boolean smoothKhist, boolean smoothPeaks){
		if(fname==null && peaks==null){return -1;}
		
		final long[][] arrays=tables.makeKhist(fname, cols, max, printHeader, printZeros, printTime, smoothKhist, gcHist, doLogScale, logWidth, logPasses, smoothRadius);
		final long[] array=arrays[0];
		final long[] gcArray=arrays[1];
		
		double avg=Tools.averageHistogram(array);
		
		if(peaks!=null){
			CallPeaks.printClass=false;
			ArrayList<String> args=new ArrayList<String>();
			if((smoothPeaks && !smoothKhist) && smoothRadius>0){//!smoothKhist because if smoothKhist is true the array will already be smoothed
				args.add("smoothradius="+smoothRadius);
				args.add("smoothprogressive=t");
			}
			CallPeaks.printPeaks(array, gcArray, peaks, overwrite, minHeight, minVolume, minWidth, 
					Tools.max(tables.filterMax()+2, minPeak), maxPeak, maxPeakCount, k, ploidy, doLogScale, logWidth, args);
		}
		return avg;
	}
	
	private void makeSketch(){
		Timer ts=new Timer();
		outstream.println("Generating sketch.");
		SketchObject.maxGenomeFraction=1;
		SketchObject.k=k;
		SketchTool sketcher=new SketchTool(sketchLength, minToDump, false, false);
		Sketch sketch=sketcher.toSketch((KmerTableSet)tables, true);
		if(sketch==null){
			errorState=true;
			System.err.println("WARNING: No sketch was produced, presumably because no kmers passed the filter criteria.");
			assert(false);
			return;
		}
		sketch.setName0(ReadWrite.stripToCore(ffSketch.name()));
		SketchTool.write(sketch, ffSketch);
		ts.stop();
		outstream.println("Sketch Time:                \t"+ts);
	}
	
	private void makeKhistAndPeaks(){
		if(THREADS>1 && (outHist!=null || outPeaks!=null) && outKmers!=null){
			Timer tout=new Timer();
			tout.start();
			Thread a=new DumpKmersThread();
			Thread b=new MakeKhistThread();
			a.start();
			b.start();
			while(a.getState()!=Thread.State.TERMINATED){
				try {
					a.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			while(b.getState()!=Thread.State.TERMINATED){
				try {
					b.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			tout.stop();
			outstream.println("Write Time:                 \t"+tout);
		}else{
			if(outHist!=null || outPeaks!=null){
				averageCount=makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, true, smoothKhist, smoothPeaks);
			}
			if(outKmers!=null){
				//			tables.dumpKmersAsText(outKmers, minToDump, true);
				tables.dumpKmersAsBytes_MT(outKmers, minToDump, maxToDump, true, null);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class DumpKmersThread extends Thread {
		
		DumpKmersThread(){}
		
		@Override
		public void run(){
			tables.dumpKmersAsBytes_MT(outKmers, minToDump, maxToDump, false, null);
		}
		
	}
	
	private class MakeKhistThread extends Thread {
		
		MakeKhistThread(){}
		
		@Override
		public void run(){
			makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, false, smoothKhist, smoothPeaks);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Hold kmers. */
	private final AbstractKmerTableSet tables;
	
	private boolean shave=false;
	private boolean rinse=false;
	private int shaveDepth=1;
	
	private double averageCount=-1;
	private long basesIn=-1;
	private long readsIn=-1;
	private int decimals=3;
	
	private long kmersRemoved=0;
	
	/** Kmer count output file */
	private String outKmers=null;
	/** Histogram output file */
	private String outHist=null;
	/** Histogram peak output file */
	private String outPeaks=null;
	
	private int smoothRadius=1;
	private boolean smoothKhist=false;
	private boolean smoothPeaks=false;
	
	private boolean errorState=false;
	
	/** Histogram columns */
	private int histColumns=2;
	/** Histogram rows */
	private int histMax=100000;
	/** Print a histogram header */
	private boolean histHeader=true;
	/** Histogram show rows with 0 count */
	private boolean histZeros=false;
	/** Add gc information to kmer histogram */
	protected boolean gcHist=false;
	
	boolean doLogScale=true;
	double logWidth=0.1;
	int logPasses=1;
	
	private long minHeight=2;
	private long minVolume=5;
	private int minWidth=3;
	private int minPeak=2;
	private int maxPeak=Integer.MAX_VALUE;
	private int maxPeakCount=12;
	
	private int ploidy=-1;
	
	private String sketchPath=null;
	private int sketchLength=10000;
	private String sketchName;
	private int sketchID;
	private final FileFormat ffSketch;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** min kmer count to dump to text */
	private int minToDump=1;
	private int maxToDump=Integer.MAX_VALUE;

	final int k;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=false;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();

	
}
