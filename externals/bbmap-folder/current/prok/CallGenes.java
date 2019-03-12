package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * This is the executable class for gene-calling.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class CallGenes {
	
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
		CallGenes x=new CallGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CallGenes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			outGff=parser.out1;
			maxReads=parser.maxReads;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffoutGff=FileFormat.testOutput(outGff, FileFormat.GFF, null, true, overwrite, append, ordered);
		ffoutAmino=FileFormat.testOutput(outAmino, FileFormat.FA, null, true, overwrite, append, ordered);

		if(ffoutGff!=null){
			assert(!ffoutGff.isSequence()) : "\nout is for gff files.  To output sequence, please use outa.";
		}
		if(ffoutAmino!=null){
			assert(!ffoutAmino.gff()) : "\nouta is for sequence data.  To output gff, please use out.";
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(b==null && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else if(a.equals("pgm") || a.equals("gm") || a.equals("model")){
				assert(b!=null);
				if(b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("default")){
					b=Data.findPath("?model.pgm");
					pgmList.add(b);
				}else{
					Tools.addFiles(b, pgmList);
				}
			}else if(b==null && new File(arg).exists() && FileFormat.isPgmFile(arg)){
				Tools.addFiles(arg, pgmList);
			}else if(a.equals("outamino") || a.equals("aminoout") || a.equals("outa") || a.equals("outaa") || a.equals("aaout") || a.equals("amino")){
				outAmino=b;
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				//ReadWrite.verbose=verbose;
				GeneCaller.verbose=verbose;
			}else if(a.equals("plus")){
				GeneModel.PROCESS_PLUS_STRAND=Tools.parseBoolean(b);
			}else if(a.equals("minus")){
				GeneModel.PROCESS_MINUS_STRAND=Tools.parseBoolean(b);
			}
			
			else if(a.equals("merge")){
				merge=Tools.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Tools.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("load16skmers")){
				GeneModel.load16Skmers=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("load23skmers")){
				GeneModel.load23Skmers=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("load5skmers")){
				GeneModel.load5Skmers=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("loadtrnakmers")){
				GeneModel.loadtRNAkmers=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("klongtrna")){
				GeneModel.kLongTRna=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("longkmers")){
				GeneModel.load16Skmers=GeneModel.load23Skmers=
						GeneModel.load5Skmers=GeneModel.loadtRNAkmers=Tools.parseBoolean(b);
			}
			
			else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}
			
			else if(a.equals("translate")){
				mode=TRANSLATE;
			}else if(a.equals("retranslate") || a.equals("detranslate")){
				mode=RETRANSLATE;
			}else if(a.equals("recode")){
				mode=RECODE;
			}
			
			else if(a.equals("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxoverlapss") || a.equals("overlapss") || a.equals("overlapsamestrand") || a.equals("moss")){
				maxOverlapSameStrand=Integer.parseInt(b);
			}else if(a.equals("maxoverlapos") || a.equals("overlapos") || a.equals("overlapoppositestrand") || a.equals("moos")){
				maxOverlapOppositeStrand=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minStartScore")){
				minStartScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minStopScore")){
				minStopScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minInnerScore") || a.equalsIgnoreCase("minKmerScore")){
				minKmerScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minOrfScore") || a.equalsIgnoreCase("minScore")){
				minOrfScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minAvgScore")){
				minAvgScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("breakLimit")){
				GeneCaller.breakLimit=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("clearcutoffs") || a.equalsIgnoreCase("clearfilters")){
				GeneCaller.breakLimit=9999;
				minOrfScore=-9999;
				minAvgScore=-9999;
				minKmerScore=-9999;
				minStopScore=-9999;
				minStartScore=-9999;
			}

			else if(a.equalsIgnoreCase("e1")){
				Orf.e1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e2")){
				Orf.e2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e3")){
				Orf.e3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("f1")){
				Orf.f1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f2")){
				Orf.f2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f3")){
				Orf.f3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("p0")){
				GeneCaller.p0=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p1")){
				GeneCaller.p1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p2")){
				GeneCaller.p2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p3")){
				GeneCaller.p3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p4")){
				GeneCaller.p4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p5")){
				GeneCaller.p5=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p6")){
				GeneCaller.p6=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("q1")){
				GeneCaller.q1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q2")){
				GeneCaller.q2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q3")){
				GeneCaller.q3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q4")){
				GeneCaller.q4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q5")){
				GeneCaller.q5=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("lookback")){
				GeneCaller.lookbackPlus=GeneCaller.lookbackMinus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackplus")){
				GeneCaller.lookbackPlus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackminus")){
				GeneCaller.lookbackMinus=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("compareto")){
				compareToGff=b;
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(pgmList.isEmpty()){
			String b=Data.findPath("?model.pgm");
			pgmList.add(b);
		}
		
		if(Shared.threads()<2){ordered=false;}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required.";
		assert(!pgmList.isEmpty()) : "At least 1 pgm file is required.";
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		pgmList=Tools.fixExtension(pgmList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGff, outAmino)){
			outstream.println((outGff==null)+", "+outGff);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outGff+", "+outAmino+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(pgmList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(outGff);
		foo.add(outAmino);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create read streams and process all data */
	void process(Timer t){
		
		GeneModel pgm=PGMTools.loadAndMerge(pgmList);
		
		if(GeneCaller.call16S || GeneCaller.call23S || GeneCaller.calltRNA || GeneCaller.call5S){
			pgm.loadLongKmers();
		}
		
		ByteStreamWriter bsw=makeBSW(ffoutGff);
		if(bsw!=null){
			bsw.forcePrint("##gff-version 3\n");
		}
		
		ConcurrentReadOutputStream ros=makeCros(ffoutAmino);
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsIn=genesOut=0;
		basesIn=bytesOut=0;
		
		for(String fname : fnaList){
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(fname);

			//Process the reads in separate threads
			spawnThreads(cris, bsw, ros, pgm);
			
			//Close the input stream
			errorState|=ReadWrite.closeStreams(cris, ros);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the output stream
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println(Tools.linesBytesOut(readsIn, basesIn, genesOut, bytesOut, 8, false));
		
		outstream.println();
		outstream.println("Gene Stops Made:      \t "+Tools.padLeft(geneStopsMade, 12));
		outstream.println("Gene Starts Made:     \t "+Tools.padLeft(geneStartsMade, 12));
		outstream.println("Gene Starts Retained: \t "+Tools.padLeft(geneStartsRetained, 12));
		outstream.println("Gene Stops Retained:  \t "+Tools.padLeft(geneStopsRetained, 12));
		outstream.println("Genes Out:            \t "+Tools.padLeft(geneStartsOut, 12));

		outstream.println();
		outstream.println("All ORF Stats:");
		outstream.println(stCds);
		
		outstream.println();
		outstream.println("Retained ORF Stats:");
		outstream.println(stCds2);
		
		outstream.println();
		outstream.println("Called ORF Stats:");
		outstream.println(stCdsPass);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		if(compareToGff!=null){
			CompareGff.main(new String[] {outGff, compareToGff});
		}
	}
	
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ByteStreamWriter bsw, ConcurrentReadOutputStream ros, GeneModel pgm){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, ros, pgm, minLen, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
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
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			genesOut+=pt.genesOutT;
			bytesOut+=pt.bytesOutT;
			
			geneStopsMade+=pt.caller.geneStopsMade;
			geneStartsMade+=pt.caller.geneStartsMade;
			geneStartsRetained+=pt.caller.geneStartsRetained;
			geneStopsRetained+=pt.caller.geneStopsRetained;
			geneStartsOut+=pt.caller.geneStartsOut;

			stCds.add(pt.caller.stCds);
			stCds2.add(pt.caller.stCds2);
			stCdsPass.add(pt.caller.stCdsPass);
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(4, 64, (Shared.threads()*2)/3) : 4);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, ConcurrentReadOutputStream ros_,
				GeneModel pgm_, final int minLen, final int tid_){
			cris=cris_;
			bsw=bsw_;
			ros=ros_;
			pgm=pgm_;
			tid=tid_;
			caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
					minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processList(ListNum<Read> ln){
			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;

//			System.err.println(reads.size());
			
			ArrayList<Orf> orfList=new ArrayList<Orf>();
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				Read r1=reads.get(idx);
				Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsInT+=r1.pairCount();
				basesInT+=initialLength1+initialLength2;
				
				if(r2!=null){
					if(merge){
						final int insert=BBMerge.findOverlapStrict(r1, r2, false);
						if(insert>0){
							r2.reverseComplement();
							r1=r1.joinRead(insert);
							r2=null;
						}
					}else if(ecco){
						BBMerge.findOverlapStrict(r1, r2, true);
					}
				}
				
				{
					//Reads are processed in this block.
					{
						ArrayList<Orf> list=processRead(r1);
						if(list!=null){orfList.addAll(list);}
					}
					if(r2!=null){
						ArrayList<Orf> list=processRead(r2);
						if(list!=null){orfList.addAll(list);}
					}
				}
			}
			
			genesOutT+=orfList.size();
			ByteBuilder bb=new ByteBuilder();
			
			if(bsw!=null){
				if(bsw.ordered){
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
					}
					bsw.add(bb, ln.id);
					bytesOutT+=bb.length();
				}else{
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
//						if(bb.length()>=32000){
//							bytesOutT+=bb.length();
//							bsw.addJob(bb);
//							bb=new ByteBuilder();
//						}
					}
					if(bb.length()>0){
						bsw.addJob(bb);
						bytesOutT+=bb.length();
					}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		ArrayList<Orf> processRead(final Read r){
			ArrayList<Orf> list=caller.callGenes(r, pgm);
			
			if(ros!=null){
				if(mode==TRANSLATE){
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						ros.add(prots, r.numericID);
					}
				}else if(mode==RETRANSLATE) {
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						ArrayList<Read> ret=detranslate(prots);
						ros.add(ret, r.numericID);
					}
				}else if(mode==RECODE) {
					if(list!=null && !list.isEmpty()){
						Read recoded=recode(r, list);
						r.mate=null;
						ArrayList<Read> rec=new ArrayList<Read>(1);
						rec.add(recoded);
						ros.add(rec, r.numericID);
					}
				}else{
					assert(false) : mode;
				}
			}
			
			return list;
		}
		
		/** Number of reads processed by this thread */
		protected long readsInT=0;
		/** Number of bases processed by this thread */
		protected long basesInT=0;
		
		/** Number of genes called by this thread */
		protected long genesOutT=0;
		/** Number of bytes written by this thread */
		protected long bytesOutT=0;
		
		protected ConcurrentReadOutputStream ros;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ByteStreamWriter bsw;
		/** Gene Model for annotation (not really needed) */
		private final GeneModel pgm;
		/** Gene Caller for annotation */
		final GeneCaller caller;
		/** Thread ID */
		final int tid;
	}
	
	public static ArrayList<Read> translate(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> prots=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==Orf.CDS){
					Read aa=translate(orf, r.bases, r.id);
					prots.add(aa);
				}
			}
			r.reverseComplement();
		}
		return prots;
	}
	
	public static Read recode(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return r;}
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==Orf.CDS){
					recode(orf, r.bases);
				}
			}
			r.reverseComplement();
		}
		return r;
	}
	
	public static ArrayList<Read> detranslate(final ArrayList<Read> prots){
		if(prots==null || prots.isEmpty()){return null;}
		ArrayList<Read> nucs=new ArrayList<Read>(prots.size());
		for(int strand=0; strand<2; strand++){
			for(Read prot : prots){
				Read nuc=detranslate(prot);
				nucs.add(nuc);
			}
		}
		return nucs;
	}
	
	public static Read translate(Orf orf, byte[] bases, String id){
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		if(orf.strand==1){orf.flip();}
		Read r=new Read(acids, null, id+"\t"+(Shared.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, Read.AAMASK);
		return r;
	}
	
	public static void recode(Orf orf, byte[] bases){
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		for(int apos=0, bpos=orf.start; apos<acids.length; apos++){
			byte aa=acids[apos];
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			for(int i=0; i<3; i++, bpos++) {
				bases[bpos]=(byte)codon.charAt(i);
			}
		}
		if(orf.strand==1){orf.flip();}
	}
	
	public static Read detranslate(Read prot){
		ByteBuilder bb=new ByteBuilder(prot.length()*3);
		for(byte aa : prot.bases){
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			bb.append(codon);
		}
		Read r=new Read(bb.array, null, prot.id, prot.numericID, 0);
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static GeneCaller makeGeneCaller(GeneModel pgm){
		GeneCaller caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
				minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
		return caller;
	}
	
	private long maxReads=-1;
	private boolean merge;
	private boolean ecco;
	
	private long readsIn=0;
	private long basesIn=0;
	private long genesOut=0;
	private long bytesOut=0;
	
	private static int minLen=60;
	private static int maxOverlapSameStrand=80;
	private static int maxOverlapOppositeStrand=120;
	
	private static float minStartScore=0.00f;
	private static float minStopScore=-0.02f;
	private static float minKmerScore=0.06f;
	private static float minOrfScore=40f; //Higher increases SNR dramatically but reduces TP rate
	private static float minAvgScore=0.10f; //Not very effective
	
	long geneStopsMade=0;
	long geneStartsMade=0;
	long geneStartsRetained=0;
	long geneStopsRetained=0;
	long geneStartsOut=0;

	ScoreTracker stCds=new ScoreTracker(Orf.CDS);
	ScoreTracker stCds2=new ScoreTracker(Orf.CDS);
	ScoreTracker stCdsPass=new ScoreTracker(Orf.CDS);
	
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> pgmList=new ArrayList<String>();
	private String outGff=null;
	private String outAmino=null;
	private String compareToGff=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffoutGff;
	private final FileFormat ffoutAmino;
	
	/** Determines how sequence is processed if it will be output */
	int mode=TRANSLATE;
	
	/** Translate nucleotides to amino acids */
	private static final int TRANSLATE=1;
	/** Translate nucleotides to amino acids,
	 * then translate back to nucleotides */
	private static final int RETRANSLATE=2;
	/** Re-encode coding regions of nucleotide
	 * sequences as a canonical codons */
	private static final int RECODE=3;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean ordered=false; //this is OK sometimes, but sometimes hangs (e.g. on RefSeq mito), possibly if a sequence produces nothing.
	
}
