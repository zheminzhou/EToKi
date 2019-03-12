package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringPair;
import var2.CallVariants2.Sample;

/**
 * @author Brian Bushnell
 * @date December 18, 2016
 *
 */
public class MergeSamples {
	
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSamples x=new MergeSamples(args);
		//x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public MergeSamples(){
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	public MergeSamples(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	/*--------------------------------------------------------------*/
	
	public void mergeSamples(ArrayList<Sample> list, ScafMap scafMap, String outVcf, String scoreHistFile){
		map=scafMap;
		ArrayList<StringPair> vcfList=new ArrayList<StringPair>(list.size());
		for(Sample s : list){vcfList.add(new StringPair(s.name, s.vcfName));}
		mergeFiles(vcfList, outVcf, scoreHistFile);
	}
	
	public void mergeFiles(ArrayList<StringPair> list, String outVcf, String scoreHistFile){
		System.err.println("Merging "+list);
		final int ways=list.size();
		ByteFile[] bfa=new ByteFile[ways];
		final boolean allowSubprocess=(ways<=4);
		for(int i=0; i<ways; i++){
			StringPair pair=list.get(i);
			FileFormat ff=FileFormat.testInput(pair.b, FileFormat.VCF, null, allowSubprocess, false);
			bfa[i]=ByteFile.makeByteFile(ff);
		}
//		System.err.println("Made byte files.");
//		assert(false) : outVcf;
//		System.err.println("Started writer.");
		
		mergeMT(outVcf, bfa);

		if(scoreHistFile!=null){
			CallVariants.writeScoreHist(scoreHistFile, scoreArray);
		}
		
//		System.err.println("Closed stream.");
	}
	
	private void mergeST(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			bswVcf=new ByteStreamWriter(outVcf, true, false, true, FileFormat.VCF);
			bswVcf.start();
		}
		
		ByteBuilder bb=new ByteBuilder(34000);
		VCFLine[] row=processRow(bfa, bb);
		while(row!=null){
//			System.err.println("Processed a line.");
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
				if(bb.length>32000){
					if(bswVcf!=null){bswVcf.print(bb);}
					bb=new ByteBuilder(34000);
				}
			}
			row=processRow(bfa, bb);
		}
//		System.err.println("Finished processing.");
		
		if(bswVcf!=null){
			if(bb.length>0){bswVcf.print(bb);}
			bswVcf.poisonAndWait();
		}
	}
	
	private void mergeMT(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			FileFormat ff=FileFormat.testOutput(outVcf, FileFormat.VCF, null, true, true, append, true);
			bswVcf=new ByteStreamWriter(ff);
			bswVcf.start();
		}
		
		ArrayList<MergeThread> alpt=spawnThreads(bswVcf);
		
		long nextID=0;
		ByteBuilder header=new ByteBuilder(34000);
		
		VCFLine[] row=processRow(bfa, header);
		while(row!=null && row[0]==null){//Header
			row=processRow(bfa, header);
		}
		if(bswVcf!=null){
			bswVcf.add(header, nextID);
			nextID++;
		}
		
		ListNum<VCFLine[]> list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
		while(row!=null){
			if(row[0]!=null){
				list.add(row);
				if(list.size()>=200){
					putList(list);
					nextID++;
					list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
				}
			}
			row=processRow(bfa, header);
		}
		if(list.size()>0){
			putList(list);
			nextID++;
		}
		
		putList(POISON_LIST);
		
		waitForFinish(alpt);
		
		if(bswVcf!=null){bswVcf.poisonAndWait();}
	}
	
	VCFLine[] processRow(ByteFile[] bfa, ByteBuilder bb){
		byte[][] lines=new byte[bfa.length][];
		for(int i=0; i<bfa.length; i++){
			byte[] line=bfa[i].nextLine();
			if(line==null){return null;}
			lines[i]=line;
		}
		
		VCFLine[] row=new VCFLine[bfa.length];
		if(lines[0][0]=='#'){
			processHeader(lines, bb);
			return row;
		}
		for(int i=0; i<lines.length; i++){
			byte[] line=lines[i];
			row[i]=new VCFLine(line);
			if(i>0){assert(row[i].pos==row[0].pos) : "\n"+row[0]+"\n"+row[i];}
		}
		return row;
	}
	
	void processHeader(byte[][] lines, ByteBuilder bb){
		String[][] matrix=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			matrix[i]=new String(lines[i]).split("=");
		}
		
		if(matrix[0][0].equals("##ploidy")){
			ploidy=Integer.parseInt(matrix[0][1]);
			bb.append("##ploidy="+ploidy+"\n");
		}else if(matrix[0][0].equals("##reads")){
			for(String[] split : matrix){
				reads+=Long.parseLong(split[1]);
			}
			bb.append("##reads="+reads+"\n");
		}else if(matrix[0][0].equals("##pairedReads")){
			for(String[] split : matrix){
				pairedReads+=Long.parseLong(split[1]);
			}
			bb.append("##pairedReads="+pairedReads+"\n");
		}else if(matrix[0][0].equals("##properlyPairedReads")){
			for(String[] split : matrix){
				properlyPairedReads+=Long.parseLong(split[1]);
			}
			properPairRate=properlyPairedReads*1.0/(Tools.max(1, reads));
			bb.append("##properlyPairedReads="+properlyPairedReads+"\n");
			bb.append("##properPairRate="+String.format(Locale.ROOT, "%.4f\n", properPairRate));
		}else if(matrix[0][0].equals("##properPairRate")){
			//do nothing
		}else if(matrix[0][0].equals("##totalQualityAvg")){
			totalQualityAvg=0;
			for(String[] split : matrix){
				totalQualityAvg+=Float.parseFloat(split[1]);
			}
			totalQualityAvg/=lines.length;
			bb.append("##totalQualityAvg="+String.format(Locale.ROOT, "%.4f\n", totalQualityAvg));
		}else if(matrix[0][0].equals("##mapqAvg")){
			mapqAvg=0;
			for(String[] split : matrix){
				mapqAvg+=Float.parseFloat(split[1]);
			}
			mapqAvg/=lines.length;
			bb.append("##mapqAvg="+String.format(Locale.ROOT, "%.2f\n", mapqAvg));
		}else if(matrix[0][0].equals("##readLengthAvg")){
			readLengthAvg=0;
			for(String[] split : matrix){
				readLengthAvg+=Float.parseFloat(split[1]);
			}
			readLengthAvg/=lines.length;
			bb.append("##readLengthAvg="+String.format(Locale.ROOT, "%.2f\n", readLengthAvg));
		}else if(matrix[0][0].startsWith("#CHROM\tPOS\t")){
			bb.append(lines[0]);
			for(int i=1; i<lines.length; i++){
				String[] split=new String(lines[i]).split("\t");
				bb.tab().append(split[split.length-1]);
			}
			bb.nl();
		}else{
			bb.append(lines[0]);
			bb.nl();
		}
	}
	
	VCFLine merge(VCFLine[] row){

//		System.err.println(row.length);
//		System.err.println(row[0]);
		
		Var sum=null;
		VCFLine best=null;
		for(VCFLine line : row){
			if(best==null || line.qual>best.qual){best=line;}
			Var v=line.toVar();
			assert(v!=null);
			if(sum==null){sum=v;}
			else{
				sum.add(v);
				sum.addCoverage(v);
			}
		}
		assert(best!=null);
		assert(sum!=null) : row.length+", "+row[0];
		
		//ByteBuilder bb, double properPairRate, double totalQualityAvg, double mapqAvg, int ploidy, ScafMap map, VarFilter filter, boolean trimWhitespace
		ByteBuilder bb=sum.toVCF(new ByteBuilder(), properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, ploidy, map, filter, trimWhitespace);
		VCFLine merged=new VCFLine(bb.toBytes());
		merged.samples.clear();
		for(VCFLine line : row){
			merged.samples.addAll(line.samples);
		}
		if(merged.qual<best.qual){
			merged.qual=best.qual;
			merged.filter=best.filter;
		}
		scoreArray[(int)merged.qual]++;
		return merged;
	}
	
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	
	final ListNum<VCFLine[]> takeList(){
		ListNum<VCFLine[]> list=null;
		while(list==null){
			try {
				list=inq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return list;
	}
	
	final void putList(ListNum<VCFLine[]> list){
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Spawn process threads */
	private ArrayList<MergeThread> spawnThreads(ByteStreamWriter bsw){
		
		//Do anything necessary prior to processing
		
		//Fill a list with MergeThreads
		ArrayList<MergeThread> alpt=new ArrayList<MergeThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new MergeThread(bsw));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(MergeThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	private void waitForFinish(ArrayList<MergeThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(MergeThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
	}
	
	private class MergeThread extends Thread {

		MergeThread(ByteStreamWriter bsw_){
			bsw=bsw_;
		}

		@Override
		public void run(){
			ListNum<VCFLine[]> list=takeList();
			while(list!=null && list!=POISON_LIST){
				processList(list);
				list=takeList();
			}
			putList(POISON_LIST);
		}

		private void processList(ListNum<VCFLine[]> list){
			ByteBuilder bb=new ByteBuilder(4096);
			for(VCFLine[] row : list){
				mergeRow(row, bb);
			}
			if(bsw!=null){bsw.add(bb, list.id);}
		}
		
		private void mergeRow(VCFLine[] row, ByteBuilder bb){
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
			}
		}
		
		private final ByteStreamWriter bsw;
		
	}
	
	/*--------------------------------------------------------------*/
	
	final ListNum<VCFLine[]> POISON_LIST=new ListNum<VCFLine[]>(null, -1);
	private final ArrayBlockingQueue<ListNum<VCFLine[]>> inq;
	private final int threads;
	
	/*--------------------------------------------------------------*/

	long readsSum;
	long pairsSum;
	int ploidy=1;
	
	double properPairRate;
	double totalQualityAvg;
	double mapqAvg;
	double readLengthAvg;
	
	long reads;
	long pairedReads;
	long properlyPairedReads;
	
	VarFilter filter;
	ScafMap map;
	boolean trimWhitespace=true;
	
	private String in1=null;
	private String out1=null;
	private String outInvalid=null;
	
	long[] scoreArray=new long[200];
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
