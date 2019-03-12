package var2;

import java.io.PrintStream;
import java.util.ArrayList;

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

/**
 * @author Brian Bushnell
 * @date January 14, 2017
 *
 */
public class FilterVCF {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FilterVCF x=new FilterVCF(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FilterVCF(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		varFilter.clear();

		boolean setSamFilter=false;
		boolean setVarFilter=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
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
			}else if(a.equals("multithreaded") || a.equals("mt")){
				multithreaded=Tools.parseBoolean(b);
			}else if(a.equals("singlethreaded") || a.equals("st")){
				multithreaded=!Tools.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("sub") || a.equals("subs")){
				Var.CALL_SUB=Tools.parseBoolean(b);
			}else if(a.equals("del") || a.equals("dels")){
				Var.CALL_DEL=Tools.parseBoolean(b);
			}else if(a.equals("ins") || a.equals("inss")){
				Var.CALL_INS=Tools.parseBoolean(b);
			}else if(a.equals("indel") || a.equals("indels")){
				Var.CALL_INS=Var.CALL_DEL=Tools.parseBoolean(b);
			}else if(a.equals("junction") || a.equals("junctions")){
				Var.CALL_JUNCTION=Tools.parseBoolean(b);
			}else if(a.equals("minscore")){
				minScore=Double.parseDouble(b);
			}
			
			else if(a.equals("clearfilters")){
				if(Tools.parseBoolean(b)){
					varFilter.clear();
					samFilter.clear();
				}
			}else if(samFilter.parse(arg, a, b)){
				setSamFilter=true;
			}else if(varFilter.parse(a, b, arg)){
				setVarFilter=true;
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(!setSamFilter){samFilter=null;}
		if(!setVarFilter){varFilter=null;}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least two input files are required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		

		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, multithreaded);
		
		if(ref!=null){ScafMap.loadReference(ref, scafMap, samFilter, true);}
//		inq=new ArrayBlockingQueue<ListNum<byte[]>>(Shared.threads()+1);
		
		//Determine how many threads may be used
		threads=Tools.min(8, Shared.threads());
	}
	
	public void loadHeaderInScafMap(){
		for(byte[] line : header){
			if(Tools.startsWith(line, "##contig=<ID=")){
				scafMap.addFromVcf(line);
			}
		}
	}
	
	public String headerToString(){
		StringBuilder sb=new StringBuilder();
		for(byte[] line : header){
			for(byte b : line){
				sb.append((char)b);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	/** Spawn process threads */
	private ArrayList<ProcessThread> spawnThreads(ByteFile bf, ByteStreamWriter bsw){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(bf, bsw, jobIDOffset));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
//				synchronized(pt){
//					try {
//						pt.wait(200);
//					} catch (InterruptedException e) {}
//				}
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			linesProcessed+=pt.linesProcessedT;
			headerLinesProcessed+=pt.headerLinesProcessedT;
			variantLinesProcessed+=pt.variantLinesProcessedT;
			variantLinesOut+=pt.variantLinesOutT;
			bytesProcessed+=pt.bytesProcessedT;

			allSuccess&=pt.success;
		}
		
		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}
	
	private void processVcfHeader(ByteFile bf, ByteStreamWriter bsw){
		byte[] line=bf.nextLine();

		ByteBuilder bb=new ByteBuilder();
		while(line!=null && (line.length==0 || line[0]=='#')){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				headerLinesProcessed++;
				bytesProcessed+=line.length;
				bb.append(line).append('\n');
				headerLinesOut++;
				header.add(line);
				if(Tools.startsWith(line, VCFFile.CHROM_POS)){
					String[] split=new String(line).split("\t");
					for(int i=9; i<split.length; i++){
						samples.add(split[i]);
					}
				}else{
					String[] split=new String(line).split("=");
					if(split.length==2){
						String a=split[0], b=split[1];
						if(a.equalsIgnoreCase("##ploidy")){
							ploidy=Integer.parseInt(b);
						}else if(a.equalsIgnoreCase("##properPairRate")){
							properPairRate=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##totalQualityAvg")){
							totalQualityAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##mapqAvg")){
							totalMapqAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##readLengthAvg")){
							readLengthAvg=(float) Double.parseDouble(b);
						}
					}
				}
			}
			line=bf.nextLine();
		}
		if(line!=null && line.length>0){bf.pushBack(line);}
		if(bb.length()>0 && bsw!=null){
			bsw.add(bb, jobIDOffset);
			jobIDOffset++;
		}
	}
	
	private void processVcfVarsST(ByteFile bf, ByteStreamWriter bsw){
		boolean varFormatOK=true;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				final boolean isHeader=(line[0]=='#');
				
				if(isHeader){
					assert(false) : "Encountered intermediate header.";
					headerLinesProcessed++;
					if(bsw!=null){bsw.println(line);}
					header.add(line);
					if(Tools.startsWith(line, VCFFile.CHROM_POS)){
						String[] split=new String(line).split("\t");
						for(int i=9; i<split.length; i++){
							samples.add(split[i]);
						}
					}
				}else{
					variantLinesProcessed++;
					VCFLine vline=new VCFLine(line);
					boolean pass=true;
					if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}
					if(pass && varFilter!=null){
						Var v=null;
						
						if(varFormatOK){
							try {
								v=vline.toVar();
							} catch (Throwable e) {
								System.err.println("WARNING: This VCF file does not support Var format.\n"
										+ "Filtering can only be done on location and quality score.\n");
								varFormatOK=false;
							}
						}

						if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
						else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
						else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
						else if(!Var.CALL_JUNCTION && vline.isJunction()){pass=false;}
						
						if(v!=null){
							pass&=varFilter.passesFilter(v, properPairRate,
									totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
						}else{
							pass&=vline.qual>=varFilter.minScore;
						}
					}
					if(pass){
						if(bsw!=null){bsw.println(line);}
						variantLinesOut++;
					}
				}
			}
			line=bf.nextLine();
		}
	}
	
	private void processVcfVarsMT (ByteFile bf, ByteStreamWriter bsw){
		ArrayList<ProcessThread> alpt=spawnThreads(bf, bsw);
//		readBytes(bf);
		waitForFinish(alpt);
	}
	
//	private void readBytes(ByteFile bf){
//		ListNum<byte[]> ln=bf.nextList();
//		while(ln!=null){
//			putBytes(ln);
//			ln=bf.nextList();
//		}
//		putBytes(POISON_BYTES);
//	}
	
	public void filter(FileFormat ff, ByteStreamWriter bsw){
		
		ByteFile bf=ByteFile.makeByteFile(ff);
		
		processVcfHeader(bf, bsw);
		
		loadHeaderInScafMap();
		assert(ScafMap.defaultScafMap().size()>0) : ScafMap.defaultScafMap()+"\n"+headerToString();

		if(multithreaded){
			processVcfVarsMT(bf, bsw);
		}else{
			processVcfVarsST(bf, bsw);
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	void process(Timer t){
		
		ByteStreamWriter bsw;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		filter(ffin1, bsw);
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Header Lines In:   \t"+headerLinesProcessed);
		outstream.println("Variant Lines In:  \t"+variantLinesProcessed);
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/

//	final void putBytes(ListNum<byte[]> list){
//		while(list!=null){
//			try {
//				inq.put(list);
//				list=null;
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//	}
//	
//	final ListNum<byte[]> takeBytes(){
//		ListNum<byte[]> list=null;
//		while(list==null){
//			try {
//				list=inq.take();
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		return list;
//	}
	
	private class ProcessThread extends Thread {
		
		ProcessThread(ByteFile bf_, ByteStreamWriter bsw_, long jobIDOffset_){
			bf=bf_;
			bsw=bsw_;
			offset=jobIDOffset_;
		}

		@Override
		public void run(){
			ListNum<byte[]> ln=bf.nextList();
//			ListNum<byte[]> ln=takeBytes();
			while(ln!=null && ln!=POISON_BYTES){
				ByteBuilder bb=new ByteBuilder(4096);
				for(byte[] line : ln){
					linesProcessedT++;
					processLine(line, bb);
				}
				if(bsw!=null){bsw.add(bb, ln.id+offset);}
				ln=bf.nextList();
//				ln=takeBytes();
			}
//			putBytes(POISON_BYTES);
			success=true;
			synchronized(this){notify();}
		}

		void processLine(byte[] line, ByteBuilder bb){
			linesProcessedT++;
			bytesProcessedT+=line.length;
			if(line.length<1){return;}
			final boolean isHeader=(line[0]=='#');
			if(isHeader){
				assert(false) : "Encountered intermediate header.";
				headerLinesProcessedT++;
				bb.append(line).append('\n');
				synchronized(header){
					header.add(line);
				}
				if(Tools.startsWith(line, VCFFile.CHROM_POS)){
					String[] split=new String(line).split("\t");
					synchronized(samples){
						for(int i=9; i<split.length; i++){
							samples.add(split[i]);
						}
					}
				}
			}else{
				variantLinesProcessedT++;
				boolean pass=true;
				
				VCFLine vline=new VCFLine(line);
				pass&=vline.qual>=minScore;
				
//				if(pass){
//					Var v=null;
//					if(varFormatOK){
//						try {
//							v=Var.fromVCF(line, scafMap, true);
//						} catch (Throwable e) {
//							System.err.println("WARNING: This VCF file does not support Var format.\n"
//									+ "Filtering can only be done on location and quality score.\n");
//							varFormatOK=false;
//						}
//					}
//
//					if(v!=null){
//						if(pass && samFilter!=null){pass&=samFilter.passesFilter(v, scafMap);}
//						if(pass && varFilter!=null){
//							if(!Var.CALL_DEL && v.type()==Var.DEL){pass=false;}
//							else if(!Var.CALL_INS && v.type()==Var.INS){pass=false;}
//							else if(!Var.CALL_SUB && v.type()==Var.SUB){pass=false;}
//							pass&=varFilter.passesFilter(v, properPairRate,
//									totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
//							assert(false);
//						}
//					}else{
//						VCFLine vline=new VCFLine(line);
//						if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}
//
//						if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
//						else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
//						else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
//
//						pass&=vline.qual>=minScore;
//					}
//				}
				
				{	
					if(pass){
						if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
						else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
						else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
						else if(!Var.CALL_JUNCTION && vline.isJunction()){pass=false;}
					}

					if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}

					if(pass && varFilter!=null){
						if(varFormatOK){
							try {
								
								final Var v;
								if(threads>1){
									//Fast but does not capture everything
									//								v=Var.fromVCF(line, scafMap, false, false);
									//Fast multithreaded but slow singlethreaded
									v=Var.fromVCF(line, scafMap, true, true);
								}else{
									//Fast singlethreaded but slow multithreaded
									//Especially with -ea
									v=vline.toVar();
								}
								
								pass&=varFilter.passesFilter(v, properPairRate,
										totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
							} catch (Throwable e) {
								System.err.println("WARNING: This VCF file does not support Var format.\n"
										+ "Filtering can only be done on location and quality score.\n"+e);
								e.printStackTrace();
								varFormatOK=false;
							}
						}
					}
				}
				
				if(pass){
					bb.append(line).append('\n');
					variantLinesOutT++;
				}
			}
		}
		
		final ByteFile bf;
		final ByteStreamWriter bsw;
		final long offset;
		boolean varFormatOK=true;

		long linesProcessedT=0;
		long headerLinesProcessedT=0;
		long variantLinesProcessedT=0;
		long variantLinesOutT=0;
		long bytesProcessedT=0;

		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/

	private long linesProcessed=0;
	private long headerLinesProcessed=0;
	private long variantLinesProcessed=0;
	private long headerLinesOut=0;
	private long variantLinesOut=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	public ArrayList<String> samples=new ArrayList<String>();
	
	SamFilter samFilter=new SamFilter();
	VarFilter varFilter=new VarFilter();
	
	/*--------------------------------------------------------------*/
	
	double minScore=0;
	
	public int ploidy=1;
	public float properPairRate=0;
	public float totalQualityAvg=30;
	public float totalMapqAvg=30;
	public float readLengthAvg=150;
	
	final int threads;
	public boolean multithreaded=false;
	private long jobIDOffset=0;
	
//	private final ArrayBlockingQueue<ListNum<byte[]>> inq;
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String ref=null;

	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	public final ScafMap scafMap=new ScafMap();
	
	/*--------------------------------------------------------------*/

	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, -1);
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
