package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Aug 23, 2013
 *
 */
public class RenameReads {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameReads x=new RenameReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RenameReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
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
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("prefix") || a.equals("p")){
				prefix=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Tools.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Tools.parseBoolean(b);
			}else if(a.equals("renamebyinsert")){
				renameByInsert=Tools.parseBoolean(b);
			}else if(a.equals("renamebytrim")){
				renameByTrim=Tools.parseBoolean(b);
			}else if(a.equals("addprefix")){
				addPrefix=Tools.parseBoolean(b);
			}else if(a.equals("prefixonly")){
				prefixOnly=Tools.parseBoolean(b);
			}else if(a.equals("underscore")){
				addUnderscore=Tools.parseBoolean(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				in1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
			
			renameByMapping=FASTQ.TAG_CUSTOM;
		}
		
		if(addUnderscore && prefix!=null && !prefix.endsWith("_")){prefix+="_";}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		if(prefix==null || prefix.length()<1){prefix="";}
		else if(!prefix.endsWith("_") && !prefixOnly){
			prefix=prefix+"_";
		}

		if(renameByInsert){
			prefix="insert=";
			FASTQ.PARSE_CUSTOM=true;
		}else if(renameByTrim){
			prefix="";
			FASTQ.PARSE_CUSTOM=true;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		if(!parser.setInterleaved){
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

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(renameByMapping){
			assert(ffout1==null || ffout1.fastq()) : "Currently renameByMapping requires fastq output.";
		}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //4567
		}
		
//		TextStreamWriter tsw=new TextStreamWriter(args[2], false, false, true);
//		tsw.start();
		
		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		long x=0;
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				final Read r2=r1.mate;
				
				if(renameByMapping){
					//Should be handled automatically, if output is fastq.
				}else if(r2!=null && (renameByInsert || renameByTrim)){
					
					r1.setMapped(true);
					r2.setMapped(true);
					x=Read.insertSizeMapped(r1, r2, false);
					if(verbose){System.err.println("True Insert: "+x);}
					if(renameByTrim){
						r1.id=r1.numericID+"_"+r1.length()+"_"+Tools.min(x, r1.length())+" 1:";
						r2.id=r2.numericID+"_"+r2.length()+"_"+Tools.min(x, r2.length())+" 2:";
					}else{
						String s=prefix+x;
						r1.id=s+" 1:"+r1.numericID;
						if(r2!=null){
							r2.id=s+" 2:"+r1.numericID;
						}
					}
					
				}else if(prefixOnly){
					r1.id=prefix;
					if(r2!=null){
						r2.id=prefix;
					}
					x++;
				}else if(addPrefix){
					r1.id=prefix+r1.id;
					if(r2!=null){
						r2.id=prefix+r2.id;
					}
					x++;
				}else{
					r1.id=prefix+x;
					if(r2!=null){
						r1.id=r1.id+" 1:";
						r2.id=prefix+x+" 2:";
					}
					x++;
				}
			}
			if(ros!=null){ros.add(reads, ln.id);}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		System.err.println("Time: "+t);
	}
	
	private PrintStream outstream=System.err;
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;

	private String prefix=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;

	private boolean overwrite=true;
	private boolean append=false;
	private boolean verbose=false;
	private long maxReads=-1;
	public boolean errorState=false;

	public boolean addUnderscore=true;
	public boolean renameByMapping=false;
	public boolean renameByInsert=false;
	public boolean renameByTrim=false;
	public boolean addPrefix=false;
	public boolean prefixOnly=false;
	
}
