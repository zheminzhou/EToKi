package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import server.ServerTools;
import shared.KillSwitch;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class SendSketch extends SketchObject {
	
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
		
		final boolean oldUnpigz=ReadWrite.USE_UNPIGZ;
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		SendSketch x=new SendSketch(args);
		
		//Run the object
		x.process(t);
		
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		Shared.setBufferLen(oldBufLen);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		assert(!x.errorState) : "This program ended in an error state.";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SendSketch(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_UNPIGZ=true;
		KILL_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		defaultParams.printFileName=true;
		defaultParams.inputVersion=Shared.BBMAP_VERSION_STRING;
		boolean setBlacklist=false;
		boolean setLocal=false;
		boolean setPrintDepth=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(a.equals("blacklist")){
				setBlacklist=true;
				parseSketchFlags(arg, a, b);
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//Do nothing
			}else if(a.equals("local")){
				local=Tools.parseBoolean(b);
				setLocal=true;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("address")){
				assert(b!=null) : "Bad parameter: "+arg;
				address=toAddress(b);
			}
			
			else if(a.equalsIgnoreCase("nt")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}else if(a.equalsIgnoreCase("silva") || a.equalsIgnoreCase("ribo")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}else if(a.equalsIgnoreCase("refseq")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}else if(a.equalsIgnoreCase("img")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}else if(a.equalsIgnoreCase("nr")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
				amino=true;
			}else if(a.equalsIgnoreCase("refseqprot") || a.equalsIgnoreCase("prokprot") 
					|| a.equalsIgnoreCase("protein") || a.equalsIgnoreCase("protien") || a.equalsIgnoreCase("prot")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
				translate=true;
			}else if(a.equalsIgnoreCase("mito")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}else if(a.equalsIgnoreCase("fungi")){
				address=toAddress(a);
				if(blacklist==null){blacklist=Blacklist.toBlacklist(a);}
			}
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxID=Integer.parseInt(b);
			}else if(a.equals("spid")){
				outSpid=Integer.parseInt(b);
			}else if(a.equals("imgid")){
				outImgID=Integer.parseInt(b);
			}else if((a.startsWith("meta_") || a.startsWith("mt_")) && b!=null){
				if(outMeta==null){outMeta=new ArrayList<String>();}
				int underscore=a.indexOf('_', 0);
				outMeta.add(a.substring(underscore+1)+":"+b);
			}
			
			else if(a.equals("outsketch") || a.equals("outs") || a.equals("sketchout") || a.equals("sketch")){
				outSketch=b;
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(b==null && in.isEmpty() && new File(arg).exists()){
				in.add(arg);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		outMeta=SketchObject.fixMeta(outMeta);
		
		if(address!=null && (address.equals(refseqAddress) || address.equals(prokProtAddress)) && !SET_AUTOSIZE_FACTOR){AUTOSIZE_FACTOR=2.0f;} //High-resolution mode
		
		while(address!=null && address.endsWith("/")){
			address=address.substring(0,  address.length()-1);
		}
		
		setFromAddress(address, setBlacklist);
		
		if(local){blacklist=null;}
		
		postParse();
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		if(!ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out, outSketch)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out+", "+outSketch+"\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(targetSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts(), defaultParams.mergePairs);
		
//		assert(false) : defaultParams.toString()+"\n"+k+", "+amino+", "+HASH_VERSION;
		if(verbose){
			if(local){outstream.println("Running in local mode.");}
			if(useWhitelist){outstream.println("Using a whitelist.");}
			if(blacklist!=null){outstream.println("Using a blacklist.");}
		}
		
		defaultParams.postParse(false);
		allowMultithreadedFastq=(in.size()==1 && Shared.threads()>2);
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
	}
	
	private void setFromAddress(String address, boolean setBlacklist){
		if(address.equals(nrAddress)){
			amino=true;
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
//			defaultParams.dbName="nr";
			assert(false) : "Need to set K.";
		}else if(address.equals(prokProtAddress)){
			translate=true;
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
//			defaultParams.dbName="nr";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.prokProtBlacklist();}
		}else if(address.equals(ntAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="nt";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.ntBlacklist();}
		}else if(address.equals(refseqAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="RefSeq";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.refseqBlacklist();}
		}else if(address.equals(silvaAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="Silva";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.silvaBlacklist();}
		}else if(address.equals(imgAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
//			defaultParams.dbName="IMG";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.imgBlacklist();}
		}else if(address.equals(mitoAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.mitoBlacklist();}
		}else if(address.equals(fungiAddress)){
			if(!setK){k=defaultK; k2=defaultK2;}
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.fungiBlacklist();}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		if(local){processLocal(t);}
		else{processRemote(t);}
	}
	
	private void processRemote(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		inSketches=tool.loadSketches_MT(defaultParams.mode, defaultParams.samplerate, defaultParams.reads, defaultParams.minEntropy, in);
		final int numLoaded=(inSketches.size());
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" in "+t);
		assert(numLoaded<=MAX_ALLOWED_SKETCHES) : "\nSendSketch is configured to send at most "+MAX_ALLOWED_SKETCHES+" to prevent overwhelming the server.\n"
				+ "If you need to compare more than that, please use CompareSketch locally instead.\n"
				+ "References can be downloaded at http://portal.nersc.gov/dna/microbial/assembly/bushnell/\n";
		t.start();
//		outstream.println(inSketches.get(0));
		
		if(numLoaded>4000){
			SEND_BUFFER_MAX_BYTES*=4;
			SEND_BUFFER_MAX_SKETCHES*=4;
		}else if(numLoaded>1000){
			SEND_BUFFER_MAX_BYTES*=2;
			SEND_BUFFER_MAX_SKETCHES*=2;
		}
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		if(defaultParams.format==DisplayParams.FORMAT_QUERY_REF_ANI || defaultParams.format==DisplayParams.FORMAT_CONSTELLATION){tsw.println(defaultParams.header());}
		
		ByteBuilder bb=new ByteBuilder();
		
		int cntr=0;
		int chunks=0;
		for(Sketch sk : inSketches){
			
			if(sk.taxID<1 || sk.taxID>=minFakeID || outTaxID>0){sk.taxID=outTaxID;}
			if(outSpid>0){sk.spid=outSpid;}
			if(outImgID>0){sk.imgID=outImgID;}
			if(outTaxName!=null){sk.setTaxName(outTaxName);}
			if(outFname!=null){sk.setFname(outFname);}
			if(outName0!=null){sk.setName0(outName0);}
			sk.setMeta(outMeta);
			
			if(bb.length==0){
				bb.append(defaultParams.toString(chunks));
				chunks++;
			}
			sk.toBytes(bb);
			cntr++;
			if(cntr>=SEND_BUFFER_MAX_SKETCHES || bb.length>SEND_BUFFER_MAX_BYTES){ //Don't allow too much data in a single transaction
				if(verbose){outstream.println("Sending:\n"+bb);}
//				outstream.println(cntr+", "+bb.length);
				
				byte[] message=bb.toBytes();
				bb.clear();
				try {
//					outstream.println("Sending to "+address+"\n"+message+"\n"); //123
					StringNum result=ServerTools.sendAndReceive(message, address);
					if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
						System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
						KillSwitch.kill();
					}
					errorState|=checkForError(result.s);
					tsw.print(result.s);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				cntr=0;
			}
		}
		
		if(bb.length>0){
			if(verbose){outstream.println("Sending:\n"+bb);}
			byte[] message=bb.toBytes();
			bb.clear();
			try {
				StringNum result=ServerTools.sendAndReceive(message, address);
				if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
					System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
					KillSwitch.kill();
				}
				errorState|=checkForError(result.s);
				tsw.print(result.s);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		tsw.println();
		tsw.poison();
		
//		outstream.println("sending "+bb.toString());
		
		if(outSketch!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(outSketch, overwrite, append, true, FileFormat.SKETCH);
			bsw.start();
			for(Sketch sk : inSketches){
				sk.toBytes(bb);
				bsw.print(bb);
				bb.clear();
			}
			bsw.poisonAndWait();
			errorState|=bsw.errorState;
		}
		
		tsw.waitForFinish();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	/** For programmatic use */
	public static String sendSketch(Sketch sk, String address, int format, int chunkNum){
		address=toAddress(address);
		ByteBuilder bb=new ByteBuilder();

		DisplayParams params=defaultParams;
		if(format>=0){
			new DisplayParams();
			params.format=format;
		}
		if(bb.length==0){bb.append(params.toString(chunkNum));}
		sk.toBytes(bb);

		byte[] message=bb.toBytes();
		try {
//			System.err.println("Sending to "+address+"\n"+new String(message)+"\n"); //123
			StringNum result=ServerTools.sendAndReceive(message, address);
			if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
				System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
				KillSwitch.kill();
			}
			return result.s;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	private static boolean checkForError(String s){
		if(s==null){return false;}
		return s.contains("java.io.IOException: Server returned HTTP response code:");
	}
	
	private void processLocal(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		
		final String message=defaultParams.toString(0);
		for(String fname : in){
			String address2=address+"/file/"+new File(fname).getAbsolutePath();
			
			if(verbose){outstream.println("Sending:\n"+message+"\nto "+address2);}
			try {
//				outstream.println("Sending to "+address2+"\n"+message+"\n"); //123
				StringNum result=ServerTools.sendAndReceive(message.getBytes(), address2);
				if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
					System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
					KillSwitch.kill();
				}
				tsw.print(result.s);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
		}
		tsw.println();
		
//		outstream.println("sending "+bb.toString());
		
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean addFiles(String a, Collection<String> list){
		int initial=list.size();
		if(a==null){return false;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){
				list.add(s);
			}
		}
		return list.size()>initial;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	private String outSketch=null;
	
	private final SketchTool tool;
	
	private ArrayList<Sketch> inSketches;

	private String address=refseqAddress;
	private boolean local=false;
	
	/*Override metadata */
	private String outTaxName=null;
	private String outFname=null;
	private String outName0=null;
	private int outTaxID=-1;
	private long outSpid=-1;
	private long outImgID=-1;
	private ArrayList<String> outMeta=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
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
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final String toAddress(String b){
		String address=b;
		if(b==null){
			//do nothing
		}else if(b.equalsIgnoreCase("nt")){
			address=ntAddress;
		}else if(b.equalsIgnoreCase("refseq")){
			address=refseqAddress;
		}else if(b.equalsIgnoreCase("silva") || b.equalsIgnoreCase("ribo")){
			address=silvaAddress;
		}else if(b.equalsIgnoreCase("img")){
			address=imgAddress;
		}else if(b.equalsIgnoreCase("refseqprot") || b.equalsIgnoreCase("prokprot") 
				|| b.equalsIgnoreCase("protein") || b.equalsIgnoreCase("protien") || b.equalsIgnoreCase("prot")){
			address=prokProtAddress;
		}else if(b.equalsIgnoreCase("refseqmito") || b.equalsIgnoreCase("mito")){
			address=mitoAddress;
		}else if(b.equalsIgnoreCase("refseqfungi") || b.equalsIgnoreCase("fungi")){
			address=fungiAddress;
		}
		return address;
	}

	public int SEND_BUFFER_MAX_BYTES=8000000;
	public int SEND_BUFFER_MAX_SKETCHES=400;
	private static final int MAX_ALLOWED_SKETCHES=100000;
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
	private static final String ntAddress="https://nt-sketch.jgi-psf.org/sketch"; //gpweb34 port 3071
	private static final String refseqAddress="https://refseq-sketch.jgi-psf.org/sketch"; //gpweb34 port 3072
	private static final String silvaAddress="https://ribo-sketch.jgi-psf.org/sketch"; //gpweb34 port 3073
	private static final String imgAddress="https://img-sketch.jgi-psf.org/sketch";
	private static final String nrAddress="https://nr-sketch.jgi-psf.org/sketch";
	private static final String prokProtAddress="https://protein-sketch.jgi.doe.gov/sketch";//http://gpweb35.nersc.gov:3074/sketch/
	private static final String mitoAddress="https://mito-sketch.jgi-psf.org/sketch";
	private static final String fungiAddress="https://fungi-sketch.jgi-psf.org/sketch";
	
}
