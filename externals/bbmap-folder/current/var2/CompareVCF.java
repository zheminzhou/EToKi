package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map.Entry;

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

/**
 * @author Brian Bushnell
 * @date January 14, 2017
 *
 */
public class CompareVCF {
	
	public static void main(String[] args){
		Timer t=new Timer();
		CompareVCF x=new CompareVCF(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public CompareVCF(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		int mode_=DIFFERENCE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("difference") || a.equals("minus") || a.equals("dif") || a.equals("diff") || a.equals("subtraction") || a.equals("subtract")){
				mode_=DIFFERENCE;
			}else if(a.equals("union") || a.equals("plus")){
				mode_=UNION;
			}else if(a.equals("intersection") || a.equals("shared")){
				mode_=INTERSECTION;
			}else if(a.equals("addsamples")){
				addSamples=Tools.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		{//Process parser fields
			in1=(parser.in1==null ? null : parser.in1.split(","));
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null || in1.length<2){throw new RuntimeException("Error - at least two input files are required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		ffin1=new FileFormat[in1.length];
		for(int i=0; i<in1.length; i++){
			ffin1[i]=FileFormat.testInput(in1[i], FileFormat.TXT, null, true, true);
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
	}
	
	public HashSet<VCFLine> getSet(FileFormat ff, HashSet<VCFLine> set){
		if(set==null){set=new HashSet<VCFLine>();}
		VCFFile vfile=new VCFFile(ff);
		samples.addAll(vfile.sampleNames);
		if(header==null){
			header=vfile.header;
			if(ScafMap.defaultScafMap()==null){
				ScafMap.setDefaultScafMap(vfile.toScafMap(null), ff.name());
			}
		}
		for(Entry<VCFLine, VCFLine> e : vfile.map.entrySet()){
			VCFLine v=e.getValue();
			if(!set.contains(v)){set.add(v);}
		}
		
		linesProcessed+=vfile.linesProcessed();
		headerLinesProcessed+=vfile.header.size();
		variantLinesProcessed+=vfile.map.size();
		bytesProcessed+=vfile.bytesProcessed();
		
		errorState|=vfile.errorState;
		return set;
	}
	
	public HashSet<VCFLine> union(){
		final HashSet<VCFLine> set=new HashSet<VCFLine>();
		for(FileFormat ff : ffin1){
			getSet(ff, set);
		}
		return set;
	}
	
	public HashSet<VCFLine> intersection(){
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.retainAll(set);}
		}
		return set0;
	}
	
	public HashSet<VCFLine> difference(){
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.removeAll(set);}
		}
		return set0;
	}
	
	ArrayList<VCFLine> toList(){
		final HashSet<VCFLine> set;
		if(mode==DIFFERENCE){
			set=difference();
		}else if(mode==UNION){
			set=union();
		}else if(mode==INTERSECTION){
			set=intersection();
		}else{
			throw new RuntimeException("Unknown mode "+mode);
		}
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(set.size());
		list.addAll(set);
		Shared.sort(list);
		return list;
	}
	
	void process(Timer t){
		
		ArrayList<VCFLine> list=toList();
		
		if(ffout1!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
			bsw.start();
			for(byte[] line : header){
				headerLinesOut++;
				bsw.println(line);
			}
			ByteBuilder bb=new ByteBuilder(33000);
			for(VCFLine line : list){
				variantLinesOut++;
				line.toText(bb);
				bb.nl();
				if(bb.length>=32000){
					bsw.print(bb);
					bb.clear();
				}
			}
			if(bb.length>0){
				bsw.print(bb);
				bb.clear();
			}
			errorState|=bsw.poisonAndWait();
		}
		
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
	
	/*--------------------------------------------------------------*/

	private long linesProcessed=0;
	private long headerLinesProcessed=0;
	private long variantLinesProcessed=0;
	private long headerLinesOut=0;
	private long variantLinesOut=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

	public ArrayList<byte[]> header=null;
	public ArrayList<String> samples=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	
	private String in1[]=null;
	private String out1=null;

	private final FileFormat ffin1[];
	private final FileFormat ffout1;
	
	public final int mode;
	
	public boolean addSamples=true;
	
	/*--------------------------------------------------------------*/
	
	public static int DIFFERENCE=0, UNION=1, INTERSECTION=2;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
