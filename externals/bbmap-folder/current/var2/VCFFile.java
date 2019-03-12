package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * @author Brian Bushnell
 * @date January 14, 2017
 *
 */
public class VCFFile {
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		VCFLine.AUTOCACHE=true;
		
		Timer t=new Timer();
		String in=null;
		long maxLines=Long.MAX_VALUE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
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
			in=parser.in1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		
		VCFFile vf=new VCFFile(in);
		t.stop();
		vf.printTime(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public VCFFile(String s){
		in1=s;
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, false);
		load();
	}
	
	public VCFFile(FileFormat ff){
		in1=ff.name();
		ffin1=ff;
		load();
	}
	
	void load(){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		
		byte[] line=bf.nextLine();
		
		while(line!=null && (line.length==0 || line[0]=='#')){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				header.add(line);
				if(Tools.startsWith(line, CHROM_POS)){
					String[] split=new String(line).split("\t");
					for(int i=9; i<split.length; i++){
						sampleNames.add(split[i]);
					}
				}
			}
			line=bf.nextLine();
		}
		if(ScafMap.defaultScafMap()==null){ScafMap.setDefaultScafMap(toScafMap(null), in1);}
		assert(ScafMap.defaultScafMap().size()>0) : ScafMap.defaultScafMap()+"\n"+headerToString();
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				final boolean isHeader=(line[0]=='#');
				
				if(isHeader){
					header.add(line);
					if(Tools.startsWith(line, CHROM_POS)){
						String[] split=new String(line).split("\t");
						for(int i=9; i<split.length; i++){
							sampleNames.add(split[i]);
						}
					}
				}else{
					VCFLine vline=new VCFLine(line);
					map.put(vline, vline);
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
	}
	
	void printTime(Timer t){
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Header Lines:      \t"+header.size());
		outstream.println("Variant Lines:     \t"+map.size());
	}
	
	public ArrayList<VCFLine> lines(){
		ArrayList<VCFLine> lines=new ArrayList<VCFLine>(map.size());
		for(Entry<VCFLine, VCFLine> e : map.entrySet()){
			lines.add(e.getValue());
		}
		return lines;
	}
	
	public ScafMap toScafMap(ScafMap sm){
		if(sm==null){sm=new ScafMap();}
		for(byte[] line : header){
			if(Tools.startsWith(line, "##contig=<ID=")){
				sm.addFromVcf(line);
			}
		}
		return sm;
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
	
	/*--------------------------------------------------------------*/
	
	public long linesProcessed() {
		return linesProcessed;
	}
	
	public long bytesProcessed() {
		return bytesProcessed;
	}
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	public ArrayList<String> sampleNames=new ArrayList<String>();
	public LinkedHashMap<VCFLine, VCFLine> map=new LinkedHashMap<VCFLine, VCFLine>();
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	
	private final FileFormat ffin1;
	
	
	/*--------------------------------------------------------------*/
	
	public static final String CHROM_POS="#CHROM\tPOS\t";
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	
}
