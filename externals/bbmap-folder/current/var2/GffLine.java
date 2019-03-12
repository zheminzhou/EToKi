package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

public class GffLine {
	
	/** Translates VCF to GFF */
	public static void main(String[] args){
		Timer t=new Timer();
		PrintStream outstream=System.err;
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
			t.outstream=outstream;
		}
		
		Parser parser=new Parser();
		String in=null;
		String out=null;
		boolean overwrite=true, append=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("in") || a.equals("vcf")){
				in=b;
			}else if(a.equals("out") || a.equals("gff")){
				out=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(in==null && b==null && i==0 && Tools.canRead(arg)){
				in=arg;
			}else if(in==null && b==null && i==1){
				out=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		translate(in, out, overwrite, append);
		t.stop("Time: \t");
	}
	
	private static void translate(String in, String out, boolean overwrite, boolean append){
		//Create output FileFormat objects
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.GFF, "gff", true, overwrite, append, false);

		//Create input FileFormat objects
		FileFormat ffin=FileFormat.testInput(in, FileFormat.VCF, "vcf", true, true);
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		ByteStreamWriter bsw=null;
		if(ffout!=null){
			bsw=new ByteStreamWriter(ffout);
			bsw.start();
		}
		
		ByteBuilder bb=new ByteBuilder(17000);
		bb.append("##gff-version   3\n");
		String header="#seqid	source	type	start	end	score	strand	phase	attributes";
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>1){
				if(line[0]=='#'){
					if(Tools.startsWith(line, "##fileformat") || Tools.startsWith(line, "##FORMAT") || 
							Tools.startsWith(line, "##INFO") || Tools.startsWith(line, "#CHROM	POS")){
						//skip
					}else{
						int i=1;
						while(i<line.length && line[i]=='#'){i++;}
						i--;
						bb.append(line, i, line.length-i);
						bb.nl();
					}
				}else{
					if(header!=null){
						bb.append(header).append('\n');
						header=null;
					}
					VCFLine vline=new VCFLine(line);
					GffLine gline=new GffLine(vline);
					gline.appendTo(bb);
					bb.nl();
				}
			}
			if(bb.length()>=16384){
				if(bsw!=null){
					bsw.print(bb);
				}
				bb.clear();
			}
		}
		if(bb.length()>0){
			if(bsw!=null){
				bsw.print(bb);
			}
			bb.clear();
		}
		bf.close();
		if(bsw!=null){bsw.poisonAndWait();}
	}
	
	//#seqid	source	type	start	end	score	strand	phase	attributes
	public GffLine(byte[] line){
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		seqid=parseSeqid ? intern(new String(line, a, b-a)) : null;
//		assert(seqid==null || seqid.equals(new String(line, a, b-a)));
//		assert(seqid!=null) : new String(line, a, b-a)+", "+a+", "+b+"\n"+line;
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		if(b==a+1 && line[a]=='.'){source=DOTS;}
		else{source=paseSource ? intern(new String(line, a, b-a)) : null;}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		if(b==a+1 && line[a]=='.'){type=DOTS;}
		else{type=parseType ? intern(new String(line, a, b-a)) : null;}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		start=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		stop=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		if(b<=a){
			//Badly formatted line; common in IMG
			return;
		}
		assert(b>a) : "Missing field 5: "+new String(line);
		if(b==a+1 && line[a]=='.'){score=-1;}
		else{score=Tools.parseInt(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		assert(b==a+1);
		strand=find(line[a], STRANDS);
//		assert(strand>0) : line[a]+", "+Arrays.toString(STRANDS)+", "+(char)line[b];
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		assert(b==a+1);
		if(line[a]=='.'){phase=-1;}
		else{phase=Tools.parseInt(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		if(b==a+1 && line[a]=='.'){attributes=DOTS;}
		else{attributes=parseAttributes ? new String(line, a, b-a) : null;}
		b++;
		a=b;
		
//		assert(strand>=0) : "\n"+this.toString()+"\n"+new String(line);
	}
	
	public GffLine(VCFLine vcf){
		seqid=vcf.scaf;
		source=DOTS;
		type="sequence_variant_obs";
		start=vcf.start()+1;
		stop=vcf.stop()+1;
		score=(float)vcf.qual;
		strand=PLUS;
		phase=-1;
		final int vtype=vcf.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]).append(' ');
		if(vtype==Var.SUB){
			bb.append(vcf.ref).append('>').append(vcf.alt);
		}else if(vtype==Var.DEL){
			bb.append("length ").append(vcf.reflen()-vcf.readlen());
		}else if(vtype==Var.INS){
			int offset=vcf.reflen();
			int length=vcf.readlen()-offset;
			bb.append(vcf.alt, offset, length);
		}else if(vtype==Var.NOCALL){
			bb.append("length ").append(vcf.reflen());
		}
		attributes=bb.toString();
		bb.clear();
	}
	
	public GffLine(Var v, double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
		seqid=v.scafName();
		source=DOTS;
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=(float)v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map);
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	public GffLine(Var v){
		seqid=v.scafName();
		source="BBTools";
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=-1;
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	public static ArrayList<GffLine> loadGffFile(String fname, String types){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.GFF, null, false, false);
		return loadGffFile(ff, types);
	}
	
	public static ArrayList<GffLine> loadGffFile(FileFormat ff, String types){

		HashSet<String> set=null;
		if(types!=null){
			String[] split=types.split(",");
			set=new HashSet<String>(split.length*2);
			for(String s : split){
				set.add(s);
			}
		}
		
		ArrayList<GffLine> list=new ArrayList<GffLine>();
		ByteFile bf=ByteFile.makeByteFile(ff);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]=='#'){
				//skip
			}else{
				GffLine gline=new GffLine(line);
				assert(gline.strand>=0) : "\n"+gline.toString()+"\n"+new String(line)+"\n";
				if(set==null || (gline.type!=null && set.contains(gline.type))){
					list.add(gline);
				}
			}
		}
		
		boolean error=bf.close();
		assert(!error) : "Problem with file "+ff.name();
		return list;
	}
	
	public static void toText(ByteBuilder bb, Var v, double properPairRate, double totalQualityAvg, 
			double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
//		assert(false);
		bb.append(v.scafName(map)).append('\t');
		bb.append('.').append('\t');
		bb.append("sequence_variant_obs").append('\t');
		bb.append(v.start+1).append('\t');
		bb.append(Tools.max(v.start+1, v.stop)).append('\t');
		bb.append(v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map), 2).append('\t');
		bb.append('+').append('\t');
		bb.append('.').append('\t');
//		System.err.println(v.typeString()+", "+v.start+", "+v.stop);
		final int vtype=v.type();
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
	}
	
	public static String toHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();
		
		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("##gff-version   3\n");
		sb.append("#BBMapVersion\t"+Shared.BBMAP_VERSION_STRING+"\n");
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append(String.format(Locale.ROOT, "#rarity\t%.5f\n", rarity));
		sb.append(String.format(Locale.ROOT, "#minAlleleFraction\t%.4f\n", minAlleleFraction));
		sb.append("#reads\t"+reads+"\n");
		sb.append("#pairedReads\t"+pairs+"\n");
		sb.append("#properlyPairedReads\t"+properPairs+"\n");
		sb.append(String.format(Locale.ROOT, "#readLengthAvg\t%.2f\n", readLengthAvg));
		sb.append(String.format(Locale.ROOT, "#properPairRate\t%.4f\n", properPairRate));
		sb.append(String.format(Locale.ROOT, "#totalQualityAvg\t%.4f\n", totalQualityAvg));
		sb.append(String.format(Locale.ROOT, "#mapqAvg\t%.2f\n", mapqAvg));
		if(ref!=null){sb.append("#reference\t"+ref+"\n");}
		
		sb.append("#seqid	source	type	start	end	score	strand	phase	attributes");
		return sb.toString();
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(seqid==null ? "." : seqid).append('\t');
		bb.append(source==null ? "." : source).append('\t');
		bb.append(type==null ? "." : type).append('\t');
		bb.append(start).append('\t');
		bb.append(stop).append('\t');
		if(score<0){bb.append('.').append('\t');}
		else{bb.append(score, 2).append('\t');}
		
		bb.append((strand>=0 ? STRANDS[strand] : (byte)'.')).append('\t');
		
		if(phase<0){bb.append('.').append('\t');}
		else{bb.append(phase).append('\t');}
		
		bb.append(attributes==null ? "." : attributes);
		return bb;
	}
	
	private static int find(byte a, byte[] array){
		for(int i=0; i<array.length; i++){
			if(array[i]==a){return i;}
		}
		return -1;
	}
	
	private static String intern(String s){
		return Data.forceIntern(s);
	}
	
	public String seqid;
	public String source;
	public String type;
	public int start;
	public int stop;
	public float score;
	public int strand;
	public int phase;
	public String attributes;

	private static final byte[] STRANDS=new byte[] {'+', '-', '?', '.'};
	public static final int PLUS=0, MINUS=1, QMARK=2, DOT=3;
	public static final String DOTS=".";

//	public static boolean parseSeqid=true;
//	public static boolean paseSource=false;
//	public static boolean parseType=true;
//	public static boolean parseScore=false;
//	public static boolean parseAttributes=false;

	public static boolean parseSeqid=true;
	public static boolean paseSource=true;
	public static boolean parseType=true;
	public static boolean parseScore=true;
	public static boolean parseAttributes=true;
	
}
