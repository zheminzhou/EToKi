package var2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import shared.Tools;
import structures.ByteBuilder;

public class VCFLine implements Comparable<VCFLine> {
	
	public VCFLine(byte[] line) {
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		scaf=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		pos=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		id=line[a]=='.' ? DOT : Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		if(b<=a+1){ref=Var.AL_MAP[line[a]];}
		else{ref=Arrays.copyOfRange(line, a, b);}
		reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		if(b<=a+1){alt=Var.AL_MAP[line[a]];}
		else{alt=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;
		
		//Trim matching suffixes for a canonical representation
		if(TRIM_TO_CANONICAL && ref.length>1 && alt.length>1 && ref.length!=alt.length){
			int suffix=0;
			for(int rpos=ref.length-1, apos=alt.length-1; rpos>0 && apos>0 && ref[rpos]==alt[apos]; rpos--, apos--){
				suffix++;
			}
			if(suffix>0){
				reflen=ref.length-suffix;
				int altlen=alt.length-suffix;
				ref=(reflen==1 ? Var.AL_MAP[ref[0]] : Arrays.copyOf(ref, reflen));
				alt=(altlen==1 ? Var.AL_MAP[alt[0]] : Arrays.copyOf(alt, altlen));
				assert(alt.length>0  && ref.length>0);
			}
		}
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		qual=(line[a]=='.' && b==a+1 ? 40 : Tools.parseDouble(line, a, b));
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		filter=Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		info=Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		int TYP=Tools.indexOfDelimited(info, "TYP=", 0, (byte)';');
		type=(TYP>=0 ? Var.typeInitialArray[info[TYP+4]] : type_old());
		assert(type>=0) : type+", "+TYP+"\n"+new String(info);
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		format=Arrays.copyOfRange(line, a, b);
		b++;
		a=b;
		
		while(b<line.length){
			while(b<line.length && line[b]!='\t'){b++;}
			if(b<=a){
				break;
			}
			byte[] sample=Arrays.copyOfRange(line, a, b);
			samples.add(sample);
			b++;
			a=b;
		}
		
		hashcode=hash();
		if(AUTOCACHE){cache();}
	}
	
	public Var toVar(){
		return makeVar(info, alt);
	}
	
	public static Var makeVar(byte[] info, byte[] alt){
		int a=0, b=0;
		
		//SN=0;STA=547693;STO=547694;TYP=SUB;
		assert(Tools.startsWith(info, "SN", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 0: "+new String(info);
		int scaf=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 1: "+new String(info);
		int start=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 2: "+new String(info);
		int stop=Tools.parseInt(info, a, b);
		b++;
		a=b;

//		assert(Tools.startsWith(info, "TYP", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 3: "+new String(info);
		final int type=Var.typeInitialArray[info[a]];
		assert(type>=0) : type+new String(info);
//		if(Tools.contains(info, SUB, a)){type=Var.SUB;}
//		else if(Tools.contains(info, DEL, a)){type=Var.DEL;}
//		else if(Tools.contains(info, INS, a)){type=Var.INS;}
//		else if(Tools.contains(info, NOCALL, a)){type=Var.NOCALL;}
//		else{assert(false) : new String(info);}
		b++;
		a=b;
		
		//R1P=20;R1M=29;R2P=25;R2M=19;
		assert(Tools.startsWith(info, "R1P", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 4: "+new String(info);
		int r1p=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 5: "+new String(info);
		int r1m=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 6: "+new String(info);
		int r2p=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 7: "+new String(info);
		int r2m=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//AD=2;DP=24;MCOV=0;PPC=0;
		assert(Tools.startsWith(info, "AD=", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 8: "+new String(info);
//		int ad=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 9: "+new String(info);
		int cov=Tools.parseInt(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 10: "+new String(info);
		int mcov=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "PPC", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 11: "+new String(info);
		int pc=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		//AF=0.0833;RAF=0.0833;LS=280;
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 12: "+new String(info);
//		double af=Tools.parseDouble(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "RAF", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 13: "+new String(info);
		double raf=Tools.parseDouble(info, a, b);
		b++;
		a=b;
		
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 14: "+new String(info);
		long ls=Tools.parseLong(info, a, b);
		b++;
		a=b;
		
		//MQS=86;MQM=43;BQS=64;BQM=32;
		assert(Tools.startsWith(info, "MQS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 15: "+new String(info);
		long mqs=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 16: "+new String(info);
		int mqm=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "BQS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 17: "+new String(info);
		long bqs=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 18: "+new String(info);
		int bqm=Tools.parseInt(info, a, b);
		b++;
		a=b;

		//EDS=18;EDM=9;IDS=1984;IDM=992;
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 19: "+new String(info);
		long eds=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 20: "+new String(info);
		int edm=Tools.parseInt(info, a, b);
		b++;
		a=b;
		
		assert(Tools.startsWith(info, "IDS", a));
		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 21: "+new String(info);
		long ids=Tools.parseLong(info, a, b);
		b++;
		a=b;

		while(b<info.length && info[b]!='='){b++;}
		a=b+1;
		while(b<info.length && info[b]!=';'){b++;}
		assert(b>a) : "Missing field 22: "+new String(info);
		int idm=Tools.parseInt(info, a, b);
		b++;
		a=b;
//		//CED=13;HMP=1;SB=0.9980;DP=24;DP4=22,0,2,0
//		assert(Tools.startsWith(info, "CED", a));
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 23: "+new String(info);
////		int ced=Tools.parseInt(info, a, b);
//		b++;
//		a=b;
//		
//		//HMP=2;AF=0.989;
//		assert(Tools.startsWith(info, "HMP", a));
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 24: "+new String(info);
////		int hmp=Tools.parseInt(info, a, b);
//		b++;
//		a=b;
//		
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 25: "+new String(info);
////		double sb=Tools.parseDouble(info, a, b);
//		b++;
//		a=b;
//
//		while(b<info.length && info[b]!='='){b++;}
//		a=b+1;
//		while(b<info.length && info[b]!=';'){b++;}
//		assert(b>a) : "Missing field 26: "+new String(info);
////		String dp4=new String(info, a, b-a);
//		b++;
//		a=b;
		
		if(type==Var.DEL || type==Var.INS){
			if(alt.length<=1){alt=Var.AL_0;}
			else if(alt.length==2){alt=Var.AL_MAP[alt[1]];}
			else{alt=Arrays.copyOfRange(alt, 1, alt.length);}
		}
		
		//GT:DP:AD:AF	1:94:93:0.989

		Var v=new Var(scaf, start, stop, alt, type);
		v.r1plus=r1p;
		v.r1minus=r1m;
		v.r2plus=r2p;
		v.r2minus=r2m;
		v.properPairCount=pc;
		v.lengthSum=ls;
		v.mapQSum=mqs;
		v.mapQMax=mqm;
		v.baseQSum=bqs;
		v.baseQMax=bqm;
		v.endDistSum=eds;
		v.endDistMax=edm;
		v.idSum=ids;
		v.idMax=idm;
		v.revisedAlleleFraction=raf;
		v.setCoverage(cov, mcov);
//		v.homopolymerCount=hmp; //derived
		
		return v;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b){
		return equals((VCFLine)b);
	}
	
	public boolean equals(VCFLine b){
		return hashcode==b.hashcode && compareTo(b)==0;
	}

	private int hash(){
		return scaf.hashCode()^Integer.rotateLeft(pos, 9)^Integer.rotateRight(pos+ref.length, 9)^Var.hash(alt);
	}
	
	@Override
	public int hashCode(){
		return hashcode;
	}
	
	public long toKey() {
		long key=Long.rotateLeft(pos, 31)^Long.rotateRight(hashcode, 10)^scaf.hashCode();
		return key&0x3FFFFFFFFFFFFFFFL;
	}
	
	@Override
	public int compareTo(VCFLine v){
		ScafMap map=ScafMap.defaultScafMap();
		assert(map!=null);
		int scafnum1=map.getScaffold(scaf).number;
		int scafnum2=map.getScaffold(v.scaf).number;
		if(scafnum1!=scafnum2){return scafnum1-scafnum2;}
		if(pos!=v.pos){return pos-v.pos;}
		final int typeA=type(), typeB=v.type();
		if(typeA!=typeB){return typeA-typeB;}
		int stop1=pos+reflen(), stop2=v.pos+reflen();
		if(stop1!=stop2){return stop1-stop2;}
		return compare(alt, v.alt);
	}
	
	public int compare(byte[] a, byte[] b){
		if(a==b){return 0;}
		if(a.length!=b.length){return b.length-a.length;}
		for(int i=0; i<a.length; i++){
			byte ca=a[i], cb=b[i];
			if(ca!=cb){return ca-cb;}
		}
		return 0;
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		return toText(bb).toString();
	}
	
	public ByteBuilder toText(ByteBuilder bb){
		bb.append(scaf).append('\t');
		bb.append(pos).append('\t');
		bb.append(id).append('\t');
		bb.append(ref).append('\t');
		bb.append(alt).append('\t');
		bb.append(qual, 2).append('\t');
		bb.append(filter).append('\t');
		bb.append(info).append('\t');
		bb.append(format);
		for(byte[] sample : samples){
			bb.tab().append(sample);
		}
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/
	
	public int reflen(){return ref.length;}
	public int readlen(){return alt.length;}
	
	public int type_old(){
		int reflen=reflen(), readlen=readlen();
		if(reflen<readlen){return Var.INS;}
		if(reflen>readlen){return Var.DEL;}
		for(byte b : alt){
			if(b!='N'){return Var.SUB;}
		}
		return Var.NOCALL;
	}
	
	public int type(){return type;}
	
	public boolean isJunction(){
		return type==Var.LJUNCT || type==Var.RJUNCT || type==Var.BJUNCT;
	}
	
	public boolean isIndel(){
		return type==Var.INS || type==Var.DEL;
	}
	
	void cache(){
//		assert(false) : AUTOCACHE;
		id=cache(id);
		if(ref.length<5){ref=cache(ref);}
		if(alt.length<5){alt=cache(alt);}
		filter=cache(filter);
		ref=cache(ref);
		format=cache(format);
	}
	
	static byte[] cache(byte[] line){
		if(line==null){return line;}
		String s=new String(line);
		byte[] old=cache.get(s);
		if(old!=null){return old;}
		if(cache.size()>20000){return line;}
		synchronized(cache){
			old=cache.get(s);
			if(old!=null){return old;}
			cache.put(s, line);
			return line;
		}
	}
	
	static byte[] cache(String s){
		if(s==null){return null;}
		byte[] old=cache.get(s);
		if(old!=null){return old;}
		synchronized(cache){
			old=cache.get(s);
			if(old!=null){return old;}
			old=s.getBytes();
			cache.put(s, old);
			return old;
		}
	}
	
	/** 0-based */
	public int start(){
		return pos-1;
	}
	
	/** 0-based */
	public int stop(){
		return Tools.max(start(), pos+reflen-2);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String scaf;
	public final int pos;
	public byte[] id;
	public byte[] ref;
	public int reflen;
	public byte[] alt;
	public double qual;
	public byte[] filter;
	public final byte[] info;
	public byte[] format;
	public final int hashcode;
	public final int type;
	public ArrayList<byte[]> samples=new ArrayList<byte[]>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static HashMap<String, byte[]> cache=new HashMap<String, byte[]>(99997);
	
	static boolean AUTOCACHE=false;
	static boolean TRIM_TO_CANONICAL=true;
	
	private static final byte[] NOCALL=cache("NOCALL");
	private static final byte[] SUB=cache("SUB");
	private static final byte[] DEL=cache("DEL");
	private static final byte[] INS=cache("INS");
	private static final byte[] LJUNCT=cache("LJUNCT");
	private static final byte[] RJUNCT=cache("RJUNCT");
	private static final byte[] BJUNCT=cache("BJUNCT");
	private static final byte[] DOT=cache(".");
	private static final byte[] PASS=cache("PASS");
	private static final byte[] FAIL=cache("FAIL");
	private static final byte[] FORMAT=cache("GT:DP:AD:AF:SC:PF");
	
	static{
		cache(Var.AL_0);
		cache(Var.AL_A);
		cache(Var.AL_C);
		cache(Var.AL_G);
		cache(Var.AL_T);
		cache(Var.AL_N);
	}

}
