package var2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import align2.QualityTools;
import dna.AminoAcid;
import fileIO.FileFormat;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Tracks data for a variation.
 * Uses half-open coordinates.
 * @author Brian Bushnell
 * @date November 4, 2016
 *
 */
public class Var implements Comparable<Var>, Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 3536617113257471595L;

	public static void main(String[] args){
		if(args.length>1){
			for(int i=1; i<args.length; i++){			
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				Parser.parseCommonStatic(arg, a, b);
			}
		}
		Timer t=new Timer();
		VarMap vmap;
		FileFormat ff=FileFormat.testInput(args[0], FileFormat.VAR, null, true, true);
		if(ff.vcf()){
			ScafMap smap=ScafMap.loadVcfHeader(ff);
			vmap=VcfLoader.loadVcf(args[0], smap, false, false);
		}else{
			vmap=VcfLoader.loadVars(args[0], null);
		}
		t.stop("Loaded "+vmap.size()+" variants.\nTime: \t");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Var(int scafnum_, int start_, int stop_, int allele_, int type_){
		this(scafnum_, start_, stop_, AL_MAP[allele_], type_);
	}

	public Var(Var v) {
		this(v.scafnum, v.start, v.stop, v.allele, v.type);
	}
	
	public Var(int scafnum_, int start_, int stop_, byte[] allele_, int type_){
		scafnum=scafnum_;
		start=start_;
		stop=stop_;
		allele=allele_;
		type=type_;
		hashcode=hash();
//		stamp=stamper.getAndIncrement();
		assert(allele.length>1 || allele==AL_0 ||
				allele==AL_A || allele==AL_C || allele==AL_G || allele==AL_T || allele==AL_N) : new String(allele_+", "+allele_.length);
		assert(start<=stop) : "\n"+Var.toBasicHeader()+"\n"+this+"\n";
//		if(type()==SUB && allele.length==1){ //TODO: 123 - mainly for testing
//			final byte call=allele[0];
//			final Scaffold scaf=ScafMap.defaultScafMap.getScaffold(scafnum);
//			final byte ref=scaf.bases[start];
////			System.err.println((char)call+"="+(char)ref+" at scaf "+scafnum+" pos "+start);
//			assert(ref!=call) : (char)call+"="+(char)ref+" at scaf "+scafnum+" pos "+start;
//		}
//		assert(false) : this.alleleCount();
	}
	
	public Var(final byte[] line, final byte delimiter){
		int a=0, b=0;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		scafnum=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		start=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		stop=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		type=typeInitialArray[line[a]];
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>=a) : "Missing field 4: "+new String(line);
		if(b==a){allele=AL_0;}
		else if(b==a+1){allele=AL_MAP[line[a]];}
		else{allele=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		r1plus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		r1minus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		r2plus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		r2minus=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 9: "+new String(line);
		properPairCount=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 10: "+new String(line);
		lengthSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 11: "+new String(line);
		mapQSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 12: "+new String(line);
		mapQMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 13: "+new String(line);
		baseQSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 14: "+new String(line);
		baseQMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 15: "+new String(line);
		endDistSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 16: "+new String(line);
		endDistMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 17: "+new String(line);
		idSum=Tools.parseLong(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 18: "+new String(line);
		idMax=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 19: "+new String(line);
		coverage=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 20: "+new String(line);
		minusCoverage=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
//		assert(b>a) : "Missing field 21: "+new String(line);
//		int contigEndDist=Tools.parseInt(line, a, b); //Unused
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		if(b>a){
//			phredScore=Tools.parseFloat(line, a, b);
			b++;
			a=b;
		}
		
		hashcode=hash();
		assert(allele.length>1 || allele==AL_0 ||
				allele==AL_A || allele==AL_C || allele==AL_G || allele==AL_T || allele==AL_N);
		assert(start<=stop) : this.toString();
		assert(type>=LJUNCT || type==type_old()) : type+", "+type_old()+", "+this.toString();
	}

	//#CHROM POS    ID        REF  ALT     QUAL
	public static Var fromVCF(byte[] line, ScafMap scafMap, boolean parseCoverage, boolean parseExtended) {
		int a=0, b=0;
		
		//CHROM
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		String scaf=new String(line, a, b-a);
		b++;
		a=b;
		
		//POS
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		int pos=Tools.parseInt(line, a, b);
		b++;
		a=b;

		//ID
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		//String id=new String(line, a, b-a);
		b++;
		a=b;
		
		//REF
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		//byte[] ref=Arrays.copyOf(line, a, b);
		int reflen=line[a]=='.' ? 0 : b-a;
		b++;
		a=b;

		//ALT
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		byte[] alt;
		if(b<=a+1){alt=AL_MAP[line[a]];}
		else{alt=Arrays.copyOfRange(line, a, b);}
		b++;
		a=b;

		//QUAL, FILTER, INFO
		
		int infoStart=b;
		
//		while(b<line.length && line[b]!='\t'){b++;}
//		assert(b>a) : "Missing field 5: "+new String(line);
//		int qual=Tools.parseInt(line, a, b);
//		b++;
//		a=b;
		
		final int start;
		final int readlen;
		if(alt.length!=reflen && alt.length>0){
			alt=Arrays.copyOfRange(alt, 1, alt.length);
			start=pos;
			if(alt.length==0){alt=AL_0;}
			else if(alt.length==1 && AL_MAP[alt[0]]!=null){alt=AL_MAP[alt[0]];}
		}else{
			start=pos-1;
		}
		readlen=alt.length;
		final int stop=start+reflen;
		assert(scaf!=null);
		assert(scafMap!=null);
		final int scafNum=scafMap.getNumber(scaf);
		assert(scafNum>=0) : scaf+"\n"+scafMap.keySet()+"\n"+scafMap.altKeySet()+"\n";
//		final Scaffold scaffold=scafMap.getScaffold(scafNum);
		
		int type=-1;
		if(parseExtended){
			type=Tools.max(type, parseVcfIntDelimited(line, "TYP=", infoStart));
			if(type<0){type=typeStartStop(start, stop, alt);}
		}else{
			type=typeStartStop(start, stop, alt);
		}
		Var v=new Var(scafNum, start, stop, alt, type);
//		System.err.println(new String(line)+"\n"+v.toString()+"\ntype="+type+"\n");
		
		if(parseCoverage){
			infoStart=Tools.indexOfDelimited(line, "R1P=", infoStart, (byte)';');
			
			//R1P=2;R1M=0;R2P=0;R2M=0;
			int r1p=Tools.max(0, parseVcfIntDelimited(line, "R1P=", infoStart));
			int r1m=Tools.max(0, parseVcfIntDelimited(line, "R1M=", infoStart));
			int r2p=Tools.max(0, parseVcfIntDelimited(line, "R2P=", infoStart));
			int r2m=Tools.max(0, parseVcfIntDelimited(line, "R2M=", infoStart));
			
			//AD=2;DP=24;MCOV=0;PPC=0;
			int cov=parseVcfIntDelimited(line, "DP=", infoStart);
			assert(cov>0) : new String(line, infoStart, line.length-infoStart);
			int mcov=parseVcfIntDelimited(line, "MCOV=", infoStart);
			
			v.coverage=cov;
			v.minusCoverage=mcov;
			v.r1plus=r1p;
			v.r1minus=r1m;
			v.r2plus=r2p;
			v.r2minus=r2m;
		}
		
		if(parseExtended){
			infoStart=Tools.indexOfDelimited(line, "PPC=", infoStart, (byte)';');
			
			//Some extended fields
			//AD=2;DP=24;MCOV=0;PPC=0;
			int pc=Tools.max(0, parseVcfIntDelimited(line, "PPC=", infoStart));
			
			//AF=0.0833;RAF=0.0833;LS=280;
			double raf=parseVcfDoubleDelimited(line, "RAF=", infoStart);
			long ls=Tools.max(0, parseVcfLongDelimited(line, "LS=", infoStart));
	
			//MQS=86;MQM=43;BQS=64;BQM=32;
			long mqs=Tools.max(0, parseVcfLongDelimited(line, "MQS=", infoStart));
			int mqm=Tools.max(0, parseVcfIntDelimited(line, "MQM=", infoStart));
			long bqs=Tools.max(0, parseVcfLongDelimited(line, "BQS=", infoStart));
			int bqm=Tools.max(0, parseVcfIntDelimited(line, "BQM=", infoStart));
			
			//EDS=18;EDM=9;IDS=1984;IDM=992;
			long eds=Tools.max(0, parseVcfLongDelimited(line, "EDS=", infoStart));
			int edm=Tools.max(0, parseVcfIntDelimited(line, "EDM=", infoStart));
			long ids=Tools.max(0, parseVcfLongDelimited(line, "IDS=", infoStart));
			int idm=Tools.max(0, parseVcfIntDelimited(line, "IDM=", infoStart));
			
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
		}
		
		return v;
	}
	
	//Parses INFO field
	private static int parseVcfIntDelimited(byte[] line, String query, int start){
		return (int)Tools.min(Integer.MAX_VALUE, parseVcfLongDelimited(line, query, start));
	}
	
	//Parses INFO field
	private static long parseVcfLongDelimited(byte[] line, String query, final int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
//		System.err.println(loc+", "+(line.length)+", "+(line.length-loc));
//		System.err.println(loc+", "+new String(line, loc, line.length-loc));
//		if(true){KillSwitch.kill();}
		if(loc<0){return -1;}
		long current=0;
		long mult=1;
		if(loc>0){
			if(line[loc+query.length()]=='-'){mult=-1; loc++;}
			for(int i=loc+query.length(); i<line.length; i++){
				final byte x=line[i];
				if(Tools.isDigit(x)){
					current=current*10+(x-'0');
				}else{
					assert(x==tab || x==colon) : x+", "+query+", "+new String(line, loc, i-loc+1);
					break;
				}
			}
		}
		return mult*current;
	}
	
	//Parses INFO field
	private static double parseVcfDoubleDelimited(byte[] line, String query, int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
		if(loc<0){return -1;}
		if(loc>0){
			loc=loc+query.length();
			int loc2=loc+1;
			while(loc2<line.length){
				byte b=line[loc2];
				if(b==tab || b==colon){break;}
				loc2++;
			}
			return Tools.parseDouble(line, loc, loc2);
		}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	//Only called by MergeSamples from VCFLine arrays.
	public void addCoverage(Var b){
		assert(this.equals(b));
		coverage+=b.coverage;
		minusCoverage+=b.minusCoverage;
	}
	
	//Called in various places such as when processing each read
	public void add(Var b){
		final int oldReads=alleleCount();
		
//		assert(oldReads>0) : this;
		assert(oldReads==0 || baseQSum/oldReads<=60) : this;
		
		assert(this.equals(b));
		r1plus+=b.r1plus;
		r1minus+=b.r1minus;
		r2plus+=b.r2plus;
		r2minus+=b.r2minus;
		properPairCount+=b.properPairCount;
		lengthSum+=b.lengthSum;
		
		mapQSum+=b.mapQSum;
		mapQMax=Tools.max(mapQMax, b.mapQMax);
		baseQSum+=b.baseQSum;
		baseQMax=Tools.max(baseQMax, b.baseQMax);

		endDistSum+=b.endDistSum;
		endDistMax=Tools.max(endDistMax, b.endDistMax);

		idSum+=b.idSum;
		idMax=Tools.max(idMax, b.idMax);

//		assert(count()>0 && count()>oldReads) : "\n"+this+"\n"+b;
		assert(alleleCount()>=oldReads) : "\n"+this+"\n"+b;
		assert(alleleCount()==oldReads+b.alleleCount()) : "\n"+this+"\n"+b;
		assert(alleleCount()==0 || baseQSum/alleleCount()<=60) : "\n"+this+"\n"+b;
	}
	
	public void add(Read r){
		SamLine sl=(SamLine)r.obj;
		final int bstart=calcBstart(r, sl);
		final int bstop=calcBstop(bstart, r);
		add(r, bstart, bstop);
	}
		
	public void add(Read r, final int bstart, final int bstop){

		final int oldReads=alleleCount();
		
		SamLine sl=(SamLine)r.obj;
		
		if(sl.strand()==0){
			if(sl.pairnum()==0){
				r1plus++;
			}else{
				r2plus++;
			}
		}else{
			if(sl.pairnum()==0){
				r1minus++;
			}else{
				r2minus++;
			}
		}
		
		lengthSum+=r.length();
		properPairCount+=(sl.properPair() ? 1 : 0);
		mapQSum+=sl.mapq;
		mapQMax=Tools.max(mapQMax, sl.mapq);
		
		int baseQ=calcBaseQ(bstart, bstop, r, sl);
		baseQSum+=baseQ;
		baseQMax=Tools.max(baseQMax, baseQ);
		
		int endDist=calcEndDist(bstart, bstop, r);
		endDistSum+=endDist;
		endDistMax=Tools.max(endDistMax, endDist);
		
		int id=(int)(1000*Read.identitySkewed(r.match, false, false, false, true));
		idSum+=id;
		idMax=Tools.max(idMax, id);
		
		assert(alleleCount()>0) : this;
		assert(alleleCount()==oldReads+1) : this;
		assert(baseQSum/alleleCount()<=60) : this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, ScafMap scafMap){
//		if(!r.containsVariants()){return null;}
		final int scafnum=scafMap.getNumber(sl.rnameS());
		return toVars(r, sl, callNs, scafnum);
	}
	
	/**
	 * @TODO This crashes on indels in the last position in the match string.
	 * @param r
	 * @param sl
	 * @param callNs
	 * @param scafnum
	 * @return A list of variants
	 */
	public static ArrayList<Var> toVars(Read r, SamLine sl, boolean callNs, final int scafnum){
		final boolean hasV=r.containsVariants();
		final boolean callSID=(CALL_DEL || CALL_INS || CALL_SUB) && hasV;
		final boolean callJ=(CALL_JUNCTION) && (hasV || r.containsClipping());
		if(!callSID && !callJ){return null;}
		
		r.toLongMatchString(false);
		if(sl.strand()==1 && !r.swapped()){
			r.reverseComplement();
			r.setSwapped(true);
		}
		
		ArrayList<Var> sidList=null, jList=null;
		if(callSID){sidList=toSubsAndIndels(r, sl, callNs, scafnum);}
		if(callJ){jList=toJunctions(r, sl, scafnum, hasV, 8);}
		
		if(sidList!=null){
			if(jList!=null){sidList.addAll(jList);}
			return sidList;
		}else{
			return jList;
		}
		

		//Note: Did not un-rcomp the read
	}
	
	/**
	 * @TODO This crashes on indels in the last position in the match string.
	 * @param r
	 * @param sl
	 * @param callNs
	 * @param scafnum
	 * @return A list of variants
	 */
	private static ArrayList<Var> toSubsAndIndels(Read r, SamLine sl, boolean callNs, final int scafnum){
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final int rpos0=sl.pos-1;
		ArrayList<Var> list=new ArrayList<Var>();

		int bstart=-1, bstop=-1;
		int rstart=-1, rstop=-1;
		int mode=-1;
		
		int mpos=0, bpos=0, rpos=rpos0;
		for(; mpos<match.length; mpos++){
			byte m=match[mpos];
			
			if(m!=mode){
				if(mode=='D'){
					bstop=bpos;
					rstop=rpos;
//					assert(false) : (char)m+", "+(char)mode+", "+rstart+", "+bstart;
					if(CALL_DEL){
						Var v=new Var(scafnum, rstart, rstop, 0, DEL);
						v.add(r, bstart, bstop);
						list.add(v);
					}
					bstart=bstop=rstart=rstop=-1;
				}else if(mode=='I'){
					bstop=bpos;
					rstop=rpos;
					int blen=bstop-bstart;
					if(CALL_INS){
						Var v;
						if(blen==1){
							v=new Var(scafnum, rstart, rstop, bases[bstart], INS);
						}else{
							v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop), INS);
						}
						v.add(r, bstart, bstop);
						list.add(v);
					}
					bstart=bstop=rstart=rstop=-1;
				}
			}
			
			if(m=='C'){
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){
				if(m=='S' || (m=='N' && callNs)){
//					System.err.println(sl.toString()+"\n"+ScafMap.defaultScafMap.getScaffold(scafnum).getSequence(sl)); //123
//					System.err.println(ScafMap.defaultScafMap.getScaffold(scafnum));
					if(CALL_SUB){
						Var v=new Var(scafnum, rpos, rpos+1, bases[bpos], SUB);
						//					System.err.println(v.toString());
						v.add(r, bpos, bpos+1);
						list.add(v);

						if(TEST_REF_VARIANTS && v.type()==SUB && v.allele.length==1){ //TODO: 123 - mainly for testing
							final byte call=v.allele[0];
							final Scaffold scaf=ScafMap.defaultScafMap().getScaffold(scafnum);
							final byte ref=scaf.bases[v.start];
							//						System.err.println((char)call+"="+(char)ref+" at scaf "+scafnum+" pos "+start);
							assert(ref!=call) : (char)call+"="+(char)ref+" at scaf "+scafnum+" pos "+v.start+"\n"
							+sl+"\n"+ScafMap.defaultScafMap().getScaffold(scafnum).getSequence(sl)+"\n";
						}
					}
					
				}
				bpos++;
				rpos++;
			}else if(m=='D'){
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				rpos++;
			}else if(m=='I'){
				if(mode!=m){
					rstart=rpos;
					bstart=bpos;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
			mode=m;
		}
		
		if(mode=='D'){
			bstop=bpos;
			rstop=rpos;
			if(CALL_DEL){
				Var v=new Var(scafnum, rstart, rstop, 0, DEL);
				v.add(r, bstart, bstop);
				list.add(v);
			}
			bstart=bstop=rstart=rstop=-1;
		}else if(mode=='I'){
			bstop=bpos;
			rstop=rpos-1;
			int blen=bstop-bstart;
			if(CALL_INS){
				Var v;
				assert(rstart<=rstop) : "\n"+rstart+", "+rstop+", "+rpos+
				"\n"+bstart+", "+bstop+", "+bpos+
				"\n"+r+"\n"+sl;
				if(blen==1){
					v=new Var(scafnum, rstart, rstop, bases[bstart], INS);
				}else{
					v=new Var(scafnum, rstart, rstop, Arrays.copyOfRange(bases, bstart, bstop), INS);
				}
				v.add(r, bstart, bstop);
				list.add(v);
			}
			bstart=bstop=rstart=rstop=-1;
		}
		
		return list;
	}
	
	/** Long-match */
	private static int countLeftClip(byte[] longmatch){
		for(int i=0; i<longmatch.length; i++){
			if(longmatch[i]!='C'){return i;}
		}
		return longmatch.length;
	}
	
	/** Long-match */
	private static int countRightClip(byte[] longmatch){
		for(int i=0, j=longmatch.length-1; i<longmatch.length; i++, j--){
			if(longmatch[j]!='C'){return i;}
		}
		return longmatch.length;
	}
	
	private static ArrayList<Var> toJunctions(Read r, SamLine sl, final int scafnum,
			boolean containsVars, final int minClip){
		final byte[] match0=r.match, match;
		final byte[] bases=r.bases;
//		final int rpos0=sl.pos-1;
		final int start, stop;
		int leftClip=countLeftClip(match0), rightClip=countRightClip(match0);
		if(leftClip==0 && rightClip==0){//try soft-clipping
			int[] rvec=new int[2];
			byte[] smatch=SoftClipper.softClipMatch(match0, minClip, false, r.start, r.stop, rvec);
			if(smatch==null){
				return null;
			}else{
				start=rvec[0];
				stop=rvec[1];
				match=smatch;
				leftClip=countLeftClip(match);
				rightClip=countRightClip(match);
			}
		}else{
			if(leftClip<minClip && rightClip<minClip){return null;}
			start=r.start;
			stop=r.stop;
			match=match0;
		}
		
		ArrayList<Var> list=new ArrayList<Var>();
		if(leftClip>=minClip){//LJUNCT
			int bpos=leftClip-1;
			byte jcall=bases[bpos];
			int jpos=start+leftClip;
			Var v=new Var(scafnum, jpos, jpos+1, jcall, LJUNCT);
			v.add(r, bpos, bpos+1);
			list.add(v);
		}
		if(rightClip>=minClip){//RJUNCT
			int bpos=bases.length-rightClip;
			byte jcall=bases[bpos];
			int jpos=stop-rightClip+1;
			Var v=new Var(scafnum, jpos, jpos+1, jcall, RJUNCT);
			v.add(r, bpos, bpos+1);
			list.add(v);
		}
		return list;
	}

	public int calcBstart(Read r, SamLine sl){
		r.toLongMatchString(false);
		byte[] match=r.match;
		final int rstart=sl.pos-1;
		final int type=type();
		
		int bstart=-1;
		
		for(int mpos=0, rpos=rstart, bpos=0; mpos<match.length; mpos++){
			byte m=match[mpos];
			if(m=='C'){
				bpos++;
			}else if(m=='m' || m=='S' || m=='N'){
				if(rpos==rstart){
					assert(type==SUB || type==NOCALL) : type+", "+bpos+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				bpos++;
				rpos++;
			}else if(m=='D'){
				if(rpos==rstart){
					assert(type==DEL) : type+", "+rpos+"\n"+new String(match);
					bstart=bpos;
					break;
				}
				rpos++;
			}else if(m=='I'){
				if(rpos==rstart && type==INS){
					bstart=bpos;
					break;
				}
				bpos++;
			}else{
				assert(false) : "Unhandled symbol "+(char)m;
			}
		}
		assert(bstart>=0);
		return bstart;
	}
	
	public int calcBstop(int bstart, Read r){
		assert(bstart>=0);
		int bstop=bstart+readlen();
		assert(bstop<=r.length());
		return bstop;
	}
	
	public int calcEndDist(int bstart, int bstop, Read r){
		int dist=Tools.min(bstart, r.length()-bstop);
		assert(dist>=0 && dist<=r.length()/2) : dist+", "+r.length()+", "+bstart+", "+bstop+"\n"+this+"\n"+
			"\n"+new String(r.match)+"\n"+r.obj+"\n";
		assert(dist<=(r.length()-readlen())/2);
		return dist;
	}
	
	public int calcBaseQ(final int bstart0, final int bstop0, Read r, SamLine sl){
		final byte[] quals=r.quality;
		if(quals==null){return Shared.FAKE_QUAL;}
		final int type=type();
		final int bstart, bstop;
		final int len=r.length();
		
		if(sl.strand()==0 || (sl.strand()==1 && r.swapped())){
			bstart=bstart0;
			bstop=bstop0;
		}else{
			bstart=len-bstop0-1;
			bstop=len-bstart0-1;
			assert(bstop-bstart==bstop0-bstart0);
		}
		
		int sum=0, avg=0;
		if(type==DEL){
			if(bstart==0){
				sum=avg=quals[0];
			}else if(bstop>=len-1){
				sum=avg=quals[len-1];
			}else{
				assert(bstop==bstart) : bstart0+", "+bstop0+", "+bstart+", "+bstop+"\n"+
						r.length()+", "+r.swapped()+", "+type()+", "+readlen()+", "+reflen()+
						"\n"+this+"\n"+new String(r.match)+"\n"+r.obj+"\n";
				
//				-1, 73, -1, 73
//				151, true, 2, 0, 1
				
				sum=quals[bstart]+quals[bstop+1];
				avg=sum/2;
			}
		}else{
			for(int i=bstart; i<bstop; i++){
				sum+=quals[i];
			}
			avg=Math.round(sum/(bstop-bstart));
		}
		return avg;
	}

	public int reflen(){
		return stop-start;
	}
	
	int readlen(){
		return (allele==null || allele.length==0 || allele[0]=='.' ? 0 : allele.length);
	}
	
	public int type(){return type;}
	
	int type_old(){
		int reflen=reflen(), readlen=readlen();
		return typeReadlenReflen(readlen, reflen, allele);
	}
	
	static int typeStartStop(int start, int stop, byte[] allele){
		final int readlen=(allele.length==0 || allele[0]=='.' ? 0 : allele.length);
		final int reflen=stop-start;
		return typeReadlenReflen(readlen, reflen, allele);
	}
	
	static int typeReadlenReflen(int readlen, int reflen, byte[] allele){
		if(reflen<readlen){return INS;}
		if(reflen>readlen){return DEL;}
//		assert(readlen>0) : start+", "+stop+", "+reflen+", "+readlen+", "+new String(allele);
//		if(reflen==0){return INS;}
//		if(readlen==0){return DEL;}
//		assert(start<=stop) : start+", "+stop;
//		assert(reflen==readlen) : reflen+", "+readlen+", "+new String(allele)+", "+start+", "+stop;
		for(byte b : allele){
			if(b!='N'){return SUB;}
		}
		return NOCALL;
	}
	
	String typeString(){
		return typeArray[type()];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Contract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b){
		return equals((Var)b);
	}
	
	public boolean equals(Var b){
		return hashcode==b.hashcode && compareTo(b)==0;
	}
	
	@Override
	public int hashCode(){
		return hashcode;
	}
	
	public long toKey() {
		long key=Long.rotateLeft(start, 30)^Long.rotateLeft(scafnum, 10)^hash(allele);
		return key&0x3FFFFFFFFFFFFFFFL;
	}
	
	@Override
	public int compareTo(Var v){
		if(scafnum!=v.scafnum){return scafnum-v.scafnum;}
		final int typeA=type(), typeB=v.type();
		int stA=start+(typeA==DEL ? -1 : 0);
		int stB=v.start+(typeB==DEL ? -1 : 0);
		if(stA!=stB){return stA-stB;}
		if(typeA!=typeB){return typeA-typeB;}
		if(stop!=v.stop){return stop-v.stop;}
		return compare(allele, v.allele);
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
		return toTextQuick(new ByteBuilder()).toString();
	}

	public ByteBuilder toTextQuick(ByteBuilder bb){
		return toText(bb, 0.99f, 30, 30, 150, 1, 2, null);
	}
	
	public static String toVarHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();
		
		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("#fileformat\tVar_"+varFormat+"\n");
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
		
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tcontigEndDist");
		sb.append("\tphredScore");
		if(extendedText){
			sb.append("\treadCount\talleleFraction\trevisedAF\tstrandRatio\tbaseqAvg\tmapqAvg\tedistAvg\tidentityAvg");
			sb.append("\tedistScore\tidentityScore\tqualityScore\tpairedScore\tbiasScore\tcoverageScore\thomopolymerScore\tscore");
		}
		return sb.toString();
	}
	
	public static String toBasicHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tcontigEndDist");
		sb.append("\tphredScore");
		return sb.toString();
	}
	
	public ByteBuilder toText(ByteBuilder bb, double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
		useIdentity=true;
		bb.append(scafnum).append('\t');
		bb.append(start).append('\t');
		bb.append(stop).append('\t');
		bb.append(typeArray[type()]).append('\t');
		for(byte b : allele){bb.append(b);}
		bb.tab();
		
		bb.append(r1plus).append('\t');
		bb.append(r1minus).append('\t');
		bb.append(r2plus).append('\t');
		bb.append(r2minus).append('\t');
		bb.append(properPairCount).append('\t');
		bb.append(lengthSum).append('\t');

		bb.append(mapQSum).append('\t');
		bb.append(mapQMax).append('\t');
		bb.append(baseQSum).append('\t');
		bb.append(baseQMax).append('\t');
		bb.append(endDistSum).append('\t');
		bb.append(endDistMax).append('\t');
		bb.append(idSum).append('\t');
		bb.append(idMax).append('\t');

		bb.append(coverage).append('\t');
		bb.append(minusCoverage).append('\t');
		
		final int scafEndDist=!doNscan ? nScan : (map==null ? start : contigEndDist(map));
		bb.append(scafEndDist).append('\t');

//		bb.append(prevBase<0 ? 'N' : (char)prevBase).append('\t');
		
		final double score=score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map);
//		bb.append(String.format(Locale.ROOT, "%.2f", toPhredScore(score))).append('\t');
		bb.append(toPhredScore(score), 2).append('\t');
		
//		if(extendedText){
//			
//			bb.append(count()).append('\t');
//			double af=alleleFraction();
//			bb.append(String.format(Locale.ROOT, "%.4f\t", af));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",revisedAlleleFraction(af, readLengthAvg)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",strandRatio()));
//			bb.append(String.format(Locale.ROOT, "%.2f\t",baseQAvg()));
//			bb.append(String.format(Locale.ROOT, "%.2f\t",mapQAvg()));
//			bb.append(String.format(Locale.ROOT, "%.2f\t",edistAvg()));
//			bb.append(String.format(Locale.ROOT, "%.2f\t",identityAvg()));
//			
//			bb.append(String.format(Locale.ROOT, "%.4f\t",edistScore()));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",identityScore()));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",qualityScore(totalQualityAvg, totalMapqAvg)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",pairedScore(properPairRate, scafEndDist)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",biasScore(properPairRate, scafEndDist)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",coverageScore(ploidy, rarity, readLengthAvg)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",homopolymerScore(map)));
//			bb.append(String.format(Locale.ROOT, "%.4f\t",score));
//		}
		
		if(extendedText){
			
			bb.append(alleleCount()).append('\t');
			final double af=alleleFraction();
			bb.append(af, 4).append('\t');
			bb.append(revisedAlleleFraction(af, readLengthAvg), 4).append('\t');
			bb.append(strandRatio(), 4).append('\t');
			bb.append(baseQAvg(), 2).append('\t');
			bb.append(mapQAvg(), 2).append('\t');
			bb.append(edistAvg(), 2).append('\t');
			bb.append(identityAvg(), 2).append('\t');
			
			bb.append(edistScore(), 4).append('\t');
			bb.append(identityScore(), 4).append('\t');
			bb.append(qualityScore(totalQualityAvg, totalMapqAvg), 4).append('\t');
			bb.append(pairedScore(properPairRate, scafEndDist), 4).append('\t');
			bb.append(biasScore(properPairRate, scafEndDist), 4).append('\t');
			bb.append(coverageScore(ploidy, rarity, readLengthAvg), 4).append('\t');
			bb.append(homopolymerScore(map), 4).append('\t');
			bb.append(score, 4).append('\t');
		}
		
		bb.length--;
		
		return bb;
	}

	public static String toVcfHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy,
			long reads, long pairs, long properPairs, long bases, String ref, ScafMap map, String sampleName, boolean trimWhitespace) {
		StringBuilder sb=new StringBuilder();
		
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##BBMapVersion="+Shared.BBMAP_VERSION_STRING+"\n");
		sb.append("##ploidy="+ploidy+"\n");
		sb.append(String.format(Locale.ROOT, "##rarity=%.5f\n", rarity));
		sb.append(String.format(Locale.ROOT, "##minallelefraction=%.5f\n", minAlleleFraction));
		sb.append("##reads="+reads+"\n");
		sb.append("##pairedReads="+pairs+"\n");
		sb.append("##properlyPairedReads="+properPairs+"\n");
		sb.append(String.format(Locale.ROOT, "##readLengthAvg=%.3f\n", (bases/Tools.max(reads, 1.0))));
		sb.append(String.format(Locale.ROOT, "##properPairRate=%.5f\n", properPairRate));
		sb.append(String.format(Locale.ROOT, "##totalQualityAvg=%.3f\n", totalQualityAvg));
		sb.append(String.format(Locale.ROOT, "##mapqAvg=%.3f\n", mapqAvg));
		if(ref!=null){sb.append("##reference="+ref+"\n");}
		
		for(Scaffold scaf : map.list){
			String name=scaf.name;
			if(trimWhitespace){name=Tools.trimWhitespace(name);}
			sb.append("##contig=<ID="+name+",length="+scaf.length+">\n");
		}
		
		{
			sb.append("##FILTER=<ID=FAIL,Description=\"Fail\">\n");
			sb.append("##FILTER=<ID=PASS,Description=\"Pass\">\n");
			
			sb.append("##INFO=<ID=SN,Number=1,Type=Integer,Description=\"Scaffold Number\">\n");
			sb.append("##INFO=<ID=STA,Number=1,Type=Integer,Description=\"Start\">\n");
			sb.append("##INFO=<ID=STO,Number=1,Type=Integer,Description=\"Stop\">\n");
			sb.append("##INFO=<ID=TYP,Number=1,Type=String,Description=\"Type\">\n");
			
			sb.append("##INFO=<ID=R1P,Number=1,Type=Integer,Description=\"Read1 Plus Count\">\n");
			sb.append("##INFO=<ID=R1M,Number=1,Type=Integer,Description=\"Read1 Minus Count\">\n");
			sb.append("##INFO=<ID=R2P,Number=1,Type=Integer,Description=\"Read2 Plus Count\">\n");
			sb.append("##INFO=<ID=R2M,Number=1,Type=Integer,Description=\"Read2 Minus Count\">\n");
			
			sb.append("##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
			sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
//			sb.append("##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Coverage\">\n");
			sb.append("##INFO=<ID=MCOV,Number=1,Type=Integer,Description=\"Minus Coverage\">\n");
			sb.append("##INFO=<ID=PPC,Number=1,Type=Integer,Description=\"Paired Count\">\n");

			sb.append("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
			sb.append("##INFO=<ID=LS,Number=1,Type=Integer,Description=\"Length Sum\">\n");

			sb.append("##INFO=<ID=MQS,Number=1,Type=Integer,Description=\"MAPQ Sum\">\n");
			sb.append("##INFO=<ID=MQM,Number=1,Type=Integer,Description=\"MAPQ Max\">\n");
			sb.append("##INFO=<ID=BQS,Number=1,Type=Integer,Description=\"Base Quality Sum\">\n");
			sb.append("##INFO=<ID=BQM,Number=1,Type=Integer,Description=\"Base Quality Max\">\n");
			
			sb.append("##INFO=<ID=EDS,Number=1,Type=Integer,Description=\"End Distance Sum\">\n");
			sb.append("##INFO=<ID=EDM,Number=1,Type=Integer,Description=\"End Distance Max\">\n");
			sb.append("##INFO=<ID=IDS,Number=1,Type=Integer,Description=\"Identity Sum\">\n");
			sb.append("##INFO=<ID=IDM,Number=1,Type=Integer,Description=\"Identity Max\">\n");
			
			sb.append("##INFO=<ID=CED,Number=1,Type=Integer,Description=\"Contig End Distance\">\n");
			sb.append("##INFO=<ID=HMP,Number=1,Type=Integer,Description=\"Homopolymer Count\">\n");
			sb.append("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n");
			sb.append("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Ref+, Ref-, Alt+, Alt-\">\n");
			
			sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
			sb.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
			sb.append("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
			sb.append("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##FORMAT=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
			sb.append("##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n");
			sb.append("##FORMAT=<ID=SC,Number=1,Type=Float,Description=\"Score\">\n");
			sb.append("##FORMAT=<ID=PF,Number=1,Type=String,Description=\"Pass Filter\">\n");
			
		}
		
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		if(sampleName!=null){sb.append('\t').append(sampleName);}
		return sb.toString();
	}
	
	//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1323_1066121
	//ChrIII_A_nidulans_FGSC_A4	133	.	T	C	204	PASS
	//	DP=27;VDB=1.614640e-01;RPB=9.947863e-01;AF1=0.5;AC1=1;DP4=5,8,6,8;MQ=49;FQ=180;PV4=1,1,0.055,1
	//	GT:PL:DP:SP:GQ	0/1:234,0,207:27:0:99
	public ByteBuilder toVCF(ByteBuilder bb, double properPairRate, double totalQualityAvg, double mapqAvg, double readLengthAvg,
			int ploidy, ScafMap map, VarFilter filter, boolean trimWhitespace){
		
		final Scaffold scaf=map.getScaffold(scafnum);
		final byte[] bases=scaf.bases;
		final int reflen=reflen(), readlen=readlen(), type=type();
		final double score=phredScore(properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, filter.rarity, ploidy, map);
		final boolean pass=(filter==null ? true :
			filter.passesFilter(this, properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, ploidy, map));
		
		bb.append(trimWhitespace ? Tools.trimWhitespace(scaf.name) : scaf.name).append('\t');
		boolean indel=(type==INS || type==DEL);
		boolean addPrevBase=true;
		final int vcfStart=start+(indel && addPrevBase ? 0 : 1);
		bb.append(vcfStart).append('\t');
		bb.append('.').append('\t');
		
		byte prevBase=(bases==null ? (byte)'N' : bases[Tools.mid(start-1, 0, bases.length-1)]);
		if(UPPER_CASE_ALLELES){prevBase=(byte) Tools.toUpperCase(prevBase);}
		
		if(addPrevBase){
			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			for(int i=0, rpos=start; i<reflen; i++, rpos++){
				bb.append(bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos]);
			}
			bb.tab();

			if(reflen==0 || allele.length<1){bb.append(prevBase);}
			bb.append(allele).append('\t');
		}else{
			if(reflen==0){
				bb.append('.');
			}else{
				for(int i=0, rpos=start; i<reflen; i++, rpos++){
					char refBase=bases==null || rpos<0 || rpos>=bases.length ? (char)'N' : (char)bases[rpos];
					if(UPPER_CASE_ALLELES){refBase=Tools.toUpperCase(refBase);}
					bb.append(refBase);
				}
			}
			bb.tab();

			if(allele.length<1){
				bb.append('.').append('\t');
			}else{
				bb.append(allele).append('\t');
			}
		}

//		bb.append(String.format(Locale.ROOT, "%.2f\t", score));
		bb.append(score, 2).append('\t');
		bb.append(pass ? "PASS\t" : "FAIL\t");

		final int scafEndDist=!doNscan ? nScan : (map==null ? start : contigEndDist(map));
		final int count=alleleCount();
		final double af=alleleFraction();
		final double raf=revisedAlleleFraction(af, readLengthAvg);
		final double strandBias=strandBiasScore(scafEndDist);

		{//INFO
			assert(Scaffold.trackStrand()==(minusCoverage>=0)) : Scaffold.trackStrand()+", "+minusCoverage;
			final int covMinus=(Scaffold.trackStrand() ? minusCoverage : coverage/2);
			final int covPlus=Tools.max(0, coverage-covMinus);
			final int refMinus=Tools.max(0, covMinus-alleleMinusCount());
			final int refPlus=Tools.max(0, covPlus-allelePlusCount());
			
			bb.append("SN=").append(scafnum).append(';');
			bb.append("STA=").append(start).append(';');
			bb.append("STO=").append(stop).append(';');
			bb.append("TYP=").append(typeArray[type()]).append(';');
			
			bb.append("R1P=").append(r1plus).append(';');
			bb.append("R1M=").append(r1minus).append(';');
			bb.append("R2P=").append(r2plus).append(';');
			bb.append("R2M=").append(r2minus).append(';');
			
			bb.append("AD=").append(count).append(';');
			bb.append("DP=").append(Tools.max(coverage, count)).append(';');
			bb.append("MCOV=").append(minusCoverage).append(';');
			bb.append("PPC=").append(properPairCount).append(';');
			
			bb.append("AF=").append(af,4).append(';');
			bb.append("RAF=").append(raf,4).append(';');
			bb.append("LS=").append(lengthSum).append(';');
			
			bb.append("MQS=").append(mapQSum).append(';');
			bb.append("MQM=").append(mapQMax).append(';');
			bb.append("BQS=").append(baseQSum).append(';');
			bb.append("BQM=").append(baseQMax).append(';');
			
			bb.append("EDS=").append(endDistSum).append(';');
			bb.append("EDM=").append(endDistMax).append(';');
			bb.append("IDS=").append(idSum).append(';');
			bb.append("IDM=").append(idMax).append(';');
			
			bb.append("CED=").append(scafEndDist).append(';');
			bb.append("HMP=").append(homopolymerCount(map)).append(';');
			bb.append("SB=").append(strandBias,4).append(';');
			
			bb.append("DP4=").append(refPlus).append(',').append(refMinus).append(',');
			bb.append(allelePlusCount()).append(',').append(alleleMinusCount()).append(';');
			
			bb.length--;
		}
		{
			bb.tab();
			bb.append("GT:DP:AD:AF:RAF:SB:SC:PF");
			bb.tab();

			bb.append(genotype(ploidy, pass));
			bb.append(':');
			bb.append(Tools.max(coverage, count));
			bb.append(':');
			bb.append(count);
			bb.append(':');
			bb.append(af,4);
			bb.append(':');

			bb.append(raf,4);
			bb.append(':');
			bb.append(strandBias,4);
			bb.append(':');
			
			bb.append(score,2);
			bb.append(':');
			bb.append(pass ? "PASS" : "FAIL");
		}
		
		return bb;
	}
	
	public int calcCopies(int ploidy){
		final double af=alleleFraction();
//		final int count=count();
		if(ploidy==1){
			return af<0.4 ? 0 : 1;
		}else if(ploidy==2){
			if(af<0.2){return 0;}
			if(af<0.8){return 1;}
			return 2;
		}
		
		int copies=(int)Math.round(ploidy*af);
		if(af>=0.5){copies=Tools.max(copies, 1);}
		return copies;
	}
	
	//TODO: Actually, I should also track the ref coverage.
	private String genotype(int ploidy, boolean pass) {
		if(!pass && noPassDotGenotype){
			if(ploidy==1){return ".";}
			else if(ploidy==2){return "./.";}
			StringBuilder sb=new StringBuilder(ploidy*2-1);
			sb.append('.');
			for(int i=1; i<ploidy; i++){
				sb.append('/').append('.');
			}
			return sb.toString();
		}
		if(ploidy==1){return "1";}
		int copies=calcCopies(ploidy);
		if(ploidy==2){
			if(copies==0){return "0/0";}
			if(copies==1){return "0/1";}
			return "1/1";
		}
		StringBuilder sb=new StringBuilder(ploidy*2);
		int refCopies=ploidy-copies;
		
		for(int i=0; i<refCopies; i++){
			sb.append(0).append('/');
		}
		for(int i=0; i<copies; i++){
			sb.append(1).append('/');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}

	private int hash(){
		return scafnum^Integer.rotateLeft(start, 9)^Integer.rotateRight(stop, 9)^hash(allele);
	}
	
	public static final int hash(byte[] a){
		int code=123456789;
		for(byte b : a){
			code=Integer.rotateLeft(code, 3)^codes[b];
		}
		return code&Integer.MAX_VALUE;
	}
	
	public int calcCoverage(ScafMap map){
		if(coverage>=0){return coverage;}
		
		Scaffold scaf=map.getScaffold(scafnum);
		coverage=scaf.calcCoverage(this);
		if(Scaffold.trackStrand()){minusCoverage=scaf.minusCoverage(this);}
		return coverage;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Scoring Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public static double toPhredScore(double score){
		if(score==0){return 0;}
		score=score*0.998;
		return 2.5*QualityTools.probErrorToPhredDouble(1-score);
	}
	
	public double phredScore(double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
		double score=score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map);
		return toPhredScore(score);
	}
	
	public double score(double properPairRate, double totalQualityAvg, double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map){
		int scafEndDist=(map==null ? start : contigEndDist(map));
		
		double cs=coverageScore(ploidy, rarity, readLengthAvg);
		if(cs==0){return 0;}
		double es=(useEdist ? edistScore() : 1);
		double qs=qualityScore(totalQualityAvg, totalMapqAvg);
		double ps=(usePairing ? pairedScore(properPairRate, scafEndDist) : 1);
		double bs=(useBias ? biasScore(properPairRate, scafEndDist) : 1);
		double is=(useIdentity ? identityScore() : 1);
		double hs=(useHomopolymer ? homopolymerScore(map) : 1);
		return Math.pow(es*qs*ps*bs*cs*is*hs, 0.2);
	}
	
	public double edistScore(){
		double lengthAvg=lengthAvg();
		double edistAvg=((edistAvg()*2+endDistMax))*0.333333333333;
		double constant=5+Tools.min(20, lengthAvg*0.1)+lengthAvg*0.01;
		double weighted=Tools.max(0.05, edistAvg-Tools.min(constant, edistAvg*0.95));
		weighted=weighted*weighted;
		return weighted/(weighted+4);
	}
	
	public double identityScore(){
		double lengthAvg=lengthAvg();
		double idAvg=0.001f*(((identityAvg()+idMax))*0.5f);
		double weighted=Tools.min(1, (idAvg*lengthAvg+(0.65f*Tools.max(1, readlen())))/lengthAvg); //Diminish impact of this var on itself
		weighted=0.75f+0.25f*weighted; //Compress the range
		return weighted;
	}
	
	public double qualityScore(double totalBaseqAvg, double totalMapqAvg){
		return baseQualityScore(totalBaseqAvg)*mapQualityScore(totalMapqAvg);
	}
	
//	public double baseQualityScore(double totalBaseqAvg){
//		double bqAvg=baseQAvg();
//		final double delta=Tools.max(0, totalBaseqAvg-bqAvg);
//	}
	
	public double baseQualityScore(double totalBaseqAvg){
		double bqAvg=baseQAvg();
		
		final double delta=totalBaseqAvg-bqAvg;
		if(delta>0){
			bqAvg=Tools.max(bqAvg*0.5, bqAvg-0.5*delta);
		}
		
		double mult=0.25;
		double thresh=12;
		if(bqAvg>thresh){
			bqAvg=bqAvg-thresh+(thresh*mult);
		}else{
			bqAvg=bqAvg*mult;
		}
		
		double baseProbAvg=1-Math.pow(10, 0-.1*bqAvg);
		double d=baseProbAvg*baseProbAvg;
		return d;
	}
	
//	public double mapQualityScore(double totalMapqAvg){
//		double mqAvg=mapQAvg();
//		double mqAvg2=(0.25f*(3*mqAvg+mapQMax));
//		final double delta=totalMapqAvg-mqAvg;
//		double score=(1-Math.pow(10, 0-.1*(mqAvg2+4)));
//		if(delta>0){
//
//		}else{
//
//		}
//	}

	public double mapQualityScore(double totalMapqAvg){
		double mqAvg=0.5f*(mapQAvg()+mapQMax);
		
		double mapProbAvg=1-Math.pow(10, 0-.1*(mqAvg+2));
		double d=mapProbAvg;
		return d;
	}
	
	public double pairedScore(double properPairRate, int scafEndDist){
		if(properPairRate<0.5){return 0.98;}
		final double count=alleleCount();
		if(count==0){return 0;}
		double rate=properPairCount/count;
		rate=rate*(count/(0.1+count));
		if(rate*1.05>=properPairRate){return Tools.max(rate, 1-0.001*properPairRate);}
		double score=((rate*1.05)/properPairRate)*0.5+0.5;
		score=Tools.max(0.1, score);
		return modifyByEndDist(score, scafEndDist);
	}
	
	public double modifyByEndDist(double x, int scafEndDist){
		if(x>=0.99 || !doNscan || scafEndDist>=nScan){return x;}
		if(scafEndDist<minEndDistForBias){return Tools.max(x, 0.98+0.02*x);}
		double delta=1-x;
		delta=delta*(scafEndDist*scafEndDist)/(nScan*nScan);
		return 1-delta;
	}
	
	public double coverageScore(int ploidy, double rarity, double readLengthAvg){
		int count=alleleCount();
		if(count==0){return 0;}
		double rawScore=count/(lowCoveragePenalty+count); //This may be too severe...
		
//		double ratio=alleleFraction();
		
		double ratio=0.98;
		if(coverage>0){
			double dif=coverage-count;
			if(dif>0){
				dif=dif-coverage*.01f-Tools.min(0.5f, coverage*.1f);
				dif=Tools.max(0.1f, dif);
			}
			ratio=(coverage-dif)/coverage;
			if(type()==SUB && revisedAlleleFraction!=-1 && revisedAlleleFraction<ratio){ratio=revisedAlleleFraction;}
			else{
				ratio=adjustForInsertionLength(ratio, readLengthAvg);
			}
			if(rarity<1 && ratio>rarity){
				double minExpected=1f/ploidy;
				if(ratio<minExpected){
					ratio=minExpected-((minExpected-ratio)*0.1);
				}
			}
		}
		
		double ratio2=Tools.min(1, ploidy*ratio);
		return rawScore*ratio2;
	}
	
	public void reviseAlleleFraction(double readLengthAvg, Scaffold scaffold, VarMap map){
		assert(type()==INS);
		final int ilen=readlen();
		if(ilen<3 || start<1 || start>=scaffold.length-2){return;}
		final byte[] bases=scaffold.bases;
		
		final double afIns=alleleFraction();
		final double rafIns=revisedAlleleFraction(afIns, readLengthAvg);
		final double revisedDif=0.55*(rafIns-afIns); //Half on the left and half on the right, on average
		final double mult=revisedDif/allele.length;
		
		for(int i=0, j=start; i<allele.length && j<scaffold.bases.length; i++, j++){
			final byte b=allele[i];
			if(b!=bases[j]){
				Var key=new Var(scaffold.number, j, j+1, b, SUB);
				Var affectedSub=map.get(key);
//				System.err.println("At pos "+j+": ref="+(char)bases[j]+", alt="+(char)b);
				if(affectedSub!=null){
					assert(key.type()==SUB);
//					System.err.println("Found "+value);
					final double subModifier=revisedDif-mult*i;
					synchronized(affectedSub){
						double afSub=affectedSub.alleleFraction();
						double rafSub=affectedSub.revisedAlleleFraction;
						double modified=afSub-subModifier;
						if(rafSub==-1){
//							System.err.println("sub="+sub+", old="+old);
							affectedSub.revisedAlleleFraction=Tools.max(afSub*0.05, modified);
						}else{
							affectedSub.revisedAlleleFraction=Tools.min(rafSub, Tools.max(afSub*0.05, modified));
						}
					}
				}
			}
		}
		
		for(int i=0, j=start-1; i<allele.length && j>=0; i++, j--){
			final byte b=allele[allele.length-1-i];
			if(b!=bases[j]){
				Var key=new Var(scaffold.number, j, j+1, b, SUB);
				Var affectedSub=map.get(key);
				if(affectedSub!=null){
					assert(key.type()==SUB);
//					System.err.println("Found "+value);
					final double subModifier=revisedDif-mult*i;
					synchronized(affectedSub){
						double afSub=affectedSub.alleleFraction();
						double rafSub=affectedSub.revisedAlleleFraction;
						double modified=afSub-subModifier;
						if(rafSub==-1){
//							System.err.println("sub="+sub+", old="+old);
							affectedSub.revisedAlleleFraction=Tools.max(afSub*0.05, modified);
						}else{
							affectedSub.revisedAlleleFraction=Tools.min(rafSub, Tools.max(afSub*0.05, modified));
						}
					}
				}
			}
		}
	}
	
	public double revisedAlleleFraction(double af, double readLengthAvg){
		if(revisedAlleleFraction!=-1){
			return revisedAlleleFraction;
		}else if(type()==INS){
			return revisedAlleleFraction=adjustForInsertionLength(af, readLengthAvg);
		}
		return af;
	}
	
	public double adjustForInsertionLength(final double ratio, final double rlen0){
		if(type()!=INS){return ratio;}
		final int ilen=readlen();
		if(ilen<2){return ratio;}
		
		final double rlen=Tools.max(ilen*1.2+6, rlen0);
		final double sites=rlen+ilen-1;
		final double goodSites=rlen-ilen*1.1-6;
		
		final double expectedFraction=goodSites/sites;
		final double revisedRatio=Tools.min(ratio/expectedFraction, 1-(1-ratio)*0.1);
		return revisedRatio;
	}
	
	public double homopolymerScore(ScafMap map){
		if(map==null){return 1;}
		
		int count=homopolymerCount(map);
//		assert(false) : count;
		if(count<2){return 1;}
		return 1f-(count*0.1f/9);
	}
	
	public int homopolymerCount(ScafMap map){
		if(map==null){return 0;}
		final byte[] bases=map.getScaffold(scafnum).bases;
		if(bases==null){return 0;}
		
		final int type=type();
		if(type==SUB){
			assert(start==stop-1) : start+", "+stop;
			final byte base=allele[0];
			int x=homopolymerCountSub(bases, start, base);
//			assert(false) : (char)base+", "+x;
			return x;
		}else if(type==INS){
			final byte base1=allele[0], base2=allele[allele.length-1];
			int i=0;
			while(i<allele.length && allele[i]==base1){i++;}
			while(i<allele.length && allele[i]==base2){i++;}
			if(i<bases.length){return 0;}
			int left=homopolymerCountLeft(bases, start, base1);
			int right=homopolymerCountRight(bases, stop+1, base2);
//			assert(false) : "INS "+(left+right+1)+" "+new String(allele)+" "+new String(bases, start-4, 8);
			return left+right+1;
		}else if(type==DEL){
			if(start<0 || start+1>=bases.length || stop<=0 || stop>=bases.length){return 0;}
			final byte base1=bases[start+1], base2=bases[stop-1];
			int pos=start+1;
			while(pos<=stop && bases[pos]==base1){pos++;}
			while(pos<=stop && bases[pos]==base2){pos++;}
			if(pos<=stop){return 0;}
			int left=homopolymerCountLeft(bases, start, base1);
			int right=homopolymerCountRight(bases, stop, base2);
//			assert(false || reflen()>10) : "DEL "+(left+right+1)+" "+new String(allele)+" "+new String(bases, start-4, stop-start+8);
			return left+right+1;
		}else{
//			assert(false) : type();
			return 0;
		}
	}
	
	public static int homopolymerCountSub(final byte[] bases, final int pos, final byte base){
//		System.err.println("A:"+pos+", "+bases.length+", "+(char)base);
		if(pos<0 || pos>=bases.length){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
//		System.err.println("B:"+pos+", "+bases.length);
		
		int count1=0;
		for(int i=pos-1, lim=Tools.max(0, pos-4); i>=lim; i--){
			if(bases[i]==base){count1++;}
			else{break;}
		}
//		System.err.println("C:"+pos+", "+bases.length+", "+count1);
		
		int count2=0;
		for(int i=pos+1, lim=Tools.min(bases.length, pos+5); i<lim; i++){
			if(bases[i]==base){count2++;}
			else{break;}
		}
//		System.err.println("D:"+pos+", "+bases.length+", "+count2);
//		System.err.println("E:"+new String(bases, pos-4, 9));
		assert(count1+count2<=8) : count1+", "+count2;
		
		return count1+count2+(count1>0 && count2>0 ? 1 : 0);
	}
	
	public static int homopolymerCountLeft(final byte[] bases, final int pos, final byte base){
		if(pos<0 || bases[pos]!=base){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
		
		int count=0;
		for(int i=pos, lim=Tools.max(0, pos-3); i>=lim; i--){
			if(bases[i]==base){count++;}
			else{break;}
		}
		return count;
	}
	
	public static int homopolymerCountRight(final byte[] bases, final int pos, final byte base){
		if(pos<0 || bases[pos]!=base){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
		
		int count=0;
		for(int i=pos, lim=Tools.min(bases.length, pos+4); i<lim; i++){
			if(bases[i]==base){count++;}
			else{break;}
		}
		return count;
	}
	
	public double biasScore(double properPairRate, int scafEndDist){
		double strandBias=strandBiasScore(scafEndDist);
		double readBias=readBiasScore(properPairRate);
		return Math.sqrt(strandBias*readBias);
	}
	
	public double strandBiasScore(int scafEndDist){
		double x=eventProb(allelePlusCount(), alleleMinusCount());
		return modifyByEndDist(x, scafEndDist);
	}
	
	public double readBiasScore(double properPairRate){
		if(properPairRate<0.5){return 0.95f;}
		
		return eventProb(r1AlleleCount(), r2AlleleCount());
	}
	
	/** Adjusted probability of a binomial event being at least this lopsided. */
	public static double eventProb(int a, int b){
		
		double allowedBias=0.75;
		double slopMult=0.95;
		
		double n=a+b;
		double k=Tools.min(a, b);
		
		double slop=n*(allowedBias*0.5);
		double dif=n-k*2;
		dif=dif-(Tools.min(slop, dif)*slopMult);
		n=k*2+dif;
//		k=n*0.5-dif;
		assert(k<=n*0.5) : a+", "+b+", "+n+", "+k+", "+slop+", "+dif;
		
		if(n>PROBLEN){
			double mult=PROBLEN/n;
			n=PROBLEN;
			k=(int)(k*mult);
		}

		int n2=(int)Math.round(n);
		int k2=Tools.min(n2/2, (int)(k+1));
		

//		if(a+b>3){
//			System.err.println(n+", "+k+", "+n2+", "+k2);
//		}
		
		double result=prob[n2][k2];
		if(result<1 || a==b || a+1==b || a==b+1){return result;}
		
		double slope=Tools.min(a, b)/(double)Tools.max(a, b);
		return (0.998+slope*0.002);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int allelePlusCount(){return r1plus+r2plus;}
	public int alleleMinusCount(){return r1minus+r2minus;}
	public int r1AlleleCount(){return r1plus+r1minus;}
	public int r2AlleleCount(){return r2plus+r2minus;}
	public int alleleCount(){return r1plus+r1minus+r2plus+r2minus;}

	public double alleleFraction(){
		int count=alleleCount();
		int cov=Tools.max(count, coverage, 1);
		return count/(double)cov;
	}
	
	public double strandRatio(){
		int plus=allelePlusCount();
		int minus=alleleMinusCount();
		if(plus==minus){return 1;}
		return (Tools.min(plus,  minus)+1)/(double)Tools.max(plus, minus);
	}
	public double baseQAvg(){return baseQSum/(double)alleleCount();}
	public double mapQAvg(){return mapQSum/(double)alleleCount();}
	public double edistAvg(){return endDistSum/(double)alleleCount();}
	public double identityAvg(){return idSum/(double)alleleCount();}
	public double lengthAvg(){return lengthSum/(double)alleleCount();}
	public double properPairRate(){return properPairCount/(double)alleleCount();}
	
	
	public void setCoverage(int coverage_, int minusCoverage_){
		coverage=coverage_;
		minusCoverage=minusCoverage_;
	}
	
	public int coverage(){
		assert(coverage>-1) : coverage+", "+this;
		return coverage;
	}
	
	public boolean hasCoverage(){
		return coverage>-1;
	}
	
	public int contigEndDist(ScafMap map){
		Scaffold scaf=map.getScaffold(scafnum);
		int len=scaf.length;
		byte[] bases=scaf.bases;
		
		int scafEndDist=Tools.max(0, Tools.min(start, len-stop));
		if(bases==null || nScan<1){return scafEndDist;}
		int limit=Tools.min(nScan, scafEndDist);
		int contigEndDist=leftContigEndDist(bases, limit);
		limit=Tools.min(limit, contigEndDist);
		contigEndDist=rightContigEndDist(bases, limit);
		return Tools.min(scafEndDist, contigEndDist);
	}
	
	public int leftContigEndDist(byte[] bases, int maxDist){
		if(start>=bases.length){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=start, lim=Tools.max(0, start-maxDist); i>=lim; i--){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=start-i-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}
	
	public int rightContigEndDist(byte[] bases, int maxDist){
		if(stop<0){return Tools.min(bases.length, maxDist+1);}
		int ns=0;
		for(int i=stop, lim=Tools.min(bases.length-1, stop+maxDist); i<=lim; i++){
			if(AminoAcid.isFullyDefined(bases[i])){
				ns=0;
			}else{
				ns++;
				if(ns>=10){
					int x=i-stop-ns+1;
					assert(x>=0);
					return x;
				}
			}
		}
		return maxDist+1;
	}
	
	public String scafName(){
		return scafName(ScafMap.defaultScafMap());
	}
	
	public String scafName(ScafMap map){
		return map.getScaffold(scafnum).name;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int scafnum;
	public final int start;
	public final int stop; //Half-open, so stop is always after start except for insertions
	public final byte[] allele;
	public final int hashcode;
	public final int type;
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/

	private int coverage=-1;
	private int minusCoverage=-1;
	
	int r1plus;
	int r1minus;
	int r2plus;
	int r2minus;
	int properPairCount;
	
	long mapQSum;
	public int mapQMax;
	
	long baseQSum;
	public int baseQMax;
	
	long endDistSum;
	public int endDistMax;
	
	long idSum;
	int idMax;
	
	long lengthSum;
	
	double revisedAlleleFraction=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean CALL_INS=true;
	public static boolean CALL_DEL=true;
	public static boolean CALL_SUB=true;
	public static boolean CALL_NOCALL=false;
	public static boolean CALL_JUNCTION=false;
	
	public static boolean extendedText=true;
	
	public static boolean noPassDotGenotype=false;
	
	public static boolean useHomopolymer=true;
	public static boolean useIdentity=true;
	public static boolean usePairing=true;
	public static boolean useBias=true;
	public static boolean useEdist=true;
	public static boolean doNscan=true;
	public static double lowCoveragePenalty=0.8;
	public static int nScan=600;
	public static int minEndDistForBias=200;
	
	final static int PROBLEN=100;
	public static final boolean UPPER_CASE_ALLELES=true;
	/** Verify that there are no variants with alt equal to ref */
	private static final boolean TEST_REF_VARIANTS=false;
	
//	static final String[] typeArray=new String[] {"NOCALL","SUB","DEL","INS"};
//	static final int NOCALL=0, SUB=1, DEL=2, INS=3;

	private static final byte colon=';';
	private static final byte tab='\t';
	public static final String[] typeArray=new String[] {"INS","NOCALL","SUB","DEL","LJUNCT","RJUNCT","BJUNCT"};
	public static final int INS=0, NOCALL=1, SUB=2, DEL=3, LJUNCT=4, RJUNCT=5, BJUNCT=6;
	public static final int VAR_TYPES=7;
	static final byte[] typeInitialArray=new byte[128];
	
	static final byte[] AL_0=new byte[0];
	static final byte[] AL_A=new byte[] {(byte)'A'};
	static final byte[] AL_C=new byte[] {(byte)'C'};
	static final byte[] AL_G=new byte[] {(byte)'G'};
	static final byte[] AL_T=new byte[] {(byte)'T'};
	static final byte[] AL_N=new byte[] {(byte)'N'};
	static final byte[][] AL_MAP=makeMap();
	static final int[] codes=makeCodes();

//	public static ScafMap scafMap;
	
	static final int[] makeCodes(){
		Random randy=new Random(1);
		int[] array=new int[256];
		for(int i=0; i<array.length; i++){
			array[i]=randy.nextInt();
		}
		return array;
	}
	
	static final byte[][] makeMap(){
		byte[][] map=new byte[128][];
		map[0]=map['.']=map['\t']=AL_0;
		map['A']=map['a']=AL_A;
		map['C']=map['c']=AL_C;
		map['G']=map['g']=AL_G;
		map['T']=map['t']=AL_T;
		map['N']=map['n']=AL_N;
		for(int i=0; i<map.length; i++){
			if(map[i]==null){map[i]=new byte[(byte)i];}
		}
		return map;
	}
	
	/** factorial[n]=n! */
	private static final double[] factorial=makeFactorialArray(PROBLEN+1);
	/** binomial[n][k] = combinations in n pick k */
	private static final double[][] binomial=makeBinomialMatrix(PROBLEN+1);
	/** prob[n][k] = probability of an event this lopsided or worse. */
	private static final double[][] prob=makeProbMatrix(PROBLEN+1);

	private static double[] makeFactorialArray(int len) {
		double[] x=new double[len];
		x[0]=1;
		for(int i=1; i<len; i++){
			x[i]=x[i-1]*i;
		}
		return x;
	}

	private static double[][] makeBinomialMatrix(int len) {
		double[][] matrix=new double[len][];
		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double nf=factorial[n];
			matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double kf=factorial[k];
				final double nmkf=factorial[n-k];
				double combinations=nf/kf;
				combinations=combinations/nmkf;
				matrix[n][k]=combinations;
			}
		}
		return matrix;
	}

	private static double[][] makeProbMatrix(int len) {
		double[][] matrix=new double[len][];
		double mult=2;
		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double[] array=matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double combinations=binomial[n][k];
				array[k]=combinations*mult;
			}
//			if(n<=12){System.err.println(Arrays.toString(array));}
			for(int k=0; k<=kmax; k++){
				array[k]=Tools.min(1, (k==0 ? 0 : array[k-1])+array[k]);
			}
//			if(n<=12){System.err.println(Arrays.toString(array));}
//			assert(array[kmax]==1) : Arrays.toString(array);
			mult*=0.5;
//			if(n<=12){System.err.println();}
		}
		return matrix;
	}
	
	static {
		Arrays.fill(typeInitialArray, (byte)-1);
		typeInitialArray['I']=INS;
		typeInitialArray['N']=NOCALL;
		typeInitialArray['S']=SUB;
		typeInitialArray['D']=DEL;
		typeInitialArray['L']=LJUNCT;
		typeInitialArray['R']=RJUNCT;
		typeInitialArray['B']=BJUNCT;
	}
	
	public static final String varFormat="1.1";
	
}
