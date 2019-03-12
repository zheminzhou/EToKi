package sketch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Tools;
import structures.AbstractBitSet;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongHashMap;
import structures.LongList;
import structures.LongPair;
import tax.ImgRecord;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch extends SketchObject implements Comparable<Sketch> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	//Array should already be hashed, sorted, unique, subtracted from Long.MAX_VALUE, then reversed.
	public Sketch(long[] array_, int[] counts_, ArrayList<String> meta_){
		this(array_, counts_, -1, -1, -1, -1, -1, -1, null, null, null, meta_);
	}
	
	public Sketch(SketchHeap heap, boolean clearFname, boolean keepCounts, ArrayList<String> meta_){
		this(heap.toSketchArray(), null, (int)heap.taxID, heap.imgID, heap.genomeSizeBases, heap.genomeSizeKmers, heap.genomeSequences,
				heap.probCorrect(), heap.taxName(), heap.name0(), heap.fname(), meta_);
		assert(counts==null);
		if(!heap.setMode && keepCounts){
			LongHashMap map=heap.map();
			counts=new int[array.length];
			for(int i=0; i<array.length; i++){
				int count=map.get(Long.MAX_VALUE-array[i]);
				assert(count>0) : array[i]+" -> "+count+"\n"+Arrays.toString(map.values())+"\n"+Arrays.toString(map.keys());
				counts[i]=count;
			}
		}
		if(heap.setMode){heap.clearSet();}
		heap.clear(clearFname);
//		System.err.println("size="+size+", genome="+this.genomeSize+", m"); : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (counts==null)+", "+heap.setMode+", "+keepCounts;
	}

	public Sketch(long[] array_, int[] counts_, int taxID_, long imgID_, long gSizeBases_, long gSizeKmers_, long gSequences_, double probCorrect_,
			String taxName_, String name0_, String fname_, ArrayList<String> meta_){
		array=array_;
		counts=counts_;
		assert(counts==null || array==null || counts.length==array.length) : (array==null ? "null" : array.length)+", "+(counts==null ? "null" : counts.length);
		taxID=taxID_;
		imgID=imgID_;
		sketchID=nextSketch.getAndIncrement();
		genomeSizeBases=gSizeBases_;
		genomeSizeKmers=gSizeKmers_;
		genomeSequences=gSequences_;
		probCorrect=probCorrect_<=0 ? 0f : (float)probCorrect_;
		meta=meta_;
		
		taxName=fix(taxName_);
		name0=fix(name0_);
		fname=fix(fname_);
		
//		if(k2>0){
//			if(useToValue2) {
//				int count=0;
//				for(long x : array){
//					count+=(int)(x&1);
//				}
//				k1Count=count;
//			}else {
//				k1Count=array.length/2;
//			}
//		}else{
//			k1Count=array.length;
//		}
		
		if(ImgRecord.imgMap!=null && imgID>=0 && taxID<0){
			ImgRecord record=ImgRecord.imgMap.get(imgID);
			if(record!=null){
				if(record.name!=null && taxName==null){taxName=record.name;}
				taxID=record.taxID;
			}
		}
	}
	
	void addMeta(String s){
		s=fixMeta(s);
		if(s==null){return;}
		if(meta==null){meta=new ArrayList<String>(1);}
		meta.add(s);
	}
	
	void setMeta(ArrayList<String> list){
		assert(meta==null);
		meta=list;
	}
	
	private static String fix(String s){
		if(s==null){return null;}
		return s.replace('\t', ' ');
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(Sketch other, int maxlen){
		final long[] a=array;
		final long[] b=other.array;
		if(maxlen<1){
			assert(false);
			maxlen=1000000;
		}
		LongList list=new LongList(Tools.min(maxlen, a.length+b.length));
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//match
				list.add(ka);
				i++;
				j++;
			}else if(ka<kb){
				list.add(ka);
				i++;
			}else{
				list.add(kb);
				j++;
			}
			if(list.size()>=maxlen){break;}
		}
		
		if(array.length==list.size()){
			for(int i=0; i<list.size; i++){
				array[i]=list.array[i];
			}
		}else{
			array=list.toArray();
		}
	}

	public void resize(int newSize) {
		if(newSize>=length()){return;}
		array=Arrays.copyOf(array, newSize);
		if(counts!=null){counts=Arrays.copyOf(counts, newSize);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Set Operations        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final Sketch intersection(Sketch sa, Sketch sb){
		Sketch shared=intersection(sa.array, sb.array, sa.counts);
		if(shared!=null){
			shared.taxID=sb.taxID;
			shared.taxName=sb.taxName;
			shared.name0=sb.name0;
			shared.fname=sb.fname;
			shared.meta=sb.meta;
			shared.imgID=sb.imgID;
			shared.spid=sb.spid;
		}
		return shared;
	}
	
	public static final Sketch intersection(long[] a, long[] b, int[] aCounts){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=new IntList();
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(aCounts!=null){
					il.add(aCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				i++;
			}else{
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null, null);
	}
	
	public static final Sketch union(Sketch sa, Sketch sb){
		Sketch shared=union(sa.array, sb.array, sa.counts, sb.counts);
		if(shared!=null){
			shared.taxID=sa.taxID;
			shared.taxName=sa.taxName;
			shared.name0=sa.name0;
			shared.fname=sa.fname;
			shared.meta=sa.meta;
			shared.imgID=sa.imgID;
			shared.spid=sa.spid;
		}
		return shared;
	}
	
	public static final Sketch union(long[] a, long[] b, int[] aCounts, int[] bCounts){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=new IntList();
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(aCounts!=null && bCounts!=null){
					il.add(aCounts[i]+bCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				ll.add(ka);
				if(aCounts!=null && bCounts!=null){
					il.add(aCounts[i]);
				}
				i++;
			}else{
				ll.add(kb);
				if(aCounts!=null && bCounts!=null){
					il.add(bCounts[i]);
				}
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null, null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Filtering          ----------------*/
	/*--------------------------------------------------------------*/

	boolean passesMeta(DisplayParams params) {
		return passesMeta(params.requiredMeta, params.bannedMeta, /*params.requiredTaxid, params.bannedTaxid,*/ params.requiredMetaAnd);
	}

	boolean passesMeta(ArrayList<String> requiredMeta, ArrayList<String> bannedMeta, /*IntList requiredTaxid, IntList bannedTaxid,*/ boolean requiredMetaAnd) {
		assert(requiredMeta!=null || bannedMeta!=null /*|| requiredTaxid!=null || bannedTaxid!=null*/);
		assert(requiredMeta==null || requiredMeta.size()>0);
		assert(bannedMeta==null || bannedMeta.size()>0);

//		if(requiredTaxid!=null && !requiredTaxid.contains(taxID)){return false;}
//		if(bannedTaxid!=null && bannedTaxid.contains(taxID)){return false;}
		
		if(requiredMeta!=null){
			if(meta==null){return false;}
			for(String tag : requiredMeta){
				if(meta.contains(tag)){
					if(!requiredMetaAnd){break;}
				}else if(requiredMetaAnd){
					return false;
				}
			}
		}
		if(bannedMeta!=null && meta!=null){
			for(String tag : bannedMeta){
				if(!meta.contains(tag)){
					return false;
				}
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/
	
	public int countMatches(Sketch other, CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){
		return countMatches(array, other.array, counts, other.counts, refHitCounts, other.taxID, buffer, present, fillPresent, taxHits, contamLevel);
	}
	
	public static final int countMatches(long[] a, long[] b, int[] aCounts, int[] bCounts, int[] refHitCounts, int bid,
			CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){

		
//		if(verbose2){System.err.println("fillPresent: "+fillPresent+", "+present);}
//		assert(fillPresent) : bid+", "+minFakeID+", "+(taxHits!=null);
		
		if(bid>0 && bid<minFakeID && taxHits!=null){
			bid=taxtree.getIdAtLevelExtended(bid, contamLevel);
		}else{
			bid=-1;
		}
		
//		assert(false) : (buffer==null)+", "+fillPresent+", "+present.cardinality();
		assert(a.length>0 && b.length>0);
		
		//Kmers hitting this reference
		int matches=0;
		
		//Kmers hitting this reference and others
		int multiMatches=0;
		
		//Kmers hitting nothing
		int noHits=0;
		
		//Kmers hitting some organism but not this reference
		int contamHits=0;
		
		//Kmers hitting something in this taxa, but not this reference
		int sameTax=0;
		
		//Kmers hitting this organism and no other taxa
		int unique2=0;
		
		//Kmers hitting only this taxa but not this organism (this count may not include everything due to early exit)
		int unique3_temp=0;
		
		//Kmers hitting multiple organisms but not this reference
		int multiContamHits=0;
		
		//Sum of query counts for shared kmers
		long depthSum=0;
		
		//Sum of query counts for shared kmers divided by ref counts for those kmers
		double depthSum2=0;//Slow, but necessary.
		
		//Number of times matching query keys occurred in all references
		long refHitSum=0;

		//Matches from k1
		int k1hits=0;
		
		//Query kmers from k1
		int k1seenQ=0;
		
		//Reference kmers from k1
		int k1seenR=0;
		int i=0, j=0;
		assert(present==null || present.capacity()==a.length);
//		assert(false) : buffer.rbs.capacity()+", "+buffer.rbs+", "+present;
		if(present!=null){
			if(fillPresent){
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					final int bit=(int)(ka&1);
					if(ka==kb){
						present.increment(i);
						matches++;
						k1hits+=bit;
						k1seenQ+=bit;
						k1seenR+=bit;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
						i++;
						j++;
					}else if(ka<kb){
						i++;
						k1seenQ+=bit;
					}else{
						j++;
						k1seenR+=bit;
					}
				}
			}else{
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					final int bit=(int)(ka&1);
					if(ka==kb){
						final int count=present.getCount(i);
						if(count>1){
							multiMatches++;
						}
						
						matches++;
						k1hits+=bit;
						k1seenQ+=bit;
						k1seenR+=bit;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique2++;}
						}
						
						i++;
						j++;
					}else if(ka<kb){
						k1seenQ+=bit;
						final int count=present.getCount(i);
						if(count>0){
							contamHits++;
							if(count>1){
								multiContamHits++;
							}
						}else{
							noHits++;
						}
						
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null){
								if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
								for(int tid : taxHitsRow){
									if(tid==bid){
										sameTax++;
										break;
									}
								}
							}
						}
						
						i++;
					}else{
						k1seenR+=bit;
						j++;
					}
				}
				
				//For the remaining query kmers, we don't know whether the reference sketch would have shared them had it been longer.
				//This section can be disabled to prevent them from being displayed.
				if(bid>0 && i<a.length-1){
					for(; i<a.length; i++){
						int[] taxHitsRow=taxHits[i];
						if(taxHitsRow!=null){
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
						}
					}
				}
			}
		}else{
			for(; i<a.length && j<b.length; ){
				final long ka=a[i], kb=b[j];
				final int bit=(int)(ka&1);
				final int count=present.getCount(i);
				if(ka==kb){
					matches++;
					k1hits+=bit;
					k1seenQ+=bit;
					k1seenR+=bit;
					if(aCounts!=null){depthSum+=aCounts[i];}
					if(refHitCounts!=null){refHitSum+=refHitCounts[i];}
					i++;
					j++;
				}else if(ka<kb){
					i++;
					k1seenQ+=bit;
				}else{
					j++;
					k1seenR+=bit;
				}
			}
		}
		
		if(k2<1){
			k1hits=matches;
			k1seenQ=i;
			k1seenR=j;
		}
		
//		if(taxHits!=null){
//			System.err.println("matches="+matches+", noHits="+noHits+", contamHits="+contamHits+", sameTax="+sameTax+", multiContamHits="+multiContamHits);
//		}

//		assert(bid<1 || unique2>=(matches-multiMatches)) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches;
//		assert(matches<1000 || multiMatches==0) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches+", "+fillPresent;
		
		if(buffer!=null){
//			System.err.println("*A) "+matches+", "+multiMatches+", "+unique2+", "+unique3_temp);
//			new Exception().printStackTrace();
//			assert(k1hits<=(matches-k1hits)) : k1hits+", "+matches;
			buffer.set(matches, multiMatches, unique2, unique2+unique3_temp, noHits, 
					contamHits, contamHits-sameTax, multiContamHits, i, j, 
					a.length, b.length, depthSum, depthSum2, refHitSum, k1hits, k1seenQ, k1seenR);
		}
		return matches;
	}
	
//	public float identity(Sketch b, float[] ret){
//		if(ret!=null){Arrays.fill(ret, 0);}
//		return identityWeighted(array, b.array, ret);
//	}
//
//	public static float identity(long[] a, long[] b){
//		int matches=countMatches(a, b);
//		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
//	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
//		System.err.println(code+", "+gSize+", "+taxID+", "+imgID+", "+name0);
		return code;
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		int x=taxName.compareTo(b.taxName);
		if(x!=0){return x;}
		if(name0!=null && b.name0!=null){return name0.compareTo(b.name0);}
		return name0!=null ? 1 : b.name0!=null ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/
	
	public ByteBuilder toHeader(){
		ByteBuilder sb=new ByteBuilder();
		return toHeader(sb);
	}
	
	public ByteBuilder toHeader(ByteBuilder bb){
		bb.append("#SZ:").append(array.length);
		bb.append("\tCD:");
		bb.append(codingArray[CODING]);
		if(deltaOut){bb.append('D');}
		if(counts!=null){bb.append('C');}
		if(aminoOrTranslate()){bb.append('M');}
		if(amino8){bb.append('8');}

		bb.append("\tK:").append(k);
		if(k2!=0){bb.append(",").append(k2);}
		if(HASH_VERSION>1){bb.append("\tH:").append(HASH_VERSION);}
		
		if(genomeSizeBases>0){bb.append("\tGS:").append(genomeSizeBases);}
		if(genomeSizeKmers>0){bb.append("\tGK:").append(genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){bb.append("\tGE:").append(ge);}
		if(genomeSequences>0){bb.append("\tGQ:"+genomeSequences);}
		if(probCorrect>0){bb.append("\tPC:"+String.format("%.4f", probCorrect));}
		if(taxID>=0){bb.append("\tID:").append(taxID);}
		if(imgID>=0){bb.append("\tIMG:").append(imgID);}
		if(spid>0){bb.append("\tSPID:").append(spid);}
		if(fname!=null){bb.append("\tFN:").append(fname);}
		if(taxName!=null){bb.append("\tNM:").append(taxName);}
		if(name0!=null){bb.append("\tNM0:").append(name0);}
		if(meta!=null){
			for(String s : meta){
				bb.append("\tMT_").append(s);
			}
		}
		return bb;
	}
	
	public ByteBuilder toBytes(){
		return toBytes(new ByteBuilder());
	}
	
	public ByteBuilder toBytes(ByteBuilder bb){
		if(CODING==A48 && deltaOut){return toBytesA48D(bb);}
		long prev=0;
		toHeader(bb);
		bb.append("\n");
		byte[] temp=null;
		if(CODING==A48){temp=new byte[12];}
		for(int i=0; i<array.length; i++){
			long key=array[i];
			int count=(counts==null ? 1 : counts[i]);
			long x=key-prev;
			if(CODING==A48){
				appendA48(x, bb, temp);
				if(count>1){
					bb.append('\t');
					appendA48(count-1, bb, temp);
				}
				bb.append('\n');
			}else if(CODING==HEX){
				bb.append(Long.toHexString(x)).append('\n');
			}else if(CODING==RAW){
				bb.append(x).append('\n');
			}else{
				assert(false);
			}
			if(deltaOut){prev=key;}
		}
		return bb;
	}
	
	//This is to make the common case fast
	private ByteBuilder toBytesA48D(ByteBuilder bb){
		assert(CODING==A48 && deltaOut);
		long prev=0;
		toHeader(bb);
		bb.append("\n");
		final byte[] temp=new byte[12];

		if(counts==null){
			for(int i=0; i<array.length; i++){
				long key=array[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, bb, temp);
					bb.append('\n');
				}
				prev=key;
			}
		}else{
			for(int i=0; i<array.length; i++){
				long key=array[i];
				int count=counts[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, bb, temp);
					if(count>1){
						bb.append('\t');
						appendA48(count-1, bb, temp);
					}
					bb.append('\n');
				}
				prev=key;
			}
		}
		return bb;
	}
	
	public static final void appendA48(long value, ByteBuilder bb, byte[] temp){
		int i=0;
//		long value=value0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			bb.append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				bb.append((char)(temp[i]+48));
			}
		}
	}
	
	public static final String toA48(long value){
		int i=0;
//		long value=value0;
		StringBuilder sb=new StringBuilder(12);
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			sb.append((char)(b+48));
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			sb.reverse();
		}
		return sb.toString();
	}
	
	@Override
	public String toString(){
		return toBytes().toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       String Parsing         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static long parseA48(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=6;
			long c=line.charAt(i);
			x|=(c-48);
		}
		return x;
	}
	
	/** Parses coverage too */
	public static long parseA48C(String line, IntList covList){
		if(line.length()==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length();
		for(; i<len; i++){
			long c=line.charAt(i);
			if(c<48){break;}
			key<<=6;
			key|=(c-48);
		}
		for(i++; i<len; i++){
			long c=line.charAt(i);
			cov<<=6;
			cov|=(c-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	public static long parseHex(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=4;
			x|=hexTable[line.charAt(i)];
		}
		if(line.charAt(0)=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Array Parsing         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static long parseA48(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	public static long parseNuc(String line){
		return parseNuc(line.getBytes());
	}
	
	/** Returns the maximal key in the sequence */
	public static long parseNuc(byte[] bases){
		if(bases.length<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		long kmer=0, rkmer=0;
		int len=0;
		
		long key=-1;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
//				long z=Tools.max(kmer, rkmer);
				final long hashcode=hash(kmer, rkmer);
				key=Tools.max(key, hashcode);
			}
		}
		return key<minHashValue ? -1 : Long.MAX_VALUE-key;
	}
	
	/** Parses coverage too */
	public static long parseA48C(byte[] line, IntList covList){
		if(line.length==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length;
		for(; i<len; i++){
			long b=line[i];
			if(b<48){break;}
			key<<=6;
			key|=(b-48);
		}
		for(i++; i<len; i++){
			long b=line[i];
			cov<<=6;
			cov|=(b-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static long parseA48(ByteBuilder bb){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(len==0){return 0;}
		long x=0;
		for(int i=0; i<len; i++){
			x<<=6;
			x|=(((long)line[i])-48);
		}
		return x;
	}
	
	public static long parseNuc(ByteBuilder bb){
		return parseNuc(bb.toBytes());
	}
	
	/** Parses coverage too */
	public static long parseA48C(ByteBuilder bb, IntList covList){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(len==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0;
		for(; i<len; i++){
			long b=line[i];
			if(b<48){break;}
			key<<=6;
			key|=(b-48);
		}
		for(i++; i<len; i++){
			long b=line[i];
			cov<<=6;
			cov|=(b-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	public static long parseHex(ByteBuilder bb){
		final int len=bb.length;
		final byte[] line=bb.array;
		if(line.length==0){return 0;}
		long x=0;
		for(int i=0; i<len; i++){
			x<<=4;
			x|=hexTable[line[i]];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long genomeSizeEstimate() {
		return array.length==0 ? 0 : Tools.min(genomeSizeKmers, genomeSizeEstimate(array[array.length-1], array.length));
	}
	
	public long genomeSizeEstimate(int minCount) {
		if(minCount<2){return genomeSizeEstimate();}
		if(length()==0){return 0;}
		long max=0;
		int num=0;
		for(int i=0; i<counts.length; i++){
			if(counts[i]>=minCount){
				max=array[i];
				num++;
			}
		}
		if(max==0){return 0;}
		long est=Tools.min(genomeSizeKmers, SketchObject.genomeSizeEstimate(max, num));
		return est;
	}
	
	public String filePrefix(){return ReadWrite.stripToCore(fname);}
	public String name(){return taxName!=null ? taxName : name0!=null ? name0 : fname;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public String fname(){return fname;}
	public int length(){return array.length;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	public void setFname(String s){
//		assert(!s.endsWith("sketch")) : s; //123
		fname=s;
	}
	
//	public float k1Fraction(){return k1Count/(float)array.length;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Assorted           ----------------*/
	/*--------------------------------------------------------------*/

	public ArrayList<LongPair> toKhist() {
		HashMap<Long, LongPair> map=new HashMap<Long, LongPair>();
		for(int i=0; i<counts.length; i++){
			int a=counts[i];
			Long key=(long)a;
			LongPair value=map.get(key);
			if(value==null){
				value=new LongPair();
				value.a=a;
				map.put(key, value);
			}
			value.b++;
		}
		ArrayList<LongPair> list=new ArrayList<LongPair>(map.size());
		list.addAll(map.values());
		Collections.sort(list);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            BitSet            ----------------*/
	/*--------------------------------------------------------------*/

	public void makeBitSets(boolean printContam, boolean index){
		assert(compareBitSet==null && indexBitSet==null);
		if(!printContam){return;}
		compareBitSet=AbstractBitSet.make(length(), bitSetBits);
		if(index){indexBitSet=AbstractBitSet.make(length(), bitSetBits);}
	}
	
	public void addToBitSet(AbstractBitSet rbs){
		compareBitSet.add(rbs);
	}
	
	public AbstractBitSet compareBitSet(){return compareBitSet;}
	
	
	public AbstractBitSet indexBitSet(){return indexBitSet;}
	
	public void mergeBitSets(){
		assert(!mergedBitSets);
		if(compareBitSet!=null && indexBitSet!=null){
			compareBitSet.setToMax(indexBitSet);
		}
		indexBitSet=null;
		mergedBitSets=true;
	}
	
	public boolean merged(){return mergedBitSets;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long[] array;
	int[] counts;
	public int taxID;
	int sketchID;//Based on loading order
	public final long genomeSequences;
	public final long genomeSizeBases;
	public final long genomeSizeKmers;
	public final float probCorrect;
//	public final int k1Count; //Number of keys made from k1 rather than k2
	private String taxName;
	private String name0;
	private String fname;
	ArrayList<String> meta;
	
	//TODO: These should move to SketchResults.
	private AbstractBitSet compareBitSet; //Used for comparison
	private AbstractBitSet indexBitSet;
	
	
	//Extended information
	public long imgID=-1;
	public long spid=-1;
//	public String seqUnitName=null;
	
	private boolean mergedBitSets=false; //TODO: Temporary for debugging
	/** Tracks the number of reference sketches sharing each kmer.
	 * Should be set to null when no longer needed. */
	private int[] refHitCounts;

	public int[] refHitCounts(){return refHitCounts;}
	public void clearRefHitCounts(){
		refHitCounts=null;
	}
	public void setRefHitCounts(int[] x) {
		refHitCounts=x;
		assert(x!=null);
	}
	
	private static AtomicInteger nextSketch=new AtomicInteger(1);
	
}
