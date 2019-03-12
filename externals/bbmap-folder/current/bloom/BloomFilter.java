package bloom;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import structures.LongList;

/**
 * Wraps a KCountArray and provides multithreaded reference loading.
 * 
 * @author Brian Bushnell
 * @date April 23, 2018
 *
 */
public class BloomFilter implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3987955563503838492L;
	
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
		
		//Create an instance of this class
		BloomFilter x=new BloomFilter(args);
		
		System.err.println(x.filter.toShortString());
		t.stop("Time: \t");
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BloomFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		KmerCountAbstract.CANONICAL=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		int k_=31;
		int bits_=1;
		int hashes_=2;
		int minConsecutiveMatches_=3;
		float memFraction=1;
		boolean rcomp_=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
				assert(k_<=31 && k_>=1);
			}else if(a.equals("hashes")){
				hashes_=Integer.parseInt(b);
				assert(hashes_<=10000 && hashes_>=1);
			}else if(a.equals("minhits")){
				minConsecutiveMatches_=Integer.parseInt(b);
				assert(minConsecutiveMatches_>=1);
			}else if(a.equals("bits")){
				bits_=Integer.parseInt(b);
			}else if(a.equals("memfraction")){
				memFraction=Float.parseFloat(b);
			}else if(a.equals("extra")){
				if(b==null){extra.clear();}
				else{
					for(String s : b.split(",")){extra.add(s);}
				}
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		k=k_;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;
		

		assert(bits==1 || bits==2 || bits==4 || bits==8 || bits==16 || bits==32) : "Bits must be a power of 2.";
		
		{//Process parser fields
			in1=parser.in1;
			in2=parser.in2;
		}
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		filterMemory=setMemory(memFraction);
		
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=load();
	}
	
	public BloomFilter(String in1_, String in2_, ArrayList<String> extra_, int k_, int bits_, int hashes_,
			int minConsecutiveMatches_, boolean rcomp_, boolean ecco_, boolean merge_, float memFraction){
		if(extra_!=null){
			for(String s : extra_){extra.add(s);}
		}

		filterMemory=setMemory(memFraction);
		
		in1=in1_;
		in2=in2_;
		k=k_;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;
		ecco=ecco_;
		merge=merge_;

		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=load();
	}

	public BloomFilter(boolean bbmapIndex_, int k_, int bits_, int hashes_, int minConsecutiveMatches_, boolean rcomp_) {
		assert(bbmapIndex_);
		filterMemory=setMemory(0.75);
		
		in1=null;
		in2=null;
		k=k_;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;

		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=loadFromIndex();
	}
	
	private static long setMemory(double mult){
		Shared.printMemory();
		
		Runtime rt=Runtime.getRuntime();
		final long mmemory=rt.maxMemory();
		final long tmemory=rt.totalMemory();
		final long fmemory=rt.freeMemory();
		final long umemory=tmemory-fmemory;
		
		double xmsRatio=Shared.xmsRatio();
		double usableMemory=Tools.max(((mmemory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), mmemory*0.45);
		double availableMemory=usableMemory-umemory;
		double filterMemory=availableMemory*mult;
		
//		System.err.println((long)(usableMemory/1000000)+", "+(long)(availableMemory/1000000)+", "+(long)(filterMemory/1000000));
		
		return (long)filterMemory;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private KCountArray7MTA load(){
		final int cbits=bits;
		final long totalBits=8*filterMemory;
		final long cells=(OVERRIDE_CELLS>0 ? OVERRIDE_CELLS : totalBits/cbits);
//		System.err.println("filterMemory="+filterMemory+", cells="+cells);
		KCountArray7MTA kca=(KCountArray7MTA)KmerCount7MTA.makeKca(in1, in2, extra==null || extra.isEmpty() ? null : extra,
				k, cbits, 0, cells, hashes, minq, rcomp, ecco, merge,
				maxReads, 1, 1, 1, 1, null, 0, Shared.AMINO_IN);
		return kca;
	}

	private KCountArray7MTA loadFromIndex(){
		KmerCountAbstract.CANONICAL=true;
		final int cbits=bits;
		final long totalBits=8*filterMemory;
		final long cells=totalBits/cbits;
		outstream.println("Filter Memory = "+String.format("%.2f GB", filterMemory/(double)(1024*1024*1024)));
		KCountArray7MTA kca=(KCountArray7MTA)KmerCount7MTA.makeKcaFromIndex(k, cbits, cells, hashes, rcomp);
		return kca;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean passes(Read r1, Read r2, final int thresh) {
		boolean pass=passes(r1, thresh);
		return pass && passes(r2, thresh);
	}
	
	public int minCount(Read r) {
		if(r==null || r.length()<k-1){return -1;}
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int min=Integer.MAX_VALUE;
		int counted=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				counted++;
				int count=getCount(kmer, rkmer);
				min=Tools.min(min, count);
				if(count==0){
//					assert(false) : counted+", "+min;
					break;
				}
			}
		}
//		assert(false) : counted+", "+min;
		return counted>0 ? min : -1;
	}
	
	public boolean hasHighCountFraction(Read r, final int thresh, final float fraction) {
		if(r==null || r.length()<k-1){return false;}
		final byte[] bases=r.bases;
		final int kmers=r.length()-k+1;
		
		final int minHigh=Math.round(fraction*kmers);
		final int maxLow=kmers-minHigh;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int counted=0;
		int low=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				counted++;
				int count=getCount(kmer, rkmer);
				if(count<thresh){
					low++;
					if(low>maxLow){return false;}
				}
			}
		}
		return true;
	}
	
	public boolean isJunk(Read r1, Read r2, int range){
		assert(bits>1);
		if(r2==null || r2.length()<k){return isJunk(r1, range);}
		if(r1.length()<k){return isJunk(r2, range);}
		if(getLeftCount(r1.bases, range)>1 || getLeftCount(r2.bases, range)>1){return false;}
//		return getRightCount(r1.bases, range)<3 && getRightCount(r2.bases, range)<3; //&& is more correct; || allows for a fuller filter. 
		return getRightCount(r1.bases, range)<3 || getRightCount(r2.bases, range)<3;
	}
	
	public boolean isJunk(Read r, int range){
		assert(bits>1);
		if(r.length()<k){return true;}
		return getLeftCount(r.bases, range)<2 && getRightCount(r.bases, range)<2;
	}
	
	private int getLeftCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		long kmer=0, rkmer=0;
		int len=0;
		final int stop=Tools.min(bases.length, k+range-1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=0; i<stop; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}
			else{
				len++;
				if(len>=k){
					counted++;
					int count=getCount(kmer, rkmer);
					min=Tools.min(min, count);
				}
			}
		}
		return counted>0 ? min : -1;
	}
	
	private int getRightCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		long kmer=0, rkmer=0;
		int len=0;
		final int start=Tools.max(0, bases.length-k-range+1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=start; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}
			else{
				len++;
				if(len>=k){
					counted++;
					int count=getCount(kmer, rkmer);
					min=Tools.min(min, count);
				}
			}
		}
		return counted>0 ? min : -1;
	}
	
	public boolean passes(Read r, final int thresh) {
		if(r==null || r.length()<k+minConsecutiveMatches-1){return true;}
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int streak=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				boolean found=contains(kmer, rkmer, thresh);
				if(found){
					streak++;
					if(streak>=minConsecutiveMatches){return false;}
				}else{streak=0;}
			}
		}
		return true;
	}
	
	public boolean matches(Read r, LongList keys, final int thresh) {
		return !passes(r, keys, thresh);
	}
	
	public boolean matchesEither(Read r1, Read r2, LongList keys, final int thresh) {
		boolean match=!passes(r1, keys, thresh);
		return match || !passes(r2, keys, thresh);
	}
	
	public boolean passes(Read r1, Read r2, LongList keys, final int thresh) {
		boolean pass=passes(r1, keys, thresh);
		return pass && passes(r2, keys, thresh);
	}
	
	public boolean passes(Read r, LongList keys, final int thresh) {
		if(r==null || r.length()<k+minConsecutiveMatches-1){return true;}
		if(minConsecutiveMatches<2){return passes(r, thresh);}
		keys.clear();
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){
				if(len>=k){keys.add(-1);}
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
				if(len>=k){keys.add(toKey(kmer, rkmer));}
			}
		}
		return passes(keys, thresh);
	}
	
	public boolean passes(final LongList keys, final int thresh) {
		assert(minConsecutiveMatches>1);
		final long[] array=keys.array;
		final int len=keys.size;
		
		for(int i=minConsecutiveMatches-1; i<len; i+=minConsecutiveMatches){
			final boolean found;
			{
				final long key=array[i];
				found=(key<0 ? false : filter.read(key)>=thresh);
			}
			if(found){
				int streak=1;
				for(int j=1; j<minConsecutiveMatches; j++){
					final long key=array[i-j];
					if(key<0 || filter.read(key)<thresh){break;}
					else{streak++;}
				}
				if(streak>=minConsecutiveMatches){return false;}
				for(int j=1; j<minConsecutiveMatches && j+i<len; j++){
					final long key=array[i+j];
					if(key<0 || filter.read(key)<thresh){break;}
					else{streak++;}
					if(streak>=minConsecutiveMatches){return false;}
				}
			}
		}
		return true;
	}
	
	public int getCount(final long kmer, final long rkmer){
		final long key=toKey(kmer, rkmer);
		return filter.read(key);
	}
	
	public int getCount(final long key){
		return filter.read(key);
	}
	
	public boolean contains(final long kmer, final long rkmer, final int thresh){
		final long key=toKey(kmer, rkmer);
		return filter.read(key)>=thresh;
	}
	
	public long toKey(final long kmer, final long rkmer){
		return (rcomp ? Tools.max(kmer, rkmer) : kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long maxReads=-1;
	boolean ecco=false;
	boolean merge=false;
	byte minq=0;

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private ArrayList<String> extra=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final KCountArray7MTA filter;
	final int k;
	final int bits;
	final int hashes;
	final int minConsecutiveMatches;
	
	final int shift;
	final int shift2;
	final long mask;
	final boolean rcomp;

//	private final long usableMemory;
	private final long filterMemory;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static long OVERRIDE_CELLS=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private transient PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
}
