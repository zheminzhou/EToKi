package bloom;

import java.io.Serializable;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 */
public abstract class KCountArray implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1590374813059942002L;

	public static KCountArray makeNew(long cells_, int cbits_, int gap_){
		return makeNew(cells_+1, cells_, cbits_, gap_, 1);
	}
	
	public static KCountArray makeNew(long keys_, long cells_, int cbits_, int gap_, int hashes_){
		return makeNew(keys_, cells_, cbits_, gap_, hashes_, null, 0);
	}
	
//	public static KCountArray makeNew(long keys_, long cells_, int cbits_, int gap_, int hashes_, boolean prefilter_){
//		if(!prefilter_){
//			return makeNew(keys_, cells_, cbits_, gap_, hashes_, 0, 0);
//		}else{
//			long totalbits=cells_*cbits_;
//			long prebits=totalbits/4;
//			long postCells=(totalbits-prebits+cbits_-1)/cbits_;
//			int prehashes=(hashes_+1)/2;
//			return makeNew(keys_, postCells, cbits_, gap_, hashes_, prebits, prehashes);
//		}
//	}
	
	public static KCountArray makeNew(long keys_, long cells_, int cbits_, int gap_, int hashes_, KCountArray prefilter, int prefilterLimit_){
//		assert(keys_>=cells_) : keys_+", "+cells_;
//		assert(cells_>1) : cells_;
		KCountArray kca;
		if(false && keys_<=cells_){
			kca=new KCountArray3(cells_, cbits_, gap_);
//			return new KCountArray4(cells_, cbits_, gap_, 2);
		}else{
//			kca=new KCountArray4(cells_, cbits_, gap_, hashes_); //Single-threaded; most accurate
//			kca=new KCountArray4MT(cells_, cbits_, gap_, hashes_); //Fast
//			kca=new KCountArray5MT(cells_, cbits_, gap_, hashes_); //Less efficient than 4MT
//			kca=new KCountArray6MT(cells_, cbits_, gap_, hashes_); //Fastest but substantial drop in accuracy

//			if(prefilter==null){
//				kca=new KCountArray7MT(cells_, cbits_, gap_, hashes_); //Like 4MT but uses primes
//			}else{
//				kca=new KCountArray8MT(cells_, cbits_, gap_, hashes_, prefilter); //Like 7MT but uses prefilter
//			}
			kca=new KCountArray7MTA(cells_, cbits_, gap_, hashes_, prefilter, prefilterLimit_); //Like 4MT but uses primes

//			if(prefilter==null){
//				kca=new KCountArray9MT(cells_, cbits_, gap_, hashes_); //Like 7MT but uses canonical kmers
//			}else{
//				kca=new KCountArray10MT(cells_, cbits_, gap_, hashes_, prefilter); //Like 8MT but uses canonical kmers
//			}
		}
		kca.initialize();
		
		return kca;
	}
	
	protected KCountArray(long cells_, int cbits_){
		this(cells_, cbits_, 0);
	}
		
	protected KCountArray(final long cells_, int cbits_, int gap_){
		gap=gap_;
		assert(cbits_<=32);
		assert(Integer.bitCount(cbits_)==1);
		assert(Long.bitCount(cells_)==1) || this.getClass()==KCountArray7MT.class : this.getClass();
		
		numArrays=64;
		arrayBits=31-Integer.numberOfLeadingZeros(numArrays);
		arrayMask=numArrays-1;
		
		while(cbits_*cells_<32*numArrays){
			assert(false) : cells_+", "+cbits_+", "+numArrays+", "+(cbits_*cells_)+"<"+(32*numArrays);
			cbits_*=2;
		} //Increases bits per cell so that at minimum each array is size 1
		
		assert(cbits_<=32);
		
		cells=cells_;
		cellBits=cbits_;
		valueMask=(cellBits==32 ? Integer.MAX_VALUE : ~((-1)<<cellBits));
		maxValue=min(Integer.MAX_VALUE, ~((-1)<<min(cellBits,31)));
		cellsPerWord=32/cellBits;
		indexShift=Integer.numberOfTrailingZeros(cellsPerWord);
		cellMask=cellsPerWord-1;
		
		if(verbose){
			System.out.println(description());
		}
	}

	protected KCountArray(final long cells_, int cbits_, int gap_, int arrays_){
		gap=gap_;
		assert(cbits_<=32);
		assert(Integer.bitCount(cbits_)==1);
		assert(Long.bitCount(cells_)==1) || this.getClass()==KCountArray7MT.class || this.getClass()==KCountArray7MTA.class || this.getClass()==KCountArray8MT.class;

		numArrays=arrays_;
		assert(Integer.bitCount(numArrays)==1) : numArrays+", "+cells_+", "+cbits_;
		arrayBits=31-Integer.numberOfLeadingZeros(numArrays);
		arrayMask=numArrays-1;

		while(cbits_*cells_<32*numArrays){
			assert(false) : cells_+", "+cbits_+", "+numArrays+", "+(cbits_*cells_)+"<"+(32*numArrays);
			cbits_*=2;
		} //Increases bits per cell so that at minimum each array is size 1

		assert(cbits_<=32) : "Why?";

		cells=cells_;
		cellBits=cbits_;
		valueMask=(cellBits==32 ? Integer.MAX_VALUE : ~((-1)<<cellBits));
		maxValue=min(Integer.MAX_VALUE, ~((-1)<<min(cellBits,31)));
		cellsPerWord=32/cellBits;
		indexShift=Integer.numberOfTrailingZeros(cellsPerWord);
		cellMask=cellsPerWord-1;

		if(verbose){
			System.out.println(description());
		}
	}

	public abstract int read(long key);
	public int read(long keys[]){throw new RuntimeException("Unimplemented.");}
	public final int read(long key, int k, boolean makeCanonical){return read(makeCanonical ? makeCanonical2(key, k) : key);}

	public abstract void write(long key, int value);

	public void increment(long key){incrementAndReturn(key, 1);}
	public void decrement(long key){decrementAndReturn(key, 1);}
	
	public final int readPrecise(long key, int k, boolean makeCanonical){
		assert(k<=32);
		int b=read(makeCanonical ? makeCanonical2(key, k) : key);
		if(b<1){return b;}
		int a=readLeft(key, k, makeCanonical);
		if(a>=b){return b;}
		int c=readRight(key, k, makeCanonical);
		if(c>=b){return b;}
		return (int)(((long)a+(long)c)/2);
//		return max(a, c);
//		int mid=Tools.min(a, b, c);
//		System.out.println("a="+a+", b="+b+", c="+c+" -> "+mid);
//		return mid;
	}
	
	public final int readPreciseMin(long key, int k, boolean makeCanonical){
		assert(k<=32);
		int b=read(makeCanonical ? makeCanonical2(key, k) : key);
		if(b<1){return b;}
		int a=readLeft(key, k, makeCanonical);
		if(a<1){return a;}
		int c=readRight(key, k, makeCanonical);
		return Tools.min(a, b, c);
	}
	
	/**
	 * @param key Kmer to evaluate
	 * @return Sum of counts of all 4 possible left-adjacent kmers
	 */
	public int readLeft(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Sum of counts of all 4 possible right-adjacent kmers
	 */
	public int readRight(long key, int k, boolean makeCanonical){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Array of counts of all 4 possible left-adjacent kmers
	 */
	public int[] readAllLeft(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}
	/**
	 * @param key Kmer to evaluate
	 * @return Array of counts of all 4 possible right-adjacent kmers
	 */
	public int[] readAllRight(final long key, final int k, boolean makeCanonical, int[] rvec){throw new RuntimeException("Unsupported.");}
	
	public void increment(long[] keys){
		synchronized(this){
			for(long key : keys){
				increment(key);
			}
		}
	}

	/** Returns incremented value. Optional method. */
	public abstract int incrementAndReturn(long key, int incr);
	
	/** Returns decremented value. Optional method. */
	public int decrementAndReturn(long key, int incr){
		throw new RuntimeException("This class "+getClass().getName()+" does not support decrementAndReturn.");
	}
	
	/** Returns unincremented value */
	public abstract int incrementAndReturnUnincremented(long key, int incr);
	
//	/** Returns unincremented value */
//	public final int incrementAndReturnUnincremented(Kmer kmer, int incr){
//		return incrementAndReturnUnincremented(kmer.xor(), incr);
//	}
	
	public int incrementAndReturnUnincremented(long[] keys, int incr){
		throw new RuntimeException("Unimplemented.");
	}
	
	public abstract long[] transformToFrequency();
	public final long[] transformToFrequency(int[][] matrix){
		long[] freq=new long[100000];
		int maxFreq=freq.length-1;

		if(cellBits!=32){
			assert(cellBits>0);
			for(int[] array : matrix){
				for(int i=0; i<array.length; i++){
					int word=array[i];
					int j=cellsPerWord;
					//				System.out.println("initial: word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					for(; word!=0; j--){
						int x=word&valueMask;
						int x2=(int)min(x, maxFreq);
						freq[x2]++;
						word=(word>>>cellBits);
						//					System.out.println("word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					}
					freq[0]+=j;
				}
			}
		}else{
			for(int[] array : matrix){
				for(int i=0; i<array.length; i++){
					int word=array[i];
					int x2=(int)min(word, maxFreq);
					freq[x2]++;
				}
			}
		}
		return freq;
	}
	
	public final long[] transformToFrequency(AtomicIntegerArray[] matrix){
		long[] freq=new long[100000];
		int maxFreq=freq.length-1;

		if(cellBits!=32){
			assert(cellBits>0);
			for(AtomicIntegerArray array : matrix){
				for(int i=0; i<array.length(); i++){
					int word=array.get(i);
					int j=cellsPerWord;
					//				System.out.println("initial: word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					for(; word!=0; j--){
						int x=word&valueMask;
						int x2=(int)min(x, maxFreq);
						freq[x2]++;
						word=(word>>>cellBits);
						//					System.out.println("word = "+word+", j = "+Integer.toHexString(j)+", cellbits="+cellBits);
					}
					freq[0]+=j;
				}
			}
		}else{
			for(AtomicIntegerArray array : matrix){
				for(int i=0; i<array.length(); i++){
					int word=array.get(i);
					int x2=(int)min(word, maxFreq);
					freq[x2]++;
				}
			}
		}
		return freq;
	}
	
	public final ByteBuilder description(){
		ByteBuilder sb=new ByteBuilder();
		long words=cells/cellsPerWord;
		int wordsPerArray=(int)(words/numArrays);
		sb.append("cells:   \t"+cells).append('\n');
		sb.append("cellBits:\t"+cellBits).append('\n');
		sb.append("valueMask:\t"+Long.toHexString(valueMask)).append('\n');
		sb.append("maxValue:\t"+maxValue).append('\n');
		sb.append("cellsPerWord:\t"+cellsPerWord).append('\n');
		sb.append("indexShift:\t"+indexShift).append('\n');
		sb.append("words:   \t"+words).append('\n');
		sb.append("wordsPerArray:\t"+wordsPerArray).append('\n');
		sb.append("numArrays:\t"+numArrays).append('\n');
		sb.append("Memory:   \t"+mem()).append('\n');
		sb.append("Usage:    \t"+String.format(Locale.ROOT, "%.3f%%",usedFraction()*100));
		return sb;
	}
	
	public final String toShortString(){
		return (gap>0 ? "gap = "+gap+"   \t " : "")+"mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+String.format(Locale.ROOT, "%.3f%%",usedFraction()*100);
	}
	
	public final String toShortString(int hashes){
		return (gap>0 ? "gap = "+gap+"   \t " : "")+("hashes = "+hashes+"   \t ")+
				"mem = "+mem()+"   \tcells = "+toKMG(cells)+"   \tused = "+String.format(Locale.ROOT, "%.3f%%",usedFraction()*100);
	}

	@Override
	public final String toString(){
		return description().toString();
	}
	
	public abstract CharSequence toContentsString();
	
	public abstract double usedFraction();
	
	public abstract double usedFraction(int mindepth);
	
	public abstract long cellsUsed(int mindepth);
	
	public final double estimateUniqueKmers(int hashes){
		double f=usedFraction();
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	public final double estimateUniqueKmers(int hashes, int mindepth){
//		assert(false) : this.getClass().getName();
		double f=usedFraction(mindepth);
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	public final double estimateUniqueKmersFromUsedFraction(int hashes, double usedFraction){
		double f=usedFraction;
		double f2=(1-Math.pow(1-f, 1.0/hashes));
		double n=(-cells)*Math.log(1-f2);
		return n;
	}
	
	public final String mem(){
		long mem=(cells*cellBits)/8;
		if(mem<(1<<20)){
			return (String.format(Locale.ROOT, "%.2f KB", mem*1d/(1<<10)));
		}else if(mem<(1<<30)){
			return (String.format(Locale.ROOT, "%.2f MB", mem*1d/(1<<20)));
		}else{
			return (String.format(Locale.ROOT, "%.2f GB", mem*1d/(1<<30)));
		}
	}
	
	public static String toKMG(long x){
		double div=1;
		String ext="";
		if(x>10000000000L){
			div=1000000000L;
			ext="B";
		}else if(x>10000000){
			div=1000000;
			ext="M";
		}else if(x>100000){
			div=1000;
			ext="K";
		}
		return String.format(Locale.ROOT, "%.2f", x/div)+ext;
	}
	
	static final AtomicIntegerArray[] allocMatrix(final int numArrays, final int wordsPerArray){
		final AtomicIntegerArray[] matrix=new AtomicIntegerArray[numArrays];
		final AllocThread[] array=new AllocThread[Tools.min(Tools.max(Shared.threads()/2, 1), numArrays)];
		final AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<array.length; i++){
			array[i]=new AllocThread(matrix, next, wordsPerArray);
		}
		for(int i=0; i<array.length; i++){array[i].start();}
		for(AllocThread at : array){
			while(at.getState()!=Thread.State.TERMINATED){
				try {
					at.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		return matrix;
	}
	
	private static class AllocThread extends Thread{
		
		AllocThread(AtomicIntegerArray[] matrix_, AtomicInteger next_, int wordsPerArray_){
			matrix=matrix_;
			next=next_;
			wordsPerArray=wordsPerArray_;
		}
		
		@Override
		public void run(){
			int x=next.getAndIncrement();
			while(x<matrix.length){
				matrix[x]=new AtomicIntegerArray(wordsPerArray);
				x=next.getAndIncrement();
			}
		}
		
		private final AtomicIntegerArray[] matrix;
		private final AtomicInteger next;
		private final int wordsPerArray;
		
	}
	
	
//	long hash(long x, int y){throw new RuntimeException("Not supported.");}
	abstract long hash(long x, int y);
	
	public static final int min(int x, int y){return x<y ? x : y;}
	public static final int max(int x, int y){return x>y ? x : y;}
	public static final long min(long x, long y){return x<y ? x : y;}
	public static final long max(long x, long y){return x>y ? x : y;}
	
	/** Any necessary initialization. */
	public void initialize(){}
	
	/** Any necessary shutdown steps. */
	public void shutdown(){}
	
	public final long cells;
	public final int cellBits;
	/** Originally this was different than valueMask in the case that valueMask was negative, but now they are the same. */
	public final int maxValue;
	public final int gap; //Set this for convenience on gapped tables to make sure you're using the right table.
	
	protected final int cellsPerWord;
	protected final int indexShift;
	protected final int cellMask;
	protected final int valueMask;
	
	protected static int minArrays=calcMinArrays();
	protected final int arrayBits;
	protected final int numArrays;
	protected final int arrayMask;
	
//	protected static final int arrayBits=6;
//	protected static final int numArrays=1<<arrayBits;
//	protected static final int arrayMask=numArrays-1;
	
	public static boolean verbose=false;
	
	private static final int calcMinArrays(){
		int x=Tools.max(Shared.threads(), 2);
		while(Integer.bitCount(x)!=1){x++;}
		return x;
	}
	
//	public static final boolean isCanonical(long key, int k){
//		assert(k>3 && k<=32);
//		short a=(short)(key&canonMask);
//		short b=AminoAcid.rcompBinaryTable[(int)((key>>(2*(k-4)))&canonMask)];
////		System.out.println("x="+Long.toBinaryString(key)+"\na="+Integer.toBinaryString(a)+"\nb="+Integer.toBinaryString(b)+"\n"+(a>=b));
////		assert(a>=b || isCanonical(AminoAcid.reverseComplementBinaryFast(key, k), k));
//		return a>=b;
//	}
	
//	public static final boolean isCanonical(long key, int k){
//		assert(k>3 && k<=32);
//		short a=(short)(key&canonMask);
//		short b=AminoAcid.rcompBinaryTable[(int)((key>>(2*(k-4)))&canonMask)];
////		System.out.println("x="+Long.toBinaryString(key)+"\na="+Integer.toBinaryString(a)+"\nb="+Integer.toBinaryString(b)+"\n"+(a>=b));
////		assert(a>=b || isCanonical(AminoAcid.reverseComplementBinaryFast(key, k), k));
//		return a>=b;
//	}
	
	public static final boolean isCanonical(long key, int k){
		assert(k>3 && k<=32);
		long b=AminoAcid.reverseComplementBinaryFast(key, k);
		return key>=b;
	}
	
	/** Assumes that the key is not canonical */
	public static final long makeCanonical(final long key, final int k){
		assert(k>3 && k<=32);
//		assert(!isCanonical(key, k));
		final long r=AminoAcid.reverseComplementBinaryFast(key, k);
		assert(r>=key);
//		assert(isCanonical(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key);
		return r;
	}
	
	
	public static final long makeCanonical2(final long key, final int k){
		assert(k>3 && k<=32);
		if(isCanonical(key, k)){return key;}
		long r=AminoAcid.reverseComplementBinaryFast(key, k);
//		assert(isCanonical(r, k)) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
//		assert(AminoAcid.reverseComplementBinaryFast(r, k)==key) : k+"\n"+Long.toBinaryString(key)+"\n"+Long.toBinaryString(r)+"\n"+Long.toBinaryString(AminoAcid.reverseComplementBinaryFast(r, k));
		return r;
	}
	
	public KCountArray prefilter(){
		throw new RuntimeException("TODO: Override");
	}
	
	public void purgeFilter(){
		throw new RuntimeException("TODO: Override");
	}
	
//	private static final short[] canonMask={0, 3, 15, 63, 255, 1023, 4095, 16383};
	private static final long canonK=4;
	static final long canonMask=(1<<(canonK*2))-1; //e.g. 255 for k=4
	
}
