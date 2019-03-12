package bloom;

import shared.Shared;

/**
 * @author Brian Bushnell
 * @date Dec 2, 2014
 *
 */
public abstract class KmerCountAbstract {

	protected static final long[] transformToFrequency(int[] count){
		long[] freq=new long[2000];
		int max=freq.length-1;
		for(int i=0; i<count.length; i++){
			int x=count[i];
			x=min(x, max);
			freq[x]++;
		}
		return freq;
	}

	protected static final long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	protected static final long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	protected static final int min(int x, int y){return x<y ? x : y;}
	protected static final int max(int x, int y){return x>y ? x : y;}
	
	public static byte minQuality=6;
	public static long readsProcessed=0;
	public static long maxReads=-1;
	public static int BUFFERLEN=500;
	
	public static float minProb=0.5f;

	public static long keysCounted=0;
	
	public static int THREADS=Shared.threads();
	public static final boolean verbose=false;
	public static boolean PREJOIN=false;
	public static boolean CANONICAL=false;
	public static boolean KEEP_DUPLICATE_KMERS=true;
	public static boolean SKETCH_MODE=false;
	public static boolean STORE_HASHED=false;
	
}
