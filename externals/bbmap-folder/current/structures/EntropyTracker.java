package structures;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;

/**
 * Tracks entropy over a sliding window.
 * @author Brian Bushnell
 * @date Oct 6, 2017
 *
 */
public class EntropyTracker {
	
	public static void main(String[] args){
		final int k=args.length>0 ? Integer.parseInt(args[0]) : 2;
		final int window=args.length>1 ? Integer.parseInt(args[1]) : 3;
		final float cutoff=args.length>2 ? Float.parseFloat(args[2]) : 0.7f;
		final boolean highPass=args.length>3 ? Tools.parseBoolean(args[3]) : true;
		
		EntropyTracker et=new EntropyTracker(k, window, false, cutoff, highPass);
		System.err.println(et);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Normal constructor.
	 * @param k_ Kmer length.
	 * @param window_ Window size in bases.
	 */
	public EntropyTracker(int k_, int window_, boolean amino_){
		this(k_, window_, amino_, -1, true);
	}
	
	/**
	 * Allows the use of passes() based on entropy.
	 * @param k_ Kmer length.
	 * @param window_ Window size in bases.
	 * @param cutoff_ Entropy cutoff, 0 (no entropy) to 1 (max entropy).
	 * @param highPass_ True passes entropy of at least cutoff; false fails.
	 */
	public EntropyTracker(int k_, int window_, boolean amino_, float cutoff_, boolean highPass_){
		
		k=k_;
		windowBases=window_;
		windowKmers=windowBases-k+1;
		amino=amino_;
		entropyCutoff=cutoff_;
		highPass=highPass_;
		
		
		assert(k>0 && k<=15 && k<windowBases) : k+", "+windowBases;
		assert(windowKmers>0 && entropyCutoff>=0 && entropyCutoff<=1) : windowKmers+", "+entropyCutoff;
		
		bitsPerBase=(amino ? 5 : 2);
		mask=(k>15 ? -1 : ~((-1)<<(bitsPerBase*k)));
		kmerSpace=(1<<(bitsPerBase*k));//Note: This should be different for amino, but decoding would be a pain
		symbolToNumber=AminoAcid.symbolToNumber(amino);
		symbolToNumber0=AminoAcid.symbolToNumber0(amino);
		
		entropy=makeEntropyArray(windowKmers);
		entropyMult=-1/Math.log(windowKmers);
		counts=new short[kmerSpace];
		countCounts=new short[windowKmers+2];
		countCounts[0]=(short)windowKmers;
		bases=new byte[windowBases];
		
//		assert(false) : amino+", "+windowBases+", "+entropyCutoff+", "+bitsPerBase+", "+
//			Integer.toBinaryString(mask)+", "+kmerSpace+", "+entropy+", "+entropyMult;
		
//		entropyDeltaPlus=makeEntropyDeltaPlus(entropy, entropyMult);
//		entropyDeltaMinus=makeEntropyDeltaMinus(entropy, entropyMult);
	}
	
	/**
	 * @param maxCount Highest possible count; equal to window size in kmers.
	 * @return Entropy array
	 */
	private static final double[] makeEntropyArray(int maxCount){
		final double[] array=new double[maxCount+2]; //+2 to handle temporary condition of an extra kmer
		final double mult=1d/maxCount;
		for(int i=1; i<array.length; i++){//First element (kmer count of 0) contains zero entropy, and yields NaN, so is skipped.
			double pk=i*mult;
			array[i]=pk*Math.log(pk);
		}
		return array;
	}
	
	/**
	 * @param entropy Entropy array
	 * @param entropyMult Multiplier applied to each array element
	 * @return entropyDeltaPlus array
	 */
	@SuppressWarnings("unused")
	private static final double[] makeEntropyDeltaPlus(double[] entropy, double entropyMult){
		final double[] array=new double[entropy.length];
		for(int i=0; i<entropy.length-1; i++){//Last element is never used.
			array[i]=(entropy[i+1]-entropy[i])*entropyMult;
		}
		return array;
	}
	
	/**
	 * @param entropy Entropy array
	 * @param entropyMult Multiplier applied to each array element
	 * @return entropyDeltaMinus array
	 */
	@SuppressWarnings("unused")
	private static final double[] makeEntropyDeltaMinus(double[] entropy, double entropyMult){
		final double[] array=new double[entropy.length];
		for(int i=1; i<entropy.length; i++){//First element is never used.
			array[i]=(entropy[i-1]-entropy[i])*entropyMult;
		}
		return array;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @return Number of unique kmers in current window. */
	public int unique(){return unique;}
	
	/** @return Number of undefined bases in current window. */
	public int ns(){return ns;}
	
	/** @return Sequence position of rightmost base in window. */
	public int rightPos(){return len-1;}
	
	/** @return Sequence position of leftmost base in window. */
	public int leftPos(){return len-windowBases;}

	/** @return Window size in bases. */
	public int windowBases() {return windowBases;}

	/** @return Entropy cutoff. */
	public float cutoff() {return entropyCutoff;}
	
	/*--------------------------------------------------------------*/
	/*----------------      Entropy Calculation     ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculate entropy in current window.
	 * @return Entropy in current window.
	 */
	public float calcEntropy(){
		if(speed==FAST){return calcEntropyFast();}
		else if(speed==MEDIUM){return calcEntropyMedium();}
		else if(speed==SLOW){return calcEntropySlow();}
		else{return calcEntropySuperSlow();}
	}
	
	/** Sub-method of calcEntropy() */
	private float calcEntropyFast(){
		final float f=(float)(currentEsum*entropyMult);
//		final float f=(float)currentEsum; //For delta arrays
		
		//Avoid potential negative numbers due to underflow
		return f>0 ? f : 0;
	}
	
	/** Sub-method of calcEntropy()
	 * Calculates entropy from countCounts using precalculated entropy array and early exit. */
	private float calcEntropyMedium(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Sum of kmer counts, used for loop early exit
		int cSum=countCounts[0];
		
		//Simpler loop with no early exit.  Slower for complex sequences, but faster for homopolymonomers.
//		for(int i=1; i<countCounts.length; i++){
		
		//Loop over all nonzero kmer counts
		for(int i=1; cSum<windowKmers; i++){
			//Number of unique kmers with count i
			final int cc=countCounts[i];
			cSum+=cc;
			
			//Entropy contribution for each unique kmer with count i
			double pklogpk=entropy[i];
			
			//Add the entropy contribution for all unique kmers with count i
			eSum+=(cc*pklogpk);
		}
		//eSum now holds negative entropy in bits.
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*entropyMult);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(String.format("%.3f", eSum)+"\t"+String.format("%.3f", entropyMult)+"\t"
					+String.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/** Sub-method of calcEntropy()
	 * Calculates entropy from counts using precalculated entropy array. */
	private float calcEntropySlow(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Loop over all nonzero kmer counts
		for(int count : counts){
			
			//Entropy contribution for this kmer's count
			double pklogpk=entropy[count];
			eSum+=pklogpk;
		}
		//eSum now holds negative entropy in bits.
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*entropyMult);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(String.format("%.3f", eSum)+"\t"+String.format("%.3f", entropyMult)+"\t"
					+String.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/**
	 * Sub-method of calcEntropy()
	 * Demonstrates explicit entropy calculation.
	 * 
	 * Definition:
	 * Entropy, or information content, can be calculated using kmer counts of a sequence.
	 * 
	 * Probability of an event (unique kmer), pk, is that kmer's count divided by the number of kmers.
	 * Entropy contribution from that kmer is -pk*log(pk).
	 * Total entropy is sum of (-pk*log(pk)) for all kmer counts.
	 * entropyMult is simply a multiplier to convert the entropy measure (in bits) to a convenient 0-1 scale,
	 * corresponding to the inverse of the maximum possible entropy.
	 */
	private float calcEntropySuperSlow(){
		//Sum of entropy contributions from each kmer
		double eSum=0;
		
		//Loop over all nonzero kmer counts
		for(int count : counts){
			//Prevent NaN and INF
			if(count>0){
				//Fraction of total represented by this kmer
				double pk=count/(double)windowKmers;

				//Entropy contribution for this kmer's count
				double npklogpk=-pk*Math.log(pk);
				eSum+=npklogpk;
			}
		}
		//eSum now holds entropy in bits.
		
		//Multiplier to convert entropy to 0-1 scale.
		double multiplier=1/Math.log(windowKmers);
		
		//Adjust entropy to 0-1 scale based on window size
		float e=(float)(eSum*multiplier);
		assert(e>=0 && e<=1) : e+", "+eSum+", "+entropyMult+"\n"+Arrays.toString(entropy)+"\n"+this;
		
		//Get rid of negative zero
		if(e<=0){e=0;}
		
		if(verbose){
			System.err.println(String.format("%.3f", eSum)+"\t"+String.format("%.3f", entropyMult)+"\t"
					+String.format("%.3f", e)+"\t"+len+"\t"+unique+"\t"+basesToString());
		}
		
		return e;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculate the average entropy of a sequence.
	 * @param bases Sequence as bytes
	 * @param allowNs True if windows containing undefined bases should be included
	 * @return Average entropy
	 */
	public float averageEntropy(byte[] bases, boolean allowNs){
		//Reset the initial state
		clear();
		
		//Position in sequence
		int i=0;
		
		//Accumulated entropy
		double sum=0;
		
		//Number of entropy measurements
		int divisor=0;
		
		//Prefill the first window
		for(final int lim=Tools.min(bases.length, windowBases); i<lim; i++){add(bases[i]);}
		
		//Calculate entropy for the first window.
		//This allows entropy to pass if it is high enough even though the sequence is shorter than window length.
		if(allowNs || ns==0){
			sum+=calcEntropy();
			divisor++;
		}
		
		//Calculate entropy for remaining windows
		for(; i<bases.length; i++){
			add(bases[i]);
			if(allowNs || ns==0){
				sum+=calcEntropy();
				divisor++;
			}
		}
		
		//Calculate the average
		double avg=(sum/(Tools.max(1, divisor)));
		return (float)avg;
	}
	
	/**
	 * Calculate entropy in the window and compare to the cutoff.
	 * If Ns are important they should be handled externally with ns().
	 * @return True if the entropy passes the cutoff.
	 */
	public boolean passes(){
		//This function should only be used if entropyCutoff is set.
		assert(entropyCutoff>=0);
		float e=calcEntropy();
		
		//XOR: highPass inverts truth of comparison.
		return highPass ^ (e<entropyCutoff);
	}
	
	/**
	 * Calculate average entropy of the sequence and compare to the cutoff.
	 * @param sequence Sequence to measure.
	 * @param allowNs True if entropy should be calculated in windows containing Ns.
	 * @return True if the average entropy passes the cutoff.
	 */
	public boolean passes(byte[] sequence, boolean allowNs){
		//This function should only be used if entropyCutoff is set.
		assert(entropyCutoff>=0);
		float e=averageEntropy(sequence, allowNs);
		
		//XOR: highPass inverts truth of comparison.
		return highPass ^ (e<entropyCutoff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Mutating Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Slide the window by adding a new base.
	 * @param b Base to add.
	 */
	public void add(final byte b){
		//Test initial state
		assert(!verify || verify()) : this;
		
		final byte oldBase=bases[pos]; //Leftmost base, about to be overwritten
		
		//Increment length
		len++;
		
		{//Add a new rightmost base
			bases[pos]=b;
			final int n=symbolToNumber0[b];
			kmer=((kmer<<bitsPerBase)|n)&mask; //Generate new rightmost kmer using the new base
			
			//Update number of Ns in current window
			if(!isFullyDefined(b)){
				ns++;
				assert(ns<=windowBases+1) : "There are more Ns than bases in the window:\n"+this;
			}
			
			if(len>=k){//Add a kmer
				//System.err.println("Adding "+kmer);
				
				final short oldCount=counts[kmer];
				
				/* Update unique kmer count */
				if(oldCount<1){
					assert(oldCount==0) : "An incoming array has negative counts: \n"+this;
					unique++;
				}
				
				/* Decrement the old countCount */
				countCounts[oldCount]--;
				
				/* countCounts[0] could be temporarily -1 at this point; for all others, min is 0. */
				assert(countCounts[oldCount]>=-1) : this;
				
				/* Increment the kmer count */
				final short newCount=counts[kmer]=(short)(oldCount+1);
				
				/* The count could at most be 1 more than the total window kmers here temporarily */
				assert(newCount<=windowKmers+1) : this;
				
				/* Increment the new countCount */
				countCounts[newCount]++;
				
				/* Update entropy */
				currentEsum=currentEsum+entropy[newCount]-entropy[oldCount];
//				currentEsum+=entropyDeltaPlus[oldCount];

				assert(!verify || unique==Tools.cardinality(counts)) : this;
				assert(!verify || (Tools.sum(countCounts)>0 && (Tools.sum(countCounts)<=windowKmers+1))) : this;
			}
		}
		
		//At this point the state is inconsistent as it may have one too many kmers.
		
		if(pos2>=0){//Remove the leftmost base
			final byte b2=bases[pos2];//This is not the leftmost base, but the base to the right of the leftmost kmer
			final int n2=symbolToNumber0[b2];
			//Generate old leftmost kmer using an internal base near the left end
			kmer2=((kmer2<<bitsPerBase)|n2)&mask;

			if(len>windowBases){//Remove a kmer, only if a base is leaving the window
				//System.err.println("Removing "+kmer2);
				
				//Update number of Ns in current window
				if(!isFullyDefined(oldBase)){
					ns--;
					assert(ns>=0) : "There are fewer than 0 Ns in the window:\n"+this;
				}
				
				assert(kmer2>=0) : "A negative kmer was observed: "+kmer2+"\n"+this;
				
				final short oldCount=counts[kmer2];
				assert(oldCount>0) : "Attempting to decrement a nonpositive count: \n"+oldCount+"\n"+this;
				
				//Decrement the old countCount
				countCounts[oldCount]--;
				
				assert(countCounts[oldCount]>=0) : "A countCount became negative: \n"+countCounts[oldCount]+"\n"+this;
				
				//Decrement the kmer count
				final short newCount=counts[kmer2]=(short)(oldCount-1);
				
				/* Increment the new countCount */
				countCounts[newCount]++;
				
				/* Update unique kmer count */
				if(newCount<1){
					assert(newCount==0): "An outgoing array has negative counts: \n"+this;
					unique--;
				}
				
				/* Update entropy */
				currentEsum=currentEsum+entropy[newCount]-entropy[oldCount];
//				currentEsum+=entropyDeltaMinus[oldCount];
				
				assert(!verify || unique==Tools.cardinality(counts)) : this;
				assert(!verify || (Tools.sum(countCounts)>=0 && (Tools.sum(countCounts)<=windowKmers))) : this;
			}
		}
		
		//Update position pointers
		//Can use modulo, but this is faster because the branch is normally skipped.
		//Ternary conditionals are also slower.
		pos++;
		pos2++;
		if(pos>=windowBases){pos=0;}
		if(pos2>=windowBases){pos2=0;}
		
		//Test final state.
		assert(!verify || verify()) : this;
	}
	
	/**
	 * Reset fields to prepare for a new sequence.
	 */
	public void clear(){
		//Reset scalars
		kmer=kmer2=len=pos=unique=ns=0;
		pos2=0-windowBases+k-1;
		currentEsum=0;
		
		//Clear mutable arrays.  Bases does not need to be cleared.
		Arrays.fill(counts, (short)0); //Time proportional to kmer space
		
		//Note - countCounts are only needed for medium speed mode.
		Arrays.fill(countCounts, (short)0); //Time proportional to window size
		
		//Sets the number of kmers with a count of zero to maximum.
		countCounts[0]=(short)windowKmers;
		
		//Verify the state was cleared
		assert(!verify || verifyClear()) : this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Verify that mutable fields were properly cleared.
	 * Throws an assertion error upon failure.
	 * @return True.
	 */
	public boolean verifyClear(){
		for(int c : counts){assert(c==0) : this;}
		for(int i=1; i<countCounts.length; i++){assert(countCounts[i]==0) : this;}
		assert(kmer==0 && kmer2==0) : this;
		assert(pos==0) : this;
		assert(len==pos) : this;
		assert(ns==0) : this;
		assert(unique==0) : this;
		assert(pos2<0) : this;
		assert(currentEsum==0) : this;
		return true;
	}
	
	/**
	 * Verify that internal state is consistent.
	 * Throws an assertion error upon failure.
	 * @return True.
	 */
	public boolean verify(){
		
		//Number of unique kmers in the window
		int existSum=0;
		
		//Total number of kmers in the window
		int countSum=0;
		
		//Number of undefined symbols in the window
		int nSum=0;
		
		//Check the kmer counts
		for(int c : counts){
			assert(c>=0 && c<=windowKmers) : "A kmer count exceeds the possible bounds.\n"+this;
			if(c>0){
				existSum++;
				countSum+=c;
			}
		}
		
		//Check the countCounts
		for(int cc : countCounts){
			assert(cc>=0 && cc<=windowKmers) : "A countCount exceeds the possible bounds.\n"+this;
		}
		
		//Count undefined symbols
		for(byte b : bases){
			if(!isFullyDefined(b)){nSum++;}
		}
		
		//Number of kmers with count 0
		final int cc0=countCounts[0];
		
		//Sum of countCounts
		final int ccSum=(int)Tools.sum(countCounts);
		
		//Sum of nonzero countCounts; should equal the number of unique kmers
		final int ccSum1=ccSum-cc0;
		
		//Do assertions
		assert(len<windowBases || ns==nSum) : this;
		assert(existSum==unique) : this;
		assert(existSum==ccSum1) : this;
		assert(existSum>=0 && existSum<=windowKmers) : this;
		assert(len<windowBases || countSum==windowKmers) : this;
		assert(ccSum==windowKmers);
		assert(pos==len%windowBases) : this;
		
		//Ensure different entropy calculation methods are consistent
		if(len>=windowKmers){
			float a=calcEntropyFast();
			float b=calcEntropyMedium();
			float c=calcEntropySlow();
			float d=calcEntropySuperSlow(); //"d" should be identical to "c".
			assert(Tools.absdif(a, b)<0.0001) : "Fast and Medium entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+this;
			assert(Tools.absdif(b, c)<0.000001) : "Medium and Slow entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+this;
			assert(c==d) : "Slow and SuperSlow entropy differ:\n"+a+"\n"+b+"\n"+c+"\n"+d+"\n"+this;
		}
		return true;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('\n');
		sb.append("kmer\t"+kmer).append('\n');
		sb.append("kmer2\t"+kmer2).append('\n');
		sb.append("pos\t"+pos).append('\n');
		sb.append("pos2\t"+pos2).append('\n');
		sb.append("len\t"+len).append('\n');
		sb.append("unique\t"+unique).append('\n');
		sb.append("ns\t"+ns).append('\n');
		sb.append('\n');
		sb.append("k\t"+k).append('\n');
		sb.append("windowBases\t"+windowBases).append('\n');
		sb.append("windowKmers\t"+windowKmers).append('\n');
		sb.append("mask\t"+mask).append('\n');
		sb.append('\n');
		sb.append("cardinality\t"+Tools.cardinality(counts)).append('\n');
		sb.append("counts\t"+Arrays.toString(counts)).append('\n');
		sb.append("ccounts\t"+Arrays.toString(countCounts)).append('\n');
		sb.append("bases\t"+basesToString()).append('\n');
		sb.append("entropy\t"+Arrays.toString(entropy)).append('\n');
		sb.append("entropyMult\t"+entropyMult).append('\n');
		sb.append("entropySum\t"+currentEsum).append('\n');
		return sb.toString();
	}
	
	/** Returns the ring buffer as a String in its correct order. */
	public String basesToString(){
		StringBuilder sb=new StringBuilder(bases.length);
		for(int i=0; i<bases.length; i++){
			byte b=(bases[(i+pos)%bases.length]);
			if(b==0){b='N';}
			sb.append((char)b);
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Current leading kmer; rightmost k bases of window */
	int kmer=0;
	/** Current trailing kmer; leftmost k-1 bases of window plus the removed base */
	int kmer2=0;
	/** Position in ring buffer to place incoming bases */
	int pos=0;
	/** Position in ring buffer to read next base for kmer2 */
	int pos2=0;
	/** Current number of processed bases.  Equal to pos without being reset at buffer wrap. */
	int len=0;
	/** Number of unique kmers in the current window */
	int unique=0;
	/** Number of undefined bases in the current window */
	int ns;

	/** Current sum of entropy from kmers in the current window */
	double currentEsum=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Mutable Arrays        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of times each kmer occurs in current window.
	 * Indexed by the kmer's numeric value.
	 * counts[0] stores the count of the kmer AAA, if k=3. */
	private final short[] counts;
	
	/** Number of instances of each number in counts.
	 * countCounts[0] stores the number of kmers with count 0.
	 * This is only needed in medium speed mode (or verify mode). */
	private final short[] countCounts;
	
	/** Ring buffer of bases in current window.
	 * Not strictly necessary, but convenient. */
	private final byte[] bases;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Kmer length for entropy calculation */
	private final int k;
	/** Window length for entropy calculation */
	private final int windowBases;
	/** Number of kmers in the window */
	private final int windowKmers;
	/** Amino acid mode */
	private final boolean amino;

	/** Bits per symbol */
	private final int bitsPerBase;
	/** Mask for sliding kmers */
	private final int mask;
	
	/** Minimum entropy to be considered "complex", on a scale of 0-1; optional */
	private final float entropyCutoff;
	/** Pass entropy values above the cutoff */
	private final boolean highPass;
	
	/** Number of possible unique kmers */
	private final int kmerSpace;
	/** A precalculated constant */
	private final double entropyMult;
	/** Array of precalculated constants */
	private final double[] entropy;

	/** Translation table yielding 0 if undefined */
	private final byte[] symbolToNumber0;
	/** Translation table yielding -1 if undefined */
	private final byte[] symbolToNumber;
	
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	//Note:  These incur fewer operations, but in testing, were not faster.
//	/** For calculating entropy running average quickly when adding a kmer.
//	 * entropyDeltaPlus[i] = (entropy[i+1]-entropy[i])*entropyMult */
//	private final double[] entropyDeltaPlus;
//	/** For calculating entropy running average quickly when removing a kmer.
//	 * entropyDeltaMinus[i] = (entropy[i-1]-entropy[i])*entropyMult */
//	private final double[] entropyDeltaMinus;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Entropy calculation speed constants.
	 * FAST is less precise for long sequences.
	 * MEDIUM is probably most precise. */
	public static final int FAST=0, MEDIUM=1, SLOW=2, SUPERSLOW=3;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Entropy calculation mode */
	public static int speed=FAST;
	
	/** Verify consistency of related data structures (slow) */
	public static boolean verify=false;
	/** Verbose output */
	public static final boolean verbose=false;
	
}
