package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.TranslateSixFrames;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.EntropyTracker;
import structures.ListNum;
import tax.TaxNode;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMakerMini extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 */
	public SketchMakerMini(SketchTool tool_, int mode_, float minEntropy_){
		
		tool=tool_;
		mode=mode_;
		
		aminoShift=AminoAcid.AMINO_SHIFT;
		if(!aminoOrTranslate()){
			shift=2*k;
			shift2=shift-2;
			mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		}else{
			shift=aminoShift*k;
			shift2=shift-aminoShift;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
		}
		
		if(AUTOSIZE && (mode==ONE_SKETCH || mode==PER_FILE)){
			heap=new SketchHeap(Tools.max(tool.stTargetSketchSize, (int)(80000*Tools.mid(1, AUTOSIZE_FACTOR, 32))), tool.minKeyOccuranceCount, tool.trackCounts);
		}else{
			heap=new SketchHeap(tool.stTargetSketchSize, tool.minKeyOccuranceCount, tool.trackCounts);
		}
		
		if(minEntropy_>0){
			eTracker=new EntropyTracker(entropyK, k, (amino || translate), minEntropy_, true);
		}else{
			eTracker=null;
		}
		
		if(translate){
			gCaller=CallGenes.makeGeneCaller(pgm);
		}else{
			gCaller=null;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public ArrayList<Sketch> toSketches(final String fname, float samplerate, long reads){
		heap.clear(false); //123
		final String simpleName;
		
		final FileFormat ffin1, ffin2;
		if(fname.indexOf('#')>=0 && FileFormat.isFastq(ReadWrite.rawExtension(fname)) && !new File(fname).exists()){
			ffin1=FileFormat.testInput(fname.replaceFirst("#", "1"), FileFormat.FASTQ, null, true, true);
			ffin2=FileFormat.testInput(fname.replaceFirst("#", "2"), FileFormat.FASTQ, null, true, true);
		}else{
			ffin1=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			ffin2=null;
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			simpleName=ffin1.simpleName();
			heap.setFname(simpleName);
			cris=ConcurrentReadInputStream.getReadInputStream(reads, true, ffin1, ffin2, null, null);
			if(samplerate!=1){cris.setSampleRate(samplerate, sampleseed);}
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		if(mode==ONE_SKETCH || mode==PER_FILE){
			 if(heap.name0()==null){heap.setName0(simpleName);}
		}
		ArrayList<Sketch> sketches=processInner(cris);
		
		errorState|=ReadWrite.closeStream(cris);
		sketchesMade++;
		return sketches;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Iterate through the reads */
	ArrayList<Sketch> processInner(ConcurrentReadInputStream cris){
		assert(heap.size()==0);
		ArrayList<Sketch> sketches=new ArrayList<Sketch>(mode==ONE_SKETCH || mode==PER_FILE ? 1 : 8);
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//As long as there is a nonempty read list...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				processReadPair(r1, r2);
				if(mode!=ONE_SKETCH && mode!=PER_FILE){
					if(heap!=null && heap.size()>0 && heap.maxLen()>=Tools.max(1, minSketchSize)){
						int size=heap.size();
						Sketch sketch=new Sketch(heap, false, tool.trackCounts, null);
						assert(sketch.array.length>0) : sketch.array.length+", "+size;
						sketches.add(sketch);
						sketchesMade++;
					}
					if(heap!=null){heap.clear(false);}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		if(mode==ONE_SKETCH || mode==PER_FILE){
			Sketch sketch=new Sketch(heap, false, tool.trackCounts, null);
			sketches.add(sketch);
			sketchesMade++;
		}
		heap.clear(true);
		return sketches;
	}

	void processReadPair(Read r1, Read r2){
		//Track the initial length for statistics
		final int initialLength1=r1.length();
		final int initialLength2=r1.mateLength();

		//Increment counters
		readsProcessed+=r1.pairCount();
		basesProcessed+=initialLength1+initialLength2;
		
		if(mode!=ONE_SKETCH && mode!=PER_FILE){
			int expectedSize=toSketchSize(initialLength1+initialLength2, -1, -1, targetSketchSize);
			if(heap==null || heap.capacity()<expectedSize){heap=new SketchHeap(expectedSize, tool.minKeyOccuranceCount, tool.trackCounts);}
		}
		
		if(tool.mergePairs && r1!=null && r2!=null){
			final int insert=BBMerge.findOverlapStrict(r1, r2, false);
			if(insert>0){
				heap.genomeSequences++;
				r2.reverseComplement();
				r1=r1.joinRead(insert);
				r2=null;
			}
		}
		
		processRead(r1);
		if(r2!=null){processRead(r2);}

		if(heap.name0()==null){
			heap.setName0(r1.id);
		}

		TaxNode tn=null;
		if(taxtree!=null && heap.taxID<0){
			try {
				tn=taxtree.parseNodeFromHeader(r1.id, true);
			} catch (Throwable e) {}
			if(tn!=null){
				heap.taxID=tn.id;
				if(heap.taxName()==null){
					heap.setTaxName(tn.name);
				}
			}
		}
		assert(heap.taxID<0 || heap.taxName()!=null) : heap.taxID+", "+heap.taxName()+", "+heap.name()+", "+tn;
	}

	public void processRead(final Read r){
		if(amino){
			processReadAmino(r);
		}else if(translate){
			processReadTranslated(r);
		}else{
			processReadNucleotide(r);
		}
	}
	
	public void processReadTranslated(final Read r){
		assert(!r.aminoacid());
		final ArrayList<Read> prots;
		if(sixframes) {
			prots=TranslateSixFrames.toFrames(r, true, false, 6);
		}else{
			ArrayList<Orf> list;
			list=gCaller.callGenes(r);
			prots=CallGenes.translate(r, list);
		}
		if(prots!=null){
			for(Read p : prots){
				processReadAmino(p);
			}
		}
	}
	
	void processReadNucleotide(final Read r){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		long kmer=0;
		long rkmer=0;
		int len=0;
		assert(!r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		if(eTracker!=null){eTracker.clear();}
		
//		assert(false) : minProb+", "+minQual+", "+(quals==null);
		
		if(quals==null || (minProb<=0 && minQual<2)){
//			System.err.println("A");
			for(int i=0; i<bases.length; i++){
//				System.err.println("B: len="+len);
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(eTracker!=null){eTracker.add(b);}
				if(x<0){len=0; rkmer=0;}else{len++;}
				
//				System.err.println("\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+"\n"
//				+AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(rkmer, k), k)+"\n"
//				+len+", "+(char)b+", "+x+", "+x2+"\n");
				
				if(len>=k){
					kmersProcessed++;
					heap.genomeSizeKmers++;
//					heap.probSum++; //Note really necessary for fasta data
					if(eTracker==null || eTracker.passes()){
//						System.err.println("Pass.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString());
						
//						assert(kmer==AminoAcid.reverseComplementBinaryFast(rkmer, k)) : 
//							"\n"+AminoAcid.kmerToString(kmer, k)+"\n"+AminoAcid.kmerToString(rkmer, k)+"\n"
//							+AminoAcid.kmerToString(AminoAcid.reverseComplementBinaryFast(rkmer, k), k)+"\n"
//							+len+", "+(char)b+", "+x+", "+x2;
						
						final long hashcode=hash(kmer, rkmer);
						//				System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
						if(hashcode>min){
							if(noBlacklist){
								heap.add(hashcode);
							}else{
								heap.checkAndAdd(hashcode);
							}
						}
					}else{
//						System.err.println("Fail.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString()+"\n"+r.toFastq()+"\n"+eTracker);
//						assert(false);
					}
				}
			}
		}else{
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];

				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);
				if(eTracker!=null){eTracker.add(b);}
				
				{//Quality-related stuff
					final byte q=quals[i];
					assert(q>=0) : Arrays.toString(quals)+"\n"+minProb+", "+minQual;
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
					if(x<0 || q<minQual){
						len=0;
						kmer=rkmer=0;
						prob=1;
					}else{
						len++;
					}
				}
				
				if(len>=k && prob>=minProb){
					kmersProcessed++;
					heap.genomeSizeKmers++;
					heap.probSum+=prob;
					if(eTracker==null || eTracker.passes()){
						final long hashcode=hash(kmer, rkmer);
						//				System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
						if(hashcode>min){
							if(noBlacklist){
								heap.add(hashcode);
							}else{
								heap.checkAndAdd(hashcode);
							}
						}
					}else{
//						System.err.println("Fail.\t"+eTracker.calcEntropy()+"\t"+eTracker.basesToString());
					}
				}
				
				//This version is slow but calculates depth better.
//				if(len>=k){
//					kmersProcessed++;
//					heap.genomeSizeKmers++;
//					final long hash=hash(kmer, rkmer);
//					//				System.err.println(kmer+"\t"+rkmer+"\t"+z+"\t"+hash);
//					if(hash>min){
//						if(prob>=minProb || (!heap.setMode && heap.contains(hash))){
//							if(noBlacklist){
//								heap.add(hash);
//							}else{
//								heap.checkAndAdd(hash);
//							}
//						}
//					}
//				}
			}
		}
//		assert(false);
	}

	void processReadAmino(final Read r){
		final byte[] bases=r.bases;
		long kmer=0;
		int len=0;
		assert(r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.acidToNumberNoStops[b];
			kmer=((kmer<<aminoShift)|x)&mask;
//			if(eTracker!=null){eTracker.add(b);}
			
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				heap.genomeSizeKmers++;
//				if(eTracker==null || eTracker.passes()){
//					assert(false) : (eTracker==null)+", "+eTracker.cutoff()+", "+eTracker.calcEntropy()+", "+r;
					long hashcode=hash(kmer, kmer);
					if(hashcode>min){
						if(noBlacklist){
							heap.add(hashcode);
						}else{
							heap.checkAndAdd(hashcode);
						}
					}
//				}
			}
		}
	}

	void processReadAmino_old_no_entropy(final Read r){
		final byte[] bases=r.bases;
		long kmer=0;
		int len=0;
		assert(r.aminoacid());
		
		final boolean noBlacklist=!(Blacklist.exists() || Whitelist.exists());
		final long min=minHashValue;
		heap.genomeSizeBases+=r.length();
		heap.genomeSequences++;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.acidToNumberNoStops[b];
			kmer=((kmer<<aminoShift)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				heap.genomeSizeKmers++;
				long hashcode=hash(kmer, kmer);
				if(hashcode>min){
					if(noBlacklist){
						heap.add(hashcode);
					}else{
						heap.checkAndAdd(hashcode);
					}
				}
			}
		}
	}
	
	public Sketch toSketch(){
		Sketch sketch=null;
		if(heap!=null && heap.size()>0){
			try {
				sketch=new Sketch(heap, false, tool.trackCounts, null);
			} catch (Throwable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			heap.clear(false);
		}
		return sketch;
	}
	
	public void add(SketchMakerMini smm){
		heap.add(smm.heap);
		readsProcessed+=smm.readsProcessed;
		basesProcessed+=smm.basesProcessed;
		kmersProcessed+=smm.kmersProcessed;
		sketchesMade+=smm.sketchesMade;
	}

	/** True only if this thread has completed successfully */
	boolean success=false;

	SketchHeap heap;

	final int aminoShift;
	final int shift;
	final int shift2;
	final long mask;
	final EntropyTracker eTracker;
	final GeneCaller gCaller;
	
	public float minEntropy() {
		// TODO Auto-generated method stub
		return eTracker==null ? -1 : eTracker.cutoff();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of bases processed */
	protected long kmersProcessed=0;
	/** Number of sketches started */
	protected long sketchesMade=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final SketchTool tool;
	final int mode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
