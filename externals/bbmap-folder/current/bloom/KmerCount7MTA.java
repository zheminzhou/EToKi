package bloom;

import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import kmer.KmerTableSet;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import sketch.SketchObject;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class KmerCount7MTA extends KmerCountAbstract {
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
//			outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>1 ? args[1] : null);
		int k=14;
		int cbits=16;
		int gap=0;
		int matrixbits=-1;
		int hashes=1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("gap")){
				gap=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Tools.parseKMG(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("canonical")){
				CANONICAL=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		int kbits=Tools.min(2*k, 62);
		if(matrixbits<0){
			matrixbits=kbits;
		}
		matrixbits=Tools.min(kbits, matrixbits);
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		
		KCountArray counts=KCountArray.makeNew(1L<<kbits, 1L<<matrixbits, cbits, gap, hashes);
		try {
			counts=count(fname1, fname2, k, cbits, gap, true, false, false, false, counts);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		counts.shutdown();
		
//		verbose=true;
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		printStatistics(counts);
		
	}

	public static void printStatistics(KCountArray counts){
		long[] freq=counts.transformToFrequency();

//		System.out.println(count+"\n");
//		System.out.println(Arrays.toString(freq)+"\n");
		
		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			int x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+String.format(Locale.ROOT, "%.3f    ",(keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+String.format(Locale.ROOT, "%.3f%%   ",(100l*x/(double)sum2))+"\t"+x);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles, int k, int cbits, 
			boolean rcomp, boolean ecco, boolean merge, boolean amino){
		return makeKca(fname1, fname2, extraFiles, k, cbits, 0, Tools.min(2*k, 35), 1, minQuality, 
				rcomp, ecco, merge, maxReads, 1, 1, 1, 2, amino);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, 
			boolean rcomp, boolean ecco, boolean merge, long maxreads, boolean amino){
		assert(matrixbits<63);
		return makeKca(fname1, fname2, extraFiles, k, cbits, gap, matrixbits, hashes, minqual, 
				rcomp, ecco, merge, maxreads, 1, 1, 1, 2, amino);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, int matrixbits, int hashes, int minqual, boolean rcomp, boolean ecco, boolean merge,
			long maxreads, int passes, int stepsize, int thresh1, int thresh2, boolean amino){
		assert(matrixbits<63);
		return makeKca(fname1, fname2, extraFiles,
				k, cbits, gap, 1L<<matrixbits, hashes, minqual, rcomp, ecco, merge, 
				maxreads, passes, stepsize, thresh1, thresh2, null, 0, amino);
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, long cells, int hashes, int minqual, boolean rcomp, boolean ecco, boolean merge,
			long maxreads, int passes, int stepsize, int thresh1, int thresh2, boolean amino){
		return makeKca(fname1, fname2, extraFiles,
				k, cbits, gap, cells, hashes, minqual, rcomp, ecco, merge, 
				maxreads, passes, stepsize, thresh1, thresh2, null, 0, amino);
	}
	
	public static KCountArray makeKca_als(ArrayList<String> fname1, ArrayList<String> fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, long cells, int hashes, int minqual, boolean rcomp, boolean ecco, boolean merge,
			long maxreads, int passes, int stepsize, int thresh1, int thresh2, KCountArray prefilter, int prefilterLimit_, boolean amino_){
		String a=null, b=null;
		ArrayList<String> list=new ArrayList<String>();
		if(fname1!=null){
			for(int i=0; i<fname1.size(); i++){
				if(i==0){a=fname1.get(i);}
				else{list.add(fname1.get(i));}
			}
		}
		if(fname2!=null){
			for(int i=0; i<fname2.size(); i++){
				if(i==0){b=fname2.get(i);}
				else{list.add(fname2.get(i));}
			}
		}
		if(extraFiles!=null){
			for(String s : extraFiles){
				list.add(s);
			}
		}
		return makeKca(a, b, list.isEmpty() ? null : list, k, cbits, gap, cells, hashes, minqual, 
				rcomp, ecco, merge, maxreads, passes, stepsize, thresh1, thresh2,
				prefilter, prefilterLimit_, amino_);
	}
	
	public static KCountArray makeKcaFromIndex(int k, int cbits, long cells, int hashes, boolean rcomp){

		final int kbits=Tools.min(2*k, 62);
		KCountArray kca=KCountArray.makeNew(1L<<kbits, cells, cbits, 0, hashes, null, 0);
		try {
			countFromIndex(k, cbits, rcomp, kca);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		kca.shutdown();
		return kca;
	}
	
	public static KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int k, int cbits, int gap, long cells, int hashes, int minqual, boolean rcomp,
			boolean ecco, boolean merge, long maxreads, int passes, int stepsize, int thresh1, int thresh2,
			KCountArray prefilter, int prefilterLimit_, final boolean amino){
		assert(!amino || !rcomp);
		final int kbits=Tools.min(2*k, 62);
//		verbose=true;
		if(verbose){System.err.println("Making kca from ("+fname1+", "+fname2+")\nk="+k+", gap="+gap+", cells="+Tools.toKMG(cells)+", cbits="+cbits);}
		
		if(fname1==null && fname2==null && extraFiles==null){
			return KCountArray.makeNew(1L<<kbits, cells, cbits, gap, hashes, prefilter, prefilterLimit_);
		}
		
		boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		long oldmax=maxReads;
		byte oldq=minQuality;
		maxReads=maxreads;
		minQuality=(byte)minqual;
		//	System.out.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<kbits, cells, cbits, gap, hashes, prefilter, prefilterLimit_);
		
//		System.out.println("a");
		{//For processing input lists
			ArrayList<String> extra2=null;
			if(fname1!=null && fname1.contains(",")){
				String[] s=fname1.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname1=s[0];
			}
			if(fname2!=null && fname2.contains(",")){
				String[] s=fname2.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname2=s[0];
			}
			if(extra2!=null){
				if(extraFiles!=null){
					for(String s : extraFiles){
						extra2.add(s);
					}
				}
				extraFiles=extra2;
			}
		}
//		System.out.println("b");
		
		if(extraFiles!=null){
			for(String s : extraFiles){
				if(fileIO.FileFormat.hasFastaExtension(s)){
					assert(!FastaReadInputStream.SPLIT_READS);
				}
			}
		}
		
//		System.out.println("c");
		
		if(passes==1){
//			System.out.println("c1");
			if(fname1!=null){
				try {
					count(fname1, fname2, k, cbits, gap, rcomp, ecco, merge, amino, kca);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						count(s, null, k, cbits, gap, rcomp, ecco, merge, amino, kca);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();

		}else{
//			System.out.println("c2");
			assert(passes>1);
			KCountArray trusted=null;
			for(int i=1; i<passes; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
				//			if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}

				try {
					count(fname1, fname2, k, cbits, rcomp, ecco, merge, kca, trusted, maxreads, thresh1, step, conservative, amino);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						try {
							count(s, null, k, cbits, rcomp, ecco, merge, kca, trusted, maxreads, thresh1, step, conservative, amino);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
				kca.shutdown();
				
				System.out.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(1L<<kbits, cells, cbits, gap, hashes, prefilter, prefilterLimit_);

			}

			try {
				count(fname1, fname2, k, cbits, rcomp, ecco, merge, kca, trusted, maxreads, thresh2, stepsize, true, amino);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					try {
						count(s, null, k, cbits, rcomp, ecco, merge, kca, trusted, maxreads, thresh2, stepsize, true, amino);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			kca.shutdown();
		}
//		System.out.println("d");
		minQuality=oldq;
		maxReads=oldmax;
		FastaReadInputStream.SPLIT_READS=oldsplit;
		
		
		return kca;
	}
	
	public static KCountArray count(String reads1, String reads2, int k, int cbits, int gap,
			boolean rcomp, boolean ecco, boolean merge, boolean amino, KCountArray counts) throws Exception{
		assert(k>=1 && (counts!=null || k<20));
		final int kbits=Tools.min(2*k, 62);
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
//		System.err.println("countFastq...  making a new cris");
		if(counts==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			counts=KCountArray.makeNew(cells, cbits, gap);
		}
		assert(gap==counts.gap);
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			ff1.preferShreds=true;
//			if(ff2!=null){ //TODO - interleaved flag
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
//		countFastq(cris, k, rcomp, count);
//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, k, rcomp, ecco, merge, amino, counts);
			cta[i].start();
		}
//		System.out.println("~1");
		for(int i=0; i<cta.length; i++){
//			System.out.println("~2");
			CountThread ct=cta[i];
			synchronized(ct){
//				System.out.println("~3");
				while(ct.getState()!=State.TERMINATED){
//					System.out.println("~4");
					try {
						ct.join(2000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
//					System.out.println("~5");
				}
			}
		}
//		System.out.println("~6");
		
		ReadWrite.closeStream(cris);
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		
		return counts;
	}
	
	public static KCountArray countFromIndex(int k, int cbits, boolean rcomp, KCountArray counts) throws Exception{
		assert(k>=1 && (counts!=null || k<20));
		final int kbits=Tools.min(2*k, 62);
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
//		System.err.println("countFastq...  making a new cris");
		if(counts==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			counts=KCountArray.makeNew(cells, cbits, 0);
		}
		
//		countFastq(cris, k, rcomp, count);
//		assert(false) : THREADS;
		final CountThread2[] cta=new CountThread2[Tools.min(Data.numChroms*THREADS_PER_CHROM, THREADS)];
		final AtomicInteger nextChrom=new AtomicInteger(0);
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread2(k, rcomp, false, counts, nextChrom);
			cta[i].start();
		}
//		System.out.println("~1");
		for(int i=0; i<cta.length; i++){
//			System.out.println("~2");
			CountThread2 ct=cta[i];
			synchronized(ct){
//				System.out.println("~3");
				while(ct.getState()!=State.TERMINATED){
//					System.out.println("~4");
					try {
						ct.join(2000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
//					System.out.println("~5");
				}
			}
		}
		
		return counts;
	}

	
	public static KCountArray count(final String reads1, final String reads2, final int k, final int cbits, 
			final boolean rcomp, final boolean ecco, boolean merge,
			KCountArray counts, final KCountArray trusted, final long maxReads, final int thresh, final int detectStepsize,
			final boolean conservative, final boolean amino)
			throws Exception{
		
		assert(!amino || !rcomp);
		assert(k>=1 && (counts!=null || k<20));
		final int kbits=Tools.min(2*k, 62);
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
//		System.out.println("k="+k+", kbits="+kbits+", mask="+Long.toHexString(mask)+", thresh="+thresh);
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
//		verbose=true;
		
		if(counts==null){
			final long cells=1L<<kbits;
			if(verbose){System.err.println("k="+k+", kbits="+kbits+", cells="+cells+", mask="+Long.toHexString(mask));}
			counts=KCountArray.makeNew(cells, cbits, 0);
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		

//		countFastq(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);

//		assert(false) : THREADS;
		CountThread[] cta=new CountThread[THREADS];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, k, rcomp, ecco, merge, amino, 
					counts, trusted, thresh, detectStepsize, conservative);
			cta[i].start();
		}
		
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		cris.close();
		if(verbose){System.err.println("Closed stream");}
		
//		System.out.println("*** after ***");
//		System.out.println("\ntrusted=\n"+trusted);
//		System.out.println("\ncount=\n"+count);
		
		return counts;
	}
	
	private static final int findOverlap(Read r1, Read r2, boolean ecc){
		if(vstrict){
			return BBMerge.findOverlapVStrict(r1, r2, ecc);
		}else{
			return BBMerge.findOverlapStrict(r1, r2, ecc);
		}
	}
	
	public static boolean vstrict=false;
	
	private static class CountThread extends Thread{
		
		CountThread(final ConcurrentReadInputStream cris_, final int k_, final boolean rcomp_,
				final boolean ecco_, boolean merge_, final boolean amino_,
				final KCountArray counts_){
			this(cris_, k_, rcomp_, ecco_, merge_, amino_, counts_, null, 2, 1, true);
		}
		
		CountThread(final ConcurrentReadInputStream cris_, final int k_, final boolean rcomp_, 
				final boolean ecco_, boolean merge_, final boolean amino_,
				final KCountArray counts_, final KCountArray trusted_, final int thresh_,
				final int detectStepsize_, final boolean conservative_){
			cris=cris_;
			k=k_;
			rcomp=rcomp_;
			ecco=ecco_;
			merge=merge_;
			counts=counts_;
			trusted=trusted_;
			thresh=thresh_;
			detectStepsize=detectStepsize_;
			conservative=conservative_;
			amino=amino_;
			assert(!amino || !rcomp);
			assert(!amino || !ecco);
		}
		
		@Override
		public void run(){
//			System.out.println("Running");
			if(trusted==null){
				count(cris, rcomp, counts);
			}else{
				count(cris, rcomp, counts, trusted, thresh, detectStepsize,  conservative);
			}
//			System.out.println("Finished: "+readsProcessedLocal);
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}
		
		private int increment(final long key, final KCountArray counts){
			if(SKETCH_MODE){
				final long code=SketchObject.hash(key);
				if(code>minHashValue){
					counts.increment(STORE_HASHED ? code : key);
					return 1;
				}else{
					return 0;
				}
			}else{
				counts.increment(key);
				return 1;
			}
		}
		
		private final void count(ConcurrentReadInputStream cris, boolean rcomp, KCountArray counts){
			assert(k>=1 && counts!=null);
//			System.out.println("Waiting for list");
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
//			System.out.println("Got list: "+(ln==null ? "null" : ln.id)+", "+(ln==null || ln.list==null ? "null" : ln.list.size()));

			long[] array=null;
			final Kmer kmer=new Kmer(k);
			if(counts.gap==0){
				final int kbits=Tools.min(2*k, 62);
				final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
				
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					//System.err.println("reads.size()="+reads.size());
					for(Read r1 : reads){
						
						Read r2=r1.mate;
						if((ecco || merge) && r1!=null && r2!=null){
							if(merge){
								final int insert=findOverlap(r1, r2, false);
								if(insert>0){
									r2.reverseComplement();
									r1=r1.joinRead(insert);
									r2=null;
								}
							}else if(ecco){
								findOverlap(r1, r2, true);
							}
						}
						readsProcessedLocal++;
						
						if(k<=maxShortKmerLength){
							array=addRead_Advanced(r1, counts, mask, array);
						}else{
							addReadBig(r1, kmer);
							addReadBig(r1.mate, kmer);
						}
//						System.out.println(r);
//						System.out.println("kmers hashed: "+keysCountedLocal);
					}
					//System.err.println("returning list");
					cris.returnList(ln);
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}else{
				final int k1=(k+1)/2;
				final int k2=k/2;
				final int kbits1=2*k1;
				final int kbits2=2*k2;
				final int gap=counts.gap;
				final long mask1=~((-1L)<<(kbits1));
				final long mask2=~((-1L)<<(kbits2));
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					//System.err.println("reads.size()="+reads.size());
					for(Read r1 : reads){
						Read r2=r1.mate;
						if((ecco || merge) && r1!=null && r2!=null){
							if(merge){
								final int insert=findOverlap(r1, r2, false);
								if(insert>0){
									r2.reverseComplement();
									r1=r1.joinRead(insert);
									r2=null;
								}
							}else if(ecco){
								findOverlap(r1, r2, true);
							}
						}
						readsProcessedLocal++;
						addReadSplit(r1, counts, k1, k2, mask1, mask2, gap, rcomp);
						if(r1.mate!=null){
							addReadSplit(r1.mate, counts, k1, k2, mask1, mask2, gap, rcomp);
						}
					}
					//System.err.println("returning list");
					cris.returnList(ln);
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		


		private final void count(final ConcurrentReadInputStream cris, final boolean rcomp,
				final KCountArray count, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			if(count.gap>0){
				countFastqSplit(cris, k, rcomp, count, trusted, thresh, detectStepsize, conservative);
				return;
			}
			assert(k>=1 && (count!=null || k<20));
			final int kbits=Tools.min(2*k, 62);
			final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			long[] array=null;
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r1 : reads){
					
					Read r2=r1.mate;
					if((ecco || merge) && r1!=null && r2!=null){
						if(merge){
							final int insert=findOverlap(r1, r2, false);
							if(insert>0){
								r2.reverseComplement();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							findOverlap(r1, r2, true);
						}
					}
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r1, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r1, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.length()));
//							System.out.println(new String(r.bases));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r1.length(); i=bs.nextClearBit(i+1)){
									r1.bases[i]='N';
									if(r1.quality!=null){r1.quality[i]=0;}
								}
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format(Locale.ROOT, "%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.length()){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							if(bs!=null){
								for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
									r2.bases[i]='N';
									if(r2.quality!=null){r2.quality[i]=0;}
								}
							}
						}
					}
					array=addRead_Advanced(r1, count, mask, array);

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		
		private void countFastqSplit(final ConcurrentReadInputStream cris, final int k, final boolean rcomp,
				final KCountArray counts, final KCountArray trusted, final int thresh, final int detectStepsize, final boolean conservative){
			assert(false) : cris.paired();
			assert(counts.gap>0);
			assert(k<=maxShortKmerLength && k>=1 && (counts!=null || k<20));
			final int kbits=2*k;
			final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
			

			final int k1=(k+1)/2;
			final int k2=k/2;
			final int kbits1=2*k1;
			final int kbits2=2*k2;
			final int gap=counts.gap;
			final long mask1=~((-1L)<<(kbits1));
			final long mask2=~((-1L)<<(kbits2));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r1 : reads){
					
					Read r2=r1.mate;
					if((ecco || merge) && r1!=null && r2!=null){
						if(merge){
							final int insert=findOverlap(r1, r2, false);
							if(insert>0){
								r2.reverseComplement();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							findOverlap(r1, r2, true);
						}
					}
					{
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r1, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r1, trusted, k, thresh, detectStepsize));
//							System.out.println("\n"+toString(bs, r.length()));
//							System.out.println(new String(r.bases));
							for(int i=bs.nextClearBit(0); i<r1.length(); i=bs.nextClearBit(i+1)){
								r1.bases[i]='N';
								r1.quality[i]=0;
							}
//							System.out.println(new String(r.bases));
//							System.out.println("used = "+String.format(Locale.ROOT, "%.3f%%",count.usedFraction()*100));
//							System.out.println("used = "+((KCountArray4)count).cellsUsed());
//							if(bs.length()<r.length()){r=null;}
						}
//						if(r!=null){addRead(r, count, k, mask, rcomp);}
						
						addReadSplit(r1, counts, k1, k2, mask1, mask2, gap, rcomp);
					}
					if(r2!=null){
						if(trusted!=null){
							BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r2, trusted, k, thresh, detectStepsize) :
								ErrorCorrect.detectTrusted(r2, trusted, k, thresh, detectStepsize));
							for(int i=bs.nextClearBit(0); i<r2.length(); i=bs.nextClearBit(i+1)){
								r2.bases[i]='N';
								r2.quality[i]=0;
							}
						}
						addReadSplit(r2, counts, k1, k2, mask1, mask2, gap, rcomp);
					}

				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		
		/**
		 * Hash a read's kmers into the KCountArray.
		 * Advanced mode processes paired reads together and sorts kmers to eliminate spurious duplicates.
		 * @param r1
		 * @param counts
		 * @param k
		 * @param mask
		 * @param rcomp
		 */
		private final long[] addRead_Advanced(Read r1, final KCountArray counts, final long mask, long[] array){
			assert(counts.gap==0) : "Gapped: TODO";
			if(PREJOIN && r1.mate!=null && r1.insert()>0){
				r1.mate.reverseComplement();
				r1=r1.joinRead();
			}
			Read r2=r1.mate;
			final int len1=Tools.max(0, r1.length()-k+1);
			final int len2=(r2==null || r2.bases==null) ? 0 : Tools.max(0, r2.length()-k+1);
			final int len=len1+len2;
			if(len<1){return array;}
			
			if(!KEEP_DUPLICATE_KMERS){
				if(array==null || array.length!=len){array=new long[len];}
				Arrays.fill(array, -1);
				fillKmerArray(r1, mask, array, 0, len1);
				if(r2!=null){fillKmerArray(r2, mask, array, len1, len);}
				
				Arrays.sort(array);
				long prev=-1;
				for(int i=0; i<array.length; i++){
					long kmer=array[i];
					if(kmer!=prev){
						keysCountedLocal+=increment(kmer, counts);
						prev=kmer;
					}
				}
			}else{
				if(len1>0){addRead(r1);}
				if(len2>0){addRead(r2);}
			}
			return array;
		}
		
		private final void addReadBig(Read r, Kmer kmer){
			if(r==null || r.bases==null){return;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			int len=0;
			
			if(bases==null || bases.length<k){return;}
			kmer.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];

				//Update kmers
				kmer.addRight(b);
				
				if(minProb>0 && quals!=null){//Update probability
					prob=prob*KmerTableSet.PROB_CORRECT[quals[i]];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*KmerTableSet.PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					prob=1;
				}else{len++;}
				
				assert(len==kmer.len);
				
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
				if(len>=k && prob>=minProb){
//					System.err.println("Incrementing xor()="+kmer.xor());
					keysCountedLocal+=increment(kmer.xor(), counts);
//					counts.incrementAndReturnUnincremented(kmer.xor(), 1);
//					keysCountedLocal++;
				}
			}
		}
		
		private final void fillKmerArray(Read r, final long mask, final long[] array, final int start, final int stop){
			if(amino){
				fillKmerArrayAmino(r, k, array, start, stop);
				return;
			}
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, k, array, start, stop);
				return;
			}
			assert(counts.gap==0);
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k+counts.gap){return;}
			
			final int passes=(rcomp ? 2 : 1);
			for(int pass=0; pass<passes; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					assert(b>=0) : r.id+", "+bases.length+", "+i+", "+b+"\n"+(bases.length<20000 ? r.toFasta().toString() : "");
					int x=AminoAcid.baseToNumber[b];
					
//					int x=AminoAcid.baseToNumber[b<0 ? 'N' : b];

					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && prob>=minProb){
							array[idx]=Tools.max(array[idx], kmer);
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void fillKmerArrayAmino(Read r, final long mask, final long[] array, final int start, final int stop){
			assert(false) : "TODO"; //TODO
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, k, array, start, stop);
				return;
			}
			assert(counts.gap==0);
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k+counts.gap){return;}
			
			for(int pass=0; pass<2; pass++){
				int len=0;
				int idx=(pass==0 ? start-k+1 : stop+k-2);
				long kmer=0;
				float prob=1;
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					assert(b>=0) : r.id+", "+bases.length+", "+i+", "+b+"\n"+(bases.length<20000 ? r.toFasta().toString() : "");
					int x=AminoAcid.baseToNumber[b];
					
//					int x=AminoAcid.baseToNumber[b<0 ? 'N' : b];

					byte q;
					if(quals==null){
						q=50;
					}else{
						q=quals[i];
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
					}

					if(x<0 || q<minQuality){
						len=0;
						kmer=0;
						prob=1;
					}else{
						kmer=((kmer<<2)|x)&mask;
						len++;
						if(len>=k && prob>=minProb){
							array[idx]=Tools.max(array[idx], kmer);
						}
					}
					if(pass==0){idx++;}else{idx--;}
				}
//				System.out.println(Arrays.toString(array));
				r.reverseComplement();
			}
		}
		
		private final void addRead(Read r){
			if(amino){
				addReadAmino(r);
				return;
			}
			assert(counts.gap==0);
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);

			final byte[] bases=r.bases;
			final byte[] quals=r.quality;

			if(bases==null || bases.length<k+counts.gap){return;}
			
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
			long kmer=0;
			long rkmer=0;
			int len=0;
			float prob=1;

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);

				final byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}

				if(x<0 || q<minQuality){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{
					len++;
					if(len>=k && prob>=minProb){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						keysCountedLocal+=increment(key, counts);
					}
				}
			}
		}
		
		private final void addReadAmino(Read r){
			final int aminoShift=AminoAcid.AMINO_SHIFT;
			assert(counts.gap==0);
			assert(k*aminoShift<64);
			assert(!PREJOIN || r.mate==null);

			final byte[] bases=r.bases;

			if(bases==null || bases.length<k+counts.gap){return;}

			final int shift=aminoShift*k;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			int len=0;

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.acidToNumber[b];
				kmer=((kmer<<aminoShift)|x)&mask;

				if(x<0){
					len=0;
					kmer=0;
				}else{
					len++;
					if(len>=k){
						keysCountedLocal+=increment(kmer, counts);
					}
				}
			}
		}
		
		private final void fillKmerArrayLong(Read r, final int k, final long[] array, final int start, final int stop){
			assert(k>maxShortKmerLength) : k;
			assert(counts.gap==0);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(array!=null);
			Kmer kmer=new Kmer(k);
			
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			kmer.clear();
			
			for(int i=0, idx=start-k+1; i<bases.length; i++, idx++){
				byte b=bases[i];
				kmer.addRight(b);
				
				byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(kmer.len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}
				
				if(!AminoAcid.isFullyDefined(b) || q<minQuality){
					kmer.clear();
					prob=1;
				}
				if(kmer.len>=k && prob>=minProb){
					array[idx]=kmer.xor();
				}
			}
		}
		
		private final void addReadSplit(Read r, final KCountArray counts, final int k1, final int k2, final long mask1, final long mask2, final int gap, boolean rcomp){
			assert(false) : "TODO";
			if(PREJOIN && r.mate!=null && r.insert()>0){
				if(verbose){System.err.println("Prejoining "+r.numericID+" at "+r.insert());}
				r.mate.reverseComplement();
				r=r.joinRead();
			}
			
			int len=0;
			int shift=k2*2;
			long kmer1=0;
			long kmer2=0;
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			assert(kmer1>=kmer2);
			
//			assert(false) : k1+", "+k2+", "+mask1+", "+mask2+", "+gap;
			
			if(verbose){System.err.println("Hashing read "+r.numericID+"; loop limits "+(k1+gap)+"-"+(bases.length));}
			for(int i=0, j=i+k1+gap; j<bases.length; i++, j++){
				int x1=AminoAcid.baseToNumber[bases[i]];
				int x2=AminoAcid.baseToNumber[bases[j]];
				
				byte q1, q2;
				if(quals==null){
					q1=50;
					q2=50;
				}else{
					q1=quals[i];
					q2=quals[j];
					prob=prob*align2.QualityTools.PROB_CORRECT[q1]*align2.QualityTools.PROB_CORRECT[q2];
					if(len>k){
						byte oldq1=quals[i-k1];
						byte oldq2=quals[j-k2];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq1]*align2.QualityTools.PROB_CORRECT_INVERSE[oldq2];
					}
				}
				
				if(x1<0 || x2<0 || q1<minQuality || q2<minQuality){
					len=0;
					kmer1=0;
					kmer2=0;
					prob=1;
				}else{
					kmer1=((kmer1<<2)|x1)&mask1;
					kmer2=((kmer2<<2)|x2)&mask2;
					len++;
					if(len>=k1 && prob>=minProb){
						
//						System.out.print("Incrementing "+Long.toHexString(kmer)+": "+count.read(kmer));
						
						long key=(kmer1<<shift)|kmer2;
//						System.err.println(Long.toHexString(key));
						
						if(verbose){System.err.println("Hashing key "+Long.toHexString(key)+" at length "+len);}
						keysCountedLocal+=increment(key, counts);
//						count.inchrement(kmer);
						
						
//						System.out.println(" -> "+count.read(kmer));
//						System.out.print("Incrementing array for "+Long.toHexString(kmer)+": "+array[(int)kmer]);
//						array[(int)kmer]++;
//						System.out.println(" -> "+array[(int)kmer]+"\n");
//						assert(array[(int)kmer]==count.read(kmer) || array[(int)kmer]>3);
					}
				}
			}
			if(rcomp){
				r.reverseComplement();
				addReadSplit(r, counts, k1, k2, mask1, mask2, gap, false);
			}
		}

		private final ConcurrentReadInputStream cris;
		private final int k;
		private final boolean rcomp;
		private final boolean ecco;
		private final boolean merge;
		private final KCountArray counts;
		private final KCountArray trusted;
		private final int thresh;
		private final int detectStepsize;
		private final boolean conservative;
		private final boolean amino;
		private long keysCountedLocal=0;
		private long readsProcessedLocal=0;
		private final long minHashValue=SketchObject.minHashValue;
	}
	
	private static class CountThread2 extends Thread{
		
		CountThread2(final int k_, final boolean rcomp_, final boolean amino_, final KCountArray counts_, AtomicInteger nextChrom_){
			k=k_;
			rcomp=rcomp_;
			counts=counts_;
			amino=amino_;
			nextChrom=nextChrom_;
			assert(!amino || !rcomp);
		}
		
		@Override
		public void run(){
			count(counts);
			
			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}
		
		private final void count(KCountArray counts){
			assert(k>=1 && counts!=null);
			final int maxCount=THREADS_PER_CHROM*Data.numChroms;
			for(int cnum=nextChrom.getAndIncrement(); cnum<maxCount; cnum=nextChrom.getAndIncrement()){
				ChromosomeArray ca=Data.getChromosome(cnum/THREADS_PER_CHROM+1);
				processChrom(ca, cnum%THREADS_PER_CHROM);
			}
		}
		
		private final void processChrom(ChromosomeArray ca, int segNum){
			assert(counts.gap==0);
			assert(k<=maxShortKmerLength);
			assert(CANONICAL);

			final byte[] bases=ca.array;
			if(bases==null || bases.length<k+counts.gap){return;}
			final int segLength=bases.length/4;
			final int start=Tools.max(0, segNum*segLength-k);
			final int stop=Tools.min(bases.length, (segNum+1)*segLength);

			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			
			long kmer=0;
			long rkmer=0;
			int len=0;

			for(int i=start; i<stop; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=(rkmer>>>2)|(x2<<shift2);

				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						counts.increment(key);
						readsProcessedLocal++;
					}
				}
			}
		}
		
		private final int k;
		private final boolean rcomp;
		private final KCountArray counts;
		private final boolean amino;
		private final AtomicInteger nextChrom;
		private long keysCountedLocal=0;
		private long readsProcessedLocal=0;
		private final long minHashValue=SketchObject.minHashValue;
	}
	
	public static int maxShortKmerLength=31;
	private static final int THREADS_PER_CHROM=4;
	
}
