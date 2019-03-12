package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.IntHashSet;
import structures.IntList;
import structures.ListNum;

/**
 * This class is designed to store kmer frequencies related to gene
 * starts, stops, and interiors.  It can be loaded from a pgm file.
 * 
 * It's possible to use multiple GeneModels; for example, one for
 * each of several GC ranges or clades.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class GeneModel {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public GeneModel(boolean fill){
		if(fill){
			fillContainers();
			fillContainerArrays();
		}
	}

	void fillContainers(){
		statsCDS.setInner(kInnerCDS, 3);
		statsCDS.setStart(kStartCDS, startFrames, startLeftOffset);
		statsCDS.setStop(kStopCDS, stopFrames, stopLeftOffset);

		for(int i=0; i<rnaContainers.length; i++){
			StatsContainer sc=rnaContainers[i];
			sc.setInner(kInnerRNA, 1);
		}

		statstRNA.setStart(kStartRNA, 14, 4);
		statstRNA.setStop(kStopRNA, 14, 6);

		stats16S.setStart(kStartRNA, 20, 7);
		stats16S.setStop(kStopRNA, 12, 16);

		stats23S.setStart(kStartRNA, 17, 3);
		stats23S.setStop(kStopRNA, 15, 12);

		stats5S.setStart(kStartRNA, 20, 5);
		stats5S.setStop(kStopRNA, 15, 5);
	}

	void fillContainerArrays(){
		rnaContainers=new StatsContainer[] {statstRNA, stats16S, stats23S, stats5S};
		allContainers=new StatsContainer[] {statsCDS, statstRNA, stats16S, stats23S, stats5S};
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean process(String genomeFname, String gffFname){
		fnames.add(ReadWrite.stripPath(genomeFname));
		FileFormat fnaFile=FileFormat.testInput(genomeFname, FileFormat.FA, null, true, true);
		FileFormat gffFile=FileFormat.testInput(gffFname, FileFormat.GFF, null, true, true);
		
		if(fnaFile==null || gffFile==null){
			errorState=true;
			return true;
		}
		filesProcessed++;
		
		ArrayList<ScafData> scafList;
		{//Scoped to save memory
			ArrayList<Read> reads=ReadInputStream.toReads(fnaFile, maxReads);
			readsProcessed+=reads.size();
			scafList=new ArrayList<ScafData>(reads.size());
			for(Read r : reads){
				basesProcessed+=r.length();
				scafList.add(new ScafData(r));
			}
		}
		{//Scoped to save memory
			ArrayList<GffLine>[] allGffLines=GffLine.loadGffFileByType(gffFile, "CDS,rRNA,tRNA");
			ArrayList<GffLine> cds=allGffLines[0];
			ArrayList<GffLine> rrna=allGffLines[1];
			ArrayList<GffLine> trna=allGffLines[2];
			genesProcessed+=cds.size();
			HashMap<String, ScafData> scafMap=makeScafMap(scafList);
			fillScafDataCDS(cds, scafMap);
			fillScafDataRNA(rrna, scafMap);
			fillScafDataRNA(trna, scafMap);
		}
		
		countBases(scafList);
		if(PROCESS_PLUS_STRAND){
			processStrand(scafList, Shared.PLUS);
		}
		if(PROCESS_MINUS_STRAND){
			for(ScafData sd : scafList){
				sd.clear();
				sd.reverseComplement();
			}
			processStrand(scafList, Shared.MINUS);
			for(ScafData sd : scafList){
				sd.clear();
				sd.reverseComplement();
			}
		}
		return false;
	}
	
	public void add(GeneModel pgm){
		for(int i=0; i<allContainers.length; i++){
			allContainers[i].add(pgm.allContainers[i]);
		}
		
		readsProcessed+=pgm.readsProcessed;
		basesProcessed+=pgm.basesProcessed;
		genesProcessed+=pgm.genesProcessed;
		filesProcessed+=pgm.filesProcessed;
		
		fnames.addAll(pgm.fnames);
		taxIds.addAll(pgm.taxIds);
		Tools.add(baseCounts, pgm.baseCounts);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public float gc(){
		long a=baseCounts[0];
		long c=baseCounts[1];
		long g=baseCounts[2];
		long t=baseCounts[3];
		return (float)((g+c)/Tools.max(1.0, a+t+g+c));
	}
	
	HashMap<String, ScafData> makeScafMap(ArrayList<ScafData> scafList){
		HashMap<String, ScafData> scafMap=new HashMap<String, ScafData>(scafList.size()*3);
		for(ScafData sd : scafList){scafMap.put(sd.name, sd);}
		for(ScafData sd : scafList){
			String name=sd.name;
			int idx=name.indexOf(' ');
			if(idx>=0){
				String prefix=name.substring(0, idx);
				if(scafMap.containsKey(prefix)){
					assert(false) : "Duplicate degenerate name: '"+name+"', '"+prefix+"'";
				}else{
					scafMap.put(prefix, sd);
				}
			}
		}
		return scafMap;
	}
	
	public void fillScafDataCDS(ArrayList<GffLine> cdsLines, HashMap<String, ScafData> scafMap){
		for(GffLine gline : cdsLines){
			ScafData sd=scafMap.get(gline.seqid);
			assert(sd!=null) : "Can't find scaffold for GffLine "+gline.seqid;
			sd.addCDS(gline);
		}
	}
	
	public void fillScafDataRNA(ArrayList<GffLine> rnaLines, HashMap<String, ScafData> scafMap){
		for(GffLine gline : rnaLines){
			ScafData sd=scafMap.get(gline.seqid);
			assert(sd!=null) : "Can't find scaffold for GffLine "+gline.seqid;
			sd.addRNA(gline);
		}
	}
	
	public void processStrand(ArrayList<ScafData> scafList, int strand){
		for(ScafData sd : scafList){
			processCDS(sd, strand);
			processRNA(sd, strand);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	private void countBases(ArrayList<ScafData> scafList){
		for(ScafData sd : scafList){
			countBases(sd.bases);
		}
	}
	
	private void countBases(byte[] bases){
		for(byte b : bases){
			int x=AminoAcid.baseToNumberACGTother[b];
			baseCounts[x]++;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Finding Codons        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void findStopCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i;//End of the stop codon
					if(isStopCodon(kmer) && !valid.get(point)){
						list.add(point);
						valid.set(point);
					}
				}
			}else{len=0;}
		}
		
		for(int i=50; i<bases.length-3; i+=2000){//Add some non-canonical sites, aka noise
			if(!valid.get(i)){
				list.add(i);
			}
		}
	}
	
	private static void findStartCodons(byte[] bases, IntList list, BitSet valid){
		final int k=3;
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int point=i-k+1;//Start of the start codon
					if(isStartCodon(kmer) && !valid.get(point)){
						list.add(point);
						valid.set(point);
					}
				}
			}else{len=0;}
		}
		
		for(int i=50; i<bases.length-3; i+=2000){//Add some non-canonical sites, aka noise
			if(!valid.get(i)){
				list.add(i);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Processing GffLines      ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void processGene(GffLine gline, ScafData sd){
		final int strand=gline.strand;
		assert(strand==sd.strand());
		final byte[] frames=sd.frames;
		int start=gline.start-1, stop=gline.stop-1;
		if(start<0 || stop>=sd.length()){return;}
		assert(start<stop);
		if(strand==Shared.MINUS){
			int x=sd.length()-start-1;
			int y=sd.length()-stop-1;
			start=y;
			stop=x;

//			String a=new String(sd.bases, start, 3);
//			String b=new String(sd.bases, stop-2, 3);
////			assert(false) : start+", "+stop+"\n"+gline+"\n"+new String(sd.bases, start, 3)+", "+new String(sd.bases, stop-2, 3);
//			outstream.println(a+", "+b+", "+start+", "+stop);
		}
		assert(start>=0) : gline.toString()+"\n"+sd.length()+"\n"+sd.name;
		markFrames(start, stop, frames, kInnerCDS);
		sd.starts.add(start);
		sd.stops.add(stop);
//		assert(gline.start!=337) : gline+"\n"+start+", "+stop;
	}
	
	private void processRnaLine(GffLine gline, ScafData sd, byte flag){
		final int strand=gline.strand;
		assert(strand==sd.strand());
		final byte[] frames=sd.frames;
		int start=gline.start-1, stop=gline.stop-1;
		if(start<0 || stop>=sd.length()){return;}
		assert(start<stop);
		if(strand==Shared.MINUS){
			int x=sd.length()-start-1;
			int y=sd.length()-stop-1;
			start=y;
			stop=x;
		}
		
		if(flag==1){
			//TODO
			statstRNA.start.processPoint(sd.bases, start, 1);
			statstRNA.stop.processPoint(sd.bases, stop, 1);
		}else if(flag==2){
			//TODO
			stats16S.start.processPoint(sd.bases, start, 1);
			stats16S.stop.processPoint(sd.bases, stop, 1);
		}else if(flag==4){
			//TODO
			stats23S.start.processPoint(sd.bases, start, 1);
			stats23S.stop.processPoint(sd.bases, stop, 1);
		}else if(flag==8){
			//TODO
			stats5S.start.processPoint(sd.bases, start, 1);
			stats5S.stop.processPoint(sd.bases, stop, 1);
		}
		
		assert(start>=0) : gline.toString()+"\n"+sd.length()+"\n"+sd.name;
		for(int i=start+kInnerRNA-1; i<stop; i++){
			frames[i]|=flag;
		}
	}
	
	/** 
	 * Each frame byte has a bit marked for valid coding frames.
	 * For example, if frames[23]=0b100, then base 23 is the last base in a kmer starting at the 3rd base in a codon.
	 * If frames[23]=0, then no coding kmer end at that location on this strand.
	 * @param start
	 * @param stop
	 * @param frames
	 * @param k
	 */
	private static void markFrames(int start, int stop, byte[] frames, int k){
		assert(start<stop) : start+", "+stop;
		for(int i=start+k-1, frameBit=(1<<((k-1)%3)), max=Tools.min(stop-3, frames.length-1); i<=max; i++){
			frames[i]=(byte)(frames[i]|frameBit);
			frameBit<<=1;
			if(frameBit>4){frameBit=1;}
		}
//		assert(false) : Arrays.toString(Arrays.copyOfRange(frames, start, start+20))+"\n"+start; //This is correct
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Counting Kmers        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processCDS(ScafData sd, int strand){
		ArrayList<GffLine> glines=sd.cdsLines[strand];
		for(GffLine gline : glines){
			assert(gline.strand==strand);
			processGene(gline, sd);
			statsCDS.addLength(gline.length());
		}
		
		statsCDS.inner.processCDSFrames(sd.bases, sd.frames);
		BitSet startSet=processEnds(sd.bases, statsCDS.start, sd.starts, 1);
		BitSet stopSet=processEnds(sd.bases, statsCDS.stop, sd.stops, 1);
//		outstream.println("Processed "+sd.starts.size+" valid starts and "+sd.stops.size+" stops.");
		sd.clear();
		findStartCodons(sd.bases, sd.starts, startSet);
		findStopCodons(sd.bases, sd.stops, stopSet);
//		outstream.println("Found "+sd.starts.size+" invalid starts and "+sd.stops.size+" stops.");
		processEnds(sd.bases, statsCDS.start, sd.starts, 0);
		processEnds(sd.bases, statsCDS.stop, sd.stops, 0);
	}
	
	private void processRNA(ScafData sd, int strand){
		sd.clear();
		ArrayList<GffLine> lines=sd.rnaLines[strand];
		for(GffLine gline : lines){
			assert(gline.strand==strand);
			final byte flag;
			if(!gline.attributes.contains("partial=true") && gline.start>1 && gline.stop<sd.length()){
				if(gline.type.equals("tRNA")){
					flag=1;
					statstRNA.addLength(gline.length());
				}else{
					assert(gline.type.equals("rRNA"));
					if(gline.attributes.contains("16S")){
						flag=2;
						stats16S.addLength(gline.length());
					}else if(gline.attributes.contains("23S")){
						flag=4;
						stats23S.addLength(gline.length());
					}else if(gline.attributes.contains("5S")){
						flag=8;
						stats5S.addLength(gline.length());
					}else{
						flag=0;
						assert(false) : gline;
					}
				}
				processRnaLine(gline, sd, flag);
			}
		}
		processRnaInner(sd);
		processRnaEnds(sd);
	}
	
	void processRnaInner(ScafData sd){
		byte[] bases=sd.bases;
		byte[] frames=sd.frames;
		final int k=kInnerRNA;//TODO: Note! This is linked to a single static variable for all RNAs.
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k){
					int vf=frames[i];
					for(int type=0; type<4; type++){
						int valid=vf&1;
						rnaContainers[type].inner.add(kmer, 0, valid);
						vf=(vf>>1);
					}
				}
			}else{len=0;}
		}
	}
	
	void processRnaEnds(ScafData sd){
		byte[] bases=sd.bases;

		final int k=stats16S.start.k;
		final int kMax=stats16S.start.kMax;
		final int mask=stats16S.start.mask;
		final long[] counts=new long[kMax];//Slow
		
		int kmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					counts[kmer]++;
				}
			}else{len=0;}
		}
		for(StatsContainer sc : rnaContainers){
			FrameStats fs=sc.start;
			for(long[] array : fs.countsFalse){
				Tools.add(array, counts);
			}
			fs=sc.stop;
			for(long[] array : fs.countsFalse){
				Tools.add(array, counts);
			}
		}
	}
	
	private static BitSet processEnds(byte[] bases, FrameStats stats, IntList list, int valid){
		BitSet points=new BitSet(bases.length);
		for(int i=0; i<list.size; i++){
			int point=list.get(i);
			stats.processPoint(bases, list.get(i), valid);
			points.set(point);
		}
		return points;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Scoring            ----------------*/
	/*--------------------------------------------------------------*/
	
//	//Assumes bases are in the correct strand
//	public float calcStartScore(int start, int stop, byte[] bases){
//		float f=scorePoint(start, bases, startStats);
////		float ss=scoreStart2(start, bases, stop, innerKmerStats);
////		if(ss>0){f=(f+0.0005f*ss);} //Does not seem to help; needs more study.
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcStopScore(int stop, byte[] bases){
//		float f=scorePoint(stop, bases, stopStats);
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcRnaStartScore(int start, int stop, byte[] bases){
//		float f=scorePoint(start, bases, rrnaStartStats);
//		return f;
//	}
//	
//	//Assumes bases are in the correct strand
//	public float calcRnaStopScore(int stop, byte[] bases){
//		float f=scorePoint(stop, bases, rrnaStopStats);
//		return f;
//	}
	
//	public static float calcKmerScore(int start, int stop, int startFrame, byte[] bases, FrameStats stats){
//
//		assert(stats.frames==3);
//		final int k=stats.k;
//		final int mask=~((-1)<<(2*k));
//		
//		int kmer=0;
//		int len=0;
//		float score=0;
//		int numKmers=0;
//		
//		for(int pos=start, currentFrame=startFrame; pos<stop; pos++){
//			final byte b=bases[pos];
//			final int x=AminoAcid.baseToNumber[b];
//
//			if(x>=0){
//				kmer=((kmer<<2)|x)&mask;
//				len++;
//				if(len>=k){
//					float prob=stats.probs[currentFrame][kmer];
//					float dif=prob-0.99f;//Prob above 1 is more likely than average
//					score+=dif;
//					numKmers++;
//				}
//			}else{
//				len=0;
//				kmer=0;
//			}
//
//			currentFrame++;
//			if(currentFrame>2){currentFrame=0;}
//		}
//		return score/Tools.max(1f, numKmers);
//	}
//	
//	/**
//	 * TODO
//	 * Evaluate the relative difference between left and right frequencies.
//	 * The purpose is to find locations where the left side looks noncoding and the right side looks coding.
//	 * Does not currently yield useful results.
//	 */
//	public static float scoreStart2(int point, byte[] bases, int stop, FrameStats stats){
//		final int k=stats.k;
//		
//		int start=point-45;
//		if(start<0 || stop>bases.length){return 0.5f;}
//
//		float left=calcKmerScore(start, Tools.min(point+k-2, bases.length), 0, bases, stats);
//		float right=calcKmerScore(point, stop-3, 0, bases, stats);
//		return right-left; //High numbers are likely to be starts; non-starts should be near 0.
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Long Kmers          ----------------*/
	/*--------------------------------------------------------------*/
	
	public synchronized void loadLongKmers(){
		assert(ssuKmers==null);
		if(load16Skmers){ssuKmers=loadLongKmers(stats16S, kLong16s, "ssu");}
		if(load23Skmers){lsuKmers=loadLongKmers(stats23S, kLong23s, "lsu");}
		if(load5Skmers){r5SKmers=loadLongKmers(stats5S, kLong5s, "5S");}
		if(loadtRNAkmers){trnaKmers=loadLongKmers(statstRNA, kLongTRna, "tRNA");}
	}
	
	private IntHashSet loadLongKmers(StatsContainer sc, int k, String prefix){
		String fname=Data.findPath("?"+prefix+"_"+k+"mers.fa.gz");
		if(!new File(fname).exists()){
			System.err.println("Can't find "+fname);
			return null;
		}
		IntHashSet set=loadLongKmers(fname, k);
		sc.kmerSet=set;
		sc.kLongLen=k;
		return set;
	}
	
	private static IntHashSet loadLongKmers(String fname, int k){//TODO: Consider making this a LongHashSet.  No reason not to...
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FA, null, false, false);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		
		IntHashSet set=new IntHashSet(1000);
		ListNum<Read> ln=cris.nextList();
		while(ln!=null && ln.size()>0){
			processList(ln, set, k);
			cris.returnList(ln);
			ln=cris.nextList();
		}
		if(ln!=null){cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());}
		ReadWrite.closeStream(cris);
		return set;
	}
	
	private static IntHashSet processList(ListNum<Read> ln, IntHashSet set, int k){
		final int mask=~((-1)<<(2*k));
		for(Read r : ln){
			final byte[] bases=r.bases;
			int kmer=0;
			int len=0;
			for(byte b : bases){
				final int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					len++;
					kmer=((kmer<<2)|num)&mask;
					if(len>=k){
						set.add(kmer);
					}
				}else{
					len=0;
				}
			}
		}
		return set;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){

		Collections.sort(fnames);
		taxIds.sort();
		
		bb.append("#BBMap "+Shared.BBMAP_VERSION_STRING+" Prokaryotic Gene Model\n");
		bb.append("#files");
		for(String fname : fnames){
			bb.tab().append(fname);
		}
		bb.nl();
		bb.append("#taxIDs");
		for(int i=0; i<taxIds.size; i++){
			bb.tab().append(taxIds.get(i));
		}
		bb.nl();
//		bb.append("#k_inner\t").append(innerKmerLength).nl();
//		bb.append("#k_end\t").append(endKmerLength).nl();
//		bb.append("#start_left_offset\t").append(startLeftOffset).nl();
//		bb.append("#start_right_offset\t").append(startRightOffset).nl();
//		bb.append("#stop_left_offset\t").append(stopLeftOffset).nl();
//		bb.append("#stop_right_offset\t").append(stopRightOffset).nl();
		bb.append("#scaffolds\t").append(readsProcessed).nl();
		bb.append("#bases\t").append(basesProcessed).nl();
		bb.append("#genes\t").append(genesProcessed).nl();
		bb.append("#GC\t").append(gc(),2).nl();
		bb.append("#ACGTN");
		for(long x : baseCounts){
			bb.tab().append(x);
		}
		bb.nl();

		for(StatsContainer sc : allContainers){sc.appendTo(bb);}
		
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/

	public StatsContainer statsCDS=new StatsContainer("CDS", Orf.CDS);
	public StatsContainer statstRNA=new StatsContainer("tRNA", Orf.tRNA);
	public StatsContainer stats16S=new StatsContainer("16S", Orf.r16S);
	public StatsContainer stats23S=new StatsContainer("23S", Orf.r23S);
	public StatsContainer stats5S=new StatsContainer("5S", Orf.r5S);

	public IntHashSet ssuKmers=null;
	public IntHashSet lsuKmers=null;
	public IntHashSet r5SKmers=null;
	public IntHashSet trnaKmers=null;
	
	StatsContainer[] rnaContainers=new StatsContainer[] {statstRNA, stats16S, stats23S, stats5S};
	StatsContainer[] allContainers=new StatsContainer[] {statsCDS, statstRNA, stats16S, stats23S, stats5S};
	
//	public final FrameStats innerKmerStats=new FrameStats("innerKmerStats", innerKmerLength, 3, 0);
//	public final FrameStats startStats=new FrameStats("startStats", endKmerLength, startFrames, startLeftOffset);
//	public final FrameStats stopStats=new FrameStats("stopStats", endKmerLength, stopFrames, stopLeftOffset);
//
//	public final FrameStats rrnaStartStats=new FrameStats("rrnaStart", 2, 16, 8);
//	public final FrameStats rrnaStopStats=new FrameStats("rrnaStop", 2, 16, 8);
//	
//	public final FrameStats trnaStats=new FrameStats("tRNA", rnaKmerLength, 1, 0);
//	public final FrameStats rrna16Sstats=new FrameStats("16S", rnaKmerLength, 1, 0);
//	public final FrameStats rrna23Sstats=new FrameStats("23S", rnaKmerLength, 1, 0);
//	public final FrameStats rrna5Sstats=new FrameStats("5S", rnaKmerLength, 1, 0);
//	public final FrameStats[] rnaKmerStats=new FrameStats[] {trnaStats, rrna16Sstats, rrna23Sstats, rrna5Sstats};
	
	/*--------------------------------------------------------------*/
	
	public ArrayList<String> fnames=new ArrayList<String>();
	public IntList taxIds=new IntList();

	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	long readsProcessed=0;
	long basesProcessed=0;
	long genesProcessed=0;
	long filesProcessed=0;
	long[] baseCounts=new long[5];
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Setters        ----------------*/
	/*--------------------------------------------------------------*/

	public void setStatics(){
		kInnerCDS=statsCDS.inner.k;
		kStartCDS=statsCDS.start.k;
		kStopCDS=statsCDS.stop.k;
		
		setStartLeftOffset(statsCDS.start.leftOffset);
		setStartRightOffset(statsCDS.start.rightOffset());
		
		setStopLeftOffset(statsCDS.stop.leftOffset);
		setStopRightOffset(statsCDS.stop.rightOffset());
		
		kInnerRNA=stats16S.inner.k;
		kStartRNA=stats16S.start.k;
		kStopRNA=stats16S.stop.k;
	}
	
	public static void setInnerK(int k){
		kInnerCDS=k;
	}
	
	public static void setStartK(int k){
		kStartCDS=k;
	}
	
	public static void setStopK(int k){
		kStopCDS=k;
	}
	
	public static void setStartLeftOffset(int x){
		startLeftOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
	}
	
	public static void setStartRightOffset(int x){
		startRightOffset=x;
		startFrames=startLeftOffset+startRightOffset+1;
//		System.err.println("startLeftOffset="+startLeftOffset+", startRightOffset="+startRightOffset+", frames="+startFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static void setStopLeftOffset(int x){
		stopLeftOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
	}
	
	public static void setStopRightOffset(int x){
		stopRightOffset=x;
		stopFrames=stopLeftOffset+stopRightOffset+1;
//		System.err.println("stopLeftOffset="+stopLeftOffset+", stopRightOffset="+stopRightOffset+", frames="+stopFrames);
//		assert(false) : endLeftOffset+", "+endRightOffset+", "+endFrames;
	}
	
	public static final boolean isStartCodon(int code){
		return code>=0 && code<=63 && isStartCodon[code];
	}
	public static final boolean isStopCodon(int code){
		return code>=0 && code<=63 && isStopCodon[code];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Class Init          ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean[] makeIsCodon(String[] codons){
		boolean[] array=new boolean[64];
		for(String s : codons){
			int x=AminoAcid.toNumber(s);
			array[x]=true;
		}
		return array;
	}
	
	public static int kInnerCDS=6;
	public static int kStartCDS=3;
	public static int kStopCDS=3;

	public static int kInnerRNA=6;
	public static int kStartRNA=3;
	public static int kStopRNA=3;
	
	public static int kLong16s=15;
	public static int kLong23s=15;
	public static int kLong5s=9;
	public static int kLongTRna=9;
	
	static int startLeftOffset(){return startLeftOffset;}
	static int startRightOffset(){return startRightOffset;}
	static int startFrames(){return startFrames;}
	
	private static int startLeftOffset=21; //21 works well for k=4
	private static int startRightOffset=8; //10 works well for k=4
	private static int startFrames=startLeftOffset+startRightOffset+1;
	
	private static int stopLeftOffset=9;
	private static int stopRightOffset=12;
	private static int stopFrames=stopLeftOffset+stopRightOffset+1;

	public static boolean PROCESS_PLUS_STRAND=true;
	public static boolean PROCESS_MINUS_STRAND=true;
	
	public static boolean load16Skmers=true;
	public static boolean load23Skmers=true;
	public static boolean load5Skmers=false;
	public static boolean loadtRNAkmers=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         More Statics         ----------------*/
	/*--------------------------------------------------------------*/
	
	//E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG[7] and one or two others (e.g., an AUU and possibly a CUG).[8][9]
	public static String[] startCodons=new String[] {"ATG", "GTG", "TTG"};
	public static String[] extendedStartCodons=new String[] {"ATG", "GTG", "TTG", "ATT", "CTG", "ATA"};
	public static String[] stopCodons=new String[] {"TAG", "TAA", "TGA"};
	public static boolean[] isStartCodon=makeIsCodon(startCodons);
	public static boolean[] isStopCodon=makeIsCodon(stopCodons);
	
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean errorState=false;
	
}
