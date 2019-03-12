package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.FloatList;
import structures.IntHashSet;
import structures.IntList;


/**
 * This class calls genes within a single thread.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class GeneCaller {
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/
	
	GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_, 
			float minStartScore_, float minStopScore_, float minInnerScore_, 
			float minOrfScore_, float minAvgScore_){
		this(minLen_, maxOverlapSameStrand_, maxOverlapOppositeStrand_, 
				minStartScore_, minStopScore_, minInnerScore_, minOrfScore_, minAvgScore_, null);
	}
	
	GeneCaller(int minLen_, int maxOverlapSameStrand_, int maxOverlapOppositeStrand_, 
			float minStartScore_, float minStopScore_, float minInnerScore_,
			float minOrfScore_, float minAvgScore_, GeneModel pgm_){
		minLen=minLen_;
		maxOverlapSameStrand=maxOverlapSameStrand_;
		maxOverlapOppositeStrand=maxOverlapOppositeStrand_;
		pgm=pgm_;
		
		minStartScore=minStartScore_;
		minStopScore=minStopScore_;
		minInnerScore=minInnerScore_;
		minOrfScore=minOrfScore_;
		minAvgScore=minAvgScore_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public ArrayList<Orf> callGenes(Read r){
		return callGenes(r, pgm);
	}
	
	public ArrayList<Orf> callGenes(Read r, GeneModel pgm_){
		pgm=pgm_;

		pgm.statsCDS.enabled=callCDS;
		pgm.statstRNA.enabled=calltRNA;
		pgm.stats16S.enabled=call16S;
		pgm.stats23S.enabled=call23S;
		pgm.stats5S.enabled=call5S;
		
		final String name=r.id;
		final byte[] bases=r.bases;

		//Lists of all longest orfs per frame
		ArrayList<Orf>[] frameLists=makeOrfs(name, bases, minLen);
		//Lists of all high-scoring orfs per frame, with potentially multiple orfs sharing stops.
		ArrayList<Orf>[] brokenLists=breakOrfs(frameLists, bases);
		
		ArrayList<Orf>[] rnaLists=null;
		if(calltRNA || call16S || call23S || call5S){
			rnaLists=makeRnas(name, bases);

			brokenLists[0].addAll(rnaLists[0]);
			brokenLists[3].addAll(rnaLists[1]);
			Collections.sort(brokenLists[0]);
			Collections.sort(brokenLists[3]);
		}
		
		boolean printAllOrfs=false;
		boolean printRnas=false;
		if(printAllOrfs){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> broken : brokenLists){
				temp.addAll(broken);
			}
			Collections.sort(temp);
			return temp;
		}
		
		if(printRnas && rnaLists!=null){
			ArrayList<Orf> temp=new ArrayList<Orf>();
			for(ArrayList<Orf> list : rnaLists){
				temp.addAll(list);
			}
			Collections.sort(temp);
			return temp;
		}
		
		stCds2.add(brokenLists);
		
		//Find the optimal path through Orfs
		ArrayList<Orf> path=findPath(brokenLists, bases);
		geneStartsOut+=path.size();
		
		stCdsPass.add(path);
		
		return path;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	ArrayList<Orf>[] makeRnas(String name, byte[] bases){
		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] array=new ArrayList[2];
		array[0]=new ArrayList<Orf>();
		array[1]=new ArrayList<Orf>();
		final float[] scores=new float[bases.length];
		final int[] kmersSeen=(pgm.lsuKmers==null && pgm.ssuKmers==null && pgm.trnaKmers==null) ? null : new int[bases.length];
		for(int strand=0; strand<2; strand++){
			for(StatsContainer sc : pgm.rnaContainers){
				if(sc.enabled){
					ArrayList<Orf> list=makeRnasForStrand(name, bases, strand, sc, scores, (sc.kmerSet==null ? null : kmersSeen));//TODO: Make this loop through all RNA types
					if(strand==1 && list!=null){
						for(Orf orf : list){
							assert(orf.strand==strand);
							orf.flip();
						}
					}
					if(list!=null){array[strand].addAll(list);}
				}
			}
			Collections.sort(array[strand]);
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return array;
	}
	
	/** 
	 * Generates lists of all max-length non-overlapping Orfs per frame.
	 * There IS overlap between frames.
	 * All Orfs come out flipped to + orientation. 
	 * */
	static ArrayList<Orf>[] makeOrfs(String name, byte[] bases, int minlen){
		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] array=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				ArrayList<Orf> list=makeOrfsForFrame(name, bases, frame, strand, minlen);
				array[frame+3*strand]=list;
				if(strand==1 && list!=null){
					for(Orf orf : list){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						orf.flip();
					}
				}
			}
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return array;
	}
	
	/**
	 * Dynamic programming phase.
	 * @param frameLists
	 * @param bases
	 * @return
	 */
	private ArrayList<Orf> findPath(ArrayList<Orf>[] frameLists, byte[] bases){
		ArrayList<Orf> all=new ArrayList<Orf>();
		for(ArrayList<Orf> list : frameLists){all.addAll(list);}
		if(all.isEmpty()){return all;}
		Collections.sort(all);
		
		for(Orf orf : all){
			orf.pathScorePlus=-999999;
			orf.pathScoreMinus=-999999;
		}
		
		int[] lastPositionScored=new int[6];
		Arrays.fill(lastPositionScored, -1);

		//Index of highest-scoring ORF in this frame, with prev on the plus strand
		int[] bestIndexPlus=new int[6];
		//Index of highest-scoring ORF in this frame, with prev on the minus strand
		int[] bestIndexMinus=new int[6];
		//Highest-scoring ORF in this frame, with prev on the plus strand
		Orf[] bestOrfPlus=new Orf[6];
		//Highest-scoring ORF in this frame, with prev on the minus strand
		Orf[] bestOrfMinus=new Orf[6];

		int[][] bestIndex=new int[][] {bestIndexPlus, bestIndexMinus};
		Orf[][] bestOrf=new Orf[][] {bestOrfPlus, bestOrfMinus};
		
		for(Orf orf : all){
			final int myListNum=3*orf.strand+orf.frame;
			calcPathScore(orf, frameLists, lastPositionScored, bestIndex);
			if(bestOrfPlus[myListNum]==null || orf.pathScorePlus>=bestOrfPlus[myListNum].pathScorePlus){
				bestOrfPlus[myListNum]=orf;
				bestIndexPlus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
			if(bestOrfMinus[myListNum]==null || orf.pathScoreMinus>=bestOrfMinus[myListNum].pathScoreMinus){
				bestOrfMinus[myListNum]=orf;
				bestIndexMinus[myListNum]=lastPositionScored[myListNum];
				assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf);
			}
		}
		
		Orf best=bestOrf[0][0];
		for(Orf[] array : bestOrf){
			for(Orf orf : array){
				if(best==null || (orf!=null && orf.pathScore()>best.pathScore())){
					best=orf;
				}
			}
		}
		ArrayList<Orf> bestPath=new ArrayList<Orf>();
		for(Orf orf=best; orf!=null; orf=orf.prev()){
			bestPath.add(orf);
		}
		Collections.sort(bestPath);
		return bestPath;
	}
	
	/**
	 * Calculate the best path to this ORF.
	 * @param orf
	 * @param frameLists
	 * @param lastPositionScored
	 * @param bestIndex
	 */
	private void calcPathScore(Orf orf, ArrayList<Orf>[] frameLists, int[] lastPositionScored, int[][] bestIndex){
		final int myListNum=3*orf.strand+orf.frame;

//		System.err.println("* "+orf);
//		System.err.println("* "+Arrays.toString(lastPositionScored));
//		System.err.println();
		
		for(int listStrand=0; listStrand<2; listStrand++){
			for(int listFrame=0; listFrame<3; listFrame++){
				int listNum=listFrame+3*listStrand;
				ArrayList<Orf> list=frameLists[listNum];
				int lastPos=lastPositionScored[listNum];
				int bestPos=bestIndex[listStrand][listNum];
				if(listStrand==0){
					calcPathScorePlus(orf, list, listStrand, lastPos, bestPos);
				}else{
					calcPathScoreMinus(orf, list, listStrand, lastPos, bestPos);
				}
			}
		}
		
//		System.err.println(myListNum+", "+Arrays.toString(lastPositionScored)+", "+frameLists[myListNum].size());
		
		lastPositionScored[myListNum]++;
		assert(frameLists[myListNum].get(lastPositionScored[myListNum])==orf) : myListNum+"\n"+orf+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum])+"\n"
			+Arrays.toString(lastPositionScored)+"\n"+frameLists[myListNum].get(lastPositionScored[myListNum]+1);
		
		//These are sanity checks to make sure that the path did not break in the middle.
		//Safe to disable.
//		assert(orf.prevPlus!=null || orf.stop<100000);
//		assert(orf.prevMinus!=null || orf.stop<100000);
//		assert(orf.pathScore>-10) : orf.pathScore+"\n"+orf+"\n"+orf.prev+"\n";
	}
	
	/**
	 * Calculate the best path to this ORF from a plus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScorePlus(final Orf orf, final ArrayList<Orf> list, final int listStrand, final int lastPos, final int bestPos){
		assert(listStrand==0);
		if(lastPos<0){
			if(orf.prevPlus==null){
				orf.pathScorePlus=orf.orfScore;
				orf.pathLengthPlus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		for(int i=lastPos, min=Tools.max(0, bestPos-lookbackPlus); i>=min || (i>0 && !found); i--){
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				float orfScore=overlap==0 ? orf.orfScore : orf.calcOrfScore(overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
				}
				
				if(overlap<1 && prevScore>0){found=true;}
				if(pathScore>=orf.pathScorePlus){
					orf.pathScorePlus=pathScore;
					orf.prevPlus=prev;
					orf.pathLengthPlus=pathLength;
//					System.err.println("Set as best");
				}
			}
			if(found && prev.stop<maxOverlap-2000 && orf.prevPlus!=null){
				System.err.println("Breaking");
				break;
			}
		}
	}
	
	/**
	 * Calculate the best path to this ORF from a minus-strand previous ORF.
	 * @param orf
	 * @param list
	 * @param listStrand
	 * @param lastPos
	 * @param bestPos
	 */
	private void calcPathScoreMinus(final Orf orf, final ArrayList<Orf> list, final int listStrand, final int lastPos, final int bestPos){
		assert(listStrand==1);
		if(lastPos<0){
			if(orf.prevMinus==null){
				orf.pathScoreMinus=orf.orfScore;
				orf.pathLengthMinus=1;
			}
			return;
		}
		if(list.isEmpty()){return;}
		
//		System.err.println("\nExamining   \t"+orf+"\nlastPos="+lastPos+", bestPos="+bestPos+", sameFrame="+sameFrame);
		boolean found=false;
		final boolean sameStrand=(orf.strand==listStrand);
		final int maxOverlap=(sameStrand ? maxOverlapSameStrand : maxOverlapOppositeStrand);
		for(int i=lastPos, min=Tools.max(0, bestPos-lookbackMinus); i>=min || (i>0 && !found); i--){
			Orf prev=list.get(i);
			assert(prev!=orf) : prev;
//			System.err.println("Comparing to \t"+prev);
			if(orf.isValidPrev(prev, maxOverlap)){
				int overlap=Tools.max(0, prev.stop-orf.start+1);
				float orfScore=overlap==0 ? orf.orfScore : orf.calcOrfScore(overlap);
				
				final float prevScore=prev.pathScore();
				final int prevLength=prev.pathLength();
				
				float pathScore;
				final int pathLength;
				if(sameStrand){
					pathLength=prevLength+1;
					pathScore=prevScore+orfScore;
					pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
				}else{
					pathLength=1;
					pathScore=prev.pathScore()+orfScore;
					pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
				}
				
				if(overlap<1 && prevScore>0){found=true;}
				if(pathScore>=orf.pathScoreMinus){
					orf.pathScoreMinus=pathScore;
					orf.prevMinus=prev;
					orf.pathLengthMinus=pathLength;
//					System.err.println("Set as best");
				}
			}
			if(found && prev.stop<maxOverlap-2000 && orf.prevMinus!=null){
				System.err.println("Breaking");
				break;
			}
		}
	}
	
	/** 
	 * Generates a list of maximal-length Orfs only (non-overlapping).
	 * All Orfs come out in native orientation (unflipped). 
	 * */
	static ArrayList<Orf> makeOrfsForFrame(String name, byte[] bases, int startFrame, int strand, int minlen){
//		assert(false) : "TODO";
		assert(minlen>=3);
		if(bases==null || bases.length<minlen){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>();
		if(!callCDS){return orfs;}
//		int mask=63;
		int code=0;
		int start=-2;
		int frame=0;
		int pos=startFrame;
		
		
		for(; pos<bases.length; pos++){
			byte b=bases[pos];
			int x=AminoAcid.baseToNumber[b];
//			code=((code<<2)|x)&mask;
			code=((code<<2)|x);
			frame++;
			if(frame==3){
				frame=0;
				if(start>=0){
					if(GeneModel.isStopCodon(code) || code<0){//NOTE: This adds a stop codon wherever there are Ns.
						int len=pos-start+1;
						if(len>=minlen){
							Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, Orf.CDS);
							orfs.add(f);
						}
						start=-1;
					}
				}else{
					if(start==-2 || (start<0 && GeneModel.isStartCodon(code))){
						start=pos-2;
					}
				}
				code=0;
			}
		}

		//Add a stop codon at the sequence end.
		if(start>=0){
			pos--;
			while(frame!=3 && frame!=-1){
				pos--;
				frame--;
			}
			int len=pos-start+1;
			if(len>=minlen){
				assert(pos<bases.length) : start+", "+pos+", "+bases.length;
				Orf f=new Orf(name, start, pos, strand, startFrame, bases, true, Orf.CDS);
				orfs.add(f);
			}
		}
		
		return orfs;
	}
	
	/** 
	 * Generates a list of maximal-length Orfs only (non-overlapping).
	 * All Orfs come out in native orientation (unflipped). 
	 * */
	ArrayList<Orf> makeRnasForStrand(String name, byte[] bases, int strand, StatsContainer sc, float[] scores, int[] kmersSeen){
		final int window=sc.lengthAvg;
		if(bases==null || bases.length<window/2){return null;}
		ArrayList<Orf> orfs=new ArrayList<Orf>(sc.type==Orf.tRNA ? 32 : 8);
		
		final FrameStats inner=sc.inner;
//		final FrameStats start=sc.start;
//		final FrameStats stop=sc.stop;
		
		final int k=inner.k;
		final int mask=inner.mask;
//		final float invLen=sc.invLengthAvg;
		final int halfWindow=window/2;
//		final int slop=Tools.max(50, window/8);
		int len=0;
		int kmer=0;
		float currentScore=0;
		float currentScoreAbs=0;
		float bias=biases[sc.type];
		
		float thresh=cutoff1[sc.type];
		float prevScore=0;
		int lastStart=0;
		
		float max=0;
		int maxPos=0;
		
		for(int pos=0; pos<bases.length; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=inner.probs[0][kmer];
					float dif=prob-bias;//Prob above 1 is more likely than average
					currentScoreAbs+=prob;
					currentScore=Tools.max(0, currentScore+dif);
				}
				
				if(currentScore>0){
					if(currentScore>max){
						max=currentScore;
						maxPos=pos;
					}
					if(prevScore<=0){
						lastStart=pos;
					}
				}else{
					int rnaLen=maxPos-lastStart;
					if(max>thresh && rnaLen>=halfWindow){
						Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
						orfs.add(orf);
						orf.orfScore=max;
						if(verbose){
							final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
							final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
							System.err.println(start2+"\t"+stop2+"\t"+max);
						}
					}
					max=0;
					lastStart=pos;
				}
//				System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
				prevScore=currentScore;
				
//				if(pos>=223000 && pos<232000){
////					System.err.println("i="+i+"\tscore="+score+"\tmax="+max+"\tmaxPos="+maxPos+"\tprevScore="+prevScore+"\tlastStart="+lastStart);
//					System.out.println(pos+"\t"+currentScore);
//				}
				
			}else{
				len=0;
				kmer=0;
			}
			scores[pos]=currentScoreAbs;
		}
		
//		System.err.println("size="+orfs.size()+", type="+Orf.typeStrings[sc.type]);
		
		
		{
			int rnaLen=maxPos-lastStart;
			if(max>thresh && rnaLen>=halfWindow){
				Orf orf=new Orf(name, lastStart, maxPos, strand, 0, bases, false, sc.type);
				orfs.add(orf);
				orf.orfScore=max;
				if(verbose){
					final int start2=(strand==0 ? lastStart : bases.length-maxPos-1);
					final int stop2=(strand==0 ? maxPos : bases.length-lastStart-1);
					System.err.println(start2+"\t"+stop2+"\t"+max);
				}
			}
		}
		
		if(kmersSeen!=null && orfs.size()>0 && sc.kmerSet!=null){
			fillKmersSeen(bases, kmersSeen, sc.kmerSet, sc.kLongLen);
		}
		
		float cutoff=cutoff2[sc.type];
		for(int i=0; i<orfs.size(); i++){
			Orf orf=orfs.get(i);
//			System.err.println(orf.orfScore);
			refineRna(orf, bases, strand, sc, scores, kmersSeen);
			if(orf.orfScore<cutoff){
				if(verbose){System.err.println("REJECT: "+orf.toStringFlipped());}
				orfs.set(i, null);
			}else{
				if(verbose){System.err.println("ACCEPT: "+orf.toStringFlipped());}
			}
		}
		Tools.condenseStrict(orfs);
		
//		assert(false);
		
//		for(int pos=0; pos<bases.length; pos++){
//			final byte b=bases[pos];
//			final int x=AminoAcid.baseToNumber[b];
//			kmer=((kmer<<2)|x)&mask;
//			
//			if(x>=0){
//				len++;
//				if(len>=k){
//					float prob=inner.probs[0][kmer];
//					float dif=prob-1.2f;//Prob above 1 is more likely than average
//					currentScore=Tools.max(0, currentScore+dif);
//					if(currentScore>0){
//						currentStreak++;
//					}else{
//						currentStreak=0;
//					}
//					if(currentScore>200 && currentStreak>1500 && currentStreak<1700){
//						Orf orf=new Orf(name, pos-currentStreak-1, pos, strand, 0, bases, false);
//						orfs.add(orf);
//						orf.orfScore=currentScore;
//						orf.startScore=start.scorePoint(orf.start, bases);
//						orf.stopScore=stop.scorePoint(orf.stop, bases);
//						currentStreak=0;
//						currentScore=0;
//					}
//				}
//			}else{
//				len=0;
//				kmer=0;
//			}
//		}
		
		return orfs;
	}
	
	void fillKmersSeen(byte[] bases, int[] kmersSeen, IntHashSet set, int k){
		final int mask=~((-1)<<(2*k));
		int kmer=0;
		int len=0;
		int seen=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final int num=AminoAcid.baseToNumber[b];
			if(num>=0){
				len++;
				kmer=((kmer<<2)|num)&mask;
				if(len>=k && set.contains(kmer)){seen++;}
			}else{
				len=0;
			}
			kmersSeen[i]=seen;
		}
	}
	
	void refineRna(Orf orf, byte[] bases, int strand, StatsContainer sc, float[] scores, int[] kmersSeen){
		if(orf==null){return;}
		if(verbose){System.err.println("REFINE: "+orf.toStringFlipped());}
		final int window=sc.lengthAvg;
		final int halfWindow=window/2;
		final int slop=Tools.max(60, 10+window/8);
		final float invWindow=sc.invLengthAvg;
		final float oldScore=orf.orfScore;
		IntList starts=new IntList();
		IntList stops=new IntList();
		FloatList startScores=new FloatList();
		FloatList stopScores=new FloatList();

		final int leftmost=Tools.max(0, orf.start-slop);
		final int rightmost=Tools.min(bases.length-1, orf.stop+slop);
		if(kmersSeen!=null){
			if(kmersSeen[leftmost]>=kmersSeen[rightmost]){
//				System.err.println("Bad: "+oldScore);
				orf.orfScore=-999;
				return;
			}else{
//				System.err.println("Good: "+oldScore);
			}
		}
		
		{//starts
			final int left=leftmost;
			final int right=Tools.min(bases.length-1, orf.stop+slop-window);
			final float thresh=cutoff3[sc.type];
			fillPoints(left, right, bases, sc.start, thresh, starts, startScores);
		}
		if(verbose){System.err.println("starts: "+starts.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(starts);}
		if(verbose){System.err.println(startScores);}

		{//stops
			final int left=Tools.max(0, orf.start-slop+window);
			final int right=rightmost;
			final float thresh=cutoff4[sc.type];
			fillPoints(left, right, bases, sc.stop, thresh, stops, stopScores);
		}
		if(verbose){System.err.println("stops: "+stops.size);}
//		if((orf.start+"").startsWith("146") || true){System.err.println(stops);}
		if(verbose){System.err.println(stopScores);}

		final float innerThresh=cutoff5[sc.type];
		
		final int minlen=window-slop;
		final int maxlen=window+slop;
		
		orf.orfScore=0;
		int lastStop=0;
		for(int i=0; i<starts.size; i++){
			final int start=starts.get(i);
			final int startSeen=kmersSeen==null ? 0 : kmersSeen[start];
			final float startScore=startScores.get(i);
//			System.err.println("start="+start);
			for(int j=lastStop; j<stops.size; j++){
				final int stop=stops.get(j);
//				System.err.println("stop="+stop);
				if(stop<start+minlen){
					lastStop=j+1;
				}else if(stop>start+window+maxlen){
//					System.err.println("broke");
					break;
				}else if(kmersSeen!=null && kmersSeen[stop]<=startSeen){//TODO: Test this
//					//skip
				}else{
					final int len=stop-start+1;
					final float stopScore=stopScores.get(j);
					final float innerScore=(scores[stop]-scores[start])/len;
//					System.err.println("innerScore="+innerScore);
					if(innerScore>=innerThresh){
						final float a=Tools.max(startScore+0.25f, 0.25f);
						final float b=Tools.max(stopScore+0.25f, 0.25f);
						final float c=innerScore*innerScore;
						final float d=(window-(2.4f*Tools.absdif(len, window)))*invWindow;
						final float score=c*d*(float)Math.sqrt(a*b);
						if(verbose && score>=0.01f*orf.orfScore){
							final int start2=(strand==0 ? start : bases.length-stop-1);
							final int stop2=(strand==0 ? stop : bases.length-start-1);
							System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
						}
						if(score>orf.orfScore){
							orf.start=start;
							orf.stop=stop;
							orf.kmerScore=innerScore*len;
//							if(verbose){
//								final int start2=(strand==0 ? start : bases.length-stop-1);
//								final int stop2=(strand==0 ? stop : bases.length-start-1);
//								System.err.println(start2+"-"+stop2+", "+startScore+", "+stopScore+", "+innerScore+", "+(score*scoreMult[sc.type])+", "+oldScore);
//							}
							orf.startScore=startScore;
							orf.stopScore=stopScore;
							orf.orfScore=score;
						}
					}else{
						assert(true);
					}
				}
			}
		}
		orf.orfScore*=scoreMult[sc.type];
	}
	
	void fillPoints(final int left, final int right, final byte[] bases, final FrameStats fs, float thresh, final IntList points, final FloatList scores){
		points.clear();
		scores.clear();
		final float minThresh=thresh;//thresh*0.05f;
//		System.err.println(left+", "+right+", "+thresh+", "+minThresh);
		while(points.size()<8 && thresh>=minThresh){
//			System.err.println("Running with thresh="+thresh);
			points.clear();
			scores.clear();
			for(int i=left; i<right; i++){
				float score=fs.scorePoint(i, bases);
//				System.err.println(i+", "+score);
				if(score>=thresh){
					points.add(i);
					scores.add(score);
				}
			}
			thresh=thresh*0.75f;
		}
	}
	
	/**
	 * Generate all possible genes from each Orf, and return them in a new set of lists.
	 * @param frameLists
	 * @param bases
	 * @return Lists of orfs.
	 */
	private ArrayList<Orf>[] breakOrfs(ArrayList<Orf>[] frameLists, byte[] bases){

		@SuppressWarnings("unchecked")
		ArrayList<Orf>[] brokenLists=new ArrayList[6];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int fnum=frame+3*strand;
				ArrayList<Orf> longest=frameLists[fnum]; //Longest Orf per stop
				ArrayList<Orf> broken=new ArrayList<Orf>(); //All high scoring Orfs, including multiple per stop, with different starts
				if(longest!=null) {
					for(Orf orf : longest){
						assert(orf.frame==frame);
						assert(orf.strand==strand);
						ArrayList<Orf> temp=breakOrf(orf, bases);
						if(temp!=null){
							broken.addAll(temp);
						}
					}
				}
				Collections.sort(broken);
				brokenLists[fnum]=broken;
			}
			//Reverse-complement bases after processing each strand
			AminoAcid.reverseComplementBasesInPlace(bases);
		}
		return brokenLists;
	}
	
	/**
	 * Generate an Orf for each possible start codon.
	 * Retain only the high-scoring ones.
	 * @param longest Longest open reading frame for a given stop.
	 * @param bases Bases, oriented for this Orf.
	 * @return List of Orfs.
	 */
	private ArrayList<Orf> breakOrf(Orf longest, byte[] bases){
		assert(longest.start<longest.stop);
		final int flipped=longest.flipped();
		if(flipped==1){longest.flip();}//Now the orf is aligned to its native strand
		
		geneStopsMade++;
		
		final FrameStats innerStats=pgm.statsCDS.inner;
		final FrameStats startStats=pgm.statsCDS.start;
		final FrameStats stopStats=pgm.statsCDS.stop;
		
		final String name=longest.scafName;
		final int start=longest.start;
		final int stop=longest.stop;
		final int strand=longest.strand;
		final int max=Tools.min(longest.stop-2, longest.stop-minLen+4);

		assert(pgm!=null) : pgm;
		assert(pgm.statsCDS!=null) : pgm;
		assert(pgm.statsCDS.inner!=null) : pgm.statsCDS;
		assert(pgm.statsCDS.inner.k>0) : pgm.statsCDS.inner;
		
		final int k=innerStats.k;
		final int mask=~((-1)<<(2*k));

		final float stopScore=stopStats.scorePoint(longest.stop, bases);
		stCds.geneStopScoreSum+=stopScore;
		stCds.geneStopScoreCount++;
		
		ArrayList<Orf> broken=new ArrayList<Orf>();
		int created=0;
		
		int codon=0;
		int kmer=0;
		int len=0;
		int numKmers=0;
		float currentScore=0;
		for(int pos=start, currentFrame=0; pos<=stop; pos++){
			final byte b=bases[pos];
			final int x=AminoAcid.baseToNumber[b];
			codon=((codon<<2)|x);
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				if(len>=k){
					float prob=innerStats.probs[currentFrame][kmer];
					float dif=prob-0.99f;//Prob above 1 is more likely than average
					currentScore+=dif;
//					outstream.println("pos="+pos+"\tdif="+String.format("%.4f", dif)+",\tscore="+String.format("%.4f", currentScore)+
//							"\tasStart="+String.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+String.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}else{
//					outstream.println("pos="+pos+"\tdif="+String.format("%.4f", 0.0)+",\tscore="+String.format("%.4f", 0.0)+
//							"\tasStart="+String.format("%.4f", pgm.calcStartScore(pos-2, bases))+"\tasStop="+String.format("%.4f", stopStats.scorePoint(pos, bases))+
//							"\tcodon="+AminoAcid.kmerToString(kmer, 3)+" frame="+(currentFrame));
				}
			}else{
				len=0;
				kmer=0;
			}
			
			currentFrame++;
//			outstream.println("pos="+pos+", codon="+AminoAcid.kmerToString(kmer, 3)+", frame="+currentFrame+", start="+start+", isStartCodon="+pgm.isStartCodon(codon));
			if(currentFrame>2){
				currentFrame=0;
				if(pos<max && created<breakLimit && (pos==start+2 || pgm.isStartCodon(codon))){
//					outstream.println(x);
					int glen=stop-pos+3;
					assert(glen>=minLen) : "glen="+glen+", minLen="+minLen+", pos="+pos+", max="+max+", start="+start;
					
					int oStart=pos-2;
					float startScore=startStats.scorePoint(oStart, bases);

					stCds.geneStartScoreSum+=startScore;
					stCds.geneStartScoreCount++;
					
					stCds.lengthSum+=(stop-(pos-2)+1);
					stCds.lengthCount++;
					
					if((startScore>=minStartScore || pos<6) /* && stopScore>=minStopScore /*|| broken.isEmpty()*/){
						Orf orf=new Orf(name, pos-2, stop, strand, longest.frame, bases, false, longest.type);
						
						geneStartsMade++;
						orf.kmerScore=currentScore;
						orf.startScore=startScore;
						orf.stopScore=stopScore;

						assert(orf.frame==longest.frame);
						assert(orf.strand==strand);

						if(strand==1){orf.flip();}
						broken.add(orf);
						created++;
					}
				}
				codon=0;
			}
		}
		
		final int size=broken.size();
		final int sizeCutoff=Tools.max(5, size/2);
		if(size<1){return broken;}
		Orf best=broken.get(0);
		Orf bestStart=broken.get(0);
		for(int i=0; i<size; i++){
			Orf orf=broken.get(i);
			
			//This fixes scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			orf.kmerScore=currentScore-orf.kmerScore;
			orf.orfScore=orf.calcOrfScore();
			if(orf.orfScore>=best.orfScore){best=orf;}
			if(orf.startScore>=bestStart.startScore){bestStart=orf;}
			
			stCds.geneInnerScoreSum+=orf.averageKmerScore();
			stCds.geneInnerScoreCount++;
		}
		
		//Sort to by score descending to eliminate low-scoring copies. 
		if(broken.size()>1){Collections.sort(broken, Feature.featureComparatorScore);}
		
		int removed=0;
		for(int i=0; i<size; i++){//Fix scores because they were generated together, from start to stop, to make this O(N) instead of O(N^2).
			Orf orf=broken.get(i);
			if(orf.averageKmerScore()<minInnerScore || orf.orfScore<minOrfScore || 
					orf.orfScore/orf.length()<minAvgScore || 
					orf.orfScore<0.5f*best.orfScore-10 || (orf.startScore<bestStart.startScore-0.55f && orf.kmerScore<best.kmerScore*1.1f && orf!=best)){
				broken.set(i, null);
				removed++;
			}else if(i>sizeCutoff){
				broken.set(i, null);
				removed++;
			}
		}
		if(removed>0){
			Tools.condenseStrict(broken);
		}
		
		geneStartsRetained+=broken.size();
		geneStopsRetained+=(broken.size()>0 ? 1 : 0);
		
		if(flipped==1){longest.flip();}
		return broken;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Current gene model.
	 * TODO: Dynamically swap this as needed for contigs with varying GC.
	 */
	GeneModel pgm;
	
	//Gene-calling cutoffs
	final int minLen;
	final int maxOverlapSameStrand;
	final int maxOverlapOppositeStrand;
	final float minStartScore;
	final float minStopScore;
	final float minInnerScore;
	final float minOrfScore;
	final float minAvgScore;

//	static float[] cutoff1=new float[] {0, 40f, 200f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 60f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.4f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 0.6f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.55f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 45f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.3f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.5f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.15f};

	//for k=4,2,2 
//	static float[] cutoff1=new float[] {0, 40f, 140f, 150f, 40f};
//	static float[] cutoff2=new float[] {0, 44f, 300f, 150f, 45f};
//	static float[] cutoff3=new float[] {0, 1.7f, 1.1f, 1.1f, 1.8f};
//	static float[] cutoff4=new float[] {0, 1.6f, 1.3f, 0.45f, 1.1f};
//	static float[] cutoff5=new float[] {0, 1.0f, 1.0f, 1.0f, 1.3f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 200f, 175f, 12f};
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};

//	//for k=5,3,3
//	static float[] cutoff1=new float[] {0, 40f, 300f, 300f, 135f};//sum score
//	static float[] cutoff2=new float[] {0, 44f, 300f, 400f, 100f};//orf score
//	static float[] cutoff3=new float[] {0, 1.7f, 1.5f, 1.5f, 3.5f};//start score
//	static float[] cutoff4=new float[] {0, 1.6f, 1.5f, 1.4f, 1.8f};//stop score
//	static float[] cutoff5=new float[] {0, 1.0f, 1.1f, 1.15f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 8f, 80f, 120f, 5f};//score mult
//	static float[] biases=new float[] {1f, 1.17f, 1.2f, 1.2f, 1.22f};//bias for sum score

	//for k=6,3,3 - this has low scores for correct 23s stops; may try k=2 or k=4 for that end 
	//Also, 5S seems to score low very for archaea
	//Order is: CDS=0, tRNA=1, r16S=2, r23S=3, r5S=4, RNA=5
	static float[] cutoff1=new float[] {0, 100f, 600f, 450f, 300f};//sum score
	static float[] cutoff2=new float[] {0, 48f, 300f, 400f, 32f};//orf score
	static float[] cutoff3=new float[] {0, 3.0f, 1.8f, 1.8f, 4.0f};//start score
	static float[] cutoff4=new float[] {0, 2.8f, 2.0f, 1.25f, 1.9f};//stop score
	static float[] cutoff5=new float[] {0, 2.8f, 1.1f, 1.15f, 1.4f};//inner score
	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f};//score mult
	static float[] biases=new float[] {1f, 1.25f, 1.2f, 1.2f, 1.22f};//bias for sum score

//	static float[] cutoff1=new float[] {0, 50f, 600f, 450f, 300f};//sum score
//	static float[] cutoff2=new float[] {0, 24f, 300f, 400f, 32f};//orf score
//	static float[] cutoff3=new float[] {0, 1.5f, 1.8f, 1.8f, 4.0f};//start score
//	static float[] cutoff4=new float[] {0, 1.4f, 2.0f, 1.25f, 1.9f};//stop score
//	static float[] cutoff5=new float[] {0, 1.4f, 1.1f, 1.15f, 1.4f};//inner score
//	static float[] scoreMult=new float[] {1f, 1f, 35f, 80f, 1.0f};//score mult
//	static float[] biases=new float[] {1f, 1.25f, 1.2f, 1.2f, 1.22f};//bias for sum score
	
	long geneStopsMade=0;
	long geneStartsMade=0;
	long geneStartsRetained=0;
	long geneStopsRetained=0;
	long geneStartsOut=0;

	ScoreTracker stCds=new ScoreTracker(Orf.CDS);
	ScoreTracker stCds2=new ScoreTracker(Orf.CDS);
	ScoreTracker stCdsPass=new ScoreTracker(Orf.CDS);
	
	public static int breakLimit=12;
	public static int lookbackPlus=70;
	public static int lookbackMinus=25;

	public static boolean callCDS=true;
	public static boolean calltRNA=true;
	public static boolean call16S=true;
	public static boolean call23S=true;
	public static boolean call5S=true;
	
//	pathScore+=p0+p1*(Tools.mid(p5*(p2+pathLength), p6*(p3-pathLength), p4));
	
	//Operon length modifiers for same strand
	public static float p0=0f;
	public static float p1=0.6f; //Sensitive - higher decreases FP and TP
	public static float p2=1.0f;//Insensitive - higher positive decreases FP and TP
	public static float p3=12f; //Higher decreases FP (substantially) and TP (slightly)
	public static float p4=-50f; //Lower decreases FP (weakly) and TP (greatly)
	public static float p5=2.0f; //Insensitive - lower increases FP and TP
	public static float p6=2f; //Greater increases FP and TP
	
//	pathScore+=q1+Tools.mid(q2*prevLength, q3+q4*prevLength, q5);
	
	//Operon length modifiers for opposite strand
	public static float q1=-35f;
	public static float q2=-1.5f; //q2 and q4 must have opposite signs
	public static float q3=-10f;
	public static float q4=1f; //(Lower [even negative] decreases FP and TP)
	public static float q5=-60f;
	
	private static PrintStream outstream=System.err;
	static boolean verbose;
	
}
