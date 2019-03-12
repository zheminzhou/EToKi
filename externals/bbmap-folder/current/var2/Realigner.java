package var2;

import align2.MSA;
import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

public class Realigner {
	
	public Realigner(){
		this(defaultMaxrows, defaultColumns, defaultPadding, defaultMsaType);
	}
	
	public Realigner(int maxrows_, int columns_, int padding_, String msaType_){
		maxrows=maxrows_;
		columns=columns_;
		padding=padding_;
		msaType=msaType_;
		msa=MSA.makeMSA(maxrows, columns+2, msaType);
	}
	
	public boolean realign(Read r, SamLine sl, final boolean unclip){
		if(!r.mapped() || sl.supplementary() || !sl.primary()){return false;}
		Scaffold scaf=map.getScaffold(sl.rnameS());
		assert(scaf!=null) : sl.rnameS();
		return realign(r, sl, scaf, unclip);
	}
	
	public boolean realign(final Read r, final SamLine sl, final Scaffold scaf, final boolean unclip){
		if(!r.mapped() || sl.supplementary() || !sl.primary()){return false;}
		int[] mSCNID=r.countMatchSymbols();
		int sumBad=mSCNID[1]+mSCNID[4]+mSCNID[5], sumIndel=mSCNID[4]+mSCNID[5];
		if(mSCNID[2]>0){//continue
		}else if(sumBad>3){//continue
		}else if(sumIndel>1 || (sumIndel>0 && mSCNID[1]>1)){//continue
		}else{return false;}

		if(mSCNID[1]<3 && mSCNID[2]==0 && mSCNID[4]<2 && mSCNID[5]<2 && sumBad<3 && sumIndel<2){return false;}
		
		if(r.length()+2>msa.maxColumns+2*padding){
			msa=MSA.makeMSA(msa.maxRows, r.length()+2+r.length()/4+2*padding, msaType);
		}
		if(r.length()+2>msa.maxRows){
			msa=MSA.makeMSA(r.length()+2+r.length()/4+2*padding, msa.maxColumns, msaType);
		}

		final int leadingClip=SamLine.countLeadingClip(sl.cigar, true, false);
		final int trailingClip=SamLine.countTrailingClip(sl.cigar, true, false);
		final int totalClip=leadingClip+trailingClip;
		final boolean clipped=totalClip>0;
		
		final int start=sl.start(true, false);
		final int stop=sl.stop(start, true, false);
		final int paddedStart=start-padding, paddedStop=stop+padding;
		final int len0=paddedStop-paddedStart+1;
		final int a=0, b=len0-1;
		if(len0>=columns){return false;}
		
		realignmentsAttempted++;
		
		SiteScore ss=null;
		final byte[] rbases=makeRbases(scaf.bases, start, stop, padding);
		byte[] qbases=r.bases;
		if(sl.strand()==Shared.MINUS){
			qbases=AminoAcid.reverseComplementBases(qbases);
		}
		
//		SiteScore oldSS=r.toSite(); //123 Slow
//		oldSS.match=r.match.clone(); //123 Slow
//		String oldSL=sl.toString(); //123 Slow
//		assert(oldSS.lengthsAgree()) : oldSS.start+", "+oldSS.stop+", "+oldSS.matchLength()+", "+oldSS.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+oldSS+"\n";
//		assert(!oldSS.matchContainsXY()) : oldSS.start+", "+oldSS.stop+", "+oldSS.matchLength()+", "+oldSS.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+oldSS+"\n";
		
		assert(!r.shortmatch()); //Otherwise convert it or change the score function.
		final int score0=msa.score(r.match);
		final int minScore=clipped ? 1 : Tools.max(1, score0-1);
		
		final int[] max;
		max=msa.fillLimited(qbases, rbases, a, b, minScore, null);
		if(max==null){return false;}

		final int[] score=msa.score(qbases, rbases, a, b, max[0], max[1], max[2], false);
		if(score==null){return false;}
		realignmentsSucceeded++;
		
		if(score[0]<=minScore || (score[0]<=score0 && !clipped)){return false;}
		realignmentsImproved++;
//		assert(false) : score[0]+", "+minScore+", "+score0;
		
		ss=new SiteScore(1, r.strand(), score[1], score[2], 1, score[0]);
		ss.setSlowScore(ss.quickScore);
		ss.score=ss.quickScore;
		ss.match=msa.traceback(qbases, rbases, a, b, score[3], score[4], score[5], false);
		assert(ss.match!=null);
		
		SiteScore oldSS2=ss.clone(); //123 Slow
		oldSS2.match=ss.match.clone(); //123 Slow
//		assert(ss.lengthsAgree()) : ss.start+", "+ss.stop+", "+ss.matchLength()+", "+ss.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+ss+"\n"
//			+oldSS.start+", "+oldSS.stop+", "+oldSS.matchLength()+", "+oldSS.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+oldSS+"\n"
//			;
//		assert(!oldSS.matchContainsXY()) : ss.start+", "+ss.stop+", "+ss.matchLength()+", "+ss.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+ss+"\n"
//			+oldSS.start+", "+oldSS.stop+", "+oldSS.matchLength()+", "+oldSS.mappedLength()+", "+scaf.length+"\n"+sl+"\n"+oldSS+"\n"
//			;
		
		//These pass
//		assert(sl.strand()==ss.strand()) : sl.strand()+", "+ss.strand()+", "+r.strand()+", "+oldSS.strand()+", "+oldSS2.strand();
//		assert(sl.strand()==r.strand());
//		assert(sl.strand()==oldSS.strand());
//		assert(sl.strand()==oldSS2.strand());
		
//		SiteScore old=r.toSite();
		
		
		//Correct for adjusted ref
		ss.start=ss.start+paddedStart;
		ss.stop=ss.stop+paddedStart;
		
		{
			int clipped2=ss.clipTipIndels(scaf.length);
			if(unclip && clipped2>0 && scaf.bases!=null && oldSS2.match[oldSS2.match.length-1]=='Y'){
				ss.unclip(qbases, scaf.bases);
//				System.err.println(sl.strand()+"\n"+new String(new String(oldSS2.match)+"\n"+new String(ss.match)));
			}
		}
		
		if(ss.matchContainsXY()){return false;}
		if(ss.matchContainsAB()){return false;}
		
		realignmentsRetained++;
		
		r.start=ss.start;
		r.stop=ss.stop;
		
//		System.err.println("old:\t"+old+"\nnew:\t"+ss+"\nold score:"+score0+"\n");
//		assert(false);
		
		r.match=ss.match;
		sl.pos=Tools.max(0, r.start)+1;
		sl.cigar=SamLine.toCigar14(r.match, r.start, r.stop, scaf.length, qbases);
		sl.optional=null;
		
//		assert((ss.start>=0) == (SamLine.countLeadingClip(sl.cigar, true, true)==0)) : ss.start+", "+ss.stop+", "+scaf.length+"\n"+sl+"\n"+ss+"\n"
//			+oldSL+"\n"+oldSS+"\n"+oldSS2+"\n"
//		;
//		assert((ss.stop<scaf.length) == (SamLine.countTrailingClip(sl.cigar, true, true)==0)) : ss.start+", "+ss.stop+", "+scaf.length+"\n"+sl+"\n"
//			+oldSL+"\n"+oldSS+"\n"+oldSS2+"\n"
//		;
		
//		assert((ss.start>=0) || (SamLine.countLeadingClip(sl.cigar, true, true)>0)) : ss.start+", "+ss.stop+", "+scaf.length+"\n"+sl+"\n"+ss+"\n"
//			+oldSL+"\n"+oldSS+"\n"+oldSS2+"\n"
//		;
//		assert((ss.stop<scaf.length) || (SamLine.countTrailingClip(sl.cigar, true, true)>0)) : ss.start+", "+ss.stop+", "+scaf.length+"\n"+sl+"\n"
//			+oldSL+"\n"+oldSS+"\n"+oldSS2+"\n"
//		;
		
		//Perform clipping
		r.match=SamLine.cigarToShortMatch_old(sl.cigar, true);
		r.setShortMatch(true);
		r.toLongMatchString(true);
		
		return true;
	}
	
	private static byte[] makeRbases(byte[] bases, int start, int stop, int padding) {
		final int start2=start-padding, stop2=stop+padding;
		byte[] out=new byte[stop2-start2+1];
		for(int opos=0, bpos=start2; opos<out.length; opos++, bpos++){
			byte b=(bpos<0 || bpos>=bases.length ? (byte)'N' : bases[bpos]);
			out[opos]=b;
		}
		return out;
	}

	long realignmentsAttempted=0;
	long realignmentsSucceeded=0;
	long realignmentsRetained=0;
	long realignmentsImproved=0;
	
	private int maxrows=602;
	private int columns=2000;
	private int padding=100;
	private String msaType;
	private MSA msa;
	
	public static int defaultMaxrows=603;
	public static int defaultColumns=2000;
	public static int defaultPadding=200;
	public static String defaultMsaType="MultiStateAligner11ts";
	public static ScafMap map=null;
	
}
