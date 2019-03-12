package jgi;

import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.SamLine;

/**
 * Filters to select only reads with substitution errors
 * for bases with quality scores in a certain interval.
 * Used for manually examining specific reads that may have
 * incorrectly calibrated quality scores, for example.
 * @author Brian Bushnell
 * @date May 5, 2015
 *
 */
public class FilterReadsWithSubs extends BBTool_ST {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		FilterReadsWithSubs bbt=new FilterReadsWithSubs(args);
		bbt.process(t);
	}
	
	public FilterReadsWithSubs(String[] args){super(args);}
	
	@Override
	public boolean parseArgument(String arg, String a, String b) {
//		System.err.println("Calling parseArgument("+arg+","+a+","+b+")");
		if(a.equals("minq")){
			minq=Tools.parseIntKMG(b);
			return true;
		}else if(a.equals("maxq")){
			maxq=Tools.parseIntKMG(b);
			return true;
		}else if(a.equals("keepperfect")){
			keepPerfect=Tools.parseBoolean(b);
			return true;
		}else if(a.equals("countindels")){
			countIndels=Tools.parseBoolean(b);
			return true;
		}
		
		//There was no match to the argument
		return false;
	}
	
	@Override
	void setDefaults(){
		SamLine.SET_FROM_OK=true;
		minq=0;
		maxq=99;
		minsubs=1;
		countIndels=true;
		keepPerfect=false;
	}
	
	@Override
	void startupSubclass() {}
	
	@Override
	void shutdownSubclass() {}
	
	@Override
	void showStatsSubclass(Timer t, long readsIn, long basesIn) {
		// TODO Auto-generated method stub

	}
	
	@Override
	boolean processReadPair(Read r1, Read r2) {
		assert(r2==null);
		final byte[] quals=r1.quality, bases=r1.bases;
		final byte[] match=(r1.match==null ? null : !r1.shortmatch() ? r1.match : Read.toLongMatchString(r1.match));
		if(match==null || quals==null || bases==null){return false;}
		
		int subs=0;
		int indels=0;
		for(int qpos=0, mpos=0, last=quals.length-1; mpos<match.length; mpos++){
			
			final byte m=match[mpos];
			final byte mprev=match[Tools.max(mpos-1, 0)];
			final byte mnext=match[Tools.min(mpos+1, match.length-1)];
			
			final byte q1=quals[qpos];
			final byte b2=bases[qpos];
			
			int sub=0, indel=0;
			if(m=='S'){
				sub=1;
			}else if(m=='I'){
				indel=1;
			}else if(m=='m'){
				if(mprev=='D' || mnext=='D'){
					indel=1;
				}
			}else if(m=='D'){
				//do nothing
			}else if(m=='C'){
				//do nothing
			}else{
				throw new RuntimeException("Bad symbol m='"+((char)m)+"'\n"+new String(match)+"\n"+new String(bases)+"\n");
			}
			subs+=sub;
			indels+=indel;
			if(q1>=minq && q1<=maxq){
				if(sub>0 || (indel>0 && countIndels)){return true;}
			}
			
			if(m!='D'){qpos++;}
		}
		return keepPerfect && subs==0 && indels==0;
	}

	public int minq, maxq, minsubs;
	public boolean countIndels;
	public boolean keepPerfect;
	
	
}
