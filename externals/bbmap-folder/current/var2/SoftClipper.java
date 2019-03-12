package var2;

import shared.Tools;
import structures.ByteBuilder;

public class SoftClipper {

	public static byte[] softClipMatch(byte[] match, int minClipLength, boolean allowMutation, 
			final int oldStart, final int oldStop, final int[] startStopRvec){

		final int matchScore=100;
		final int subScore=-200;
		final int subScore2=-100;
		final int insScore=-200;
		final int delScore=-200;
		final int delScore2=-10;
		final int clipScore=-1;
		final int nScore=1;

		int insCount=0;
		int delCount=0;
		
		long score=0;
		long maxScore=0;
		int maxPos=-1;
		int maxStart=-1;
		int currentStart=-1;
		byte current='?';
		
		for(int mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
//			long prevScore=score;
			
			if(m=='m' || m=='N' || m=='R'){
				if(score==0){currentStart=mpos;}
				
				score=score+(m=='m' ? matchScore : nScore);
				
				if(score>maxScore){
					maxScore=score;
					maxPos=mpos;
					maxStart=currentStart;
				}
			}else{
				if(m=='S' || m=='s'){
					score=score+(m==current ? subScore2 : subScore);
				}else if(m=='D'){
					score=score+(m==current ? delScore2 : delScore);
					delCount++;
				}else if(m=='I' || m=='X' || m=='Y'){
					score=score+insScore;
					insCount++;
				}else if(m=='C'){
					score=score+clipScore;
				}
				score=Tools.max(0, score);
			}
			current=m;
		}
		
		if(maxScore<1){return match;}
		final int leftClipM=maxStart;
		final int rightClipM=(match.length-maxPos-1);
		int leftClip=0, rightClip=0;
		for(int i=0; i<match.length; i++){
			byte m=match[i];
			if(i<maxStart){
				leftClip+=(m=='D' ? 0 : 1);
			}else if(i>maxPos){
				rightClip+=(m=='D' ? 0 : 1);
			}
		}
		if(leftClip<minClipLength && rightClip<minClipLength){return match;}
		int start=oldStart, stop=oldStop;
		if(delCount==0){
			final byte[] array=allowMutation ? match : match.clone();
			for(int i=0; i<leftClip; i++){array[i]='C';}
			for(int i=0, j=array.length-1; i<rightClip; i++, j--){array[j]='C';}
			startStopRvec[0]=start;
			startStopRvec[1]=stop;
			return array;
		}
		
		ByteBuilder bb=new ByteBuilder(match.length);
		if(leftClip>=minClipLength){
			for(int mpos=0, processed=0; mpos<match.length; mpos++){
				byte m=match[mpos];
				if(mpos>=leftClipM){
					bb.append(m);
				}else{
					if(m=='D'){
						start++;
					}else if(m=='I'){
						start--;
						bb.append('C');
						processed++;
					}else{
						bb.append('C');
						processed++;
					}
				}
			}
		}else{
			bb.append(match);
		}
		if(rightClip>=minClipLength){
			bb.reverseInPlace();
			byte[] temp=bb.toBytes();
			bb.clear();
			for(int mpos=0, processed=0; mpos<temp.length; mpos++){
				byte m=temp[mpos];
				if(mpos>=rightClipM){
					bb.append(m);
				}else{
					if(m=='D'){
						stop--;
					}else if(m=='I'){
						stop++;
						bb.append('C');
						processed++;
					}else{
						bb.append('C');
						processed++;
					}
				}
			}
			bb.reverseInPlace();
		}
		startStopRvec[0]=start;
		startStopRvec[1]=stop;
		return bb.toBytes();
	}
	
}
