package prok;

import java.util.ArrayList;
import java.util.Locale;

public class ScoreTracker {
	
	public ScoreTracker(int type_){
		type=type_;
	}
	
	public void add(ScoreTracker st){
		geneStartScoreSum+=st.geneStartScoreSum;
		geneStopScoreSum+=st.geneStopScoreSum;
		geneInnerScoreSum+=st.geneInnerScoreSum;
		lengthSum+=st.lengthSum;
		
		geneStartScoreCount+=st.geneStartScoreCount;
		geneStopScoreCount+=st.geneStopScoreCount;
		geneInnerScoreCount+=st.geneInnerScoreCount;
		lengthCount+=st.lengthCount;
	}
	
	public void add(ArrayList<Orf>[] array){
		for(ArrayList<Orf> list : array){add(list);}
	}
	
	public void add(ArrayList<Orf> list){
		if(list==null){return;}
		for(Orf orf : list){
			if(orf.type==type){add(orf);}
		}
	}
	
	public void add(Orf orf){
		if(orf==null || orf.type!=type){return;}
		geneStartScoreSum+=orf.startScore;
		geneStopScoreSum+=orf.stopScore;
		geneInnerScoreSum+=orf.averageKmerScore();
		lengthSum+=orf.length();
		
		geneStartScoreCount++;
		geneStopScoreCount++;
		geneInnerScoreCount++;
		lengthCount++;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("Start Score:          \t "+String.format(Locale.ROOT, "%.4f\n", geneStartScoreSum/geneStartScoreCount));
		sb.append("Stop Score:           \t "+String.format(Locale.ROOT, "%.4f\n", geneStopScoreSum/geneStopScoreCount));
		sb.append("Inner Score:          \t "+String.format(Locale.ROOT, "%.4f\n", geneInnerScoreSum/geneInnerScoreCount));
		sb.append("Length:               \t "+String.format(Locale.ROOT, "%.2f", lengthSum/(double)lengthCount));
		return sb.toString();
	}
	
	long geneStartScoreCount=0;
	long geneStopScoreCount=0;
	long geneInnerScoreCount=0;
	long lengthCount=0;
	
	double geneStartScoreSum=0;
	double geneStopScoreSum=0;
	double geneInnerScoreSum=0;
	long lengthSum=0;
	
	final int type;
	
}
