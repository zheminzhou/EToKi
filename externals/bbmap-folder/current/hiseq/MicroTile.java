package hiseq;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

public class MicroTile {

	public MicroTile(){}

	public MicroTile(int lane_, int tile_, int x1_, int x2_, int y1_, int y2_){
		lane=lane_;
		tile=tile_;
		x1=x1_;
		x2=x2_;
		y1=y1_;
		y2=y2_;
	}
	
	void process(){
		if(tracker!=null){tracker.process();}
	}
	
	public boolean contains(int x, int y){
		return x>=x1 && x<=x2 && y>=y1 && y<=y2;
	}
	
	@Override
	public String toString(){
		return lane+", "+tile+", "+x1+", "+x2+", "+y1+", "+y2;
	}
	
	public double averageQuality(){
		return readCount==0 ? 0 : qualitySum/readCount;
	}
	
	public double percentErrorFree(){
		return readCount==0 ? 0 : errorFreeSum/readCount;
	}
	
	public double hitPercent(){
		long count=kmerCount();
		return count==0 ? 0 : hits*100.0/count;
	}
	
	public double uniquePercent(){
		long count=kmerCount();
		return count==0 ? 0 : misses*100.0/count;
	}
	
	public double avgG(){
		return tracker==null ? 0 : tracker.avg('G');
	}
	
	public double maxG(){
		return tracker==null ? 0 : tracker.max('G');
	}

	public long kmerCount(){return hits+misses;}
	
	public void add(MicroTile mt) {
		hits+=mt.hits;
		misses+=mt.misses;
		readCount+=mt.readCount;
		qualitySum+=mt.qualitySum;
		errorFreeSum+=mt.errorFreeSum;

		for(int i=0; i<acgtn.length; i++){
			acgtn[i]+=mt.acgtn[i];
		}
		homoPolyGCount+=mt.homoPolyGCount;
		homoPolyGSum+=mt.homoPolyGSum;
		if(TRACK_CYCLES){
			tracker.add(mt.tracker);
		}
	}
	
	public void add(Read r){
		if(r==null){return;}
		final int len=r.length();
		if(len<1){return;}
		
		readCount++;
		qualitySum+=r.avgQualityByProbabilityDouble(true, len);
		errorFreeSum+=100*r.probabilityErrorFree(true, len);
		
		final byte[] bases=r.bases;
		int maxPolyG=0, currentPolyG=0;
		for(int i=0; i<len; i++){
			byte b=bases[i];
			byte x=AminoAcid.baseToNumberACGTN[b];
			acgtn[x]++;
			if(b=='G'){
				currentPolyG++;
				maxPolyG=Tools.max(currentPolyG, maxPolyG);
			}else{
				currentPolyG=0;
			}
		}
		homoPolyGCount+=(maxPolyG>=MIN_POLY_G ? 1 : 0);
		homoPolyGSum+=maxPolyG;
		if(TRACK_CYCLES){
			tracker.add(r);
		}
	}
	
	public void toText(ByteBuilder bb){
		bb.append(lane).append('\t');
		bb.append(tile).append('\t');
		bb.append(x1).append('\t');
		bb.append(x2).append('\t');
		bb.append(y1).append('\t');
		bb.append(y2).append('\t');
		bb.append(readCount).append('\t');
		
		bb.append(uniquePercent(), 3).append('\t');
		bb.append(averageQuality(), 3).append('\t');
		bb.append(percentErrorFree(), 3).append('\t');
		bb.append(discard);
		
		for(int i=0; i<5; i++){
			bb.tab().append(acgtn[i]);
		}
		bb.tab().append(homoPolyGCount);
		bb.tab().append(homoPolyGSum);
		
		bb.nl();
	}
	
	public long hits;
	public long misses;
	public long readCount;
	public double qualitySum;
	public double errorFreeSum;
	
	public int discard=0;
	
	public int lane;
	public int tile;
	public int x1, x2;
	public int y1, y2;
	
	//TODO: These fields are not currently parsed.
	public long[] acgtn=new long[5];
	public long homoPolyGCount;
	public long homoPolyGSum;
	
	public final CycleTracker tracker=TRACK_CYCLES ? new CycleTracker() : null;

	public static int MIN_POLY_G=25;
	public static boolean TRACK_CYCLES=false;
	
}
