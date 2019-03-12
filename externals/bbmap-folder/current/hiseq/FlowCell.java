package hiseq;

import java.util.ArrayList;

import fileIO.ByteFile;
import shared.Tools;
import structures.IntList;

public class FlowCell {
	
	public FlowCell(String fname){
//		TextFile tf=new TextFile(fname);
//		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
//			if(line.charAt(0)=='#'){
//
//			}else{
//				String[] split=line.split("\t");
//				int lane=Integer.parseInt(line)
//			}
//		}
		ByteFile bf=ByteFile.makeByteFile(fname, false);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]=='#'){
				String[] split=new String(line).split("\t");
//				System.err.println(new String(line));
				if(Tools.startsWith(line, "#xSize")){
					Tile.xSize=Integer.parseInt(split[1]);
				}else if(Tools.startsWith(line, "#ySize")){
					Tile.ySize=Integer.parseInt(split[1]);
				}
				
				else if(Tools.startsWith(line, "#avgQuality")){
					avgQuality=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#avgUnique")){
					avgUnique=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#avgErrorFree")){
					avgErrorFree=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#avgG")){
					avgG=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#stdQuality")){
					stdQuality=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#stdUnique")){
					stdUnique=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#stdErrorFree")){
					stdErrorFree=Double.parseDouble(split[1]);
				}else if(Tools.startsWith(line, "#stdG")){
					stdG=Double.parseDouble(split[1]);
				}
				
				else if(Tools.startsWith(line, "#reads")){
					readsProcessed=Long.parseLong(split[1]);
				}else if(Tools.startsWith(line, "#avgReads")){
					avgReads=Double.parseDouble(split[1]);
				}
				
			}else{
				int lane=0, tile=0, x1=0, y1=0, x2=0, y2=0, reads=0;
				int discard=line[line.length-1]-'0';
				
				int i=0;
				for(; i<line.length && line[i]!='\t'; i++){
					lane=lane*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					tile=tile*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					x1=x1*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					x2=x2*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					y1=y1*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					y2=y2*10+line[i]-'0';
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					reads=reads*10+line[i]-'0';
				}
				i++;
				
				
				double unique=0, quality=0, error=0;
				for(; i<line.length && line[i]!='\t'; i++){
					if(line[i]!='.'){
						unique=unique*10+line[i]-'0';
					}
				}
				i++;

				for(; i<line.length && line[i]!='\t'; i++){
					if(line[i]!='.'){
						quality=quality*10+line[i]-'0';
					}
				}
				i++;
				
				for(; i<line.length && line[i]!='\t'; i++){
					if(line[i]!='.'){
						error=error*10+line[i]-'0';
					}
				}
				i++;

				unique*=0.001;
				quality*=0.001;
				error*=0.001;
				
				MicroTile mt=getMicroTile(lane, tile, x1, y1);
				assert(mt.x1==x1 && mt.x2==x2 && mt.y1==y1 && mt.y2==y2) : "Micro-tile size seems to be different:\n"
						+ "xsize="+Tile.xSize+", ysize="+Tile.ySize;
				mt.readCount=reads;
				mt.discard=discard;

				mt.misses=(int)(unique*reads*0.01);
				mt.hits=reads-mt.misses;
				mt.qualitySum=reads*quality;
				mt.errorFreeSum=reads*error;
				
//				assert(false) : mt.percentErrorFree()+", "+mt.averageQuality()+", "+mt.uniquePercent();
//				assert(false) : reads+", "+unique+", "+quality+", "+error+"\n"
//					+mt.misses+", "+mt.hits+", "+mt.qualitySum+", "+mt.errorFreeSum+"\n"+new String(line);
			}
		}
	}
	
	public FlowCell(){}
	
	public MicroTile getMicroTile(String id) {//This method is NOT threadsafe
		ihp.parse(id);
		return getMicroTile(ihp);
	}
	
	public MicroTile getMicroTile(IlluminaHeaderParser ihp){
		return getLane(ihp.lane).getMicroTile(ihp.tile, ihp.x, ihp.y);
	}
	
	public MicroTile getMicroTile(int lane, int tile, int x, int y){
		return getLane(lane).getMicroTile(tile, x, y);
	}
	
	public Lane getLane(int lane){
		while(lanes.size()<=lane){lanes.add(new Lane(lanes.size()));}
		return lanes.get(lane);
	}

	public ArrayList<MicroTile> toList() {
		ArrayList<MicroTile> mtList=new ArrayList<MicroTile>();
		for(Lane lane : lanes){
			if(lane!=null){
				for(Tile tile : lane.tiles){
					if(tile!=null){
						for(ArrayList<MicroTile> ylist : tile.xlist){
							if(ylist!=null){
								for(MicroTile mt : ylist){
									if(mt!=null){
										mtList.add(mt);
									}
								}
							}
						}
					}
				}
			}
		}
		return mtList;
	}
	
	public ArrayList<MicroTile> calcStats(){
		ArrayList<MicroTile> mtList=toList();
		readsProcessed=0;
		for(MicroTile mt : mtList){
			mt.process();
			readsProcessed+=mt.readCount;
		}
		avgReads=readsProcessed*1.0/Tools.max(1, mtList.size());
		minCountToUse=(long)Tools.min(500, avgReads*0.25f);
		int toKeep=0;
		for(MicroTile mt : mtList){
			if(mt.readCount>=minCountToUse){toKeep++;}
		}
		
		IntList avgQualityList=new IntList(toKeep);
		IntList avgUniqueList=new IntList(toKeep);
		IntList avgErrorFreeList=new IntList(toKeep);
		IntList avgGList=new IntList(toKeep);
		
		for(MicroTile mt : mtList){
			if(mt!=null && mt.readCount>=minCountToUse){
				avgQualityList.add((int)(1000*mt.averageQuality()));
				avgUniqueList.add((int)(1000*mt.uniquePercent()));
				avgErrorFreeList.add((int)(1000*mt.percentErrorFree()));
				avgGList.add((int)(100000*mt.avgG()));
			}
		}
		
		int[] avgQualityArray=avgQualityList.toArray();
		int[] avgUniqueArray=avgUniqueList.toArray();
		int[] avgErrorFreeArray=avgErrorFreeList.toArray();
		int[] avgGArray=avgGList.toArray();
		
		avgQuality=Tools.averageDouble(avgQualityArray)*0.001;
		avgUnique=Tools.averageDouble(avgUniqueArray)*0.001;
		avgErrorFree=Tools.averageDouble(avgErrorFreeArray)*0.001;
		avgG=Tools.averageDouble(avgGArray)*0.00001;
		
		stdQuality=Tools.standardDeviation(avgQualityArray)*0.001;
		stdUnique=Tools.standardDeviation(avgUniqueArray)*0.001;
		stdErrorFree=Tools.standardDeviation(avgErrorFreeArray)*0.001;
		stdG=Tools.standardDeviation(avgGArray)*0.00001;
		
		return mtList;
	}
	
	public FlowCell widen(int target){
		if(readsProcessed<1){
			System.err.println("Warning: Zero reads processed.");
			return this;
		}
		if(readsProcessed<target){
			return this;
		}
		FlowCell fc=this;
		while(fc.avgReads<target){
			FlowCell fc2=fc.widen();
			fc2.calcStats();
			if(fc2.avgReads<=fc.avgReads){
				unwiden();
				return fc;
			}
			fc=fc2;
		}
		return fc;
	}
	
	public void unwiden(){
		if(Tile.xSize>Tile.ySize){Tile.ySize/=2;}
		else{Tile.xSize/=2;}
	}
	
	public FlowCell widen(){
		if(Tile.xSize>=Tile.ySize){Tile.ySize*=2;}
		else{Tile.xSize*=2;}
		System.err.println("Widening to "+Tile.xSize+"x"+Tile.ySize);
		ArrayList<MicroTile> list=toList();
		FlowCell fc=new FlowCell();
		for(MicroTile mt : list){
			MicroTile mt2=fc.getMicroTile(mt.lane, mt.tile, mt.x1, mt.y1);
			mt2.add(mt);
		}
		return fc;
	}
	
	public ArrayList<Lane> lanes=new ArrayList<Lane>();
	
	long readsProcessed;
	
	public double avgReads;
	public double minCountToUse;
	
	public double avgQuality;
	public double avgUnique;
	public double avgErrorFree;
	public double avgG;
	
	public double stdQuality;
	public double stdUnique;
	public double stdErrorFree;
	public double stdG;
	
	private IlluminaHeaderParser ihp=new IlluminaHeaderParser();
	
}
