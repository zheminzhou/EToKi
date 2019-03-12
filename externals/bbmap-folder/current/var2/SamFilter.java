package var2;

import java.util.HashSet;

import fileIO.ReadWrite;
import shared.Tools;
import stream.SamLine;

public class SamFilter {
	
	public boolean parse(String arg, String a, String b){

		if(a.equals("min") || a.equals("minpos")){
			minPos=Tools.parseIntKMG(b);
			assert(minPos<=maxPos) : "minPos>maxPos";
		}else if(a.equals("max") || a.equals("maxpos")){
			maxPos=Tools.parseIntKMG(b);
			assert(minPos<=maxPos) : "minPos>maxPos";
		}else if(a.equals("minreadmapq") || a.equals("minsammapq") || a.equals("minmapq")){
			minMapq=Tools.parseIntKMG(b);
		}else if(a.equals("maxreadmapq") || a.equals("maxsammapq") || a.equals("maxmapq")){
			maxMapq=Tools.parseIntKMG(b);
		}else if(a.equals("mapped")){
			includeMapped=Tools.parseBoolean(b);
		}else if(a.equals("unmapped")){
			includeUnmapped=Tools.parseBoolean(b);
		}else if(a.equals("secondary") || a.equals("nonprimary")){
			includeNonPrimary=Tools.parseBoolean(b);
		}else if(a.equals("supplimentary")){
			includeSupplimentary=Tools.parseBoolean(b);
		}else if(a.equals("duplicate") || a.equals("duplicates")){
			includeDuplicate=Tools.parseBoolean(b);
		}else if(a.equals("qfail") || a.equals("samqfail")){
			includeQfail=Tools.parseBoolean(b);
		}else if(a.equals("lengthzero")){
			includeLengthZero=Tools.parseBoolean(b);
		}else if(a.equals("invert")){
			invert=Tools.parseBoolean(b);
		}else if(a.equals("contigs")){
			addContig(b);
		}else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Filters            ----------------*/
	/*--------------------------------------------------------------*/
	
	void addContig(String s){
		if(s==null){return;}
		if(s.indexOf(',')>=0){
			for(String s2 : s.split(",")){
				addContig(s2);
			}
		}
		if(contigs==null){contigs=new HashSet<String>();}
		contigs.add(s);
		if(s.indexOf('_')>=0){addContig(s.replace('_', ' '));}
		String[] split=s.split("\\s+");
		if(split.length>0 && !split[0].equals(s)){contigs.add(split[0]);}
	}
	
	public boolean passesFilter(SamLine sl){
		if(sl==null){return false;}
		return invert^matchesFilter(sl);
	}
	
	boolean matchesFilter(SamLine sl){
		if(sl==null){return false;}
		if(!includeLengthZero && sl.length()<1){return false;}
		
		if(!sl.mapped()){return includeUnmapped;}
		else if(!includeMapped){return false;}

		if(!includeNonPrimary && !sl.primary()){return false;}
		if(!includeSupplimentary && sl.supplementary()){return false;}
		if(!includeDuplicate && sl.duplicate()){return false;}

		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=sl.start(true, false);
			final int stop=sl.stop(start, true, false);
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}

		if(minMapq>Integer.MIN_VALUE || maxMapq<Integer.MAX_VALUE){
			if(sl.mapq>maxMapq || sl.mapq<minMapq){return false;}
		}
		
		if(contigs!=null){
			String rname=sl.rnameS();
			if(rname==null){return false;}
			return contigs.contains(rname);
		}
		
		return true;
	}
	
	public boolean passesFilter(VCFLine vl){
		if(vl==null){return false;}
		return invert^matchesFilter(vl);
	}
	
	boolean matchesFilter(VCFLine vl){
		if(vl==null){return false;}
		
		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=vl.pos-1;
			final int stop=start+(Tools.max(0, vl.reflen-1));
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}
		
		if(contigs!=null){
			String rname=vl.scaf;
			if(rname==null){return false;}
			return contigs.contains(rname);
		}
		
		return true;
	}
	
	public boolean passesFilter(Var v, ScafMap map){
		if(v==null){return false;}
		return invert^matchesFilter(v, map);
	}
	
	boolean matchesFilter(Var v, ScafMap map){
		if(v==null){return false;}
		
		if(minPos>Integer.MIN_VALUE || maxPos<Integer.MAX_VALUE){
			final int start=v.start;
			final int stop=v.stop;
			if(!Tools.overlap(start, stop, minPos, maxPos)){return false;}
		}
		
		if(contigs!=null){
			String rname=map.getName(v.scafnum);
			if(rname==null){return false;}
			return contigs.contains(rname);
		}
		
		return true;
	}
	
	public boolean passesFilter(String name){
		if(name==null){return false;}
		return invert^matchesFilter(name);
	}
	
	boolean matchesFilter(String name){
		if(name==null){return false;}
		if(contigs!=null){
			return contigs.contains(name);
		}
		return true;
	}
	
	public void clear() {
		minMapq=Integer.MIN_VALUE;
		maxMapq=Integer.MAX_VALUE;
	}
	
	public int minPos=Integer.MIN_VALUE;
	public int maxPos=Integer.MAX_VALUE;
	public int minMapq=Integer.MIN_VALUE;
	public int maxMapq=Integer.MAX_VALUE;
	public boolean includeUnmapped=true;
	public boolean includeMapped=true;
	public boolean includeSupplimentary=true;
	public boolean includeQfail=false;
	public boolean includeDuplicate=true;
	public boolean includeNonPrimary=false;
	public boolean includeLengthZero=false;
	public HashSet<String> contigs=null;
	public boolean invert=false;
	
	public void setSamtoolsFilter(){
		ReadWrite.SAMTOOLS_IGNORE_FLAG=0;
		if(!includeUnmapped){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;}
		if(!includeNonPrimary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!includeSupplimentary){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SUPPLIMENTARY;}
		if(!includeQfail){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_QFAIL;}
		if(!includeDuplicate){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_QFAIL;}
	}
	
	/*--------------------------------------------------------------*/
	
}
