package sketch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import fileIO.TextStreamWriter;
import json.JsonObject;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;

public class SketchResults extends SketchObject {
	
	SketchResults(Sketch s){
		sketch=s;
	}
	
	SketchResults(Sketch s, ArrayList<Sketch> sketchList_, int[][] taxHits_){
		sketch=s;
		refSketchList=sketchList_;
		taxHits=taxHits_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addMap(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, CompareBuffer buffer) {

		if(map.isEmpty()){return;}
		list=addToList(map, params, list);
		
		if((true || params.needContamCounts())){
			recompare(buffer, params);
		}
		if(params.taxFilter!=null){
			int removed=0;
			for(int i=0; i<list.size(); i++){
				Comparison c=list.get(i);
				if((c.taxID()>0 && !params.taxFilter.passesFilter(c.taxID())) || (c.name()!=null && !params.taxFilter.passesFilterByNameOnly(c.name()))){
					list.set(i, null);
					removed++;
				}
			}
			if(removed>0){
				Tools.condenseStrict(list);
			}
		}
	}
	
	public void recompare(CompareBuffer buffer, DisplayParams params){
//		assert(makeIndex || !AUTOSIZE);
		
		assert(!sketch.merged());
		sketch.mergeBitSets();
		
//		System.err.println(sketch.compareBitSet());
//		assert(false) : sketch.compareBitSet().getClass();
		
		for(Comparison c : list){
			c.recompare(buffer, taxHits, params.contamLevel());
		}
		Collections.sort(list, params.comparator);
		Collections.reverse(list);
	}
	
	private static ArrayList<Comparison> addToList(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, ArrayList<Comparison> old){

//		System.err.println(map.size());
//		System.err.println(map.keySet());
		
		ArrayList<Comparison> al=(old==null ? new ArrayList<Comparison>(map.size()) : old);
		for(Entry<Integer, Comparison> e : map.entrySet()){
			al.add(e.getValue());
		}
		Shared.sort(al, params.comparator);
		Collections.reverse(al);
		while(al.size()>params.maxRecords*2+10){
			al.remove(al.size()-1);
		}
		return al;
	}
	
	int filterMeta(DisplayParams params) {
		if(refSketchList==null || refSketchList.isEmpty() || !params.hasMetaFilters()){return 0;}
		return filterMeta(params.requiredMeta, params.bannedMeta, params.requiredMetaAnd);
	}
	
	/**
	 * @param requiredMeta Required metadata tags
	 * @param bannedMeta Banned metadata tags
	 * @param requiredMetaAnd True if all required tags are needed, false if any is sufficient
	 * @return Number of remaining sketches
	 */
	private int filterMeta(ArrayList<String> requiredMeta, ArrayList<String> bannedMeta, boolean requiredMetaAnd) {
		if(refSketchList==null || refSketchList.isEmpty()){return 0;}
		int removed=0;
		for(int i=0; i<refSketchList.size(); i++){
			Sketch ref=refSketchList.get(i);
			if(!ref.passesMeta(requiredMeta, bannedMeta, requiredMetaAnd)){
				refSketchList.set(i, null);
				removed++;
			}
		}
		if(removed>0){Tools.condenseStrict(refSketchList);}
		return refSketchList.size();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int primaryTax(int level){
		//I have no idea how to implement this...
		IntHashMap map=new IntHashMap();
		assert(false);
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String recordBreak="\n"; //"\n\n"
	
	void writeResults(DisplayParams params, TextStreamWriter tsw){
		ByteBuilder sb=toText(params);
		tsw.print(sb);
	}

	public ByteBuilder toText(DisplayParams params){
		assert(params.postParsed);
		if(params.json()){
			JsonObject j=params.toJson(this);
			return j.toText();
		}
		final ByteBuilder sb=params.queryHeader(sketch);
		if(params.format==DisplayParams.FORMAT_QUERY_REF_ANI || params.format==DisplayParams.FORMAT_CONSTELLATION){
			if(list==null || list.isEmpty()){return sb;}
			int idx=0;
			int prevTaxID=0;
			for(Comparison c : list){
				params.formatComparison(c, sb, prevTaxID);
				prevTaxID=c.taxID();
				idx++;
				if(idx>=params.maxRecords){break;}
			}
		}else{
			sb.append(recordBreak);

			if(list==null || list.isEmpty()){
				sb.append("No hits.\n");
			}else{
				if(params.format==DisplayParams.FORMAT_MULTICOLUMN){sb.append(params.header()).append('\n');}
				int idx=0;
				int prevTaxID=0;
				for(Comparison c : list){
					params.formatComparison(c, sb, prevTaxID);
					prevTaxID=c.taxID();
					idx++;
					if(idx>=params.maxRecords){break;}
				}
			}
		}
		return sb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final Sketch sketch;
	public ArrayList<Sketch> refSketchList;
	public int[][] taxHits;
	public ArrayList<Comparison> list;
	public int totalRecords=0;
	
}
