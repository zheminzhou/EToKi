package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class VarMap {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	VarMap(ScafMap scafMap_){
		this(scafMap_, -1, -1, -1, -1, -1);
	}

	VarMap(ScafMap scafMap_, int ploidy_, double pairingRate_, double totalQualityAvg_, double mapqAvg_, double readLengthAvg_){
		scafMap=scafMap_;
		ploidy=ploidy_;
		properPairRate=pairingRate_;
		totalQualityAvg=totalQualityAvg_;
		totalMapqAvg=mapqAvg_;
		readLengthAvg=readLengthAvg_;
		maps=new ConcurrentHashMap[WAYS];
		for(int i=0; i<WAYS; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
//	public static VarMap loadVars(String fname, ScafMap scafMap){
//		final ByteFile bf=ByteFile.makeByteFile(fname, true);
//		final VarMap varMap=new VarMap(scafMap);
//		final byte delimiter='\t';
//		int ploidy=-1;
//		double pairingRate=-1;
//		double mapqAvg=-1;
//		double totalQualityAvg=-1;
//		double readLengthAvg=-1;
//		byte[] line=bf.nextLine();
//		while(line!=null && line.length>0){
//			if(line[0]!='#'){
//				Var v=new Var(line, delimiter);
//				varMap.addUnsynchronized(v);
//			}else{
//				String[] split=new String(line).split("\t");
//				String a=split[0], b=(split.length>1 ? split[1] : null);
//				assert(split.length>1) : new String(line);
//				if(a.equalsIgnoreCase("#ploidy")){
//					ploidy=Integer.parseInt(b);
//				}else if(a.equalsIgnoreCase("#pairingRate")){
//					pairingRate=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#totalQualityAvg")){
//					totalQualityAvg=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#mapqAvg")){
//					mapqAvg=Double.parseDouble(b);
//				}else if(a.equalsIgnoreCase("#readLengthAvg")){
//					readLengthAvg=Double.parseDouble(b);
//				}
//			}
//			line=bf.nextLine();
//		}
//		bf.close();
//		varMap.ploidy=ploidy;
//		varMap.properPairRate=(double)pairingRate;
//		varMap.totalQualityAvg=(double)totalQualityAvg;
//		varMap.totalMapqAvg=(double)mapqAvg;
//		varMap.readLengthAvg=(double)readLengthAvg;
//		return varMap;
//	}
//	
//	//#CHROM POS    ID        REF  ALT     QUAL
//	public static VarMap loadVcf(String fname, ScafMap scafMap){
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		VarMap varMap=new VarMap(scafMap);
//		byte[] line=bf.nextLine();
//		while(line!=null && line.length>0){
//			if(line[0]!='#'){
//				Var v;
//				try {
//					v = Var.fromVCF(line, scafMap);
//					varMap.addUnsynchronized(v);
//				} catch (Exception e) {
//					System.err.println("Unable to parse VCF line: '"+new String(line)+"'");
////					throw new RuntimeException(e);
//				}
//			}else{
//				String[] split=new String(line).split("=");
//				if(split.length==2){
//					String a=split[0], b=split[1];
//					if(a.equalsIgnoreCase("##ploidy")){
//						varMap.ploidy=Integer.parseInt(b);
//					}else if(a.equalsIgnoreCase("##properPairRate")){
//						varMap.properPairRate= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##totalQualityAvg")){
//						varMap.totalQualityAvg= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##mapqAvg")){
//						varMap.totalMapqAvg= Double.parseDouble(b);
//					}else if(a.equalsIgnoreCase("##readLengthAvg")){
//						varMap.readLengthAvg= Double.parseDouble(b);
//					}
//				}
//			}
//			line=bf.nextLine();
//		}
//		bf.close();
//		return varMap;
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	
	public boolean containsKey(Var v) {
		return get(v)!=null;
	}
	
	Var get(final Var v){
		final int way=v.start&MASK;
		return maps[way].get(v);
	}
	
	public long size(){
		long size=0;
		for(int i=0; i<maps.length; i++){size+=maps[i].size();}
		return size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Adders            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int add(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		synchronized(map){
			Var old=map.get(v);
			if(old==null){
				map.put(v, v);
				return 1;
			}
			else{
				synchronized(old){
					old.add(v);
				}
			}
		}
		return 0;
	}
	
	int addUnsynchronized(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		Var old=map.get(v);
		if(old==null){
			map.put(v, v);
			return 1;
		}
		old.add(v);
		return 0;
	}
	
	int dumpVars(HashMap<Var, Var> mapT){
		int added=0;
		@SuppressWarnings("unchecked")
		ArrayList<Var>[] absent=new ArrayList[WAYS];
		for(int i=0; i<WAYS; i++){
			absent[i]=new ArrayList<Var>();
		}
		for(Entry<Var, Var> e : mapT.entrySet()){
			Var v=e.getValue();
			final int way=v.start&MASK;
			ConcurrentHashMap<Var, Var> map=maps[way];
			Var old=map.get(v);
			if(old==null){absent[way].add(v);}
			else{
				synchronized(old){
					old.add(v);
				}
			}
		}
		
		mapT.clear();
		for(int way=0; way<WAYS; way++){
			ConcurrentHashMap<Var, Var> map=maps[way];
			ArrayList<Var> list=absent[way];
			synchronized(map){
				for(Var v : list){
					Var old=get(v);
					if(old==null){
						map.put(v, v);
						added++;
					}
					else{
						synchronized(old){
							old.add(v);
						}
					}
				}
			}
		}
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/

	public long[] processVariantsST(VarFilter filter, long[] scoreArray, long[] ploidyArray) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			long[] types2=processVariants(map, filter, scoreArray, ploidyArray, false);
			types2=processVariants(map, filter, ploidyArray, scoreArray, true);
			types2=processVariants(map, filter, ploidyArray, scoreArray, false);
			Tools.add(types, types2);
		}
		return types;
	}
	
	public long[] addSharedVariantsST(VarFilter filter, VarMap sharedVarMap) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[Var.VAR_TYPES];
		for(int i=0; i<maps.length; i++){
			long[] types2=addSharedVariants(maps[i], sharedVarMap.maps[i]);
			Tools.add(types, types2);
		}
		return types;
	}
	
	public long[] processVariantsMT(VarFilter filter, long[] scoreArray, long[] ploidyArray) {
		processVariantsMT_inner(filter, null, null, false);
		processVariantsMT_inner(filter, null, null, true);
		return processVariantsMT_inner(filter, scoreArray, ploidyArray, false);
	}
	
	private long[] processVariantsMT_inner(VarFilter filter, long[] scoreArray, long[] ploidyArray, boolean processInsertions) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(WAYS);
		for(int i=0; i<WAYS; i++){
			ProcessThread pt=new ProcessThread(maps[i], filter, scoreArray!=null, ploidyArray!=null, processInsertions);
			alpt.add(pt);
			pt.start();
		}
		
		long[] types=new long[Var.VAR_TYPES];
		boolean success=true;
		for(ProcessThread pt : alpt){
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			if(pt.types!=null){
				Tools.add(types, pt.types);
			}
			if(scoreArray!=null){Tools.add(scoreArray, pt.scoreArray);}
			if(ploidyArray!=null){Tools.add(ploidyArray, pt.ploidyArray);}
			success&=pt.success;
		}
		
		//Track whether any threads failed
//		if(!success){errorState=true;}
		
		return types;
	}
	
	private class ProcessThread extends Thread {
		
		ProcessThread(Map<Var, Var> map_, VarFilter filter_, boolean trackScores, boolean trackPloidy, boolean processInsertions_){
			map=map_;
			filter=filter_;
			scoreArray=(trackScores ? new long[200] : null);
			ploidyArray=(trackPloidy ? new long[ploidy+1] : null);
			processInsertions=processInsertions_;
		}
		
		@Override
		public void run(){
			types=processVariants(map, filter, scoreArray, ploidyArray, processInsertions);
			success=true;
		}
		
		final VarFilter filter;
		final Map<Var, Var> map;
		long[] types;
		final long[] scoreArray;
		final long[] ploidyArray;
		boolean processInsertions;
		boolean success=false;
	}

	private long[] processVariants(Map<Var, Var> map, VarFilter filter, long[] scoreArray, long[] ploidyArray, boolean processInsertions) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		Iterator<Entry<Var, Var>> iterator=map.entrySet().iterator();
		long[] types=new long[Var.VAR_TYPES];
		while(iterator.hasNext()){
			Entry<Var, Var> entry=iterator.next();
			final Var v=entry.getValue();
			
			if(processInsertions){
				assert(readLengthAvg>0);
				if(v.type()==Var.INS){
					synchronized(v){
						v.reviseAlleleFraction(readLengthAvg, scafMap.getScaffold(v.scafnum), this);
					}
				}
			}else{
				boolean pass=filter.passesFast(v);
				if(pass){
					v.calcCoverage(scafMap);
					pass=filter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
				}
				if(pass){
					types[v.type()]++;
					if(scoreArray!=null){scoreArray[(int)v.phredScore(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap)]++;}
					if(ploidyArray!=null){ploidyArray[v.calcCopies(ploidy)]++;}
				}else{
					iterator.remove();
				}
			}
		}
		return types;
	}

	private long[] addSharedVariants(Map<Var, Var> map, Map<Var, Var> sharedMap) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		for(Var v : sharedMap.keySet()){
			if(!map.containsKey(v)){
				Var v2=new Var(v);
				map.put(v2, v2);
			}
		}
		
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : sharedMap.keySet()){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	public Var[] toArray(boolean sort) {
		Var[] array=new Var[(int)size()];
		int i=0;
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Entry<Var, Var> e : map.entrySet()){
				array[i]=e.getValue();
				i++;
			}
		}
		if(sort){Shared.sort(array);}
		return array;
	}
	
	public long[] calcCoverage(ScafMap scafMap) {
		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Entry<Var, Var> e : map.entrySet()){
				Var v=e.getValue();
				v.calcCoverage(scafMap);
				types[v.type()]++;
			}
		}
		return types;
	}
	
	public long[] countTypes() {
		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Entry<Var, Var> e : map.entrySet()){
				types[e.getValue().type()]++;
			}
		}
		return types;
	}
	
//	public void writeVarFile_(FileFormat ff, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref){
//		Var[] array=toArray(true);
//		ByteStreamWriter bsw=new ByteStreamWriter(ff);
//		bsw.start();
//		ByteBuilder bb=new ByteBuilder(33000);
//		bb.append(Var.toVarHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
//				ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
//		for(Var v : array){
//			v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);//TODO: Track depth
//			bb.nl();
//			if(bb.length()>16384){
//				bsw.print(bb);
//				bb.clear();
//			}
//		}
//		if(bb.length()>0){
//			bsw.print(bb);
//			bb.clear();
//		}
//		bsw.poisonAndWait();
//	}
//	
//	public void writeVcfFile_(String fname, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref, String sampleName, boolean trimWhitespace){
//		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, true, false, false);
//		writeVcfFile(ff, filter, reads, pairs, properPairs, bases, ref, sampleName, trimWhitespace);
//	}
//	
//	public void writeVcfFile_(FileFormat ff, VarFilter filter, long reads, long pairs, long properPairs, long bases, String ref, String sampleName, boolean trimWhitespace){
//		Var[] array=toArray(true);
//		ByteStreamWriter bsw=new ByteStreamWriter(ff);
//		bsw.start();
//		ByteBuilder bb=new ByteBuilder(33000);
//		bb.append(Var.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
//				ploidy, reads, pairs, properPairs, bases, ref, scafMap, sampleName, trimWhitespace)).append('\n');
//		for(Var v : array){
//			v.toVCF(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, trimWhitespace);
//			bb.nl();
//			if(bb.length()>16384){
//				bsw.print(bb);
//				bb.clear();
//			}
//		}
//		if(bb.length()>0){
//			bsw.print(bb);
//			bb.clear();
//		}
//		bsw.poisonAndWait();
//	}
	
	
	public void clear() {
		properPairRate=-1;
		pairedInSequencingRate=-1;
		totalQualityAvg=-1;
		totalMapqAvg=-1;
		readLengthAvg=-1;
		for(int i=0; i<maps.length; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
	@Override
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var v : map.keySet()){
				v.toTextQuick(sb);
				sb.nl();
			}
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public int ploidy=-1;
	public double properPairRate=-1;
	public double pairedInSequencingRate=-1;
	public double totalQualityAvg=-1;
	public double totalMapqAvg=-1;
	public double readLengthAvg=-1;
	public final ScafMap scafMap;
	final ConcurrentHashMap<Var, Var>[] maps;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final int WAYS=4;
	public static final int MASK=WAYS-1;
	
}
