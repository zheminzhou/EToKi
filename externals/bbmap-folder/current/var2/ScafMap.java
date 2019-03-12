package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;

public class ScafMap {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	public ScafMap(){}

	/*--------------------------------------------------------------*/

	public static ScafMap loadSamHeader(String fname){return loadSamHeader(fname, null);}
	public static ScafMap loadSamHeader(FileFormat ff){return loadSamHeader(ff, null);}
	
	public static ScafMap loadSamHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, false, false);
		return loadSamHeader(ff, scafMap);
	}
	
	public static ScafMap loadSamHeader(FileFormat ff, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(ff);
		if(scafMap==null){scafMap=new ScafMap();}
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			//assert(false) : new String(line);
			if(Tools.startsWith(line, "@SQ\t")){
				scafMap.add(line);
			}else if(line[0]!='@'){
				break;
			}
			line=bf.nextLine();
		}
		bf.close();

		return scafMap;
	}

	/*--------------------------------------------------------------*/

	public static ScafMap loadVcfHeader(String fname){return loadVcfHeader(fname, null);}
	public static ScafMap loadVcfHeader(FileFormat ff){return loadVcfHeader(ff, null);}
	
	public static ScafMap loadVcfHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.VCF, null, false, false);
		return loadVcfHeader(ff, scafMap);
	}
	
	public static ScafMap loadVcfHeader(FileFormat ff, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(ff);
		if(scafMap==null){scafMap=new ScafMap();}
		byte[] line=bf.nextLine();
//		System.err.println(new String(line));
//		assert(line!=null) : bf.name();
//		assert(line.length>0) : bf.name();
		while(line!=null && line.length>0){
			//assert(false) : new String(line);
			if(Tools.startsWith(line, "##contig=<ID=")){
				scafMap.addFromVcf(line);
			}else if(line[0]!='#'){
				break;
			}
			line=bf.nextLine();
//			System.err.println(new String(line));
		}
		bf.close();
//		assert(false);
		return scafMap;
	}

	/*--------------------------------------------------------------*/

	public static ScafMap loadReference(String fname, boolean makeDefault){return loadReference(fname, null, null, makeDefault);}
	public static ScafMap loadReference(FileFormat ff, boolean makeDefault){return loadReference(ff, null, null, makeDefault);}
	
	public static ScafMap loadReference(String fname, ScafMap scafMap, SamFilter samFilter, boolean makeDefault){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return loadReference(ff, scafMap, samFilter, makeDefault);
	}
	
	public static ScafMap loadReference(FileFormat ff, ScafMap map, SamFilter samFilter, boolean makeDefault){
		if(defaultScafMapFile!=null && defaultScafMapFile.equals(ff.name())){return defaultScafMap;}
		ArrayList<Read> reads=ReadInputStream.toReads(ff, -1);
		if(map==null){map=new ScafMap();}
		for(Read r : reads){
			if(samFilter==null || samFilter.passesFilter(r.id)){
				map.addScaffold(r);
			}
		}
		if(makeDefault){setDefaultScafMap(map, ff.name());}
		return map;
	}

	public static ScafMap defaultScafMap(){return defaultScafMap;}
	public static void setDefaultScafMap(ScafMap map, String fname){
		SamLine.RNAME_AS_BYTES=false;
		assert(fname==null || defaultScafMapFile==null || !fname.equals(defaultScafMapFile)) : fname;
		defaultScafMap=map;
		defaultScafMapFile=fname;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Inherited           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){
		map.clear();
		alt.clear();
		list.clear();
	}
	
	public int size(){return list.size();}
	
//	public boolean containsKey(String s){
//		return getScaffold(s)!=null;
//	}
	
	public Set<String> keySet() {
		return map.keySet();
	}
	
	public Set<String> altKeySet() {
		return alt.keySet();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Adders             ----------------*/
	/*--------------------------------------------------------------*/
	
	public Scaffold addScaffold(Read r){
		Scaffold scaf=map.get(r.id);
		if(scaf==null){
			scaf=new Scaffold(r.id, size(), r.length());
			add(scaf);
		}
		scaf.bases=r.bases;
		assert(scaf.bases.length==scaf.length);
		return scaf;
	}
	
	public Scaffold add(byte[] line){
		Scaffold scaf=new Scaffold(line, size());
		Scaffold old=map.get(scaf.name);
		if(old!=null){return old;}
		return add(scaf);
	}
	
	//##contig=<ID=Escherichia_coli,length=4639675>
	public Scaffold addFromVcf(byte[] line){
		int comma=Tools.indexOf(line, (byte)',');
		String name=new String(line, 13, comma-13);
		int length=Tools.parseInt(line, comma+8, line.length-1);
		Scaffold scaf=new Scaffold(name, size(), length);
		Scaffold old=map.get(scaf.name);
		if(old!=null){return old;}
		return add(scaf);
	}
	
	public Scaffold add(String s, int len){
		Scaffold scaf=map.get(s);
		if(scaf!=null){return scaf;}
		scaf=new Scaffold(s, size(), len);
		return add(scaf);
	}
	
	private Scaffold add(Scaffold scaf){
		assert(!map.containsKey(scaf.name));
		assert(size()==scaf.number);
		
		list.add(scaf);
		map.put(scaf.name, scaf);
		String s=scaf.name;
		if(TRIM_WHITESPACE_ALSO){
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					String s2=s.substring(0, i);
					boolean b=alt.containsKey(s2);
					assert(!b);
					if(!b){
						alt.put(s2, scaf);
					}
					break;
				}
			}
		}
		return scaf;
	}
	
	public void addCoverage(SamLine sl){
		if(!sl.mapped()){return;}
		Scaffold scaf=getScaffold(sl);
		scaf.add(sl);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int getNumber(String s){
		Scaffold value=getScaffold(s);
		return value==null ? -1 : value.number;
	}
	
	public String getName(int number){
		return number>=size() ? null : list.get(number).name;
	}
	
	public int getLength(int number){
		return number>=size() ? 0 : list.get(number).length;
	}
	
	public Scaffold getScaffold(int number){
		return number>=size() ? null : list.get(number);
	}
	
	public Scaffold getScaffold(String s){
		Scaffold value=map.get(s);
		if(value==null){value=alt.get(s);}
		if(value==null && TRIM_WHITESPACE_ALSO){
			int index=-1;
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					index=i;
					break;
				}
			}
			if(index>0){
				String sub=s.substring(0, index);
				value=alt.get(sub);
				if(value==null){value=map.get(sub);}
			}
			assert(value!=null) : "Scaffold not present in reference: "+s+"\n"+index+"\n"+s.substring(0, Tools.max(1, index))
				+"\n"+keySet()+"\n"+altKeySet()+"\n";
		}
		assert(value!=null) : s+"\n"+keySet()+"\n"+altKeySet()+"\n";
		return value;
	}
	
	public Scaffold getScaffold(SamLine sl){
		return getScaffold(sl.rnameS());
	}
	
	public int getCoverage(Var v){
		Scaffold scaf=getScaffold(v.scafnum);
		return scaf.calcCoverage(v);
	}
	
	public long lengthSum() {
		long sum=0;
		for(Scaffold scaf : list){
			sum+=scaf.length;
		}
		return sum;
	}
	
	public void clearCoverage() {
		for(Scaffold scaf : list){scaf.clearCoverage();}
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		for(Scaffold sc : list){sb.append(sc).append('\n');}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	final ArrayList<Scaffold> list=new ArrayList<Scaffold>();
	final HashMap<String, Scaffold> map=new HashMap<String, Scaffold>();
	private final HashMap<String, Scaffold> alt=new HashMap<String, Scaffold>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ScafMap defaultScafMap=null;
	private static String defaultScafMapFile=null;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static boolean TRIM_WHITESPACE_ALSO=true;
	
}
