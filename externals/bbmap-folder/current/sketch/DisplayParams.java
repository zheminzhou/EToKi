package sketch;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;

import json.JsonObject;
import shared.Colors;
import shared.Tools;
import structures.ByteBuilder;
import tax.PrintTaxonomy;
import tax.TaxFilter;
import tax.TaxNode;
import tax.TaxTree;

public class DisplayParams implements Cloneable {
	
	@Override
	public DisplayParams clone(){
		try {
			DisplayParams copy=(DisplayParams)super.clone();
			if(taxFilter!=null){
				copy.taxFilter=taxFilter.deepCopy();
			}
			return copy;
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	public DisplayParams parseDoubleHeader(String s){
		if(!s.startsWith("##")){return this;}
		StringBuilder sb=new StringBuilder();
		for(int i=2; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='\n'){break;}
			sb.append(c);
		}
		return parseDoubleHeaderLine(sb.toString());
	}
	
	public DisplayParams parseDoubleHeaderLine(String line) {
		if(line.startsWith("##")){line=line.substring(2);}
		else{assert(!line.startsWith("#")) : line;}
		if(line.length()<1){return this;}
		
		DisplayParams params=this.clone();
		
		String[] args=line.split(" ");
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;} //Normally handled by PreParser, but not in this case.
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			boolean x=params.parse(arg, a, b);
//			assert(x) : "Unknown parameter "+arg+"\n"+line;
			if(!x){System.err.println("Warning: Unknown parameter "+arg);}
		}
		if(SketchObject.verbose2){System.err.println("Made it to post-parse.  taxFilter="+params.taxFilter);}
		params.postParse(true);
		if(SketchObject.verbose2){System.err.println("Passed post-parse.  taxFilter="+params.taxFilter);}
		
		return params;
	}
	
	public boolean parse(String arg, String a, String b){
	
		if(a.equals("chunk")){
			chunkNum=Integer.parseInt(b);
		}else if(a.equals("minhits")  || a.equals("hits")){
			minHits=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minwkid") || a.equalsIgnoreCase("wkid")){
			minWKID=Float.parseFloat(b);
			if(minWKID>1){minWKID/=100;}
			assert(minWKID<=1) : "minWKID should between 0 and 1";
		}else if(a.equalsIgnoreCase("minid") || a.equalsIgnoreCase("id") || a.equalsIgnoreCase("minani") || a.equalsIgnoreCase("ani")){
			minANI=Float.parseFloat(b);
			if(minANI>1){minANI/=100;}
			assert(minANI<=1) : "minANI should between 0 and 1";
			if(minANI>0){
				minWKID=(float)Tools.max(minWKID, Comparison.aniToWkid(minANI, 32));//Lowest possible minWKID for this ANI
			}
		}else if(a.equals("minbases")){
			minBases=Integer.parseInt(b);
		}else if(a.equals("minsizeratio")){
			minSizeRatio=Float.parseFloat(b);
//			assert(minSizeRatio>=0f && minSizeRatio<=1.0f) : "\nminSizeRatio must be between 0 and 1, inclusive.\n";
			if(minSizeRatio>1){minSizeRatio=1f/minSizeRatio;}
		}else if(a.equals("records") || a.equals("maxrecords") || a.equals("results")){
			maxRecords=Integer.parseInt(b);
			assert(maxRecords>=1) : "Max records must be at least 1.";
		}else if(a.equals("format")){
			assert(b!=null) : "Invalid format: "+arg;
			if(Tools.isDigit(b.charAt(0))){
				format=Integer.parseInt(b);
			}else if(b.equalsIgnoreCase("json")){
				format=FORMAT_JSON;
			}else if(b.equalsIgnoreCase("constellation")){
				format=FORMAT_CONSTELLATION;
			}else{
				assert(false) : "Invalid format: "+arg;
			}
		}else if(a.equalsIgnoreCase("json")){
			if(Tools.parseBoolean(b)){
				format=FORMAT_JSON;
			}else{
				if(format==FORMAT_JSON){format=default_format;}
			}
		}else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
			taxLevel=TaxTree.parseLevel(b);//TODO: Change to extended
		}
		
		else if(a.equalsIgnoreCase("printtax") || a.equalsIgnoreCase("printtaxa")){
			printTax=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printoriginalname") || a.equalsIgnoreCase("printseqname") || a.equalsIgnoreCase("printname0") || a.equals("pn0")){
			printOriginalName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printfilename") || a.equalsIgnoreCase("printfname")){
			printFileName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printimg")){
			printImg=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printcompleteness") || a.equalsIgnoreCase("completeness") || a.equalsIgnoreCase("printcomplt")){
			printCompleteness=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printani") || a.equalsIgnoreCase("ani")){
			printAni=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printscore") || a.equalsIgnoreCase("score")){
			printScore=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printevalue") || a.equalsIgnoreCase("evalue")){
			printEValue=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("trackcounts")){
			trackCounts=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth") || a.equalsIgnoreCase("depth")){
			printDepth=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth2") || a.equalsIgnoreCase("depth2")){
			printDepth2=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("actualdepth") || a.equalsIgnoreCase("printactualdepth")){
			printActualDepth=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printvolume") || a.equalsIgnoreCase("volume")){
			printVolume=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printavgrefhits") || a.equalsIgnoreCase("printrefhits") || a.equalsIgnoreCase("avgrefhits") || a.equalsIgnoreCase("refhits")){
			printRefHits=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sortByDepth")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.depthComparator;}
		}else if(a.equalsIgnoreCase("sortByDepth2")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.depth2Comparator;}
		}else if(a.equalsIgnoreCase("sortByVolume")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.volumeComparator;}
		}else if(a.equalsIgnoreCase("sortByScore")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.scoreComparator;}
		}
		
		else if(a.equalsIgnoreCase("printUMatches") || a.equalsIgnoreCase("printUHits") || a.equalsIgnoreCase("printUnique")){
			printUnique=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches2") || a.equalsIgnoreCase("printUnique2") || a.equalsIgnoreCase("unique2")){
			printUnique2=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches3") || a.equalsIgnoreCase("printUnique3") || a.equalsIgnoreCase("unique3")){
			printUnique3=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUContam")){
			printUContam=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printNoHit")){
			printNoHit=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits") || a.equalsIgnoreCase("contam") || a.equalsIgnoreCase("printcontam")){
			printContam=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits2") || a.equalsIgnoreCase("contam2") || a.equalsIgnoreCase("printcontam2")){
			if(b==null || b.length()<1){
				printContam2=true;
			}else if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}else{
				printContam2=Tools.parseBoolean(b);
			}
		}else if(a.equalsIgnoreCase("contamLevel")){
			if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}
		}
		
		else if(a.equalsIgnoreCase("reportAniOnly") || a.equalsIgnoreCase("AniOnly")){
			reportAniOnly=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printMatches")){
			printMatches=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printLength")){
			printLength=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxID")){
			printTaxID=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSize")){
			printGSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("gSizeKMG")){
			gSizeKMG=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGKmers")){
			printGKmers=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxName")){
			printTaxName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSeqs")){
			printGSeqs=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGBases")){
			printGBases=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("minEntropy") || a.equalsIgnoreCase("entropy") || a.equalsIgnoreCase("efilter")){
			minEntropy=Float.parseFloat(b);
		}
		
		else if(a.equalsIgnoreCase("printColors") || a.equalsIgnoreCase("colors") || a.equalsIgnoreCase("color")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(b==null || b.length()<1){
				printColors=true;
			}else if(b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
				printColors=true;
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")){
				printColors=false;
			}else{
				printColors=true;
				if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
					colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				}else{
					colorLevel=TaxTree.stringToLevelExtended(b);
				}
			}
			setColors=true;
//			System.err.println("Parsed "+arg); //123
		}else if(a.equalsIgnoreCase("colorLevel")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
			}else{
				colorLevel=TaxTree.stringToLevelExtended(b);
			}
//			System.err.println("Parsed "+arg); //123
		}
		
		else if(a.equalsIgnoreCase("printRefDivisor") || a.equalsIgnoreCase("printRDiv")){
			printRefDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQueryDivisor") || a.equalsIgnoreCase("printQDiv")){
			printQueryDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printRefSize") || a.equalsIgnoreCase("printRSize")){
			printRefSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQuerySize") || a.equalsIgnoreCase("printQSize")){
			printQuerySize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printContamHits") || a.equalsIgnoreCase("printCHits")){
			printContamHits=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printIntersection") || a.equalsIgnoreCase("intersection") || a.equalsIgnoreCase("intersect")){
			printIntersection=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("mergePairs") || a.equalsIgnoreCase("merge")){
			mergePairs=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printAll")){
			if(Tools.parseBoolean(b)){
				setPrintAll();
			}
		}
		
		else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
		}else if(a.equals("reads")){
			reads=Tools.parseKMG(b);
		}else if(a.equals("mode") || a.equalsIgnoreCase("single") || a.equalsIgnoreCase("singlesketch") || a.equalsIgnoreCase("onesketch")
				|| a.equalsIgnoreCase("persequence") || a.equalsIgnoreCase("sequence") || a.equalsIgnoreCase("pertaxa") 
				|| a.equalsIgnoreCase("perheader") || a.equalsIgnoreCase("perfile")){
			mode=SketchObject.parseMode(arg, a, b);
		}
		
		//For format 3
		else if(a.equalsIgnoreCase("useTaxidName") || a.equalsIgnoreCase("useTaxidAsName")){
			useTaxidName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useImgName") || a.equalsIgnoreCase("useImgAsName")){
			useImgName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useTaxName") || a.equalsIgnoreCase("useTaxAsName")){
			useTaxName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useFilePrefixName") || a.equalsIgnoreCase("useFilePrefixAsName")){
			useFilePrefixName=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("taxfilter") || a.equalsIgnoreCase("taxfilterset")){
			if(b==null){taxFilter=null;}
			else{
				if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
				taxFilter.clearSet();
				taxFilter.makeSet();
				taxFilter.addNumbers(b, false);
//				System.err.println("A:\t"+this);
			}
		}else if(a.equalsIgnoreCase("taxfilterlevel")){//TODO:  Change to extended
			int temp=TaxTree.parseLevel(b);
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setLevel(temp, false);
//			System.err.println("B:\t"+this);
		}else if(a.equalsIgnoreCase("taxfilterinclude") || a.equalsIgnoreCase("taxfilterexclude")){
			boolean temp=Tools.parseBoolean(b);
			if(a.equalsIgnoreCase("taxfilterexclude")){temp=!temp;}
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setInclude(temp);
//			System.err.println("C:\t"+this);
		}else if(a.equalsIgnoreCase("taxfilterstring")){
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setContainsString(b);
//			System.err.println("D:\t"+this);
		}
		
		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount") || a.equalsIgnoreCase("mincount") || a.equalsIgnoreCase("minKeyOccuranceCount")){
			minKeyOccuranceCount=Tools.max(1, Integer.parseInt(b));
		}
		
		//TODO:  Eventually remove support for "amino" and "k" and just support "hamino" and "hk"
		//This stands for "header amino" and "header k".
		
		//Parameters for compatibility verification
		else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("hk")){
//			System.err.println("A: k="+k+", k2="+k2+", arg="+arg);
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "\nBad argument "+arg+"\n"+b+"\n";
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
				if(k==k2){k2=0;}
//				System.err.println("B: k="+k+", k2="+k2+", split="+Arrays.toString(split));
			}else{
				k=Integer.parseInt(b);
//				System.err.println("C: k="+k+", k2="+k2);
			}
		}else if(a.equalsIgnoreCase("hashversion") || a.equalsIgnoreCase("hv")){
			hashVersion=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("amino") || a.equalsIgnoreCase("hamino")){
			amino=Tools.parseBoolean(b);
			if(amino){translate=false;}
		}else if(a.equalsIgnoreCase("translate")){
			translate=Tools.parseBoolean(b);
			if(translate){amino=false;}
		}else if(a.equalsIgnoreCase("sixframes")){
			sixframes=Tools.parseBoolean(b);
			if(sixframes){amino=false; translate=true;}
		}
		
		else if(a.equalsIgnoreCase("requiredmeta") || a.equalsIgnoreCase("rmeta")){
			if(b==null){requiredMeta=null;}
			else{
				String[] split2=b.split(",");
				requiredMeta=new ArrayList<String>(split2.length);
				for(String mt : split2){
					assert(mt.indexOf(':')>=0) : "Metadata tags must contain ':' symbol: "+mt;
					requiredMeta.add(mt);
				}
			}
		}else if(a.equalsIgnoreCase("bannedmeta") || a.equalsIgnoreCase("bmeta")){
			if(b==null){bannedMeta=null;}
			else{
				String[] split2=b.split(",");
				bannedMeta=new ArrayList<String>(split2.length);
				for(String mt : split2){
					assert(mt.indexOf(':')>=0) : "Metadata tags must contain ':' symbol: "+mt;
					bannedMeta.add(mt);
				}
			}
		}
		
//		else if(a.equalsIgnoreCase("requiredtaxid") || a.equalsIgnoreCase("rtaxid")){
//			if(b==null){requiredTaxid=null;}
//			else{
//				String[] split2=b.split(",");
//				requiredTaxid=new IntList(split2.length);
//				for(String mt : split2){
//					requiredTaxid.add(Integer.parseInt(mt));
//				}
//				if(requiredTaxid.isEmpty()){requiredTaxid=null;}
//			}
//		}else if(a.equalsIgnoreCase("bannedtaxid") || a.equalsIgnoreCase("btaxid")){
//			if(b==null){bannedTaxid=null;}
//			else{
//				String[] split2=b.split(",");
//				bannedTaxid=new IntList(split2.length);
//				for(String mt : split2){
//					bannedTaxid.add(Integer.parseInt(mt));
//				}
//				if(bannedTaxid.isEmpty()){bannedTaxid=null;}
//			}
//		}
		
		else if(a.equalsIgnoreCase("requiredmetaand") || a.equalsIgnoreCase("rmetaand")){
			requiredMetaAnd=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("requiredmetaor") || a.equalsIgnoreCase("rmetaor")){
			requiredMetaAnd=!Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("bbversion")){
			inputVersion=b;
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public void postParse(boolean requireTree){
		assert(!postParsed);
		synchronized(this){
			if(postParsed){return;}
			if(taxFilter!=null){
				if(taxFilter.size()==0 && taxFilter.containsString()==null){
					System.err.println("Eliminating empty TaxFilter.");
					taxFilter=null;
				}
			}
			if(taxFilter!=null && requireTree){
				assert(SketchObject.taxtree!=null) : "No taxtree loaded.";
				taxFilter.setTree(SketchObject.taxtree);
				taxFilter.promote();
			}
			postParsed=true;
		}
	}
	
	@Override
	public String toString(){
		return toString(-1);
	}
	
	public String toString(int chunkNum){
		StringBuilder sb=new StringBuilder();
		sb.append("##");
		sb.append("hits=").append(minHits);
		if(chunkNum>=0){sb.append(" chunk=").append(chunkNum);}
		sb.append(" wkid=").append(String.format(Locale.ROOT, "%.5f",minWKID));
		if(minANI>0){sb.append(" id=").append(String.format(Locale.ROOT, "%.5f",minANI));}
		if(minBases>0){sb.append(" minbases=").append(minBases);}
		if(minSizeRatio>0){sb.append(" minsizeratio=").append(String.format(Locale.ROOT, "%.5f",minSizeRatio));}
		sb.append(" records=").append(maxRecords);
		sb.append(" format=").append(format);
		sb.append(" level=").append(taxLevel);
		if(inputVersion!=null){sb.append(" bbversion=").append(inputVersion);}
		
		if(k!=SketchObject.defaultK || k2!=0 || k!=SketchObject.k || k2!=SketchObject.k2){
			assert(k>0 && k2>=0 && k2<k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			assert(SketchObject.k>0 && SketchObject.k2>=0 && SketchObject.k2<SketchObject.k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			sb.append(" hk=").append(SketchObject.k).append(',').append(SketchObject.k2);
		}
		if(SketchObject.amino){sb.append(" hamino=").append(SketchObject.amino);} //TODO: This conflicts with Parser flag
		if(SketchObject.translate){sb.append(" translate=").append(SketchObject.translate);}
		if(SketchObject.sixframes){sb.append(" sixframes=").append(SketchObject.sixframes);}
		if(SketchObject.HASH_VERSION>1){sb.append(" hashversion=").append(SketchObject.HASH_VERSION);}
		
		if(true || printTax!=default_printTax){sb.append(" printTax=").append(printTax);}
		if(true || printOriginalName!=default_printOriginalName){sb.append(" pn0=").append(printOriginalName);}
		if(true || printFileName!=default_printFileName){sb.append(" printfname=").append(printFileName);}
		if(true || printImg!=default_printImg){sb.append(" printImg=").append(printImg);}
		if(true || printAni!=default_printAni){sb.append(" printAni=").append(printAni);}
		if(true || printCompleteness!=default_printCompleteness){sb.append(" printCompleteness=").append(printCompleteness);}

		if(true || printUnique!=default_printUnique){sb.append(" printUMatches=").append(printUnique);}
		if(true || printUnique2!=default_printUnique2){sb.append(" printUnique2=").append(printUnique2);}
		if(true || printUnique3!=default_printUnique3){sb.append(" printUnique3=").append(printUnique3);}
		if(true || printUContam!=default_printUContam){sb.append(" printUContam=").append(printUContam);}
		if(true || printNoHit!=default_printNoHit){sb.append(" printNoHit=").append(printNoHit);}
		if(true || printContam!=default_printContam){sb.append(" contam=").append(printContam);}
		if(true){sb.append(" contam2=").append(printContam2 ? TaxTree.extendedToLevel(contamLevel)+"" : "f");}

		if(true || printScore!=default_printScore){sb.append(" printScore=").append(printScore);}
		if(true || printEValue!=default_printEValue){sb.append(" printEValue=").append(printEValue);}
		
		if(true || printDepth!=default_printDepth){sb.append(" printDepth=").append(printDepth);}
		if(true || printDepth2!=default_printDepth2){sb.append(" printDepth2=").append(printDepth2);}
		if(true || printActualDepth!=default_printActualDepth){sb.append(" printActualDepth=").append(printActualDepth);}
		if(true || printVolume!=default_printVolume){sb.append(" printVolume=").append(printVolume);}
		if(true || printRefHits!=default_printRefHits){sb.append(" printRefHits=").append(printRefHits);}
		
		if(true || printMatches!=default_printMatches){sb.append(" printMatches=").append(printMatches);}
		if(true || printLength!=default_printLength){sb.append(" printLength=").append(printLength);}
		if(true || printTaxID!=default_printTaxID){sb.append(" printTaxID=").append(printTaxID);}
		if(true || printGSize!=default_printGSize){sb.append(" printGSize=").append(printGSize);}
		if(true || gSizeKMG!=default_gSizeKMG){sb.append(" gSizeKMG=").append(gSizeKMG);}
		if(true || printGKmers!=default_printGKmers){sb.append(" printGKmers=").append(printGKmers);}
		if(true || printTaxName!=default_printTaxName){sb.append(" printTaxName=").append(printTaxName);}
		if(true || printGSeqs!=default_printGSeqs){sb.append(" printGSeqs=").append(printGSeqs);}
		if(true || printGBases!=default_printGBases){sb.append(" printGBases=").append(printGBases);}
		if(true || minEntropy!=default_minEntropy){sb.append(" minEntropy=").append(String.format("%.4f", minEntropy));}
		if(comparator!=Comparison.scoreComparator){sb.append(" ").append(comparator.toString());}
		
		if(taxFilter!=null){
			sb.append(" taxfilterlevel=").append(taxFilter.taxLevel());
			sb.append(" taxfilterinclude=").append(taxFilter.include());
			if(taxFilter.containsString()!=null){
				sb.append(" taxfilterstring=").append(taxFilter.containsString());
			}
			Integer[] temp=taxFilter.taxSet();
			sb.append(" taxfilterset=");
			if(temp!=null){
				for(int i=0; i<temp.length; i++){
					if(i>0){sb.append(',');}
					sb.append(temp[i]);
				}
			}
		}
		
		if(useTaxidName){sb.append(" useTaxidName=").append(useTaxidName);}
		if(useImgName){sb.append(" useImgName=").append(useImgName);}
		if(useTaxName){sb.append(" useTaxName=").append(useTaxName);}
		
		if(true){sb.append(" colors=").append(printColors ? TaxTree.extendedToLevel(colorLevel)+"" : "f");}
		
		if(minKeyOccuranceCount!=default_minKeyOccuranceCount){sb.append(" minKeyOccuranceCount=").append(minKeyOccuranceCount);}
		
//		if(printColors && colorLevel!=default_colorLevel){sb.append(" colorLevel=").append(TaxTree.extendedToLevel(colorLevel));}
		

		if(printRefDivisor){sb.append(" printRefDivisor=").append(printRefDivisor);}
		if(printQueryDivisor){sb.append(" printQueryDivisor=").append(printQueryDivisor);}
		if(printRefSize){sb.append(" printRefSize=").append(printRefSize);}
		if(printQuerySize){sb.append(" printQuerySize=").append(printQuerySize);}
		if(printContamHits){sb.append(" printContamHits=").append(printContamHits);}
		if(printIntersection){sb.append(" printIntersection=").append(printIntersection);}
		if(mergePairs){sb.append(" mergePairs=").append(mergePairs);}
		
		if(reads>-1){sb.append(" reads=").append(reads);}
		if(mode!=default_mode){sb.append(" mode=").append(mode);}
		if(samplerate!=default_samplerate){sb.append(" samplerate=").append(String.format(Locale.ROOT, "%.4f",samplerate));}

		if(!requiredMetaAnd){sb.append(" requiredmetaand="+requiredMetaAnd);}
		if(requiredMeta!=null && !requiredMeta.isEmpty()){
			sb.append(" rmeta=");
			for(String s : requiredMeta){
				sb.append(s);
				sb.append(',');
			}
			sb.setLength(sb.length()-1);
		}
		if(bannedMeta!=null && !bannedMeta.isEmpty()){
			sb.append(" bmeta=");
			for(String s : bannedMeta){
				sb.append(s);
				sb.append(',');
			}
			sb.setLength(sb.length()-1);
		}
//		if(requiredTaxid!=null && !requiredTaxid.isEmpty()){
//			sb.append(" rtaxid=");
//			for(int i=0; i<requiredTaxid.size; i++){
//				sb.append(requiredTaxid.get(i));
//				sb.append(',');
//			}
//			sb.setLength(sb.length()-1);
//		}
//		if(bannedTaxid!=null && !bannedTaxid.isEmpty()){
//			sb.append(" btaxid=");
//			for(int i=0; i<bannedTaxid.size; i++){
//				sb.append(bannedTaxid.get(i));
//				sb.append(',');
//			}
//			sb.setLength(sb.length()-1);
//		}
		
		sb.append('\n');
		return sb.toString();
	}
	
	public boolean compatible(){
		return SketchObject.k==k && SketchObject.k2==k2 && SketchObject.aminoOrTranslate()==aminoOrTranslate() && hashVersion==SketchObject.HASH_VERSION;
	}
	
	public void setPrintAll(){
		printTax=true;
		printOriginalName=true;
		printImg=true;
		printAni=true;
		printCompleteness=true;
		printScore=true;
		printEValue=true;
		printDepth=true;
		printDepth2=true;
		printVolume=true;
		printRefHits=true;
		
		printMatches=true;
		printLength=true;
		printTaxID=true;
		printGSize=true;
		printGKmers=true;
		printTaxName=true;
		printGSeqs=true;
		printGBases=true;
		
//		printColors=true;

		printUnique=true;
		printUnique2=true;
		printUnique3=true;
		printUContam=true;
		printNoHit=true;
		printContam=true;
		printContam2=true;
		
		printRefDivisor=true;
		printQueryDivisor=true;
		printRefSize=true;
		printQuerySize=true;
		printContamHits=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             JSON             ----------------*/
	/*--------------------------------------------------------------*/
	
	public JsonObject toJson(SketchResults sr){
		JsonObject j=toJson(sr.sketch);
		if(sr.list!=null){
			for(Comparison c : sr.list){
				JsonObject jc=toJson(c);
				j.add(c.name(), jc);
			}
		}
		return j;
	}

	public JsonObject toJson(Sketch sk){
		assert(format==FORMAT_JSON);
		
		JsonObject j=new JsonObject();
		j.add("Name", sk.name());
		if(dbName!=null){j.add("DB", dbName);}
		j.add("SketchLen", sk.length());
		
		j.add("Seqs", sk.genomeSequences);
		j.add("Bases", sk.genomeSizeBases);
		j.add("gSize", sk.genomeSizeEstimate());
		if(sk.probCorrect<1 && sk.probCorrect>0){j.add("Quality", sk.probCorrect);}
		if(sk.counts!=null){
			double d=Tools.averageDouble(sk.counts);
			j.add("AvgCount", d);
			j.add("Depth", d);
		}
		
		if(sk.imgID>0){j.add("IMG", sk.imgID);}
		if(sk.spid>0){j.add("spid", sk.spid);}
		if(sk.taxID>0 && sk.taxID<SketchObject.minFakeID){j.add("TaxID", sk.taxID);}

		if(printFileName && sk.fname()!=null && !sk.fname().equals(sk.name())){j.add("File", sk.fname());}
		if(printOriginalName && sk.name0()!=null && !sk.name0().equals(sk.name())){j.add("SeqName", sk.name0());}
		
		if(sk.meta!=null){
			for(String st : sk.meta){
				int colon=st.indexOf(':');
				j.add(st.substring(0,  colon), st.substring(colon+1));
			}
		}
		return j;
	}
	
	public JsonObject toJson(Comparison c){
		final int tid=c.taxID;
		
		JsonObject j=new JsonObject();
		
		//Text fields
		if(printTaxName){j.add("taxName", c.taxName()==null ? "." : c.taxName());}
		if(printOriginalName){j.add("seqName", c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				j.add("taxonomy", SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical));
			}else{
				j.add("taxonomy", (Object)null);
			}
		}
		
		j.add("WKID", 100*c.wkid());
		j.add("KID", 100*c.kid());
		
		//Primary fields
		if(printAni){j.add((aminoOrTranslate() ? "AAI" : "ANI"), 100*c.ani());}
		if(printCompleteness){j.add("Complt", 100*c.completeness());}
		if(printContam){j.add("Contam", 100*c.contamFraction());}
		if(printContam2){j.add("Contam2", 100*c.contam2Fraction());}
		if(printUContam){j.add("uContam", 100*c.uContamFraction());}
		if(printScore){j.add("Score", c.score());}
		if(printEValue){j.add("E-Val", String.format("%5.2e", c.eValue()));}
		
		if(printDepth){j.add("Depth", c.depth(printActualDepth));}
		if(printDepth2){j.add("Depth2", c.depth2(printActualDepth));}
		if(printVolume){j.add("Volume", c.volume()+0.001);}
		if(printRefHits){j.add("RefHits", c.avgRefHits());}
		
		if(printMatches){j.add("Matches", c.hits());}
		if(printUnique){j.add("Unique", c.uHits());}
		if(printUnique2){j.add("Unique2", c.unique2());}
		if(printUnique3){j.add("Unique3", c.unique3());}
		if(printNoHit){j.add("noHit", c.noHits());}
		if(printLength){j.add("Length", c.maxDivisor());}
		if(printTaxID){j.add("TaxID", tid>=SketchObject.minFakeID ? -1 : tid);}
		if(printImg){j.add("ImgID", String.format(Locale.ROOT, "\t%d", c.imgID()));}
		if(printGBases){j.add("gBases", c.genomeSizeBases());}
		if(printGKmers){j.add("gKmers", c.genomeSizeKmers());}
		if(printGSize){j.add("gSize", c.genomeSizeEstimate());}
		if(printGSeqs){j.add("gSeqs", c.genomeSequences());}
		
		//Raw fields
		if(printRefDivisor){j.add("rDiv", c.refDivisor());}
		if(printQueryDivisor){j.add("qDiv", c.queryDivisor());}
		if(printRefSize){j.add("rSize", c.refSize());}
		if(printQuerySize){j.add("qSize", c.querySize());}
		if(printContamHits){j.add("cHits", c.contamHits());}
		
		if(printIntersection){
			Sketch intersection=Sketch.intersection(c.a, c.b);
			j.add("intersection", intersection.toString());
		}
		
		return j;
	}
	
	public boolean json(){return format==FORMAT_JSON;}
	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	ByteBuilder queryHeader(Sketch sk){
		ByteBuilder sb=new ByteBuilder();
		if(format>2){return sb;}
		
		String color=toColor(sk.taxID);
		if(color!=null){sb.append(color);}
		
		sb.append("\nQuery: ").append(sk.name()==null ? "." : sk.name());
		if(dbName!=null){sb.append("\tDB: ").append(dbName);}
		sb.append("\tSketchLen: ").append(sk.length());
		sb.append("\tSeqs: ").append(sk.genomeSequences).append(' ');
		sb.append("\t"+(aminoOrTranslate() ? "SeqLen" : "Bases")+": ").append(sk.genomeSizeBases);
		sb.append("\tgSize: ").append(sk.genomeSizeEstimate());
		if(sk.probCorrect<1 && sk.probCorrect>0){sb.append("\tQuality: ").append(sk.probCorrect, 4);}
		if(sk.counts!=null){
			double d=Tools.averageDouble(sk.counts);
			sb.append("\tAvgCount: ").append(d, 3);
			sb.append("\tDepth: ").append(Tools.observedToActualCoverage(d), 3);
		}

		if(sk.imgID>0){sb.append("\tIMG: ").append(sk.imgID);}
		if(sk.spid>0){sb.append("\tspid: ").append(sk.spid);}
		if(sk.taxID>0 && sk.taxID<SketchObject.minFakeID){sb.append("\tTaxID: ").append(sk.taxID);}

		if(printFileName && sk.fname()!=null && !sk.fname().equals(sk.name())){sb.append("\tFile: "+sk.fname());}
		if(printOriginalName && sk.name0()!=null && !sk.name0().equals(sk.name())){sb.append("\tSeqName: "+sk.name0());}
		
		if(sk.meta!=null){
			for(String st : sk.meta){
				sb.append("\t").append(st.replaceFirst(":", ": "));
			}
		}
		
		if(color!=null){sb.append(Colors.RESET);}
		
		return sb;
	}
	
	int toColorTid(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return 0;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		return tn==null || tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E) ? 0 : tn.id;
	}
	
	String toColor(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return null;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		if(tn==null){
			return null;
		}else{
			if(tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E)){return Colors.WHITE;}
			else{
//				System.err.println("*"+tn.id+", "+tn.id%Colors.colorArray.length);
				return Colors.colorArray[tn.id%Colors.colorArray.length];
			}
		}
	}
	
	String header(){
		if(format==FORMAT_JSON){return null;}
		final String ani=(aminoOrTranslate() ? "AAI" : "ANI");
		if(format==FORMAT_QUERY_REF_ANI || format==FORMAT_CONSTELLATION){
			if(reportAniOnly){return "#Query\tRef\t"+ani;}
			if(format==FORMAT_QUERY_REF_ANI){return "#Query\tRef\t"+ani+"\tQSize\tRefSize\tQBases\tRBases"+(printTaxID ? "\tQTaxID\tRTaxID" : "");}
			if(format==FORMAT_CONSTELLATION){return "#Query\tRef\tKID\tWKID\t"+ani+"\tCmplt\tQSize\tRefSize\tQBases\tRefBases";}
		}
		
		StringBuilder sb=new StringBuilder();
		
		//Numeric fields
		if(true){sb.append("WKID");}
		if(true){sb.append("\tKID");}
		if(printAni){sb.append("\t"+ani);}
		if(printCompleteness){sb.append("\tComplt");}
		if(printContam){sb.append("\tContam");}
		if(printContam2){sb.append("\tContam2");}
		if(printUContam){sb.append("\tuContam");}
		if(printScore){sb.append("\tScore");}
		if(printEValue){sb.append("\tE-Val\t");}
		
		if(printDepth){sb.append("\tDepth");}
		if(printDepth2){sb.append("\tDepth2");}
		if(printVolume){sb.append("\tVolume");}
		if(printRefHits){sb.append("\tRefHits");}
		if(printMatches){sb.append("\tMatches");}
		if(printUnique){sb.append("\tUnique");}
		if(printUnique2){sb.append("\tUnique2");}
		if(printUnique3){sb.append("\tUnique3");}
		if(printNoHit){sb.append("\tnoHit");}
		if(printLength){sb.append("\tLength");}
		if(printTaxID){sb.append("\tTaxID");}
		if(printImg){sb.append("\tImgID    ");}
		if(printGBases){sb.append("\tgBases");}
		if(printGKmers){sb.append("\tgKmers");}
		if(printGSize){sb.append("\tgSize");}
		if(printGSeqs){sb.append("\tgSeqs");}
		
		
		//Raw fields
		if(printRefDivisor){sb.append("\trDiv");}
		if(printQueryDivisor){sb.append("\tqDiv");}
		if(printRefSize){sb.append("\trSize");}
		if(printQuerySize){sb.append("\tqSize");}
		if(printContamHits){sb.append("\tcHits");}
		
		//Text fields
		if(printTaxName){sb.append("\ttaxName");}
		if(printOriginalName){sb.append("\tseqName");}
		if(printTax && SketchObject.taxtree!=null){sb.append("\ttaxonomy");}
		
		return sb.toString();
	}
	
	void formatComparisonColumnwise(Comparison c, ByteBuilder sb, int prevTid){
		final int tid=c.taxID;
		boolean reset=false;
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
//		sb.append(String.format(Locale.ROOT, "%.2f%%\t%.2f%%", 100*c.idMinDivisor(), 100*c.idMaxDivisor()));
		sb.append(100*c.wkid(), 2).append('%').append('\t');
		sb.append(100*c.kid(), 2).append('%');
		
//		if(printAni){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.ani()));}
//		if(printCompleteness){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.completeness()));}
//		if(printContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.contamFraction()));}
//		if(printContam2){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.contam2Fraction()));}
//		if(printUContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.uContamFraction()));}
		
		if(printAni){sb.append('\t').append(100*c.ani(), 2).append('%');}
		if(printCompleteness){sb.append('\t').append(100*c.completeness(), 2).append('%');}
		if(printContam){sb.append('\t').append(100*c.contamFraction(), 2).append('%');}
		if(printContam2){sb.append('\t').append(100*c.contam2Fraction(), 2).append('%');}
		if(printUContam){sb.append('\t').append(100*c.uContamFraction(), 2).append('%');}
		if(printScore){sb.append('\t').append(c.scoreS());}
		if(printEValue){sb.append('\t').append(String.format("%5.2e", c.eValue()));}
		
		if(printDepth){sb.append('\t').append(c.depthS(printActualDepth));}
		if(printDepth2){sb.append('\t').append(c.depth2S(printActualDepth));}
		if(printVolume){sb.append('\t').append(c.volumeS());}
		if(printRefHits){sb.append('\t').append(c.avgRefHitsS());}
		
		if(printMatches){sb.append('\t').append(c.hits());}
		if(printUnique){sb.append('\t').append(c.uHits());}
		if(printUnique2){sb.append('\t').append(c.unique2());}
		if(printUnique3){sb.append('\t').append(c.unique3());}
		if(printNoHit){sb.append('\t').append(c.noHits());}
		if(printLength){sb.append('\t').append( c.maxDivisor());}
		if(printTaxID){sb.append('\t').append(tid>=SketchObject.minFakeID ? -1 : tid);}
		if(printImg){sb.append('\t').append(c.imgID());}
		if(printGBases){sb.append('\t').append(c.genomeSizeBases());}
		if(printGKmers){sb.append('\t').append(c.genomeSizeKmers());}
		if(printGSize){
			long size=c.genomeSizeEstimate();
			if(gSizeKMG){
				sb.append('\t').append(toKMG(size));
			}else{
				sb.append('\t').append(size);
			}
		}
		if(printGSeqs){sb.append('\t').append(c.genomeSequences());}
		
		//Raw fields
		if(printRefDivisor){sb.append('\t').append(c.refDivisor());}
		if(printQueryDivisor){sb.append('\t').append(c.queryDivisor());}
		if(printRefSize){sb.append('\t').append(c.refSize());}
		if(printQuerySize){sb.append('\t').append(c.querySize());}
		if(printContamHits){sb.append('\t').append(c.contamHits());}
		
		//Text fields
		if(printTaxName){sb.append('\t').append(c.taxName()==null ? "." : c.taxName());}
		if(printOriginalName){sb.append('\t').append(c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			sb.append('\t');
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				sb.append(SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical));
			}else{
				sb.append('.');
			}
		}
		if(printTaxName && !printOriginalName && c.taxName()==null && c.name0()!=null){sb.append('\t').append(c.name0());} //Extra column
		
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
		
		if(printIntersection){
			Sketch intersection=Sketch.intersection(c.a, c.b);
			sb.append(intersection.toString());
			sb.append('\n');
		}
		
	}
	
	String toKMG(long value){
		if(value<10000000L){return Long.toString(value);}
		value+=5;
		if(value<1000000000L){return value/1000L+"K";}
		if(value<1000000000000L){return value/1000000L+"M";}
		if(value<1000000000000000L){return value/1000000000L+"G";}
		return value/1000000000000L+"T";
	}
	
	void formatComparison3Column(Comparison c, ByteBuilder sb, int prevTid){
		Sketch query=c.a;
		final long sea=Tools.max(1, c.a.genomeSizeEstimate());
		final long seb=Tools.max(1, c.b.genomeSizeEstimate());
		final long ba=Tools.max(1, c.a.genomeSizeBases);
		final long bb=Tools.max(1, c.b.genomeSizeBases);
		final String qName=format==FORMAT_CONSTELLATION ? (useFilePrefixName ? query.filePrefix() : ""+query.sketchID) : useTaxidName ? ""+query.taxID :
			useImgName ?  ""+query.imgID : useTaxName ? query.taxName() : query.name();
		final String rName=format==FORMAT_CONSTELLATION ? (useFilePrefixName ? c.b.filePrefix() : ""+c.b.sketchID) : useTaxidName ? ""+c.taxID() :
			useImgName ?  ""+c.imgID() : useTaxName ? c.taxName() : c.name();
		final int tid=c.taxID;
		boolean reset=false;
		
		sb.append(qName).append('\t');
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}

//		sb.append(rName).append(String.format(Locale.ROOT, "\t%.2f\t%.3f", 100*c.ani(), sea/(float)seb));
//		sb.append(rName).append(String.format(Locale.ROOT, "\t%.2f\t%d\t%d\t%d", 100*c.ani(), sea, seb, ba));
		
		//"#Query\tRef\tKID\tWKID\tANI\tCmplt\tQSize\tRefSize\tQBases\tRefBases";

		float kid=100*c.kid();
		float wkid=100*c.wkid();
		float ani=100*c.ani();
		float complt=100*c.completeness();
		
		sb.append(rName).append('\t');
		if(reportAniOnly){
			sb.append(ani, 2).append('\t');
		}else if(format==FORMAT_CONSTELLATION){
			sb.append(kid, 2).append('\t');
			sb.append(wkid, 2).append('\t');
			sb.append(ani, 2).append('\t');
			sb.append(complt, 2).append('\t');
			sb.append(sea).append('\t');
			sb.append(seb).append('\t');
//			sb.append(ba).append('\t');
//			sb.append(bb).append('\t');
		}else{
			sb.append(ani, 2).append('\t');
			sb.append(sea).append('\t');
			sb.append(seb).append('\t');
			sb.append(ba).append('\t');
			sb.append(bb).append('\t');
			if(printTaxID){sb.append(c.a.taxID).append('\t');}
			if(printTaxID){sb.append(c.b.taxID).append('\t');}
		}
		sb.setLength(sb.length()-1);
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
		
//		System.err.println(sb);
	}
	
	void formatComparison(Comparison c, ByteBuilder sb, int prevTaxID){
		if(format==FORMAT_MULTICOLUMN){
			formatComparisonColumnwise(c, sb, prevTaxID);
			return;
		}else if(format==FORMAT_QUERY_REF_ANI || format==FORMAT_CONSTELLATION){
			formatComparison3Column(c, sb, prevTaxID);
			return;
		}
		String complt=(printCompleteness ? String.format(Locale.ROOT, "\tcomplt %.2f%%%%", 100*c.completeness()) : "");
		String contam=(printContam ? String.format(Locale.ROOT, "\tcontam %.2f%%%%", 100*c.contamFraction()) : "");
//		String score=(printScore ? String.format(Locale.ROOT, "\tscore %.2f", c.score2()) : "");
		String score=(printScore ? "\tscore "+c.scoreS() : "");
		String depth=(printDepth ? "\tdepth "+c.depthS(printActualDepth) : "");
		String depth2=(printDepth2 ? "\tdepth2 "+c.depth2S(printActualDepth) : "");
		String volume=(printVolume ? "\tvolume "+c.volumeS() : "");
		String ccs=complt+contam+score;
		
		if(format==FORMAT_OLD){
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d",
					100*c.wkid(), 100*c.kid(), c.hits(), c.minDivisor())+"\ttaxID "+c.taxID()+
					(printImg ? "\timgID "+c.imgID() : "")+"\tgKmers "+c.genomeSizeKmers()+"\t"+
					(c.taxName()==null ? "." : c.taxName())+
					((printOriginalName || (c.taxName()==null && c.name0()!=null)) ? "\t"+(c.name0()==null ? "." : c.name0()) : "")+"\n");
			if(printTax && SketchObject.taxtree!=null){
				if(c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
					TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
					if(tn!=null){
						PrintTaxonomy.printTaxonomy(tn, sb, SketchObject.taxtree, TaxTree.DOMAIN, SketchObject.skipNonCanonical);
					}
				}
				sb.append('\n');
			}
		}else{
			ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
			if(SketchObject.taxtree!=null && c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
				TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
				while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
					tnl.add(tn);
					tn=SketchObject.taxtree.getNode(tn.pid);
				}
			}
			
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d\t",
					100*c.wkid(), 100*c.kid(), c.hits(), c.minDivisor()));
			sb.append("\ttaxID ").append(c.taxID()).append('\t');
			if(printImg){sb.append("\timgID ").append(c.imgID()).append('\t');}
			sb.append(c.taxName()).append('\t');
			if(printOriginalName || (c.taxName()==null && c.name0()!=null)){sb.append(c.name0()).append('\t');}
			
			if(printTax){
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
			}
			sb.append('\n');
			
			tnl.clear();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//These are shared with SketchObject
	//They do not affect anything and are just for the server to validate remote settings.
	private int hashVersion=SketchObject.HASH_VERSION;
	private int k=SketchObject.k;
	private int k2=SketchObject.k2;
	boolean amino=SketchObject.amino;
	boolean translate=SketchObject.translate;
	boolean sixframes=SketchObject.sixframes;
	private boolean aminoOrTranslate(){return amino | translate;}
	
	boolean postParsed=false;
	
	boolean amino(){return amino;}
	
	//These are unique
	public int maxRecords=default_maxRecords;
	public float minANI=0;
	public int minBases=0;
	public float minSizeRatio=0;
	public float minWKID=default_minWKID;
	public int format=default_format;
	
	/** For tracking unique SendSketch queries */
	public int chunkNum=-1;
	public int minHits=default_minHits;
	public int taxLevel=default_taxLevel;
	public int mode=default_mode;
	public float samplerate=default_samplerate;
	public long reads=default_reads;
	public int minKeyOccuranceCount=default_minKeyOccuranceCount;
	public String inputVersion=null;
	
	public String dbName=null;

	boolean hasMetaFilters(){return requiredMeta!=null || bannedMeta!=null/* || requiredTaxid!=null || bannedTaxid!=null*/;}
	
	boolean requiredMetaAnd=true;
	ArrayList<String> requiredMeta=null;
	ArrayList<String> bannedMeta=null;
//	IntList requiredTaxid=null;
//	IntList bannedTaxid=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Columns        ----------------*/
	/*--------------------------------------------------------------*/
	
	//For format 2
	public boolean printTax=default_printTax;
	public boolean printOriginalName=default_printOriginalName;
	public boolean printFileName=default_printFileName;
	public boolean printImg=default_printImg;
	public boolean printAni=default_printAni;
	public boolean printCompleteness=default_printCompleteness;
	public boolean printScore=default_printScore;
	public boolean printEValue=default_printEValue;

	private boolean trackCounts=default_trackCounts;
	public boolean printDepth=default_printDepth;
	public boolean printDepth2=default_printDepth2;
	public boolean printActualDepth=default_printActualDepth;
	public boolean printVolume=default_printVolume;
	public boolean printRefHits=default_printRefHits;
	
	public boolean printLength=default_printLength;
	public boolean printTaxID=default_printTaxID;
	public boolean printGSize=default_printGSize;
	public boolean gSizeKMG=default_gSizeKMG;
	public boolean printGKmers=default_printGKmers;
	public boolean printTaxName=default_printTaxName;
	public boolean printGSeqs=default_printGSeqs;
	public boolean printGBases=default_printGBases;
	
	public float minEntropy=default_minEntropy;

	public boolean printUnique=default_printUnique;
	public boolean printUnique2=default_printUnique2;
	public boolean printUnique3=default_printUnique3;
	public boolean printUContam=default_printUContam;
	public boolean printNoHit=default_printNoHit;

	public boolean printColors=default_printColors;
	public boolean setColors=false;
	public int colorLevel=default_colorLevel;
	
	/** TODO: Note this is conflated between printing %contam and calculating things based on contam hits. */
	public boolean printContam=default_printContam;
	public boolean printContam2=default_printContam2;
	private int contamLevel=default_contamLevel;
	
	/** Raw fields */
	public boolean printMatches=default_printMatches;
	
	public boolean printRefDivisor=false;
	public boolean printQueryDivisor=false;
	public boolean printRefSize=false;
	public boolean printQuerySize=false;
	public boolean printContamHits=false;
	
	public boolean mergePairs=false;
	public boolean printIntersection=false;
	
	//For format 3 or 5
	public boolean useTaxidName=false;
	public boolean useImgName=false;
	public boolean useTaxName=false;
	public boolean useFilePrefixName=false;
	public boolean reportAniOnly=false;
	
	public TaxFilter taxFilter=null;

	/** Make sure the settings are consistent, for CompareSketch.
	 * This is not yet complete. */
	public boolean checkValid(){
		if(printUnique2 || printUnique3){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		if(printContam2){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		return true;
	}
	
	public boolean trackCounts() {
		return trackCounts || printDepth || printDepth2 || printVolume || comparator!=Comparison.scoreComparator; //|| minKeyOccuranceCount>1;
	}
	
	public boolean needContamCounts() {
		return printContam || printContam2 || printContamHits || printUnique || printUnique2 || printUnique3 || printUContam || printNoHit; // || true
	}
	
	public boolean needIndex(){
		return printContam2 || printUnique2 || printUnique3;
	}

	public int contamLevel() {
		return needIndex() ? contamLevel : -1;
	}
	
	public int compare(Comparison a, Comparison b){
		return comparator.compare(a, b);
	}
	
	public Comparator<Comparison> comparator=Comparison.scoreComparator;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int FORMAT_OLD=0, FORMAT_MULTICOLUMN=2, FORMAT_QUERY_REF_ANI=3, FORMAT_JSON=4, FORMAT_CONSTELLATION=5;
	
	public static final int default_maxRecords=20;
	public static final float default_minWKID=0.0001f;
	public static final int default_format=FORMAT_MULTICOLUMN;
	public static final boolean default_printTax=false;
	public static final boolean default_printOriginalName=false;
	public static final boolean default_printFileName=false;
	public static final boolean default_printImg=false;
	public static final boolean default_printAni=true;
	public static final boolean default_printCompleteness=true;
	public static final boolean default_printScore=false;
	public static final boolean default_printEValue=false;
	
	public static final boolean default_trackCounts=false;
	public static final boolean default_printDepth=false;
	public static final boolean default_printDepth2=false;
	public static final boolean default_printActualDepth=true;
	public static final boolean default_printVolume=false;
	public static final boolean default_printRefHits=false;

	public static final boolean default_printContam=true;
	public static final boolean default_printContam2=false;
	
	public static final boolean default_printMatches=true;
	public static final boolean default_printLength=false;
	public static final boolean default_printTaxID=true;
	public static final boolean default_printGSize=true;
	public static final boolean default_gSizeKMG=true;
	public static final boolean default_printGKmers=false;
	public static final boolean default_printTaxName=true;
	public static final boolean default_printGSeqs=true;
	public static final boolean default_printGBases=false;

	public static final float default_minEntropy=0.66f;
	public static final float default_minEntropy_amino=0.70f;

	public static final boolean default_printUnique=true;
	public static final boolean default_printUnique2=false;
	public static final boolean default_printUnique3=false;
	public static final boolean default_printUContam=false;
	public static final boolean default_printNoHit=false;

	public static final boolean default_printColors=true;
	public static final int default_colorLevel=TaxTree.FAMILY_E;

	public static final int default_taxLevel=TaxTree.SPECIES;
	public static final int default_contamLevel=TaxTree.GENUS_E;
	
	public static final int default_mode=SketchObject.ONE_SKETCH;
	
	public static final int default_minHits=3;
	public static final float default_samplerate=1;
	public static final long default_reads=-1;
	public static final int default_minKeyOccuranceCount=1;
	
}
