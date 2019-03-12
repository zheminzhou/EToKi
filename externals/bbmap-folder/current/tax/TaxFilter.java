package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;

import fileIO.TextFile;
import shared.KillSwitch;
import shared.PreParser;
import shared.Tools;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 30, 2015
 *
 */
public class TaxFilter {
	
	public static void main(String[] args){
		String regex=args[0];
		String s=args[1];
		Pattern regexPattern=Pattern.compile(regex);
		boolean b=regexPattern.matcher(s).matches();
		System.err.println(regex);
		System.err.println(s);
		System.err.println(b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public static TaxFilter makeFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		String names=null;
		String ids=null;

		String giTableFile=null;
		String taxTreeFile=null;
		String accessionFile=null;

		int taxLevelE=-1;
		int reqLevels=0;
		boolean include=false;
		boolean promote=true;
		String regex=null;
		String contains=null;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("table") || a.equals("gi")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("level") || a.equals("taxlevel")){
				taxLevelE=TaxTree.parseLevelExtended(b);
//				System.err.println("Set taxLevelE = "+TaxTree.levelToStringExtended(taxLevelE)); //123
			}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.toLowerCase().split(",");
				reqLevels=0;
				for(String s : split2){
					int level=TaxTree.parseLevel(s);
					reqLevels|=(1<<level);
				}
			}else if(a.equals("name") || a.equals("names")){
				names=b;
			}else if(a.equals("regex")){
				regex=b;
			}else if(a.equals("contains")){
				contains=b;
			}else if(a.equals("printnodesadded") || a.equals("printnodes")){
				printNodesAdded=Tools.parseBoolean(b);
			}else if(a.equals("promote")){
				promote=Tools.parseBoolean(b);
			}else if(a.equals("include")){
				include=Tools.parseBoolean(b);
			}else if(a.equals("exclude")){
				include=!Tools.parseBoolean(b);
			}else if(a.equals("requirepresent")){
				REQUIRE_PRESENT=Tools.parseBoolean(b);
				TaxTree.SHOW_WARNINGS=REQUIRE_PRESENT;
			}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("taxids")){
				ids=b;
			}
		}
		
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		TaxFilter filter=new TaxFilter(giTableFile, taxTreeFile, accessionFile, taxLevelE, reqLevels, include, promote, null, regex, contains);
		filter.addNames(names);
		filter.addNumbers(ids, true);
		return filter;
	}
	
	/**
	 * Constructor.
	 * @param tree_
	 */
	public TaxFilter(TaxTree tree_){
		tree=tree_;
		taxLevelE=-1;
		reqLevels=0;
		include=true;
		promote=true;
		taxSet=null;
		regex=null;
		regexPattern=null;
		containsString=null;
	}
	
	/**
	 * Constructor.
	 */
	public TaxFilter(TaxTree tree_, int taxLevelE_, int reqLevels_, boolean include_, boolean promote_, HashSet<Integer> taxSet_, String regex_, String contains_){
		tree=tree_;
		taxLevelE=taxLevelE_;
		reqLevels=reqLevels_;
		include=include_;
		promote=promote_;
		taxSet=(taxSet_==null ? new HashSet<Integer>() : taxSet_);
		regex=regex_;
		regexPattern=(regex==null ? null : Pattern.compile(regex));
		containsString=contains_;
	}
	
	/**
	 * Constructor.
	 */
	public TaxFilter(String giTableFile, String taxTreeFile, String accessionFile, int taxLevelE_, int reqLevels_,
			boolean include_, boolean promote_, HashSet<Integer> taxSet_, String regex_, String contains_){
		taxLevelE=taxLevelE_;
		reqLevels=reqLevels_;
		include=include_;
		promote=promote_;
		taxSet=(taxSet_==null ? new HashSet<Integer>() : taxSet_);

		tree=loadTree(taxTreeFile);
		loadGiTable(giTableFile);
		loadAccession(accessionFile);
		regex=regex_;
		regexPattern=(regex==null ? null : Pattern.compile(regex));
		containsString=contains_;
	}
	
	/** Alter taxSet and taxLevel to ensure they intersect with the file contents. */
	public void reviseByBestEffort(String fname){
		HashSet<Integer> desired=new HashSet<Integer>();
		int currentLevelE=taxLevelE;
		for(int i : taxSet){
			int x=tree.getIdAtLevelExtended(i, currentLevelE);
			desired.add(x);
		}
		HashSet<Integer> present=new HashSet<Integer>();
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith(">")){
				TaxNode tn=tree.parseNodeFromHeader(line.substring(1), true);
				if(tn!=null){
					present.add(tree.getIdAtLevelExtended(tn.id, currentLevelE));
				}
			}
		}
		
		while(currentLevelE<TaxTree.LIFE_E){
			//Intersect
			for(Integer i : desired){
				if(present.contains(i)){
					if(currentLevelE!=taxLevelE){System.err.println("Widened filter from "+
							TaxTree.levelToStringExtended(taxLevelE)+" to "+TaxTree.levelToStringExtended(currentLevelE));}
					taxLevelE=currentLevelE;
					taxSet=desired;
					return;
				}
			}
			currentLevelE++;
			
			HashSet<Integer> desired2=new HashSet<Integer>();
			for(int i : desired){
				desired2.add(tree.getIdAtLevelExtended(i, currentLevelE));
			}
			desired=desired2;
			
			HashSet<Integer> present2=new HashSet<Integer>();
			for(int i : present){
				present2.add(tree.getIdAtLevelExtended(i, currentLevelE));
			}
			present=present2;
		}
	}
	
	public static boolean validArgument(String a){
		if(a.equals("table") || a.equals("gi")){
		}else if(a.equals("tree") || a.equals("taxtree")){
		}else if(a.equals("accession")){
		}else if(a.equals("taxpath")){
		}else if(a.equals("level") || a.equals("taxlevel")){
		}else if(a.equals("name") || a.equals("names")){
		}else if(a.equals("regex")){
		}else if(a.equals("contains")){
		}else if(a.equals("besteffort")){
		}else if(a.equals("include")){
		}else if(a.equals("promote")){
		}else if(a.equals("exclude")){
		}else if(a.equals("printnodesadded") || a.equals("printnodes")){
		}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("taxids")){
		}else if(a.equals("requirepresent")){
		}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
		}else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	void loadAccession(String accessionFile){
		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			assert(tree!=null);
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
	}
	
	static void loadGiTable(String fname){
		if(fname==null){return;}
		if(PRINT_STUFF){outstream.println("Loading gi table.");}
		GiToNcbi.initialize(fname);
	}
	
	static TaxTree loadTree(String fname){
		if(fname==null){return null;}
		TaxTree tt=TaxTree.loadTaxTree(fname, PRINT_STUFF ? outstream : null, true, false);
		assert(tt.nameMap!=null);
		return tt;
	}
	
	public void addNames(String names){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addName(name);
		}
	}
	
	public boolean addName(String name){
		{
			TaxNode tn=tree.parseNodeFromHeader(name, true);
			if(tn!=null){return addNode(tn);}
		}
		List<TaxNode> list=tree.getNodesByNameExtended(name);
		boolean success=false;
		assert(list!=null) : "Could not find a node for '"+name+"'";
		for(TaxNode tn : list){
			success=addNode(tn)|success;
		}
		return success;
	}
	
	public void addNumbers(String numbers, boolean promote){
		if(numbers==null){return;}
		String[] array=numbers.split(",");
		for(String s : array){
			if(Tools.isDigit(s.charAt(0))){
				final int x=Integer.parseInt(s);
				addNumber(x, promote);
			}else if(new File(s).exists()){
				TextFile tf=new TextFile(s);
				for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
					addNumbers(line, promote);
				}
			}else{
				assert(tree!=null) : "Using organism names requires a taxonomic tree loaded; please use a numeric NCBI taxID.";
				final int x=tree.parseNameToTaxid(s);
				assert(x>0) : "Can't find a tax node for "+s;
				addNumber(x, promote);
			}
		}
	}
	
	public boolean addNumber(int taxID, boolean promote){
		if(promote){
			TaxNode tn=tree.getNode(taxID);
			assert(tn!=null) : "Could not find a node for '"+taxID+"'";
			return addNode(tn);
		}else{
			return taxSet.add(taxID);
		}
	}
	
	public boolean addNode(TaxNode tn){
		if(tn==null){return false;}
		taxSet.add(tn.id);
		if(printNodesAdded){System.err.println("Added node "+tn);}//123
		while(tn.id!=tn.pid && tn.levelExtended<taxLevelE){
			tn=tree.getNode(tn.pid);
			if(tn.levelExtended<=taxLevelE){
				if(printNodesAdded){System.err.println("Added node "+tn);}//123
				taxSet.add(tn.id);
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean passesFilter(final Read r){
		return passesFilter(r.id);
	}
	
	public boolean passesFilterByNameOnly(final String name){
		if(regexPattern!=null){
			boolean b=matchesRegex(name);
			if(b!=include){return false;}
		}
		if(containsString!=null){
			boolean b=containsString(name);
			if(b!=include){return false;}
		}
		return true;
	}
	
	public boolean passesFilter(final String name){
		if(!passesFilterByNameOnly(name)){return false;}
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		TaxNode tn=tree.parseNodeFromHeader(name, true);
		if(tn==null){tn=tree.getNodeByName(name);}
//		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node for '"+name+"'";
		
		if(REQUIRE_PRESENT && tn==null){
			KillSwitch.kill("ERROR: Could not find node for '"+name+"'"
					+ "\nTo bypass this error, add the flag 'requirepresent=f'");
		}
		
//		assert(false) : passesFilter(tn);
		return passesFilter(tn);
	}
	
	public boolean passesFilter(final int id){
//		if((taxSet==null || taxSet.isEmpty()) && reqLevels==0){return !include;}
		if((taxSet==null || taxSet.isEmpty()) && reqLevels==0){return true;}
		TaxNode tn=tree.getNode(id);
//		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node number "+id;
		
		if(REQUIRE_PRESENT && tn==null){
//			assert(false) : ("ERROR: Could not find node for "+id);
			KillSwitch.kill("ERROR: Could not find node for "+id);
		}
		
		return passesFilter(tn);
	}
	
	boolean passesFilter(final TaxNode tn0){
		TaxNode tn=tn0;
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		if(tn==null){
			assert(!REQUIRE_PRESENT) : "Null TaxNode.";
			return !include && reqLevels==0;
		}
		boolean found=taxSet.contains(tn.id);
//		System.err.println("found="+found+", node="+tn);
		int levels=1<<tn.level;
		while((!found || (levels&reqLevels)!=reqLevels) && tn.id!=tn.pid){
			tn=tree.getNode(tn.pid);
			levels|=(1<<tn.level);
			if(promote){found=found||taxSet.contains(tn.id);}
//			System.err.println("found="+found+", node="+tn);
		}
//		assert(false) : levels+", "+reqLevels+", "+tn0+", "+tree.getAncestors(tn0.pid);
		return include==found && (levels&reqLevels)==reqLevels;
	}
	
	boolean matchesRegex(String s){
		return regexPattern.matcher(s).matches();
	}
	
	boolean containsString(String s){
		return s.toLowerCase().contains(containsString);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Clone             ----------------*/
	/*--------------------------------------------------------------*/

	public TaxFilter deepCopy() {
		TaxFilter copy=null;
		try {
			copy=(TaxFilter) this.clone();
			if(taxSet!=null){
				copy.taxSet=new HashSet<Integer>();
				copy.taxSet.addAll(taxSet);
			}
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return copy;
	}
	
	public void clearSet(){
		taxSet=null;
	}
	
	public void makeSet(){
		taxSet=new HashSet<Integer>();
	}
	
	public void setInclude(boolean b){
		include=b;
	}
	
	public void setLevel(final int newLevel, boolean promote){
		final int newLevelE=TaxTree.levelToExtended(newLevel);
		assert(newLevelE>=taxLevelE || newLevelE<1 || taxSet==null && taxSet.isEmpty()) : "taxLevel may only be increased when the set is non-empty.";
		taxLevelE=newLevelE;
		if(promote){promote();}
	}
	
	public void promote(){
		if(taxSet!=null && !taxSet.isEmpty() && taxLevelE>0){
			ArrayList<Integer> list=new ArrayList<Integer>(taxSet.size());
			list.addAll(taxSet);
			taxSet.clear();
			for(Integer i : list){addNumber(i, true);}
		}
	}
	
	public int size(){return (taxSet==null ? 0 : taxSet.size());}

	public int taxLevel(){return TaxTree.extendedToLevel(taxLevelE);}
	public Integer[] taxSet(){
		return (taxSet==null || taxSet.isEmpty()) ? null : taxSet.toArray(new Integer[0]);
	}
	public boolean include(){return include;}
	public void setTree(TaxTree tree_){tree=tree_;}
	public TaxTree tree(){return tree;}
	public void setContainsString(String s){containsString=s;}
	public String containsString(){return containsString;}
	@Override
	public String toString(){return ""+taxSet;}
	
	public boolean isEmpty() {return taxSet.isEmpty();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private TaxTree tree;
	
	/** Level at which to filter */
	private int taxLevelE;
	
	/** Branch must contain ancestors at these levels (bitflag) */
	private final int reqLevels;
	
	/** Set of numeric NCBI TaxIDs */
	private HashSet<Integer> taxSet;

	private boolean include;
	private boolean promote;
	
	private String regex;
	private final Pattern regexPattern;
	private String containsString;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private static PrintStream outstream=System.err;
	
	/** Print loading messages */
	static boolean PRINT_STUFF=true;
	
	public static boolean REQUIRE_PRESENT=true;
	public static boolean printNodesAdded=true;

}
