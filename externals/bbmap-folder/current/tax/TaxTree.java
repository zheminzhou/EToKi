package tax;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntLongHashMap;

/**
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxTree implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8478466459993996423L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		assert(args.length>=4) : "TaxTree syntax:\ntaxtree.sh names.dmp nodes.dmp merged.dmp tree.taxtree.gz\n";
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=(Shared.threads()>2 ? 11 : 9);
		ReadWrite.PIGZ_BLOCKSIZE=256;
		ReadWrite.PIGZ_ITERATIONS=60;
		
		Timer t=new Timer();
		TaxTree tree=new TaxTree(args);
		t.stop();
		
		outstream.println("Retained "+tree.nodeCount+" nodes:");
		
		for(int i=tree.treeLevelsExtended.length-1; i>=0; i--){
			outstream.print(tree.nodesPerLevelExtended[i]+"\t"+taxaNamesExtended[i]);
			if(verbose){
				int lim=10;
				for(int j=0; j<lim && j<tree.treeLevelsExtended[i].length; j++){
					TaxNode n=tree.treeLevelsExtended[i][j];
					outstream.print("\n"+n+" -> "+tree.nodes[n.pid]);
				}
				for(int j=tree.treeLevelsExtended[i].length-lim; j<tree.treeLevelsExtended[i].length; j++){
					if(j>=lim){
						TaxNode n=tree.treeLevelsExtended[i][j];
						outstream.print("\n"+n+" -> "+tree.nodes[n.pid]);
					}
				}
			}
			outstream.println();
		}
		
		
		outstream.println();
		outstream.println("Time: \t"+t);
		
		if(args.length>2){//Write a tree
			ReadWrite.write(tree, args[3], true);
		}
	}
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i>3){
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	private void assignStrains(){

		outstream.println("Assigning strains.");
		int strains=0, substrains=0;
		TaxNode bacteria=getNode(BACTERIA_ID); //Can't do a name lookup since the names are not hashed yet
		TaxNode archaea=getNode(ARCHAEA_ID);
		assert(bacteria.name.equalsIgnoreCase("Bacteria"));
		assert(archaea.name.equalsIgnoreCase("Archaea"));

		ArrayList<TaxNode> bactList=new ArrayList<TaxNode>();
		ArrayList<TaxNode> archList=new ArrayList<TaxNode>();
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK && tn.minParentLevelExtended<=SPECIES_E){
				if(descendsFrom(tn, bacteria)){
					bactList.add(tn);
				}else if(descendsFrom(tn, archaea)){
					archList.add(tn);
				}
			}
		}
		
		ArrayList<TaxNode> prokList=new ArrayList<TaxNode>(bactList.size()+archList.size());
		prokList.addAll(bactList);
		prokList.addAll(archList);
		
		for(TaxNode tn : prokList){
			if(tn.maxDescendantLevelIncludingSelf()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==SPECIES_E || parent.levelExtended==SUBSPECIES_E){
					tn.levelExtended=STRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(STRAIN_E);
					strains++;
				}
			}
		}
		
//		outstream.println("Assigned "+strains+" strains.");
		for(TaxNode tn : prokList){
			if(tn.maxDescendantLevelIncludingSelf()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==STRAIN_E){
					tn.levelExtended=SUBSTRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(SUBSTRAIN_E);
					substrains++;
				}
			}
		}
//		outstream.println("Assigned "+substrains+" substrains.");
	}
	
	private void assignStrainsOld(){

		outstream.println("Assigning strains.");
		int strains=0, substrains=0;
		TaxNode bacteria=getNode(BACTERIA_ID); //Can't do a name lookup since the names are not hashed
		assert(bacteria.name.equalsIgnoreCase("Bacteria"));
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==SPECIES_E && commonAncestor(parent, bacteria)==bacteria){
//					nodesPerLevelExtended[STRAIN_E]++;
//					nodesPerLevelExtended[tn.levelExtended]--;
					tn.levelExtended=STRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(STRAIN_E);
					strains++;
				}
			}
		}
//		outstream.println("Assigned "+strains+" strains.");
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==STRAIN_E && commonAncestor(parent, bacteria)==bacteria){
//					nodesPerLevelExtended[SUBSTRAIN_E]++;
//					nodesPerLevelExtended[tn.levelExtended]--;
					tn.levelExtended=SUBSTRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(SUBSTRAIN_E);
					substrains++;
				}
			}
		}
//		outstream.println("Assigned "+substrains+" substrains.");
	}
	
	private TaxTree(String[] args){
		this(args[0], args[1], args[2], args);
	}
	
	private TaxTree(String namesFile, String nodesFile, String mergedFile, String[] args){
		
		if(args!=null) {
			Parser parser=parse(args);
		}
		
		nodes=getNames(namesFile);
		getNodes(nodesFile, nodes);
		
		mergedMap=getMerged(mergedFile);
		
		countChildren();
		outstream.println("Counted children.");
		int rounds=percolate();
		outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		
		if(assignStrains){
			assignStrains();
			rounds=percolate();
			outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		}
		
		if(simplify){
			if(verbose){outstream.println("Simplifying.");}
			int removed=simplify(nodes);
			if(verbose){outstream.println("Removed "+removed+" nodes.");}
			rounds=percolate();
			outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		}
		int errors=test(nodes);
//		assert(errors==0); //Not possible since the tree is wrong.
		if(errors>0) {
			System.err.println("Found "+errors+" errors in tree.");
		}
		
		for(TaxNode n : nodes){
			if(n!=null){
				nodesPerLevel[n.level]++;
				nodesPerLevelExtended[n.levelExtended]++;
			}
		}
		
//		for(int i=0; i<nodesPerLevel.length; i++){
//			treeLevels[i]=new TaxNode[nodesPerLevel[i]];
//		}
		for(int i=0; i<nodesPerLevelExtended.length; i++){
			treeLevelsExtended[i]=new TaxNode[nodesPerLevelExtended[i]];
		}
		
//		{
//			int[] temp=new int[nodesPerLevel.length];
//			for(TaxNode n : nodes){
//				if(n!=null){
//					int level=n.level;
//					treeLevels[level][temp[level]]=n;
//					temp[level]++;
//				}
//			}
//		}
		
		{
			int[] temp=new int[nodesPerLevelExtended.length];
			for(TaxNode n : nodes){
				if(n!=null){
					int level=n.levelExtended;
					treeLevelsExtended[level][temp[level]]=n;
					temp[level]++;
				}
			}
		}
		nodeCount=(int)Tools.sum(nodesPerLevelExtended);
		
	}

	public static final TaxTree loadTaxTree(String taxTreeFile, PrintStream outstream, boolean hashNames, boolean hashDotFormat){
		if(taxTreeFile==null){return null;}
		return loadTaxTree(taxTreeFile, null, null, null, outstream, hashNames, hashDotFormat);
	}
	
	public static final TaxTree loadTaxTree(String taxTreeFile, String taxNameFile, String taxNodeFile, String taxMergedFile, PrintStream outstream, 
			boolean hashNames, boolean hashDotFormat){
		if(taxTreeFile!=null || taxNodeFile==null){
			TaxTree tree=sharedTree(taxTreeFile, hashNames, hashDotFormat, outstream);
			if(tree!=null){return tree;}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=defaultTreeFile();}
		assert(taxTreeFile!=null || (taxNameFile!=null && taxNodeFile!=null)) : "Must specify both taxname and taxnode files.";
		Timer t=new Timer();
		if(outstream!=null){outstream.print("\nLoading tax tree; ");}
		final TaxTree tree;
		if(taxTreeFile!=null){
			tree=ReadWrite.read(TaxTree.class, taxTreeFile, true);
		}else{
			tree=new TaxTree(taxNameFile, taxNodeFile, taxMergedFile, null);
		}
		t.stop();
		if(hashNames){
			if(outstream!=null){outstream.println("Hashing names.");}
			tree.hashNames(hashDotFormat);
		}
		if(outstream!=null){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		if(ALLOW_SHARED_TREE){sharedTree=tree;}
		return tree;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static TaxNode[] getNames(String fname){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(200000);
		int max=0;
		
		TextFile tf=new TextFile(fname, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.contains("scientific name")){
				String[] split=delimiter.split(s, 3);
				assert(split.length==3) : s;
				int id=Integer.parseInt(split[0]);
				String name=split[1];
				if(id==1 && name.equalsIgnoreCase("root")){name="Life";}
				max=Tools.max(max, id);
				list.add(new TaxNode(id, name));
			}
		}
		
		TaxNode[] nodes=new TaxNode[max+1];
		for(TaxNode n : list){
			assert(nodes[n.id]==null || nodes[n.id].equals(n)) : nodes[n.id]+" -> "+n;
			nodes[n.id]=n;
		}
		
		return nodes;
	}
	
	private void countChildren(){
		for(TaxNode child : nodes){
			if(child!=null && child.pid!=child.id){
				TaxNode parent=getNode(child.pid);
				if(parent!=child){parent.numChildren++;}
			}
		}
	}
	
	//TODO - This could be finished in 2 passes using the childTable.
	private int percolate(){
		boolean changed=true;
		int rounds=0;
		while(changed){
			changed=false;
			rounds++;
			for(TaxNode child : nodes){
				if(child!=null && child.pid!=child.id){
					TaxNode parent=getNode(child.pid);
					changed=(child.discussWithParent(parent) | changed);
				}
			}
			
			if(!changed){break;}
			changed=false;
			rounds++;
			for(int i=nodes.length-1; i>=0; i--){
				TaxNode child=nodes[i];
				if(child!=null && child.pid!=child.id){
					TaxNode parent=getNode(child.pid);
					changed=(child.discussWithParent(parent) | changed);
				}
			}
		}
		return rounds;
	}
	
	public synchronized void hashNames(boolean genusDotSpecies){
		assert(nameMap==null);
		assert(nameMapLower==null);
		final int size=((int)Tools.mid(2, (nodes.length+(genusDotSpecies ? nodesPerLevelExtended[SPECIES_E] : 0))*1.5, Shared.MAX_ARRAY_LEN));
		nameMap=new HashMap<String, ArrayList<TaxNode>>(size);
		nameMapLower=new HashMap<String, ArrayList<TaxNode>>(size);
		
		//Hash the names, both lowercase and uppercase
		for(TaxNode n : nodes){
			if(n!=null){
				String name=n.name;
				if(name.indexOf('_')>=0){
					name=name.replace('_', ' ').trim();
				}
				if(name!=null && !name.equals("environmental samples")){
					{
						ArrayList<TaxNode> list=nameMap.get(name);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMap.put(name, list);
						}
						list.add(n);
					}
					{
						String lc=name.toLowerCase();
						ArrayList<TaxNode> list=nameMapLower.get(lc);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMapLower.put(lc, list);
						}
						list.add(n);
					}
				}
			}
		}
		
		//Hash G.species versions of the names, both lowercase and uppercase
		if(genusDotSpecies){
			ByteBuilder bb=new ByteBuilder(64);
			for(TaxNode n : nodes){
				if(n!=null && n.levelExtended==SPECIES_E){
					String name=n.name;
					if(name.indexOf('_')>=0){
						name=name.replace('_', ' ').trim();
					}
					if(name!=null && !name.equals("environmental samples")){
						final String dotFormat=dotFormat(name, bb);
						if(dotFormat!=null){
							{
								ArrayList<TaxNode> list=nameMap.get(dotFormat);
								if(list==null){
									list=new ArrayList<TaxNode>();
									nameMap.put(dotFormat, list);
								}
								list.add(n);
							}
							{
								String lc=dotFormat.toLowerCase();
								ArrayList<TaxNode> list=nameMapLower.get(lc);
								if(list==null){
									list=new ArrayList<TaxNode>();
									nameMapLower.put(lc, list);
								}
								list.add(n);
							}
						}
					}
				}
			}
		}
	}
	
	/** Transform e.g. "Homo sapiens" to "H.sapiens" */
	private static String dotFormat(String name, ByteBuilder buffer){
		if(name==null || name.indexOf('.')>=0){return null;}
		final int firstSpace=name.indexOf(' ');
		if(firstSpace<0 || firstSpace>=name.length()-1){return null;}
		final int lastSpace=name.lastIndexOf(' ');
		if(firstSpace!=lastSpace){return null;}
		final String a=name.substring(0, firstSpace);
		final String b=name.substring(lastSpace+1);
		final char ca=a.charAt(0);
		final char cb=b.charAt(0);
		if(!Tools.isUpperCase(ca) || !Tools.isLowerCase(cb)){return null;}
		if(buffer==null){buffer=new ByteBuilder(2+b.length());}
		else{buffer.clear();}
		buffer.append(ca).append('.').append(b);
		return buffer.toString();
	}
	
	public synchronized void hashChildren(){
		assert(childMap==null);
		int nodesWithChildren=0;
		for(TaxNode tn : nodes){
			if(tn!=null && tn.numChildren>0){nodesWithChildren++;}
		}
		childMap=new HashMap<TaxNode, ArrayList<TaxNode>>((int)Tools.mid(2, nodesWithChildren*1.5, Shared.MAX_ARRAY_LEN));
		for(TaxNode tn : nodes){
			if(tn!=null){
				if(tn.numChildren>0){
					childMap.put(tn, new ArrayList<TaxNode>(tn.numChildren));
				}
			}
		}
		for(TaxNode tn : nodes){
			if(tn!=null){
				if(tn.id!=tn.pid){
					ArrayList<TaxNode> list=childMap.get(getNode(tn.pid));
					if(list!=null){list.add(tn);}
				}
			}
		}
	}
	
	public ArrayList<TaxNode> getChildren(TaxNode parent){
		if(parent.numChildren<1){return null;}
		if(childMap!=null){return childMap.get(parent);}
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(parent.numChildren);
		for(TaxNode tn : nodes){
			if(tn!=null && tn.id!=tn.pid && tn.pid==parent.id){
				list.add(tn);
			}
		}
		return list;
	}
	
	private static TaxNode[] getNodes(String fname, TaxNode[] nodes){
		
		int max=0;
		
		LinkedHashMap<String, int[]> oddNames=new LinkedHashMap<String, int[]>();
		
		TextFile tf=new TextFile(fname, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			String[] split=delimiter.split(s, 4);
			assert(split.length==4) : s;
			int id=-1, pid=-1, level=-1, levelExtended=-1;
			
			id=Integer.parseInt(split[0]);
			try {
				pid=Integer.parseInt(split[1]);
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("Bad line: "+s+"\n"+Arrays.toString(split));
			}
			boolean alt=false;
			{
				String key=split[2];
				Integer obj0=levelMap.get(key);
				Integer obj=levelMapExtended.get(key);
				assert(obj!=null) : "No level found for "+key+"; line="+Arrays.toString(split);
				
				if(obj0==null){
					obj0=altLevelMap.get(key);
					alt=true;
				}
				if(obj0!=null){
					level=obj0;
					levelExtended=obj;
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}
				}else{
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}else{
						int[] count=oddNames.get(key);
						if(count==null){
							count=new int[1];
							oddNames.put(key, count);
						}
						count[0]++;
					}
				}
			}
			max=Tools.max(max, id);
			TaxNode n=nodes[id];
			assert(n!=null && n.pid<0) : n+" -> "+s;
			n.pid=pid;
			n.level=level;
			n.levelExtended=levelExtended;
			n.setOriginalLevel(levelExtended);
			n.setCanonical(!alt);
			assert(n.canonical()==n.isSimple() || n.levelExtended==NO_RANK_E) : n.canonical()+", "+n.isSimple()+", "+n.level+", "+n.levelExtended+"\n"+n.toString()+"\n";
		}
		
		if(oddNames.size()>0){
			outstream.println("Found "+oddNames.size()+" unknown taxonomic levels:");
			if(verbose){
				for(String s : oddNames.keySet()){
					outstream.println(oddNames.get(s)[0]+"\t"+s);
				}
			}
		}
		
		return nodes;
	}
	
	private static IntHashMap getMerged(String mergedFile) {
		if(mergedFile==null){return null;}
		String[] lines=TextFile.toStringLines(mergedFile);
		if(lines.length<1){return null;}
		IntHashMap map=new IntHashMap((int)(lines.length*1.3));
		for(String line : lines){
			String[] split=delimiterTab.split(line);
			int a=Integer.parseInt(split[0]);
			int b=Integer.parseInt(split[2]);
			map.put(a, b);
		}
		return map;
	}
	
	private int simplify(TaxNode nodes[]){
		
		int failed=test(nodes);
		
		int removed=0;
		int reassigned=0;
		
		if(reassign){
			boolean changed=true;
			int changedCount=0;
			while(changed){
				changed=false;
				for(int i=0; i<nodes.length; i++){
					TaxNode n=nodes[i];
					if(n!=null && n.levelExtended<1){
						int pid=n.pid;
						TaxNode parent=nodes[pid];
						assert(parent!=null) : n;
						if(n.maxDescendantLevelIncludingSelf()<SUBSPECIES_E &&
								parent.minAncestorLevelIncludingSelf()<=SPECIES_E && parent.minAncestorLevelIncludingSelf()>=SUBSPECIES_E){
//						if(parent.levelExtended==SPECIES_E || parent.levelExtended==SUBSPECIES_E){
							changed=true;
							n.levelExtended=SUBSPECIES_E;
							n.level=SUBSPECIES;
							changedCount++;
						}
					}
				}
			}
			System.err.println("Assigned levels to "+changedCount+" unranked nodes.");
		}
		
		
		if(skipNorank){//Skip nodes with unknown taxa
			if(verbose){outstream.println("A0");}
			
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					int pid=n.pid;
					TaxNode parent=nodes[pid];
					assert(parent!=null) : n;
					assert(parent!=n || pid==1) : n+", "+pid;
					while(parent.levelExtended<1 && n.levelExtended>parent.levelExtended){
						//System.err.println("Reassigned from "+parent);
						assert(parent.id!=parent.pid);
						parent=nodes[parent.pid];
						n.pid=parent.id;
						reassigned++;
					}
				}
			}
			
			for(int i=0; i<nodes.length; i++){
				if(nodes[i]!=null && nodes[i].levelExtended<0){
					System.err.println("Removed "+nodes[i]);
					nodes[i]=null;
					removed++;
				}
			}
			if(verbose){outstream.println("Skipped "+reassigned+" unranked parents, removed "+removed+" invalid nodes.");}
		}
		
		if(inferRankLimit>0){//Infer level for unset nodes (from "no rank")
			if(verbose){outstream.println("A");}
			int changed=1;
			while(changed>0){
				changed=0;
				for(final TaxNode n : nodes){
					if(n!=null){
						if(n.levelExtended==0){
							TaxNode parent=nodes[n.pid];
							if(n!=parent && parent.levelExtended>0 && parent.levelExtended<=inferRankLimit+1){
								n.levelExtended=Tools.max(1, parent.levelExtended-1);
								assert(n.levelExtended>0 && n.levelExtended<=parent.levelExtended && n.levelExtended<=inferRankLimit);
								changed++;
							}
						}
					}
				}
				if(verbose){outstream.println("changed: "+changed);}
			}
			
//			outstream.println("B");
//			for(TaxNode n : nodes){
//				if(n!=null && n.level==0){
//					n.level=-1;
//				}
//			}
		}
		
		failed=test(nodes);
		
//		if(reassign){//Skip nodes with duplicate taxa
//			if(verbose){outstream.println("D");}
//			int changed=1;
//			while(changed>0){
//				changed=0;
//				for(final TaxNode n : nodes){
//					if(n!=null){
//						TaxNode parent=nodes[n.pid];
//						TaxNode grandparent=nodes[parent.pid];
//						assert(n.level<=parent.level || parent.level<1 || !parent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//						assert(parent.level<=grandparent.level || grandparent.level<1 || !grandparent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//
//						while(parent!=grandparent && (parent.level<0 || (parent.level==grandparent.level && !parent.canonical()) ||
//								n.level>parent.level || (n.level==parent.level))){
//							parent=grandparent;
//							grandparent=nodes[parent.pid];
//							n.pid=parent.id;
//							reassigned++;
//							changed++;
//						}
//					}
//				}
//				if(verbose){outstream.println("changed: "+changed);}
//			}
//			if(verbose){outstream.println("E");}
//			for(int i=0; i<nodes.length; i++){
//				if(nodes[i]!=null && nodes[i].level<0){
//					nodes[i]=null;
//					removed++;
//				}
//			}
//		}
		
		failed=test(nodes);

		if(verbose){outstream.println("F");}
		{//Ensure the tree is now clean
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					TaxNode parent=nodes[n.pid];
					TaxNode grandparent=nodes[parent.pid];
					assert(n==parent || n.levelExtended<=parent.levelExtended || !n.canonical() || n.levelExtended<1 || parent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
					assert(parent==grandparent || parent.levelExtended<=grandparent.levelExtended || !parent.canonical() || parent.levelExtended<1 || grandparent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
				}
			}
		}
		
//		if(verbose){System.err.println("Reassignments: "+reassigned);}
		
		return removed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int test(TaxNode[] nodes){
		int failed=0;
		for(final TaxNode n : nodes){
			if(n!=null){
				TaxNode parent=nodes[n.pid];
				try {
					assert(n==parent || n.level<=parent.level || parent.level<1 || !parent.canonical()) : 
						"\n"+n+" -> "+parent+", level="+n.level+", plevel="+parent.level+", pcanon="+parent.canonical()+"\n"
								+ "levelE="+n.levelExtended+", plevelE="+parent.levelExtended;
					assert(n==parent || n.levelExtended<=parent.levelExtended || parent.levelExtended<1) : n+" -> "+parent;
//				assert(n==parent || n.level<parent.level || parent.level<1 || !n.canonical() || !parent.canonical()) : n+" -> "+parent;
					if(n!=parent && n.level>parent.level && parent.level>=1 && n.canonical() && parent.canonical()){
						if(verbose){outstream.println("Error: "+n+" -> "+parent);}
						failed++;
					}else if(n!=parent && parent.levelExtended>=1 && n.levelExtended>=parent.levelExtended){
//					if(verbose){outstream.println("Error: "+n+" -> "+parent);}
//					failed++;
					}
					assert(n!=parent || n.id<=1) : n;
				} catch (Throwable e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					failed++;
				}
			}
		}
		if(verbose || failed>0){outstream.println(failed+" nodes failed.");}
		return failed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Print Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public String toSemicolon(final TaxNode tn0, boolean skipNonCanonical){
		StringBuilder sb=new StringBuilder();
		if(tn0==null){return "Not found";}
		String semi="";
		ArrayList<TaxNode> list=toAncestors(tn0, skipNonCanonical);
		boolean addTaxLevel=true;
		for(int i=list.size()-1; i>=0; i--){
			sb.append(semi);
			TaxNode tn=list.get(i);
			if(tn.id!=LIFE_ID || list.size()==1){
				if(addTaxLevel && tn.canonical() && !tn.levelChanged() && tn.isSimple()){
					sb.append(tn.levelToStringShort()).append(':');
				}
				sb.append(tn.name);
				semi=";";
			}
		}
		return sb.toString();
	}
	
	public IntList toAncestorIds(final TaxNode tn0, boolean skipNonCanonical){
		if(tn0==null){return null;}
		IntList list=new IntList(8);
		
		TaxNode tn=tn0;
		while(tn!=null){
			if(!skipNonCanonical || tn.isSimple()){
				if(tn.id!=CELLULAR_ORGANISMS_ID || tn==tn0){list.add(tn.id);}
			}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		if(list.isEmpty()){list.add(tn0.id);}
		return list;
	}
	
	public ArrayList<TaxNode> toAncestors(final TaxNode tn0, boolean skipNonCanonical){
		if(tn0==null){return null;}
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(8);
		
		TaxNode tn=tn0;
		while(tn!=null){
			if(!skipNonCanonical || tn.isSimple()){
				if(tn.id!=CELLULAR_ORGANISMS_ID || tn==tn0){list.add(tn);}
			}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		if(list.isEmpty()){list.add(tn0);}
		return list;
	}
	
	public String toDir(TaxNode node, String root){
		StringBuilder sb=new StringBuilder();
		if(root==null){root="";}
		sb.append(root);
		if(root.length()>0 && !root.endsWith("/")){sb.append('/');}
		IntList list=toAncestorIds(node, false);
		list.reverse();
		assert(list.get(0)==1) : list + "," +getNode(list.get(0));
		for(int i=0; i<list.size(); i++){
			sb.append(list.get(i));
			sb.append('/');
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public static int getID(String s){return GiToNcbi.getID(s);}
	
	public static int getID(byte[] s){return GiToNcbi.getID(s);}
	
	/** Return the ancestor with taxonomic level at least minLevel */
	public TaxNode getNode(String s, int minLevelExtended){
		TaxNode tn=parseNodeFromHeader(s, true);
		while(tn!=null && tn.levelExtended<minLevelExtended && tn.pid!=tn.id){
			tn=getNode(tn.pid);
		}
		return tn;
	}
	
	public boolean descendsFrom(final int child, final int parent){
		TaxNode cn=getNode(child), pn=getNode(parent);
		assert(cn!=null) : "Invalid taxID: "+child;
		assert(pn!=null) : "Invalid taxID: "+parent;
		return descendsFrom(cn, pn);
	}
	
	public boolean descendsFrom(TaxNode child, TaxNode parent){
		assert(child!=null && parent!=null) : "Null parameters.";
		if(child==null || parent==null){return false;}
		
		while(child!=parent && child.levelExtended<=parent.levelExtended && child.id!=child.pid){
			child=getNode(child.pid);
		}
		return child==parent;
	}
	
	public int commonAncestor(final int a, final int b){
		TaxNode an=getNode(a), bn=getNode(b);
		assert(an!=null) : "Invalid taxID: "+a;
		assert(bn!=null) : "Invalid taxID: "+b;
		TaxNode cn=commonAncestor(an, bn);
		assert(cn!=null) : "No common ancestor: "+an+", "+bn;
		if(cn==null){return -1;}
		return cn.id;
	}
	
	public TaxNode commonAncestor(TaxNode a, TaxNode b){
		assert(a!=null && b!=null) : "Null parameters.";
		if(a==null){return b;}
		if(b==null){return a;}
		
		while(a!=b){
			if(a.levelExtended<b.levelExtended){
				a=getNode(a.pid);
			}else{
				b=getNode(b.pid);
			}
		}
		return a;
	}
	
	public TaxNode highestAncestor(TaxNode a){
		assert(a!=null);
		while(a.id!=a.pid){a=getNode(a.pid);}
		return a;
	}
	
	public TaxNode parseNodeFromHeader(String header, boolean bestEffort){
		if(header==null || header.length()<2){return null;}
		if(header.charAt(0)=='>'){header=header.substring(1);}
		TaxNode tn;
		if(!SILVA_MODE){
			final char delimiter=ncbiHeaderDelimiter(header);
			if(delimiter==' '){
				tn=getNodeNewStyle(header);
			}else{
				tn=getNodeOldStyle(header, delimiter);
				if(tn==null && delimiter=='|'){
//					System.err.println("A: "+header);
					int id=-1;
					String[] split=header.split("\\|");
					if(AccessionToTaxid.LOADED()){
						for(int i=0; i<split.length && id<0; i++){//Try accessions first
							if(AccessionToTaxid.isValidAccession(split[i])){
								id=AccessionToTaxid.get(split[i]);
							}
						}
					}
					for(int i=0; i<split.length && id<0; i++){//Then names
						id=parseNameToTaxid(split[i]);
					}
//					System.err.println("E: "+id);
					if(id>=0){tn=getNode(id);}
//					System.err.println("F: "+tn);
				}
			}
		}else{
			tn=getNodeSilva(header, bestEffort);
		}
		return tn;
	}
	
	public static char ncbiHeaderDelimiter(String header){
		for(int i=0; i<header.length(); i++){
			final char c=header.charAt(i);
			if(c=='|' || c=='~'){
				assert(i>0) : "i="+i+"; malformatted header '"+header+"'";
				return c;
			}else if(Character.isWhitespace(c)){
				return ' ';
			}
		}
		return ' ';
	}
	
	TaxNode getNodeSilva(String s, boolean bestEffort){
		if(s==null){return null;}
		if(s.length()>=5 && s.startsWith("tid") && (s.charAt(3)=='|' || s.charAt(3)=='~') && Tools.isDigit(s.charAt(4))){
			return getNodeOldStyle(s, s.charAt(3));
		}
		String[] split=s.split(";");
		
		int number=-1;
		for(int i=split.length-1; number<0 && i>=0; i--){
			String last=split[i];
			int paren=last.indexOf('(');
			if(paren>=0){last=last.substring(0, paren);}
			last=last.trim();
			
			if(!last.startsWith("uncultured") && !last.startsWith("unidentified")){
				number=parseNameToTaxid(last);
			}
			
			if(number>=0){return getNode(number);}
			else if(!bestEffort){break;}
		}
		return null;
	}
	
	private TaxNode getNodeOldStyle(final String s, char delimiter){
		{
			int index=s.indexOf(delimiter);
			if(index<0){
				delimiter='~';
				index=s.indexOf(delimiter);
				if(index<0){
					delimiter='_';
					index=s.indexOf(delimiter);
				}
			}
			int number=-1;
			
			Throwable e=null;
			
			if(index==2 && s.length()>3 && s.startsWith("gi") && Tools.isDigit(s.charAt(3))){
//				System.err.println("Parsing gi number.");
				
				if(GiToNcbi.isInitialized()){
					try {
						number=GiToNcbi.parseGiToNcbi(s, delimiter);
					} catch (Throwable e2) {
						e=e2;
					}
				}else{
					assert(!CRASH_IF_NO_GI_TABLE) : "To use gi numbers, you must load a gi table.\n"+s;
				}
//				if(number!=-1){System.err.println("number="+number);}
			}else if(index==4 && s.length()>5 && s.startsWith("ncbi") && Tools.isDigit(s.charAt(5))){
//				System.err.println("Parsing ncbi number.");
				number=GiToNcbi.parseNcbiNumber(s, delimiter);
			}else if(index==3 && s.length()>4 && s.startsWith("tid") && Tools.isDigit(s.charAt(4))){
//				System.err.println("Parsing ncbi number.");
				number=GiToNcbi.parseNcbiNumber(s, delimiter);
			}else if(index==3 && s.length()>4 && s.startsWith("img") && Tools.isDigit(s.charAt(4))){
//				System.err.println("Parsing ncbi number.");
				long img=parseDelimitedNumber(s, delimiter);
				ImgRecord record=imgMap.get(img);
				number=(record==null ? -1 : record.taxID);
			}
			
			if(number<0 && index>=0 && (delimiter=='|' || delimiter=='~')){
				String[] split=(delimiter=='|' ? delimiterPipe.split(s) : delimiterTilde.split(s));
				if(AccessionToTaxid.LOADED()){
					number=parseAccessionToTaxid(split);
				}
				if(number<0){
					number=parseHeaderNameToTaxid(split);
				}
			}
			
			if(number<0 && e!=null){
				assert(false) : e;
				throw new RuntimeException(e);
			}
			
			//TaxServer code could go here...
			
			if(number>=0){return getNode(number);}
		}
		if(verbose){System.err.println("Can't process name "+s);}
		if(Tools.isDigit(s.charAt(0)) && s.length()<=9){
			try {
				return getNode(Integer.parseInt(s));
			} catch (NumberFormatException e) {
				//ignore
			}
		}
		return null;
	}
	
	/** Parse a delimited number from a header, or return -1 if formatted incorrectly. */
	static long parseDelimitedNumber(String s, char delimiter){
		if(s==null){return -1;}
		int i=0;
		while(i<s.length() && s.charAt(i)!=delimiter){i++;}
		i++;
		if(i>=s.length() || !Tools.isDigit(s.charAt(i))){return -1;}
		
		long number=0;
		while(i<s.length()){
			char c=s.charAt(i);
			if(c==delimiter || c==' ' || c=='\t'){break;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			number=(number*10)+(c-'0');
			i++;
		}
		return number;
	}
	
	private TaxNode getNodeNewStyle(final String s){
		
		int space=s.indexOf(' ');
		int number=-1;
		
		if(AccessionToTaxid.LOADED()){
			if(space>0){
				number=AccessionToTaxid.get(s.substring(0, space));
			}else{
				number=AccessionToTaxid.get(s);
			}
		}
		
		if(number<0 && Tools.isDigit(s.charAt(0)) && s.length()<=9 && space<0){
			try {
				return getNode(Integer.parseInt(s));
			} catch (NumberFormatException e) {
				//ignore
			}
		}
		
		if(number<0 && space>0){
			number=parseNameToTaxid(s.substring(space+1));
		}
		
		if(number>-1){return getNode(number);}
		if(space<0 && s.indexOf('_')>0){
			return getNodeNewStyle(s.replace('_', ' '));
		}
		return null;
	}
	
	public int parseAccessionToTaxid(String[] split){
		if(split.length<4){
			return -1;
		}
		int ncbi=AccessionToTaxid.get(split[3]);
		return ncbi;
	}
	
	public int parseHeaderNameToTaxid(String[] split){
		if(split.length<5){
			return -1;
		}
		return parseNameToTaxid(split[4]);
	}
	
	public int parseNameToTaxid(String name){
//		assert(false) : name+", "+(nameMap==null)+", "+(nameMap==null ? 0 : nameMap.size());
		List<TaxNode> list=null;
		
		list=getNodesByNameExtended(name);
		
		if(list==null || list.size()>1){return -1;}
		return list.get(0).id;
	}
	
	public List<TaxNode> getNodesByNameExtended(String name){
		List<TaxNode> list=null;
		
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		name=name.replaceAll("_", " ").trim();
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		String[] split2=name.split(" ");
		
		if(split2.length>7){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6]+" "+split2[7];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>6){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>5){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>4){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4];
			list=getNodesByName(term);
//			System.err.println("5:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>3){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3];
			list=getNodesByName(term);
//			System.err.println("4:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>2){
			String term=split2[0]+" "+split2[1]+" "+split2[2];
			list=getNodesByName(term);
//			System.err.println("3:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>1){
			String term=split2[0]+" "+split2[1];
			list=getNodesByName(term);
//			System.err.println("2:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>0){
			String term=split2[0];
			list=getNodesByName(term);
//			System.err.println("1:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		return null;
	}
	
//	public TaxNode getNode(byte[] s){
//		if(Tools.indexOf(s, (byte)'|')>=0){return getNode(GiToNcbi.getID(s));}
//
//		{
//			int index=Tools.indexOf(s, (byte)'|');
//			if(index<0){index=Tools.indexOf(s, (byte)'_');}
//			int number=-1;
//			if(index==2 && s.length>3 && Tools.startsWith(s, "gi") && Tools.isDigit(s[3])){
////				System.err.println("Parsing gi number.");
//				number=GiToNcbi.parseGiToNcbi(s);
//			}else if(index==4 && s.length>5 && Tools.startsWith(s, "ncbi") && Tools.isDigit(s[5])){
////				System.err.println("Parsing ncbi number.");
//				number=GiToNcbi.getID(s);
//			}
//			if(number>=0){return getNode(number);}
//		}
//		if(verbose){System.err.println("Can't process name "+new String(s));}
//		if(Tools.isDigit(s[0]) && s.length<=9){
//			try {
//				return getNode(Tools.parseInt(s, 0, s.length));
//			} catch (NumberFormatException e) {
//				//ignore
//			}
//		}
//		return null;
//	}
	
	public int getParentID(int id){
		assert(id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return -1;}
		TaxNode tn=nodes[id];
		if(tn==null && mergedMap!=null){tn=getNode(mergedMap.get(id), true);}
		return tn==null ? -1 : tn.pid;
	}
	
	public TaxNode getNode(int id){
		assert(id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return null;}
		TaxNode tn=nodes[id];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}
	
	public TaxNode getNode(int id, boolean skipAssertion){
		assert(skipAssertion || id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return null;}
		TaxNode tn=nodes[id];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel){
		return getNodeAtLevel(id, minLevel, DOMAIN);
	}
	
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE){
		return getNodeAtLevelExtended(id, minLevelE, DOMAIN_E);
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel, int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
		return getNodeAtLevelExtended(id, minLevelExtended, maxLevelExtended);
	}
	
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE, int maxLevelE){
		TaxNode tn=getNode(id);
		while(tn!=null && tn.pid!=tn.id && tn.levelExtended<minLevelE){
			TaxNode temp=getNode(tn.pid);
			if(temp==null || temp.levelExtended>maxLevelE){break;}
			tn=temp;
		}
		return tn;
	}
	
	public int getIdAtLevelExtended(int taxID, int taxLevelExtended){
		if(taxLevelExtended<0){return taxID;}
		TaxNode tn=getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<taxLevelExtended){
			tn=getNode(tn.pid);
			if(tn.levelExtended>taxLevelExtended){break;}
			taxID=tn.id;
		}
		return taxID;
	}

	public TaxNode getNodeByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		if(list==null || list.isEmpty()){return null;}
		if(list.size()==1){return list.get(0);}
		assert(false) : "Found multiple nodes for '"+s+"':\n"+list+"\n";
		TaxNode a=list.get(0);
		for(int i=1; i<list.size(); i++){
			TaxNode b=list.get(i);
			//Keep the most specific node
//			if(a==null || (b!=null && b.minAncestorLevelIncludingSelf()<a.minAncestorLevelIncludingSelf())){//not necessary
			if(b.minAncestorLevelIncludingSelf()<a.minAncestorLevelIncludingSelf()){
				a=b;
			}
		}
		return a;
	}
	public List<TaxNode> getNodesByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		return list;
	}
	private List<TaxNode> getNodesByName(String s, boolean lowercase){
		if(s==null){return null;}
		if(s.indexOf('_')>=0){s=s.replace('_', ' ');}
		if(lowercase){s=s.toLowerCase();}
//		System.err.println("Searching for "+s);
		final HashMap<String, ArrayList<TaxNode>> map=(lowercase ? nameMapLower : nameMap);
		assert(map!=null) : "Tax names were not hashed.";
		ArrayList<TaxNode> list=map.get(s);
		if(list!=null){return list;}
//		System.err.println("No matches for '"+s+"'");
		
//		assert(false) : nameMap.containsKey(s)+", "+nameMapLower.containsKey(s);
		
		if(s.indexOf('_')<0 && s.indexOf(' ')<0){return null;}
		String[] split=delimiter2.split(lowercase ? s.toLowerCase() : s, 8);
//		System.err.println("Array: "+Arrays.toString(split));
		list=map.get(split[split.length-1]);
		if(list==null){return list;}
//		System.err.println(list==null ? "No matches for "+split[split.length-1] : "Found list( "+list.size()+")");
		
		int matchCount=0;
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){matchCount++;}
		}
		if(matchCount==list.size()){return list;}
		if(matchCount<1){return null;}
		ArrayList<TaxNode> hits=new ArrayList<TaxNode>(matchCount);
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){hits.add(tn);}
		}
		return hits;
	}
	public ArrayList<TaxNode> getAncestors(int id){
		TaxNode current=getNode(id);
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		while(current!=null && current.pid!=current.id){//ignores root
			list.add(current);
			current=getNode(current.pid);
		}
		//optionally add root here
		return list;
	}
	
	public void increment(IntList ids, IntList counts, boolean sync){
		
		ids.sort();
		ids.getUniqueCounts(counts);
		
		if(!sync){
			for(int i=0; i<ids.size; i++){
				int id=ids.get(i);
				int count=counts.get(i);
				incrementRaw(id, count);
			}
		}else{
			synchronized(this){
				for(int i=0; i<ids.size; i++){
					int id=ids.get(i);
					int count=counts.get(i);
					incrementRaw(id, count);
				}
			}
		}
	}
	
	public void incrementRaw(int id, long amt){
		assert(id>=0 && id<nodes.length) : "TaxID "+id+" is out of range."+(id<0 ? "" : "  Possibly the taxonomy data needs to be updated.");
		assert(nodes[id]!=null) : "No node for TaxID "+id+"; possibly the taxonomy data needs to be updated.";
		nodes[id].incrementRaw(amt);
	}
	
	public void percolateUp(){
		for(int i=0; i<treeLevelsExtended.length; i++){
			percolateUp(i);
		}
	}
	
	public void percolateUp(final int fromLevel){
		final TaxNode[] stratum=treeLevelsExtended[fromLevel];
		for(final TaxNode n : stratum){
			n.incrementSum(n.countRaw);
			TaxNode parent=nodes[n.pid];
			if(n!=parent){
				parent.incrementSum(n.countSum);
			}
		}
	}
	
	/** Add this amount to the node and all its ancestors. */
	public void percolateUp(TaxNode node, long amt){
		if(amt==0){return;}
		if(verbose){System.err.println("percolateUp("+amt+") node: "+node);}
		while(node.id!=node.pid){
			node.incrementSum(amt);
			node=nodes[node.pid];
		}
		node.incrementSum(amt);
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit){
		return gatherNodesAtLeastLimit(limit, 0, nodesPerLevelExtended.length-1);
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit, final int minLevel, final int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
//		assert(false) : limit+", "+minLevel+", "+maxLevel+", "+minLevelExtended+", "+maxLevelExtended;
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		for(int i=minLevelExtended; i<nodesPerLevelExtended.length && i<=maxLevelExtended; i++){
			list.addAll(gatherNodesAtLeastLimitExtended(i, limit));
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimitExtended(final int fromLevelExtended, final long limit){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		final TaxNode[] stratum=treeLevelsExtended[fromLevelExtended];
		for(final TaxNode n : stratum){
			if(n.countSum>=limit){
				list.add(n);
				TaxNode parent=nodes[n.pid];
				if(n!=parent){
					percolateUp(parent, -n.countSum);//123 This was negative for some reason
				}
			}
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Static Initializers      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxaNames.length; i++){
			map.put(taxaNames[i], i);
			map.put(taxaNames[i].toUpperCase(), i);
		}
		return map;
	}

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeLevelMapExtended() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxaNamesExtended.length; i++){
			map.put(taxaNamesExtended[i], i);
			map.put(taxaNamesExtended[i].toUpperCase(), i);
		}
		return map;
	}

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeAltLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxaNames.length; i++){
			map.put(taxaNames[i], i);
			map.put(taxaNames[i].toUpperCase(), i);
		}
		
		//Add synonyms
//		map.put("subfamily", map.get("family"));
//		map.put("tribe", map.get("family"));
//		map.put("varietas", map.get("subspecies"));
//		map.put("subgenus", map.get("genus"));
//		map.put("forma", map.get("subspecies"));
//		map.put("species group", map.get("genus"));
//		map.put("species subgroup", map.get("genus"));
//		map.put("cohort", map.get("class"));
//		map.put("subclass", map.get("class"));
//		map.put("infraorder", map.get("order"));
//		map.put("superorder", map.get("class"));
//		map.put("subphylum", map.get("phylum"));
//		map.put("infraclass", map.get("class"));
//		map.put("superkingdom", map.get("division"));
//		map.put("parvorder", map.get("order"));
//		map.put("superclass", map.get("phylum"));
//		map.put("superphylum", map.get("kingdom"));
//		map.put("subkingdom", map.get("kingdom"));
//		map.put("superfamily", map.get("order"));
//		map.put("superkingdom", map.get("domain"));
//		map.put("suborder", map.get("order"));
//		map.put("subtribe", map.get("family"));
		
		for(String[] array : taxaNamesExtendedMatrix){
			String head=array[array.length-1];
			Integer value=map.get(head);
			assert(value!=null) : head;
			for(String key : array){
				if(key!=head){
					assert(!map.containsKey(key)) : "Map already contains key "+key+": "+Arrays.toString(array);
					map.put(key, value);
					map.put(key.toUpperCase(), value);
				}
			}
		}
		
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Size             ----------------*/
	/*--------------------------------------------------------------*/
	
	public long toSize(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMap==null){return -1L;}
		final long x=refseqSizeMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	public long toSizeC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMapC==null){return -1L;}
		final long x=refseqSizeMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	public int toSeqs(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMap==null){return -1;}
		final int x=refseqSeqMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	public long toSeqsC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMapC==null){return -1L;}
		final long x=refseqSeqMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	public int toNodes(TaxNode tn){
		if(tn==null){return 0;}
		if(nodeMapC==null){return -1;}
		final int x=nodeMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	public void loadSizeFile(String fname){
		if(fname==null){return;}
		assert(refseqSizeMap==null);
		refseqSizeMap=new IntLongHashMap();
		refseqSizeMapC=new IntLongHashMap();
		refseqSeqMap=new IntHashMap();
		refseqSeqMapC=new IntLongHashMap();
		nodeMapC=new IntHashMap();
		
		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		final byte delimiter='\t';
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]!='#'){
				int a=0, b=0;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 0: "+new String(line);
				int tid=Tools.parseInt(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 1: "+new String(line);
				long size=Tools.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 2: "+new String(line);
				long csize=Tools.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 3: "+new String(line);
				int seqs=Tools.parseInt(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 4: "+new String(line);
				long cseqs=Tools.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 5: "+new String(line);
				int cnodes=Tools.parseInt(line, a, b);
				b++;
				a=b;

				if(refseqSizeMap!=null && size>0){refseqSizeMap.put(tid, size);}
				if(refseqSizeMapC!=null && csize>0){refseqSizeMapC.put(tid, csize);}
				if(refseqSeqMap!=null && seqs>0){refseqSeqMap.put(tid, seqs);}
				if(refseqSeqMapC!=null && cseqs>0){refseqSeqMapC.put(tid, cseqs);}
				if(nodeMapC!=null && cnodes>0){nodeMapC.put(tid, cnodes);}
			}
		}
		bf.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             IMG              ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int imgToNcbi(long img){
		ImgRecord ir=imgMap.get(img);
//		assert(false) : "\n"+img+"\n"+imgMap.get(img)+"\n"+562+"\n"+imgMap.get(562)+"\n"+imgMap.size()+"\n"+IMGHQ+"\n"+defaultImgFile()+"\n";
		return ir==null ? -1 : ir.taxID;
	}
	
	public TaxNode imgToNcbiNode(long img){
		int tid=imgToNcbi(img);
		return tid<1 ? null : getNode(tid);
	}
	
//	public static int loadIMGOld(String fname, boolean storeName, PrintStream outstream){
//		assert(imgMap==null);
//		if(fname==null){return 0;}
//		ImgRecord2.storeName=storeName;
//		if(outstream!=null){System.err.println("Loading IMG.");}
//		Timer t=new Timer(outstream, false);
//		ImgRecord2[] array=ImgRecord2.toArray(fname);
//		int x=loadIMG(array);
//		t.stopAndPrint();
//		return x;
//	}
	
	public static int loadIMG(String fname, boolean storeName, PrintStream outstream){
		assert(imgMap==null);
		if(fname==null){return 0;}
		ImgRecord.storeName=storeName;
		if(outstream!=null){System.err.println("Loading IMG.");}
		Timer t=new Timer(outstream, false);
		ImgRecord[] array=ImgRecord.toArray(fname, IMG_HQ);
		int x=loadIMG(array);
		t.stopAndPrint();
		return x;
	}
	
	public static int loadIMG(ImgRecord[] array){
		assert(imgMap==null);
		imgMap=new HashMap<Long, ImgRecord>((int)(array.length*1.5));
		for(ImgRecord record : array){
			imgMap.put(record.imgID, record);
		}
		return imgMap.size();
	}
	
//	/** Look up a TaxNode from the ncbi TaxID */
//	TaxNode getTaxNodeSafe(int tid){
//		TaxNode tn=null;
//		try {
//			tn=getNode(tid);
//		} catch (Throwable e) {
//			if(verbose){e.printStackTrace();}
//		}
//		return tn;
//	}
	
	@Deprecated
	public static int parseLevel(String b){
		final int level;
		if(b==null){level=-1;}
		else if(Tools.isNumeric(b.charAt(0))){
			level=Integer.parseInt(b);
		}else{
			level=stringToLevel(b.toLowerCase());
		}
		return level;
	}
	
	public static int parseLevelExtended(String b){
		final int level;
		if(b==null){level=-1;}
		else if(Tools.isNumeric(b.charAt(0))){
			level=levelToExtended(Integer.parseInt(b));
		}else{
			level=stringToLevelExtended(b.toLowerCase());
		}
		return level;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final TaxNode[] nodes;
	public final int[] nodesPerLevel=new int[taxaNames.length];
	public final int[] nodesPerLevelExtended=new int[taxaNamesExtended.length];
	public final int nodeCount;
	public final IntHashMap mergedMap;

//	public final TaxNode[][] treeLevels=new TaxNode[taxaNames.length][];
	public final TaxNode[][] treeLevelsExtended=new TaxNode[taxaNamesExtended.length][];
	
	HashMap<String, ArrayList<TaxNode>> nameMap;
	HashMap<String, ArrayList<TaxNode>> nameMapLower;
	HashMap<TaxNode, ArrayList<TaxNode>> childMap;
//	HashMap<TaxNode, Long> sizeMap;//TODO - will probably break serialization
	public HashMap<String, ArrayList<TaxNode>> nameMap(){return nameMap;}
	
	public int minValidTaxa=0;
	
	public boolean simplify=true;
	public boolean reassign=true;
	public boolean skipNorank=false;
	public int inferRankLimit=0;//levelMap.get("species");
	
	//Node Statistics
	private IntLongHashMap refseqSizeMap;
	private IntLongHashMap refseqSizeMapC;
	private IntHashMap refseqSeqMap;
	private IntLongHashMap refseqSeqMapC;
	private IntHashMap nodeMapC;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean assignStrains=true;
	public static boolean SILVA_MODE=false;
	public static boolean CRASH_IF_NO_GI_TABLE=true;

	public static boolean verbose=false;
	public static boolean SHOW_WARNINGS=false;
	
	private static HashMap<Long, ImgRecord> imgMap;

	public static boolean ALLOW_SHARED_TREE=true;
	private static TaxTree sharedTree;
	
	//TODO: Check proper-construction of double-checked synchronize
	private static TaxTree sharedTree(String fname, boolean hashNames, boolean hashDotFormat, PrintStream outstream) {
		if(!ALLOW_SHARED_TREE){return null;}
		if(sharedTree==null && fname!=null){
			if("auto".equalsIgnoreCase(fname)){fname=defaultTreeFile();}
			synchronized(TaxTree.class){
				if(sharedTree==null){
					if(outstream!=null){outstream.println("Loading tax tree.");}
					Timer t=new Timer(outstream, false);
					setSharedTree(ReadWrite.read(TaxTree.class, fname, true), hashNames, hashDotFormat);
					t.stopAndPrint();
				}
			}
		}
		if(hashNames && sharedTree.nameMap==null){
			synchronized(sharedTree){
				if(sharedTree.nameMap==null){
					if(outstream!=null){outstream.println("Hashing names.");}
					Timer t=new Timer(outstream, false);
					sharedTree.hashNames(hashDotFormat);
					t.stopAndPrint();
				}
			}
		}
		return sharedTree;
	}
	
	private static synchronized void setSharedTree(TaxTree tree, boolean hashNames, boolean hashDotFormat){
		assert(ALLOW_SHARED_TREE);
		assert(sharedTree==null);
		sharedTree=tree;
		if(hashNames && sharedTree.nameMap==null){
			synchronized(sharedTree){
				if(sharedTree.nameMap==null){
					sharedTree.hashNames(hashDotFormat);
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int stringToLevel(String s){return altLevelMap.get(s);}
	public static final boolean levelMapExtendedContains(String s){return levelMapExtended.containsKey(s);}
	public static final int stringToLevelExtended(String s){return levelMapExtended.get(s);}
	public static final String levelToString(int x){return taxaNames[x];}
	public static final String levelToStringExtended(int x){return taxaNamesExtended[x];}
	public static final String levelToStringShort(int x){return taxaNamesShort[x];}
	
	private static final String[] taxaNames=new String[] {
		"no rank", "subspecies", "species", "genus",
		"family", "order", "class", "phylum",
		"kingdom", "superkingdom", "domain", "life"
	};
	
	private static final String[][] taxaNamesExtendedMatrix=new String[][] {
		{"no rank"},
		{"substrain", "strain", "forma", "varietas", "subspecies"},
		{"species"},
		{"species subgroup", "species group", "subgenus", "genus"},
		{"subtribe", "tribe", "subfamily", "family"},
		{"superfamily", "parvorder", "infraorder", "suborder", "order"},
		{"superorder", "subcohort", "cohort", "infraclass", "subclass", "class"},
		{"superclass", "subphylum", "phylum"},
		{"superphylum", "subkingdom", "kingdom"},
		{"superkingdom"},
		{"domain"},
		{"life"}
	};
	
	private static final String[] taxaNamesExtended=makeNamesExtended();
	
	private static final String[] makeNamesExtended(){
		ArrayList<String> list=new ArrayList<String>();
		for(String[] s : taxaNamesExtendedMatrix){
			for(String ss : s){
				list.add(ss);
			}
		}
		return list.toArray(new String[0]);
	}
	
	private static final String[] taxaNamesShort=new String[] {
			"nr", "ss", "s", "g",
			"f", "o", "c", "p",
			"k", "sk", "d", "l"
	};
	
	public static final int NO_RANK=0, SUBSPECIES=1, SPECIES=2, GENUS=3,
			FAMILY=4, ORDER=5, CLASS=6, PHYLUM=7, KINGDOM=8, SUPERKINGDOM=9, DOMAIN=10, LIFE=11;
	
	public static final int LIFE_ID=1;
	public static final int CELLULAR_ORGANISMS_ID=131567;
	public static final int BACTERIA_ID=2; //Is this safe?  Who knows...
	public static final int ARCHAEA_ID=2157;
	public static final int EUKARYOTA_ID=2759;
	public static final int METAZOA_ID=33208, ANIMALIA_ID=33208;
	public static final int VIRIDIPLANTAE_ID=33090, PLANTAE_ID=33090;
	public static final int FUNGI_ID=4751;
	public static final int VIRUSES_ID=10239;
	public static final int VIROIDS_ID=12884;
	
	private static final HashMap<String, Integer> levelMap=makeLevelMap();
	private static final HashMap<String, Integer> levelMapExtended=makeLevelMapExtended();
	private static final HashMap<String, Integer> altLevelMap=makeAltLevelMap();
	
	public static final int NO_RANK_E=NO_RANK,
			SUBSTRAIN_E=stringToLevelExtended("substrain"), STRAIN_E=stringToLevelExtended("strain"),
			SUBSPECIES_E=stringToLevelExtended("subspecies"),
			SPECIES_E=stringToLevelExtended("species"), GENUS_E=stringToLevelExtended("genus"),
			FAMILY_E=stringToLevelExtended("family"), ORDER_E=stringToLevelExtended("order"),
			CLASS_E=stringToLevelExtended("class"), PHYLUM_E=stringToLevelExtended("phylum"),
			KINGDOM_E=stringToLevelExtended("kingdom"), SUPERKINGDOM_E=stringToLevelExtended("superkingdom"),
			DOMAIN_E=stringToLevelExtended("domain"), LIFE_E=stringToLevelExtended("life");

	private static final int[] levelToExtended=new int[] {
			NO_RANK_E, SUBSPECIES_E, SPECIES_E, GENUS_E, FAMILY_E,
			ORDER_E, CLASS_E, PHYLUM_E, KINGDOM_E, SUPERKINGDOM_E, DOMAIN_E, LIFE_E
		};
	
	private static final int[] extendedToLevel=makeExtendedToLevel();
	
	private static int[] makeExtendedToLevel(){
		int len=0;
		for(String[] array : taxaNamesExtendedMatrix){
			len+=array.length;
		}
		int[] ret=new int[len];
		
		int pos=0;
		for(int level=0; level<taxaNamesExtendedMatrix.length; level++){
			String[] array=taxaNamesExtendedMatrix[level];
			for(int i=0; i<array.length; i++){
				ret[pos]=level;
				pos++;
			}
		}
		return ret;
	}
	
	public static final int levelToExtended(int level){
		return level<0 ? level : levelToExtended[level];
	}
	
	public static final int extendedToLevel(int extended){
		return extended<0 ? -1 : extendedToLevel[extended];
	}
	
	private static final Pattern delimiterTab = Pattern.compile("\t");
	private static final Pattern delimiter = Pattern.compile("\t\\|\t");
	private static final Pattern delimiterPipe = Pattern.compile("\\|");
//	private static final Pattern delimiterBang = Pattern.compile("\\!");
	private static final Pattern delimiterTilde = Pattern.compile("\\~");
	private static final Pattern delimiter2 = Pattern.compile("[\\s_]+");

	public static String TAX_PATH="/global/projectb/sandbox/gaag/bbtools/tax/latest";

	public static final String defaultTableFile(){return defaultTableFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultTreeFile(){return defaultTreeFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultAccessionFile(){return defaultAccessionFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultPatternFile(){return defaultPatternFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultImgFile(){return defaultImgFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultSizeFile(){return defaultSizeFile.replaceAll("TAX_PATH", TAX_PATH);} 
	
	public static boolean IMG_HQ=false;
	private static final String defaultImgFile="TAX_PATH/imgDump.txt";
	private static final String defaultTableFile="TAX_PATH/gitable.int1d.gz";
	private static final String defaultTreeFile="TAX_PATH/tree.taxtree.gz";
	private static final String defaultPatternFile="TAX_PATH/patterns.txt";
	private static final String defaultSizeFile="TAX_PATH/taxsize.tsv.gz";
//	private static final String defaultAccessionFile2="TAX_PATH/shrunk.dead_nucl.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.dead_prot.accession2taxid.gz,TAX_PATH/shrunk.dead_wgs.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.nucl_est.accession2taxid.gz,TAX_PATH/shrunk.nucl_gb.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.nucl_gss.accession2taxid.gz,TAX_PATH/shrunk.nucl_wgs.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.pdb.accession2taxid.gz,TAX_PATH/shrunk.prot.accession2taxid.gz";

	private static final String defaultAccessionFile="TAX_PATH/shrunk.prot.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_wgs.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_gb.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_prot.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_est.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_wgs.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_gss.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_nucl.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.pdb.accession2taxid.gz";
	
	private static PrintStream outstream=System.out;
	
}
