package structures;
import java.util.Locale;

import dna.Data;
import driver.Translator2;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;


public class CoverageArray2 extends CoverageArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 8242586595591123194L;
	
	public static void main(String[] args){
		runSpeedTest(args);
		
//		translateGenomeBuild(args);
	}
	
	public static void runSpeedTest(String[] args){
		
		long time1=System.nanoTime();
		
		CoverageArray2 ca=(CoverageArray2)read(args[1]);
		ca.chromosome=Byte.parseByte(args[0]);
		long time2=System.nanoTime();
		
//		int dot=args[1].lastIndexOf(".");
//		String outfile=args[1].substring(0,dot)+".ca";
		
		args[1]=args[1].replace('\\', '/');
		int slash=args[1].lastIndexOf('/');
		String outfile;
		if(slash<1){
			outfile="coverage-chr"+ca.chromosome+"-build"+Data.GENOME_BUILD+".ca";
		}else{
			outfile=args[1].substring(0,slash+1)+"coverage-chr"+ca.chromosome+"-build"+Data.GENOME_BUILD+".ca";
		}
		
		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+String.format(Locale.ROOT, "%.3f seconds", (time2-time1)/1000000000d));

		long time3=System.nanoTime();
		ReadWrite.write(ca, outfile, false);
		ca=null;
		System.gc();
		ca=(CoverageArray2)read(outfile);
		long time4=System.nanoTime();
		
		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+String.format(Locale.ROOT, "%.3f seconds", (time4-time3)/1000000000d));
		
		
	}
	
	@Deprecated
	/** Legacy human code */
	public static void translateGenomeBuild(String[] args){

		Timer t=new Timer();

		int inBuild=Integer.parseInt(args[0]);
		int outBuild=Integer.parseInt(args[1]);
		String root=args[2];
		
		translateGenomeBuild(inBuild, outBuild, root);
		
		t.stop();
		System.out.println("Time:\t"+t);
		
	}
	
	@Deprecated
	/** Legacy human code */
	public static void translateGenomeBuild(int inBuild, int outBuild, String root){
		root=root.replace('\\', '/');
		if(!root.endsWith("/")){root+="/";}
		
		CoverageArray2[] out=new CoverageArray2[27];
		
		for(int chrom=1; chrom<out.length; chrom++){
			out[chrom]=new CoverageArray2(chrom, 500);
		}
		
		final byte PLUS=Shared.PLUS;
		
		for(int chrom=1; chrom<=25; chrom++){
			String infile=root+"coverage-chr"+chrom+"-build"+inBuild+".ca.zip";
			CoverageArray2 ca1=ReadWrite.read(CoverageArray2.class, infile, true);
			for(int loc1=ca1.minIndex; loc1<=ca1.maxIndex; loc1++){
				char cov=(char)ca1.get(loc1);
				int[] xform=Translator2.translate(inBuild, outBuild, chrom, PLUS, loc1);
				if(xform!=null){
					int chrom2=(int)xform[0];
					int loc2=xform[2];
					out[chrom2].set(loc2, cov);
				}
			}
			ca1=null;
			System.out.println("Read "+infile);
		}
		
		for(int chrom=1; chrom<=25; chrom++){
			String outfile=root+"coverage-chr"+chrom+"-build"+outBuild+".ca.zip";
			out[chrom].resize(out[chrom].maxIndex+1);
			ReadWrite.write(out[chrom], outfile, false);
			out[chrom]=null;
			System.out.println("Wrote "+outfile);
		}
		
	}
	
//	public CoverageArray2(){
//		this((int)-1);
//	}
//
//	public CoverageArray2(int chrom){
//		this(chrom, 1<<24);
//	}
	
	public CoverageArray2(int chrom, int initialLen){
		super(chrom);
		array=KillSwitch.allocChar1D(initialLen);
	}
	
	/**
	 * @param loc
	 * @param amt
	 */
	@Override
	public void increment(int loc, int amt) {
		set(loc, get(loc)+amt);
	}
	
	/**
	 * @param loc
	 */
	@Override
	public void increment(int loc) {
		set(loc, get(loc)+1);
	}

	@Override
	public void incrementRange(int min, int max, int amt) {
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}else if(max<0){max=-1;}
		for(int i=min; i<=max; i++){
			int val=array[i]+amt;
			if(val>Character.MAX_VALUE){
				val=Character.MAX_VALUE;
				 if(!OVERFLOWED){
					 System.err.println("Note: Coverage capped at "+(int)(Character.MAX_VALUE)+"; please use the flag 32bit for higher values.");
					 OVERFLOWED=true;
				 }
			}
			array[i]=(char)val;
		}
		minIndex=min(min, minIndex);
		maxIndex=max(max, maxIndex);
	}
	
	
	@Override
	public void set(int loc, int val){
		
		if(loc>=array.length){//Increase size
			int newlen=1+(7*max(array.length, loc))/4;
			assert(newlen>loc);
			resize(newlen);
			assert(array.length==newlen);
		}else if(loc<0){
//			minIndex=min(0, minIndex);
//			maxIndex=max(0, maxIndex);
			return;
		}
		
		if(val>Character.MAX_VALUE && !OVERFLOWED){
			System.err.println("Note: Coverage capped at "+(int)(Character.MAX_VALUE)+"; please use the flag 32bit for higher values.");
			OVERFLOWED=true;
		}
		array[loc]=(val>Character.MAX_VALUE ? Character.MAX_VALUE : (char)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}
	
	@Override
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}
	
	@Override
	public void resize(int newlen){
//		System.err.println("Resized CoverageArray "+chromosome+" to "+newlen);
		char[] temp=KillSwitch.allocChar1D(newlen);
		int lim=min(array.length, newlen);
		assert(lim>maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<=maxIndex; i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array[i]);
		}
		sb.append(']');
		return sb.toString();
	}
	
	public char[] array;
	@Override
	public int length(){return maxIndex-minIndex+1;}
	@Override
	public int arrayLength(){return array.length;}
	
	private static boolean OVERFLOWED=false;
	/**
	 * 
	 */
//	private static final long serialVersionUID = -7493066925636540386L;
	
}
