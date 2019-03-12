package shared;

import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import dna.Data;

public class Shared {
	
	public static void main(String[] args){
		COMMAND_LINE=args;
		mainClass=Shared.class;
		assert(false) : fullCommandline();
	}

	public static int LOGICAL_PROCESSORS=CALC_LOGICAL_PROCESSORS();
	private static int THREADS=setThreads(-1);
	
	private static int READ_BUFFER_NUM_BUFFERS=setBuffers();
	private static int READ_BUFFER_LENGTH=200;
	private static long READ_BUFFER_MAX_DATA=400000;
	
	/** Temporary, for testing; should be made non-global */
	public static boolean AMINO_IN=false;
	
	//TODO:  For some reason, it seems as though GAPBUFFER must equal exactly 1/2 of GAPLEN.  Not good; 1/4 would be far better.
	
	public static final int GAPBUFFER=64; //TODO:  Seems to break less than 64, for some reason
	public static final int GAPBUFFER2=2*GAPBUFFER;
	public static final int GAPLEN=128; //TODO: May break when over 128
	public static final int MINGAP=GAPBUFFER2+GAPLEN;
	public static final int GAPCOST=Tools.max(1, GAPLEN/64);
	public static final byte GAPC='-';

	public static String BBMAP_VERSION_STRING="38.42";
	public static String BBMAP_VERSION_NAME="Gravel Coast";
	
	public static boolean TRIM_READ_COMMENTS=false;
	public static boolean TRIM_RNAME=false; //For mapped sam reads
	
	public static boolean USE_JNI=false;//Data.GENEPOOL;
	public static boolean USE_MPI=false;
	public static boolean MPI_KEEP_ALL=true;
	/** Use ConcurrentReadInputStreamMPI instead of D */
	public static boolean USE_CRISMPI=true;
	public static int MPI_RANK=0;
	public static int MPI_NUM_RANKS=1;
	
	public static int FASTA_WRAP=70;
	public static byte FAKE_QUAL=30;
	
	public static boolean FIX_EXTENSIONS=true;
	
	/** True if assertions are enabled. */
	private static boolean EA=false;
	
	public static boolean EA(){return EA;}

	public static String BBMAP_CLASS=null;
	public static Class<?> mainClass=null;
	public static String[] COMMAND_LINE=null;
	
	public static final byte PLUS=0;
	public static final byte MINUS=1;
	/** Index with strand number */
	public static final String[] strandCodes={"+", "-", "?"};
	public static final char[] strandCodes2={'+', '-', '?'};
	
	public static List<String> JVM_ARGS(){
		return ManagementFactory.getRuntimeMXBean().getInputArguments();
	}
	
	public static String fullCommandline(){
		StringBuilder sb=new StringBuilder();
		sb.append("java ");
		for(String s : JVM_ARGS()) {
			sb.append(s).append(' ');
		}
		sb.append("-cp "+System.getProperty("java.class.path")+" ");
		sb.append(mainClass.getCanonicalName()).append(' ');
		for(String s : COMMAND_LINE) {
			sb.append(s).append(' ');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}
	
	/** Directory in which to write temp files */
	private static String TMPDIR=getTmpdir();
//	static{assert(false) : "TMPDIR="+TMPDIR;}
	
	private static String getTmpdir(){
		String s=System.getenv("SLURM_TMP");
		if(s==null){s=System.getenv("TMPDIR");}
		if(s!=null){s=(s+"/").replaceAll("//", "/").replaceAll("\\\\", "/");}
		return s;
	}
	
	public static String tmpdir(){return TMPDIR;}
	
	public static String setTmpdir(String s){
		if(s==null){TMPDIR=null;}
		else{
			s=s.replaceAll("\\\\", "/");
			if(!s.endsWith("/")){s=s+"/";}
			TMPDIR=s.replaceAll("//", "/");
		}
		return TMPDIR;
	}
	
	/** Anomaly probably resolved as of v.20.1
	 * This variable should be TRUE for normal users and FALSE for me. */
	public static boolean anomaly=!(System.getProperty("user.dir")+"").contains("/bushnell/") && !Data.WINDOWS;
	
	public static final char[] getTLCB(int len){
		char[] buffer=TLCB.get();
		if(buffer==null || buffer.length<len){
			buffer=new char[len];
			if(len<1000000){TLCB.set(buffer);}
		}
		return buffer;
	}
	private static final ThreadLocal<char[]> TLCB=new ThreadLocal<char[]>();
	
	/*--------------------------------------------------------------*/
	/*----------------           Threads            ----------------*/
	/*--------------------------------------------------------------*/

	public static int capThreads(int t) {
		assert(THREADS>0) : THREADS;
		final int old=THREADS;
		THREADS=Tools.mid(1, t, old);
		assert(THREADS>0) : THREADS;
		return old;
	}
	
	public static int setThreads(String x){
		int y=LOGICAL_PROCESSORS;
		if(x!=null && !x.equalsIgnoreCase("auto")){
			y=Integer.parseInt(x);
		}
		return setThreads(y);
	}
	
	public static int setThreads(int x){
		if(x>0){
			THREADS=x;
		}else{
			THREADS=Tools.max(1, LOGICAL_PROCESSORS);
		}
		setBuffers();
		assert(THREADS>0) : THREADS;
		return THREADS;
	}
	
	public static int threads(){
		assert(THREADS>0) : THREADS;
		return THREADS;
	}
	
	public static int CALC_LOGICAL_PROCESSORS(){
		final int procs=Tools.max(1, Runtime.getRuntime().availableProcessors());
		int slots=procs;
		Map<String,String> env=System.getenv();
		String s=env.get("NSLOTS");//Genepool
		boolean success=false;
		if(s!=null){
			int x=slots;
			try {
				x=Tools.max(1, Integer.parseInt(s));
				success=true;
			} catch (NumberFormatException e) {
				//ignore
			}
			if(x<=16){slots=x;}
		}
		if(!success){
			s=env.get("SLURM_CPUS_ON_NODE");//All SLURM systems
			if(s!=null){
				int x=slots;
				try {
					x=Tools.max(1, Integer.parseInt(s));
					success=true;
				} catch (NumberFormatException e) {
					//ignore
				}
				slots=x;
			}
		}

//		if(slots>8 && (slots*2==procs || (slots==16 && procs==40))){return procs;}//hyperthreading
		return Tools.min(slots, procs);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Buffers            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int capBuffers(int num){
		return setBuffers(Tools.min(num, READ_BUFFER_NUM_BUFFERS));
	}
	
	public static int READ_BUFFER_NUM_BUFFERS() {
		return READ_BUFFER_NUM_BUFFERS;
	}
	
	public static int setBuffers(){
		return setBuffersFromThreads(THREADS);
	}
	
	public static int setBuffersFromThreads(int threads){
		return setBuffers(Tools.max(4, (threads*3)/2));
	}
	
	public static int setBuffers(int num){
//		assert(READ_BUFFER_NUM_BUFFERS==0 || READ_BUFFER_NUM_BUFFERS==num) : READ_BUFFER_NUM_BUFFERS+" -> "+num; //TODO: 123
		num=Tools.max(2, num);
		return READ_BUFFER_NUM_BUFFERS=num;
	}
	
	public static int numBuffers(){
		return READ_BUFFER_NUM_BUFFERS;
	}
	
	public static int bufferLen(){
		return READ_BUFFER_LENGTH;
	}
	
	public static long bufferData(){
		return READ_BUFFER_MAX_DATA;
	}
	
	public static void capBufferLen(int x){
		if(x!=READ_BUFFER_LENGTH){setBufferLen(Tools.min(x, READ_BUFFER_LENGTH));}
	}
	
	public static int setBufferLen(int x){
//		assert(false) : READ_BUFFER_LENGTH+" -> "+x; //TODO: 123
		assert(x>0);
		return READ_BUFFER_LENGTH=x;
	}
	
	public static long setBufferData(long x){
//		assert(false) : READ_BUFFER_MAX_DATA+" -> "+x; //TODO: 123
		assert(x>0);
		return READ_BUFFER_MAX_DATA=x;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Memory            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean LOW_MEMORY=false;
	
	/** Ratio of -Xms to -Xmx parameters */
	public static final double xmsRatio(){
		Runtime rt=Runtime.getRuntime();
		return rt.totalMemory()*1.0/rt.maxMemory();
	}
	
	public static long memAvailable(int readThreads){
		long usableMemory;
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
			usableMemory=(long)Tools.max(((memory-48000000-(Tools.max(readThreads, 4)*400000))*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
		}
		return usableMemory;
	}
	
	public static long memTotal(){
		Runtime rt=Runtime.getRuntime();
		return rt.maxMemory();
	}
	
	public static long memFree(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
	}
	
	public static long memAvailable(){
		Runtime rt=Runtime.getRuntime();
		return rt.maxMemory()-rt.totalMemory()+rt.freeMemory();
	}
	
	/** An estimate, for preallocation */
	public static long memAvailableAdvanced(){
		Runtime rt=Runtime.getRuntime();
		final long mmemory=rt.maxMemory();
		final long tmemory=rt.totalMemory();
		final long fmemory=rt.freeMemory();
		final long umemory=tmemory-fmemory;

		double xmsRatio=Shared.xmsRatio();
		double usableMemory=Tools.max(((mmemory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), mmemory*0.45);
		double availableMemory=usableMemory-umemory;

		return (long)availableMemory;
	}
	
	public static long memUsed(){
		Runtime rt=Runtime.getRuntime();
		return rt.maxMemory()-rt.freeMemory();
	}
	
	/** Print statistics about current memory use and availability */
	public static final void printMemory(){
		try{
			if(GC_BEFORE_PRINT_MEMORY){
				System.gc();
				System.gc();
			}
			Runtime rt=Runtime.getRuntime();
			long mmemory=rt.maxMemory()/1000000;
			long tmemory=rt.totalMemory()/1000000;
			long fmemory=rt.freeMemory()/1000000;
			long umemory=tmemory-fmemory;
			System.err.println("Memory: "+"max="+mmemory+"m, total="+tmemory+"m, "+"free="+fmemory+"m, used="+umemory+"m");
		}catch(Throwable t){}
	}
	
	public static final Random threadLocalRandom(){
		Random randy;
		try {
			randy=ThreadLocalRandom.current();
		} catch (Throwable e) {//In case the JVM does not support ThreadLocalRandom;
			randy=new Random();
		}
		return randy;
	}
	
	/** Do garbage collection prior to printing memory usage */
	public static boolean GC_BEFORE_PRINT_MEMORY=false;
	
	/** For Matt N. */
	public static String comment;
	
	/*--------------------------------------------------------------*/
	/*----------------            Java 8            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final void sort(int[] array){sort(array, 0, array.length);}
	public static final void sort(int[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}

	public static final void sort(long[] array){sort(array, 0, array.length);}
	public static final void sort(long[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final void sort(float[] array){sort(array, 0, array.length);}
	public static final void sort(float[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}

	public static final <T extends Comparable<? super T>> void sort(T[] array){sort(array, 0, array.length);}
	public static final <T extends Comparable<? super T>> void sort(T[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to); //Supported in pre-Java 8.
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final <T extends Comparable<? super T>> void sort(ArrayList<T> list){
		try {
			if(!parallelSort || list.size()<=parallelSortLength || THREADS<2){
				Collections.sort(list); //Supported in pre-Java 8.
				return;
			}
			
			{//If this block causes compile errors, just replace the whole function body with "Shared.sort(list, comparator);"
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Comparable[0]);
				list.clear();
				Arrays.parallelSort(array);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	public static final <T> void sort(ArrayList<T> list, Comparator<? super T> comparator){
		try {
			if(!parallelSort){
				Collections.sort(list, comparator); //Supported in pre-Java 8.
				return;
			}
			
			{//If this block causes compile errors, just replace the whole function body with "Shared.sort(list, comparator);"
				if(list.size()<=parallelSortLength || THREADS<2){
					list.sort(comparator);
					return;
				}
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Object[0]);
				list.clear();
				Arrays.parallelSort(array, comparator);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}

	public static void closeStream(PrintStream outstream) {
		if(outstream!=null){
			synchronized(outstream){
				if(outstream!=System.err && outstream!=System.out){
					outstream.close();
				}
			}
		}
	}
	
	public static final long MAX_ARRAY_LEN=Integer.MAX_VALUE-20;
	
	public static final int parallelSortLength=10000;
	public static boolean disableParallelSort=false;
	public static boolean parallelSort=testParallelSort();
	
	public static double javaVersion=parseJavaVersion();
	
	private static double parseJavaVersion(){
		String s=System.getProperty("java.version");
		if(s==null){return 1.6;}
		int dots=0;
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<s.length() && dots<2; i++){
			char c=s.charAt(i);
			if(c=='.'){dots++;}
			else if(!Tools.isDigit(c)){break;}
			if(dots>1){break;}
			sb.append(c);
		}
		return Double.parseDouble(sb.toString());
	}
	
	public static void setParallelSort(boolean x){
		if(x){
			disableParallelSort=false;
			parallelSort=testParallelSort();
		}else{
			disableParallelSort=true;
			parallelSort=false;
		}
	}
	
	/** This tests to see if the java version supports parallel sort, which was not true prior to 1.8. */
	private static boolean testParallelSort(){
		Method m=null;
		try {
			m=Arrays.class.getMethod("parallelSort", new Class[] {Object[].class, Comparator.class});
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		} catch (Throwable t) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}
		return m!=null;
	}
	
	static{
		assert(EA=true);//Sets the EA flag
		KillSwitch.addBallast();
	}
	
}
