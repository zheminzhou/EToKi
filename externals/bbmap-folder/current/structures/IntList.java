package structures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * Efficient equivalent of ArrayList for int 
 * @author Brian Bushnell
 * @date Sep 20, 2014
 *
 */
public final class IntList{
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** For benchmarking.
	 * Accepts an optional argument for list length. */
	public static void main(String[] args){
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		benchmark(length);
	}
	
	/** For benchmarking */
	private static void benchmark(final int length){
		Timer t=new Timer();

		System.gc();
		
		{
			System.err.println("\nIntList:");
			Shared.printMemory();
			t.start();
			IntList list=new IntList();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nArrayList:");
			Shared.printMemory();
			t.start();
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nLinkedList:");
			Shared.printMemory();
			t.start();
			LinkedList<Integer> list=new LinkedList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		

		
		{
			System.err.println("\nIntList:");
			Shared.printMemory();
			t.start();
			IntList list=new IntList();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			list.shuffle();
			t.stop("Shuf Time:  \t");
			t.start();
			list.sort();
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nArrayList:");
			Shared.printMemory();
			t.start();
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			Collections.shuffle(list);
			t.stop("Shuf Time:  \t");
			t.start();
			Collections.sort(list);
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nLinkedList:");
			Shared.printMemory();
			t.start();
			LinkedList<Integer> list=new LinkedList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			Collections.shuffle(list);
			t.stop("Shuf Time:  \t");
			t.start();
			Collections.sort(list);
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public IntList(){this(256);}
	
	public IntList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocInt1D(initial);
	}
	
	public IntList copy() {
		IntList copy=new IntList(size);
		copy.addAll(this);
		return copy;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){size=0;}
	
	public final void set(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(int value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public void subtractFrom(int value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	public final void add(int x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Bulk Operations       ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addAll(IntList counts) {
		final int[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	public void shuffle() {
		if(size<2){return;}
		Random randy=Shared.threadLocalRandom();
		for(int i=0; i<size; i++){
			int j=randy.nextInt(size);
			int temp=array[i];
			array[i]=array[j];
			array[j]=temp;
		}
	}
	
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Resizing           ----------------*/
	/*--------------------------------------------------------------*/
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	public final void setSize(final int size2) {
		if(size2<array.length){resize(size2);}
		size=size2;
	}
	
	public final void shrink(){
		if(size==array.length){return;}
		array=KillSwitch.copyOf(array, size);
	}
	
	public final void shrinkToUnique(){
		//Assumes sorted.
		if(size<=0){
			shrink();
			return;
		}
		
		int unique=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]!=array[i-1]){unique++;}
		}
		if(unique==array.length){return;}
		int[] alt=KillSwitch.allocInt1D(unique);
		
		alt[0]=array[0];
		for(int i=1, j=1; j<unique; i++){
			if(array[i]!=array[i-1]){
				alt[j]=array[i];
				j++;
			}
		}
		
		array=alt;
		size=alt.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	public int lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	//Slow; for validation
	//Does not assume sorted
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	public boolean contains(int x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	public int[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	/** Assumes this is sorted.
	 * Reduces the list to a set of unique values;
	 * stores their counts in a second list. */
	public void getUniqueCounts(IntList counts) {
		counts.size=0;
		if(size<=0){return;}

		int unique=1;
		int count=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]==array[i-1]){
				count++;
			}else{
				array[unique]=array[i];
				unique++;
				counts.add(count);
				count=1;
			}
		}
		if(count>0){
			counts.add(count);
		}
		size=unique;
		assert(counts.size==size);
	}
	
	public boolean sorted(){//slow
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	public int size() {
		return size;
	}
	
	public boolean isEmpty() {
		return size<1;
	}
	
	public int capacity() {
		return array.length;
	}
	
	public int freeSpace() {
		return array.length-size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=0){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
				sb.append(comma+array[i]);
				comma=", ";
		}
		sb.append(']');
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Holds values. */
	public int[] array;
	/** Number of values in the primary array. */
	public int size=0;
	
}
