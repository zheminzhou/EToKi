package structures;

import java.util.ArrayList;
import java.util.LinkedList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;



public final class FloatList{
	
	public static void main(String[] args){
		Timer t=new Timer();
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		
		System.gc();
		
		{
			System.err.println("\nFloatList:");
			Shared.printMemory();
			t.start();
			FloatList list=new FloatList();
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
			ArrayList<Float> list=new ArrayList<Float>();
			for(int i=0; i<length; i++){
				list.add((float)i);
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
			LinkedList<Float> list=new LinkedList<Float>();
			for(int i=0; i<length; i++){
				list.add((float)i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
	}
	
	public FloatList(){this(256);}
	
	public FloatList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocFloat1D(initial);
	}
	
	public FloatList copy() {
		FloatList copy=new FloatList(size);
		copy.addAll(this);
		return copy;
	}
	
	public void clear(){size=0;}
	
	public final void set(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(float value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public void subtractFrom(float value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	public final float get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	public float lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	public final void add(float x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	//Slow; for validation
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	public void addAll(FloatList counts) {
		final float[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	public boolean contains(float x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	public final void setSize(final int size2) {
		if(size2<array.length){resize(size2);}
		size=size2;
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
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
		float[] alt=KillSwitch.allocFloat1D(unique);
		
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
	
	public float[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
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
	
	/** Assumes this is sorted.
	 * Reduces the list to a set of unique values;
	 * stores their counts in a second list. */
	public void getUniqueCounts(FloatList counts) {
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
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	public boolean sorted(){
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
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public float[] array;
	public int size=0;
	
}
