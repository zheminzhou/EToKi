package structures;

import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;



public final class LongList{
	
	public LongList(){this(256);}
	
	public LongList(int initial){
		assert(initial>0);
		array=KillSwitch.allocLong1D(initial);
	}
	
	public void clear(){
		size=0;
	}
	
	public final void set(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(long value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	public final void incrementBy(LongList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	public final void incrementBy(long[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	public final void append(LongList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	public final void append(long[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	public final void add(long x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
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
	
	public final double stdev(){
		if(size<2){return 0;}
		double sum=sum();
		double avg=sum/size;
		double sumdev2=0;
		for(int i=0; i<size; i++){
			long x=array[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/size);
	}
	
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	/** Assumes list is sorted */
	public final long median(){
		if(size<1){return 0;}
		int idx=percentile(0.5);
		return array[idx];
	}
	
	/** Assumes list is sorted */
	public final long mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		long prev=array[0];
		long best=prev;
		for(int i=0; i<size; i++){
			long x=array[i];
			if(x==prev){streak++;}
			else{
				if(streak>bestStreak){
					bestStreak=streak;
					best=prev;
				}
				streak=1;
				prev=x;
			}
		}
		if(streak>bestStreak){
			bestStreak=streak;
			best=prev;
		}
		return best;
	}
	
	public int percentile(double fraction){
		if(size<2){return size-1;}
		assert(sorted());
		double target=(sum()*fraction);
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return size-1;
	}
	
	//TODO: This could be done in-place.
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
		long[] alt=KillSwitch.allocLong1D(unique);
		
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
	
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		for(int i=0; i<x.length; i++){
			x[i]=array[i];
		}
		return x;
	}
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	public void sortSerial() {
		if(size>1){Arrays.sort(array, 0, size);}
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
	
	public int capacity() {
		return array.length;
	}
	
	public int freeSpace() {
		return array.length-size;
	}
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public long[] array;
	/** Highest occupied index plus 1, i.e., lowest unoccupied index */
	public int size=0;
	
}
