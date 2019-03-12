package structures;
import java.util.ArrayList;
import java.util.Arrays;


public class Range implements Comparable<Range>{
	
	/** A numeric range, assuming 0-based, base-centered numbering. */
	public Range(int aa, int bb){
		
		assert(aa<=bb) : aa+">"+bb;
		a=aa;
		b=bb;
		length=b-a+1;
	}
	
	public static Range toRange(String s){
		String[] s2=s.replace("[","").replace("]","").replace("(","").replace(")","").replace(",","").split("-");
		
		int a, b;
		if(s2.length==1){
			a=b=Integer.parseInt(s2[0]);
		}else{
			a=Integer.parseInt(s2[0]);
			b=Integer.parseInt(s2[1]);
		}
		return new Range(a, b);
	}
	
	@Override
	public int compareTo(Range other) {
		if(a<other.a){return -1;}
		if(a>other.a){return 1;}
		
		if(b<other.b){return -1;}
		if(b>other.b){return 1;}
		
		return 0;
	}
	
	public boolean includes(int p){
		return p>=a && p<=b;
	}
	
	public boolean intersects(int p1, int p2){
		return overlap(a, b, p1, p2);
	}
	
	public boolean includes(int p1, int p2){
		assert(p1<=p2);
		return p1>=a && p2<=b;
	}
	
	public boolean intersects(Range other){
		return intersects(other.a, other.b);
	}
	
	public boolean touches(Range other){
		if(intersects(other.a, other.b)){return true;}
		return b==other.a-1 || a==other.b+1;
	}
	
	public boolean includes(Range other){
		return includes(other.a, other.b);
	}
	
	@Override
	public boolean equals(Object other){
		return equals((Range)other);
	}
	
	public Range merge(Range other){
		assert(touches(other));
		Range r=new Range(min(a, other.a), max(b, other.b));
		
		assert(r.includes(this));
		assert(r.includes(other));
		assert(r.length<=length+other.length);
		return r;
	}
	
	public boolean equals(Range other){
		return a==other.a && b==other.b;
	}
	
	@Override
	public int hashCode(){
		return new Long(Long.rotateLeft(a, 16)^b).hashCode();
	}
	
	@Override
	public String toString(){
		return "("+a+(a==b ? "" : (" - "+b))+")";
	}
	
	
	public static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	
	public static Range[] toRanges(int[] ...arrays){
		int len=0;
		int[] combined=null;

		if(arrays.length==1){
			combined=arrays[0];
			len=combined.length;
		}else{
			for(int i=0; i<arrays.length; i++){
				len+=arrays[i].length;
			}
			combined=new int[len];
			for(int i=0, index=0; i<arrays.length; i++){
				for(int j=0; j<arrays[i].length; j++){
					combined[index]=arrays[i][j];
					index++;
				}
			}
			Arrays.sort(combined);
		}
		
		ArrayList<Range> list=new ArrayList<Range>(16);
		int start=combined[0], last=combined[0];
		
//		System.out.println(Arrays.toString(combined));
		
		for(int i=0; i<len; i++){
			int x=combined[i];
			if(x>last+1){
				list.add(new Range(start, last));
				start=last=x;
			}else{
				last=x;
			}
		}
		list.add(new Range(start, last));
		return list.toArray(new Range[list.size()]);
	}
	
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public final int a;
	public final int b;
	public final int length;

	public Object obj1=null;
	public Object obj2=null;
}
