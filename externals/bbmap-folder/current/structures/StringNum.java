package structures;

public class StringNum implements Comparable<StringNum> {

	public StringNum(String s_, long n_){
		s=s_;
		n=n_;
	}

	public long increment(){
		return (n=n+1);
	}
	
	public long increment(long x){
		return (n=n+x);
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(StringNum o) {
		if(n<o.n){return -1;}
		if(n>o.n){return 1;}
		return s.compareTo(o.s);
	}

	@Override
	public String toString(){
		return s+"\t"+n;
	}
	
	public boolean equals(StringNum other){
		if(other==null){return false;}
		if(s==other.s){return true;}
		if(s==null || other.s==null){return false;}
		return compareTo(other)==0;
	}
	
	/*--------------------------------------------------------------*/

	public final String s;
	public long n;

}
