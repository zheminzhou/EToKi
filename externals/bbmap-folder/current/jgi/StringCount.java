package jgi;

/**
 * Object holding a String and numbers, for tracking the number of read and base hits per scaffold.
 * Used by BBDuk and Seal.
 */
class StringCount implements Comparable<StringCount>{
	
	public StringCount(String name_, int len_, long reads_, long bases_){
		name=name_;
		length=len_;
		reads=reads_;
		bases=bases_;
	}
	@Override
	public final int compareTo(StringCount o){
		if(bases!=o.bases){return o.bases>bases ? 1 : -1;}
		if(reads!=o.reads){return o.reads>reads ? 1 : -1;}
		return name.compareTo(o.name);
	}
	public final boolean equals(StringCount o){
		return compareTo(o)==0;
	}
	@Override
	public final int hashCode(){
		return name.hashCode();
	}
	@Override
	public final String toString(){
		return name+"\t"+length+"\t"+reads+"\t"+bases;
	}
	
	/*--------------------------------------------------------------*/
	
	public final String name;
	public final int length;
	public final long reads, bases;
}