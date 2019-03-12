package kmer;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * @author Brian Bushnell
 * @date May 14, 2015
 *
 */
public class AtomicShortArray {
	
	public AtomicShortArray(int length_){
		assert(length_>=0);
		length=length_;
		intArray=new AtomicIntegerArray((length+1)/2);
		assert(false) : "TODO";
	}
	
//	public short set(int position, short value){
//		in
//		intArray
//	}
	
	private AtomicIntegerArray intArray;
	private final int length;
	
}
