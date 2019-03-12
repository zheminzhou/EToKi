package structures;

import shared.Tools;
import stream.Read;

public class Quantizer {
	
	public static boolean parse(String arg, String a, String b){
		if(a.equals("quantize")){
			if(b!=null && b.equalsIgnoreCase("sticky")){
				STICKY=true;
				return true;
			}
		}else if(a.equals("quantizesticky")){
			STICKY=Tools.parseBoolean(b);
			return true;
		}
		
		if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
			return Tools.parseBoolean(b);
		}
		return setArray(b);
	}
	
	private static boolean setArray(String s){
		final byte[] array;
		if(s.charAt(0)=='/'){
			int quant=Integer.parseInt(s.substring(1));
			assert(quant>0 && quant<128);
			if(quant==1){return false;}
			ByteBuilder bb=new ByteBuilder();
			for(int i=0, max=Read.MAX_CALLED_QUALITY(); i<=max; i+=quant){
				bb.append((byte)i);
			}
			array=bb.toBytes();
		}else{
			array=Tools.parseByteArray(s, ",");
		}
		setArray(array);
		return true;
	}
	
	private static void setArray(byte[] a){
		quantizeArray=a;
		qualityRemapArray=makeQualityRemapArray(quantizeArray);
	}
	
	public static void quantize(Read r1, Read r2){
		quantize(r1.quality);
		if(r2!=null){quantize(r2.quality);}
	}
	
	public static void quantize(byte[] quals){
		if(quals==null){return;}
		byte prev=0;
		for(int i=0; i<quals.length; i++){
			final byte qOld=quals[i];
			byte q=qualityRemapArray[qOld];
			if(STICKY && q!=prev && prev>0 && q>0 && Tools.absdif(qOld, prev)<=Tools.absdif(qOld, q)){q=prev;}
			quals[i]=q;
			prev=q;
		}
	}
	
	private static final byte[] makeQualityRemapArray(byte[] quantizeArray) {
		byte[] array=new byte[128];
		for(int i=0; i<array.length; i++){
			byte q=0;
			for(byte x : quantizeArray){
				if((i>0 && q==0 && x>0) || Tools.absdif(x, i)<=Tools.absdif(q, i)){q=x;}
			}
			array[i]=q;
		}
		return array;
	}
	
//	private static byte[] quantizeArray={0, 8, 13, 22, 27, 32, 37}; //Old
	private static byte[] quantizeArray={0, 14, 21, 27, 32, 36};
	private static byte[] qualityRemapArray=makeQualityRemapArray(quantizeArray);
	private static boolean STICKY=true;
	
}
