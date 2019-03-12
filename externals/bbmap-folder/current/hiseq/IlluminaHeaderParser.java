package hiseq;

import shared.KillSwitch;
import shared.Tools;
import stream.Read;

/**
 * Parses an Illumina header to gather positional information.
 * @author Brian Bushnell
 * @date Aug 22, 2018
 *
 */
public class IlluminaHeaderParser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Expected Format       ----------------*/
	/*--------------------------------------------------------------*/
	
	//@VP2-06:112:H7LNDMCVY:2:2437:14181:20134 (Novaseq)
	//2402:6:1101:6337:2237/1
	//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
	//HISEQ07:419:HBFNEADXX:1:1101:1238:2072
	//A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT
	
	//	@HWI-Mxxxx or @Mxxxx - MiSeq
	//	@HWUSI - GAIIx
	//	@HWI-Dxxxx - HiSeq 2000/2500
	//	@Kxxxx - HiSeq 3000(?)/4000
	//	@Nxxxx - NextSeq 500/550
	//
	//	AAXX = Genome Analyzer 
	//	BCXX = HiSeq v1.5 
	//	ACXX = HiSeq High-Output v3 
	//	ANXX = HiSeq High-Output v4 
	//	ADXX = HiSeq RR v1 
	//	AMXX, BCXX =HiSeq RR v2 
	//	ALXX = HiSeqX 
	//	BGXX, AGXX = High-Output NextSeq 
	//	AFXX = Mid-Output NextSeq 
	//	5 letter/number = MiSeq
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void parse(Read r){
		parse(r.id);
	}
	
	public void parse(String id_) {
		reset(id_);
		
		try {
			commentSeparator=findCommentSeparator();
			if(PARSE_COORDINATES){parseCoordinates();}
			if(PARSE_COMMENT){parseComment();}
		} catch (Throwable e) {
			System.err.println("Trouble parsing header "+id_);
			KillSwitch.throwableKill(e);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse lane, tile, x, and y coordinates */
	private void parseCoordinates(){
		pos=commentSeparator;
		goBackSeveralColons(4);
		lane=parseInt();
		if(!Tools.isDigit(id.charAt(pos))){//Hiseq 3000?  I'm not really sure what the header looks like for this block
			while(pos<limit && id.charAt(pos)!=':'){pos++;}
			pos++;
			lane=parseInt();
		}

		tile=parseInt();
		x=parseInt();
		y=parseInt();
	}
	
	/** Parse the comment field for pair number, chastity filter, and barcode */
	private void parseComment(){
		pos=commentSeparator+1;
		pairCode=parseChar();
		chastityCode=parseChar();
		controlBits=parseInt();
		barcode=id.substring(pos);
	}
	
	/** Clear all fields and point to a new header */
	private void reset(String id_){
		id=id_;
		limit=(id==null ? -1 : id.length());
		pos=-1;
		commentSeparator=-1;

		lane=-1;
		tile=-1;
		x=-1;
		y=-1;
		
		pairCode='?';
		barcode=null;
		chastityCode='?';
		controlBits=-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Find the first instance of the comment separator (' ' or '/')
	 * @return Position of the separator, or header length if not found
	 */
	private int findCommentSeparator(){
		for(int i=0; i<limit; i++){
			char c=id.charAt(i);
			if(c==' ' || c== '/'){return i;}
		}
		return limit;
	}
	
	/**
	 * Decrement pos until target colons have been encountered
	 * @modifies pos
	 * @param target number of colons to backtrack
	 */
	private void goBackSeveralColons(int target){
		for(int colons=0; pos>=0; pos--){
			if(id.charAt(pos)==':'){
				colons++;
				if(colons==target){break;}
			}
		}
		pos++;
	}
	
	/**
	 * Parse an integer from the current location in the string, 
	 * and advance the position to the beginning of the next number.
	 * @modifies pos
	 * @return The parsed number
	 */
	private int parseInt(){
		int current=0;
		assert(Tools.isDigit(id.charAt(pos))) : id;
		while(pos<limit && Tools.isDigit(id.charAt(pos))){
			current=current*10+(id.charAt(pos)-'0');
			pos++;
		}
		pos++;
		return current;
	}
	
	/**
	 * Parse a character from the current location in the string, 
	 * and advance the position to the next code character.
	 * @modifies pos
	 * @return The parsed char
	 */
	private char parseChar(){
		char c=id.charAt(pos);
		assert(c!=':');
		pos++;
		assert(pos>=limit || id.charAt(pos)==':');
		pos++;
		assert(pos>=limit || id.charAt(pos)!=':');
		return c;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Flowcell lane number */
	public int lane;
	/** Tile number */
	public int tile;
	/** X-coordinate within tile (pixels) */
	public int x;
	/** Y-coordinate within tile (pixels) */
	public int y;
	
	/** '1' for read 1, '2' for read 2 */
	public char pairCode;
	/** Read barcode */
	public String barcode;
	/** 'Y' for fail, 'N' for pass */
	public char chastityCode;
	/** A number for control bits set */
	public int controlBits;
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Read header */
	private String id;
	/** Length of header */
	private int limit;
	/** Current position, typically the start of the next token */
	private int pos;
	/** Position after coordinates, but before pairnum; typically space or slash */
	private int commentSeparator;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse lane, tile, x, and y coordinates */
	public static boolean PARSE_COORDINATES=true;
	/** Parse the comment field for pair number, chastity filter, and barcode */
	public static boolean PARSE_COMMENT=false;
	
}
