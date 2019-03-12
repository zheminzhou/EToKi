package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.GffLine;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

public class VcfWriter {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	public VcfWriter(VarMap varMap_, VarFilter filter_, long reads_,
			long pairs_, long properPairs_, long bases_, String ref_, boolean trimWhitespace_, String sampleName_){
		
		varMap=varMap_;
		filter=filter_;
		trimWhitespace=trimWhitespace_;
		
		reads=reads_;
		pairs=pairs_;
		properPairs=properPairs_;
		bases=bases_;
		ref=ref_;
		sampleName=sampleName_;
		
		ploidy=varMap.ploidy;
		properPairRate=varMap.properPairRate;
		pairedInSequencingRate=varMap.pairedInSequencingRate;
		totalQualityAvg=varMap.totalQualityAvg;
		totalMapqAvg=varMap.totalMapqAvg;
		readLengthAvg=varMap.readLengthAvg;
		scafMap=varMap.scafMap;
		
		threads=Tools.max(1, Shared.threads());
		
		array=varMap.toArray(true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean writeVcf(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVcfFile(fname);
		return vw.errorState;
	}
	
	public static boolean writeVar(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVarFile(fname);
		return vw.errorState;
	}
	
	public static boolean writeGff(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeGffFile(fname);
		return vw.errorState;
	}

	public void writeVcfFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VCF, "vcf", true, true, false, true);
		writeFile(ff, VCFMODE);
	}

	public void writeVarFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VAR, "var", true, true, false, true);
		writeFile(ff, VARMODE);
	}

	public void writeGffFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.GFF, "gff", true, true, false, true);
		writeFile(ff, GFFMODE);
	}

	public void writeVcfFile(final FileFormat ff){
		assert(ff.vcf());
		writeFile(ff, VCFMODE);
	}

	public void writeVarFile(final FileFormat ff){
		assert(ff.var());
		writeFile(ff, VARMODE);
	}

	public void writeGffFile(final FileFormat ff){
		assert(ff.gff());
		writeFile(ff, GFFMODE);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeFile(final FileFormat ff, final int mode){
		assert(ff.ordered());
		final ArrayBlockingQueue<ListNum<Var>> inq=new ArrayBlockingQueue<ListNum<Var>>(threads+1);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		
		writeHeader(bsw, mode);
		
		ArrayList<ProcessThread> alpt=spawnThreads(bsw, inq, mode);
		
		makeLists(inq);
		
		waitForFinish(alpt);
		
		errorState=bsw.poisonAndWait()|errorState;
	}
	
	private void writeHeader(ByteStreamWriter bsw, int mode){
		ByteBuilder bb=new ByteBuilder(1000);
		if(mode==VCFMODE){
			bb.append(Var.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref, scafMap, sampleName, trimWhitespace)).append('\n');
		}else if(mode==VARMODE){
			bb.append(Var.toVarHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else if(mode==GFFMODE){
			bb.append(GffLine.toHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else{assert(false);}
		bsw.add(bb, 0);
	}

	
	/** Spawn process threads */
	private ArrayList<ProcessThread> spawnThreads(ByteStreamWriter bsw, final ArrayBlockingQueue<ListNum<Var>> inq, int mode){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, bsw, inq, mode));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			if(verbose){outstream.println("Waiting for thread "+pt.tid);}
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			allSuccess&=pt.success;
		}
		
		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}
	
	void makeLists(ArrayBlockingQueue<ListNum<Var>> inq){
		ArrayList<Var> list=new ArrayList<Var>(LIST_SIZE);
		long nextJobID=1;
		for(Var v : array){
			list.add(v);
			if(list.size()>=LIST_SIZE){
				putVars(new ListNum<Var>(list, nextJobID), inq);
				nextJobID++;
				list=new ArrayList<Var>(LIST_SIZE);
			}
		}
		if(list.size()>0){
			putVars(new ListNum<Var>(list, nextJobID), inq);
			nextJobID++;
			list=null;
		}
		if(verbose){outstream.println("tid "+0+" done making var lists.");}
		putVars(POISON_VARS, inq);
		if(verbose){outstream.println("tid "+0+" done poisoning.");}
	}
	
	final void putVars(ListNum<Var> list, final ArrayBlockingQueue<ListNum<Var>> inq){
		if(verbose){outstream.println("tid "+0+" putting vlist "+list.id+", size "+list.size());}
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose){outstream.println("tid "+0+" done putting vlist");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread{
		
		public ProcessThread(int tid_, ByteStreamWriter bsw_, ArrayBlockingQueue<ListNum<Var>> inq_, int mode_){
			tid=tid_;
			bsw=bsw_;
			inq=inq_;
			mode=mode_;
		}
		
		@Override
		public void run(){
			makeBytes();
		}
		
		void makeBytes(){
			if(verbose){outstream.println("tid "+tid+" started makeBytes.");}
			
			ListNum<Var> list=takeVars();
			final ByteBuilder bb=new ByteBuilder();
			while(list!=null && list.size()>0 && list!=POISON_VARS){
				if(mode==VCFMODE){
					for(Var v : list){
						v.toVCF(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, trimWhitespace);
						bb.nl();
					}
				}else if(mode==VARMODE){
					for(Var v : list){
						v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);//TODO: Track depth
						bb.nl();
					}
				}else if(mode==GFFMODE){
					for(Var v : list){
						GffLine.toText(bb, v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);
						bb.nl();
					}
				}
				bsw.add(new ByteBuilder(bb.toBytes()), list.id);
				list=takeVars();
				
				bb.clear();
				bb.shrinkTo(SHRINK_SIZE);
			}
			if(list==POISON_VARS){
				putVars(POISON_VARS, inq);
			}
			if(verbose){outstream.println("tid "+tid+" done making bytes.");}
		}
		
		final ListNum<Var> takeVars(){
			if(verbose){outstream.println("tid "+tid+" taking vlist");}
			ListNum<Var> list=null;
			while(list==null){
				try {
					list=inq.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" took vlist "+list.id+", size "+list.size());}
			return list;
		}
		
		final ArrayBlockingQueue<ListNum<Var>> inq;
		final int tid;
		final ByteStreamWriter bsw;
		final int mode;
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	final Var[] array;
	final int ploidy;
	final double properPairRate;
	final double pairedInSequencingRate;
	final double totalQualityAvg;
	final double totalMapqAvg;
	final double readLengthAvg;
	final ScafMap scafMap;
	final VarMap varMap;
	final VarFilter filter;
	final boolean trimWhitespace;
	
	final String sampleName;
	final long reads;
	final long pairs;
	final long properPairs;
	final long bases;
	final String ref;
	
	final int threads;
	
	boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int VARMODE=0, VCFMODE=1, GFFMODE=2;
	private static final ListNum<Var> POISON_VARS=new ListNum<Var>(null, -1);
	
	private static boolean verbose=false;

	private static final int LIST_SIZE=200;
	private static final int SHRINK_SIZE=1000*LIST_SIZE;
	
	/** Print status messages to this output stream */
	protected static PrintStream outstream=System.err;
	
}
