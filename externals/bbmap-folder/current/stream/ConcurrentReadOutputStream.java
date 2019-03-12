package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import shared.Shared;

/**
 * @author Brian Bushnell
 * @date Jan 26, 2015
 *
 */
public abstract class ConcurrentReadOutputStream {
	
	/*--------------------------------------------------------------*/
	/*----------------           Factory            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, int maxSize, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, null, null, null, maxSize, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, int maxSize, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, null, null, maxSize, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int maxSize, CharSequence header, boolean useSharedHeader){
		return getStream(ff1, ff2, qf1, qf2, maxSize, header, useSharedHeader, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}
	
	public static ConcurrentReadOutputStream getStream(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
			int maxSize, CharSequence header, boolean useSharedHeader, final boolean mpi, final boolean keepAll){
		if(mpi){
			final int rank=Shared.MPI_RANK;
			final ConcurrentReadOutputStream cros0;
			if(rank==0){
				cros0=new ConcurrentGenericReadOutputStream(ff1, ff2, qf1, qf2, maxSize, header, useSharedHeader);
			}else{
				cros0=null;
			}
			final ConcurrentReadOutputStream crosD;
			if(Shared.USE_CRISMPI){
				assert(false) : "To support MPI, uncomment this.";
				crosD=null;
//				crosD=new ConcurrentReadOutputStreamMPI(cros0, rank==0);
			}else{
				crosD=new ConcurrentReadOutputStreamD(cros0, rank==0);
			}
			return crosD;
		}else{
			return new ConcurrentGenericReadOutputStream(ff1, ff2, qf1, qf2, maxSize, header, useSharedHeader);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	ConcurrentReadOutputStream(FileFormat ff1_, FileFormat ff2_){
		ff1=ff1_;
		ff2=ff2_;
		ordered=(ff1==null ? true : ff1.ordered());
	}
	
	public abstract void start();
	
	public final boolean started(){return started;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public abstract void add(ArrayList<Read> list, long listnum);
	
	public abstract void close();
	
	public abstract void join();
	
	public abstract void resetNextListID();
	
	public abstract String fname();
	
	/** Return true if this stream has detected an error */
	public abstract boolean errorState();

	public abstract boolean finishedSuccessfully();
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public long basesWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.basesWritten();}
		if(rsw2!=null){x+=rsw2.basesWritten();}
		return x;
	}
	
	public long readsWritten(){
		long x=0;
		ReadStreamWriter rsw1=getRS1();
		ReadStreamWriter rsw2=getRS2();
		if(rsw1!=null){x+=rsw1.readsWritten();}
		if(rsw2!=null){x+=rsw2.readsWritten();}
		return x;
	}
	
	public abstract ReadStreamWriter getRS1();
	public abstract ReadStreamWriter getRS2();
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final FileFormat ff1, ff2;
	public final boolean ordered;
	
	boolean errorState=false;
	boolean finishedSuccessfully=false;
	boolean started=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;
	
}
