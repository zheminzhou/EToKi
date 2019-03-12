package assemble;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Thread for exploring connectivity graph between contigs.
 * @author Brian Bushnell
 * @date July 12, 2018
 *
 */
public abstract class AbstractProcessContigThread extends Thread {

	AbstractProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
		contigs=contigs_;
		next=next_;
	}
	
	@Override
	public void run(){
		processContigs(contigs);
	}

	public final void processContigs(ArrayList<Contig> contigs){
		for(int cnum=next.getAndIncrement(); cnum<contigs.size(); cnum=next.getAndIncrement()){
			Contig c=contigs.get(cnum);
			processContigLeft(c, leftCounts, rightCounts, extraCounts);
			processContigRight(c, leftCounts, rightCounts, extraCounts);
		}
	}

	abstract void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts);

	abstract void processContigRight(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts);

	final int[] leftCounts=new int[4];
	final int[] rightCounts=new int[4];
	final int[] extraCounts=new int[4];

	final ArrayList<Contig> contigs;
	final AtomicInteger next;

	int lastLength=-1;
	int lastTarget=-1;
	int lastExitCondition=-1;
	int lastOrientation=-1;
	long edgesMadeT=0;

}
