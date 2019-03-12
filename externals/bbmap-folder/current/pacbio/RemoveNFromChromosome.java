package pacbio;

import java.io.File;

import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * @author Brian Bushnell
 * @date Jul 19, 2012
 *
 */
public class RemoveNFromChromosome {
	
	public static void main(String[] args){
		int ingenome=Integer.parseInt(args[0]);
		int outgenome=Integer.parseInt(args[1]);
		int padding=Integer.parseInt(args[2]);
		
		String outRoot=Data.ROOT_GENOME+outgenome+"/";
		File f=new File(outRoot);
		if(!f.exists()){
			f.mkdirs();
		}
		
		Data.setGenome(ingenome);
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			ChromosomeArray cha=Data.getChromosome(chrom);
			Data.unload(chrom, true);
			ChromosomeArray chb=new ChromosomeArray(chrom, Shared.PLUS, 0, cha.countDefinedBases()+2*padding+1);
			chb.maxIndex=-1;
			for(int i=0; i<padding; i++){
				chb.set(i, 'N');
			}
			for(int i=0; i<cha.maxIndex; i++){
				byte b=cha.get(i);
				if(b!='N'){
					chb.set(chb.maxIndex+1, b);
				}
			}
			for(int i=0; i<padding; i++){
				chb.set(chb.maxIndex+1, 'N');
			}
			ReadWrite.write(chb, outRoot+"chr"+chrom+Data.chromExtension(), false);
			
		}
		
		FastaToChromArrays2.writeInfo(outgenome, Data.numChroms, Data.name, Data.genomeSource, false, false);
		
	}
	
}
