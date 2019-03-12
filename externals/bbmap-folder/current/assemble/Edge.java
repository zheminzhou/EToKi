package assemble;

public class Edge {
	
	public Edge(int origin_, int destination_, int length_, int orientation_){
		origin=origin_;
		destination=destination_;
		length=length_;
		orientation=orientation_;
	}
	
	int origin;
	int destination;
	int length;
	int orientation; //0 left kmer, 1 left rkmer, 2 right kmer, 3 right rkmer (of dest)
	
}
