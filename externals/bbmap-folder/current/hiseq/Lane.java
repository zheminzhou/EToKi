package hiseq;

import java.util.ArrayList;

public class Lane {
	
	public Lane(int lane_){
		lane=lane_;
	}
	
	public MicroTile getMicroTile(int tile, int x, int y){
		return getTile(tile).get(x, y);
	}
	
	public Tile getTile(int index){
		while(tiles.size()<=index){tiles.add(null);}
		Tile t=tiles.get(index);
		if(t==null){
			t=new Tile(lane, index);
			tiles.set(index, t);
		}
		return t;
	}
	
	public ArrayList<Tile> tiles=new ArrayList<Tile>();
	
	public int lane;
	
}
