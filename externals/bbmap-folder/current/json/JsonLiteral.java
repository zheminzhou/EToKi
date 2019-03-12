package json;

import java.util.Locale;

public class JsonLiteral {
	
	public JsonLiteral(String s_){
		s=s_;
	}
	
	public JsonLiteral(double value, int decimals){
		s=String.format(Locale.ROOT, "%."+decimals+"f", value);
	}
	
	@Override
	public String toString(){return s;}
	
	private final String s;
	
}
