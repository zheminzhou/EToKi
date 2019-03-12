package json;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import structures.ByteBuilder;

public class JsonObject {

	public static void main(String[] args){
		JsonObject bob=new JsonObject("name", "bob");
		JsonObject joe=new JsonObject("name", "joe");
		JsonObject sue=new JsonObject("name", "sue");
		JsonObject dan=new JsonObject("name", "dan");
		bob.add("joe", joe);
		bob.add("sue", sue);
		joe.add("dan", dan);
		bob.add("a",1);
		bob.add("b",2);
		bob.add("c","3");
		bob.add("a","4");
		dan.add("e",5);
		dan.add("f","6");
		sue.add("g","7");

		System.out.println("dan:\n"+dan);
		System.out.println("sue:\n"+sue);
		System.out.println("joe:\n"+joe);
		System.out.println("bob:\n"+bob);
		
		ArrayList<JsonObject> list=new ArrayList<JsonObject>();
		list.add(joe);
		list.add(sue);
		list.add(dan);
		System.out.println("list:\n"+toString(list));
	}
	
	public JsonObject(){}
	
	public JsonObject(String key, Object value){
		add(key, value);
	}
	
//	public JsonObject(String name_){
//		name=name_;
//	}
//	
//	public JsonObject(String name_, String key, Object value){
//		name=name_;
//		add(key, value);
//	}

	/** Does not add quotes for strings */
	public void addLiteral(String key0, double value, int decimals){
		omap.put(key0, new JsonLiteral(value, decimals));
	}

	/** Does not add quotes for strings.
	 * This method should be used with caution as it can produce incorrectly formatted files. */
	public void addLiteral(String key0, String value){
		omap.put(key0, new JsonLiteral(value));
	}

	public void add(String key0, Object value){
		if(value!=null && value.getClass()==JsonObject.class){
			add(key0, (JsonObject)value);
			return;
		}
		int x=2;
		String key=key0;
		if(omap==null){omap=new LinkedHashMap<String, Object>(8);}
		while(omap.containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		omap.put(key, value);
	}

	public void add(final String key0, JsonObject value){
		int x=2;
		String key=key0;
		if(jmap==null){jmap=new LinkedHashMap<String, JsonObject>(8);}
		while(jmap.containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		jmap.put(key, value);
	}
	
	public static String toString(ArrayList<JsonObject> list) {
		ByteBuilder sb=new ByteBuilder();
		int commas=list.size()-1;
		for(JsonObject j : list){
			j.append(0, sb, false);
			if(commas>0){
				sb.append(",\n");
			}
			commas--;
		}
		return sb.toString();
	}
	
	public ByteBuilder toText(){
		return toText(null, 0, false);
	}
	
	public ByteBuilder toText(ByteBuilder sb, int level, boolean inArray){
		if(sb==null){sb=new ByteBuilder();}
		append(level, sb, inArray);
		return sb;
	}
	
	public String toString(String name){
		ByteBuilder sb=new ByteBuilder();
		sb.append('{').append('\n');
		for(int i=0; i<padmult; i++){sb.append(' ');}
		sb.append('"').append(name).append('"').append(':').append(' ');
		toText(sb, 1, false);
		sb.append('\n').append('}');
		return sb.toString();
	}
	
	public static String toString(Object[] array){
		ByteBuilder sb=new ByteBuilder();
		appendArray(sb, array, 0);
		return sb.toString();
	}
	
	@Override
	public String toString(){
		return toText(null, 0, false).toString();
	}
	
	public void append(int level, ByteBuilder sb, boolean inArray){
		int pad=padmult*level;
		int pad2=padmult*(level+1);
		
		sb.append('{');
		if(!inArray){sb.append('\n');}
		
		int commas=(omap==null ? 0 : omap.size())+(jmap==null ? 0 : jmap.size())-1;
		
		if(omap!=null){
			for(Entry<String, Object> e : omap.entrySet()){
				String key=e.getKey();
				Object value=e.getValue();
				if(!inArray){for(int i=0; i<pad2; i++){sb.append(' ');}}
				
				appendEntry(sb, key, value, level, inArray);

				if(commas>0){sb.append(',');}
				if(!inArray){sb.append('\n');}
				commas--;
			}
		}
		
		if(jmap!=null){
			for(Entry<String, JsonObject> e : jmap.entrySet()){
				String key=e.getKey();
				JsonObject value=e.getValue();
				if(!inArray){for(int i=0; i<pad2; i++){sb.append(' ');}}
				appendKey(sb, key);
				
				value.append(level+(inArray ? 0 : 1), sb, inArray);
				if(commas>0){sb.append(',');}
				if(!inArray){sb.append('\n');}
				commas--;
			}
		}
		
		if(!inArray){for(int i=0; i<pad; i++){sb.append(' ');}}
		sb.append('}');
	}
	
	private static void appendEntry(ByteBuilder sb, String key, Object value, int level, boolean inArray){
		appendKey(sb, key);
		appendValue(sb, value, level, inArray);
	}
	
	private static void appendKey(ByteBuilder sb, String key){
		sb.append('"').append(key).append("\": ");
	}
	
	private static void appendValue(ByteBuilder sb, Object value, int level, boolean inArray){
		final Class<?> c=(value==null ? null : value.getClass());
		if(c==null || value==null){
			sb.append("null");
		}else if(c==String.class){
			sb.append("\"").append(value.toString()).append('"');
		}else if(c==Double.class && decimals>=0){
			sb.append(((Double)value).doubleValue(), decimals);
		}else if(c==Float.class && decimals>=0){
			sb.append(((Float)value).floatValue(), decimals);
		}else if(c==JsonObject.class){
			((JsonObject)value).append(level+(inArray ? 0 : 1), sb, inArray);
		}else if(c.isArray()){
			appendArray(sb, (Object[])value, level);
		}else if(c==Boolean.class || value instanceof Number){//long, int, boolean
			sb.append(value.toString());
		}else if(value instanceof Collection){
			appendCollection(sb, (Collection<?>)value, level);
		}else{//Default behavior for unhandled classes
			sb.append("\"").append(value.toString()).append('"');
		}
	}
	
	private static void appendArray(ByteBuilder sb, Object[] array, int level){
		int commas=(array==null ? 0 : array.length)-1;
		sb.append('[');
		if(array!=null){
			for(Object value : array){
				appendValue(sb, value, level, noNewlinesInArrays);
				if(commas>0){sb.append(',').append(' ');}
				commas--;
			}
		}
		sb.append(']');
	}
	
	private static void appendCollection(ByteBuilder sb, Collection<?> stuff, int level){
		int commas=(stuff==null ? 0 : stuff.size())-1;
		sb.append('[');
		if(stuff!=null){
			for(Object value : stuff){
				appendValue(sb, value, level, noNewlinesInArrays);
				if(commas>0){sb.append(',').append(' ');}
				commas--;
			}
		}
		sb.append(']');
	}

	public String getString(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==String.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (String)o;
	}

	public Long getLong(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Long.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Long)o;
	}

//	public Double getDouble(String key){
//		if(smap==null){return null;}
//		Object o=smap.get(key);
//		if(o==null){return null;}
//		assert(o.getClass()==Double.class) : "Wrong class: "+o.getClass()+"\n"+o;
//		return (Double)o;
//	}

	public Double getDouble(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		if(o.getClass()==Long.class){
			return ((Long)o).doubleValue();
		}
		assert(o.getClass()==Double.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Double)o;
	}

	public Number getNumber(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		Class<?> c=o.getClass();
		assert(c==Double.class || c==Long.class) : "Wrong class: "+c+"\n"+o;
		return (Number)o;
	}

	public Object[] getArray(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Object[].class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Object[])o;
	}

	public JsonObject getJson(String key){
		if(jmap==null){return null;}
		return jmap.get(key);
	}
	
	public int jmapSize(){return jmap==null ? 0 : jmap.size();}
	public int omapSize(){return omap==null ? 0 : omap.size();}
	
//	public String name;
	public LinkedHashMap<String, Object> omap;
	public LinkedHashMap<String, JsonObject> jmap;

	private static int decimals=-1;
	private static String decimalFormat="%."+decimals+"f";
	public static synchronized void setDecimals(int d){
		if(d!=decimals){
			d=decimals;
			decimalFormat="%."+decimals+"f";
		}
	}
	
	public static final int padmult=3;
	public static boolean noNewlinesInArrays=true;
	
}
