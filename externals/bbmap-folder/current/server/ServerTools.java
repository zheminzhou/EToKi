package server;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.Map.Entry;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;

import fileIO.ReadWrite;
import shared.Tools;
import structures.ByteBuilder;
import structures.StringNum;

public class ServerTools {
	
	public static void main(String[] args){
		String address=args[0];
		int rounds=1;
		String message="";
		if(args.length>1){rounds=Integer.parseInt(args[1]);}
		if(args.length>2){message=args[2];}
		byte[] messageBytes=message.getBytes();
		
		long[] times=new long[rounds];
		StringNum response=null;
		long prevTime=System.nanoTime();
		for(int i=0; i<rounds; i++){
			response=sendAndReceive(messageBytes, address);
			long currentTime=System.nanoTime();
			times[i]=currentTime-prevTime;
			prevTime=currentTime;
			System.out.println(times[i]);
		}
		
		System.out.println(response.s);
		
		Arrays.sort(times);
		long sum=Tools.sum(times);
		System.out.println("Avg:    \t"+sum/1000000.0+" ms");
		System.out.println("QPS:    \t"+(rounds*1000000000/sum)+" ms");
		System.out.println("Median: \t"+(times[rounds/2]/1000000.0)+" ms");
		
	}
	
	public static ByteBuilder readPage(String address, boolean convert){
    	if(convert){address=PercentEncoding.commonSymbolToCode(address);}
//    	assert(false) : address;
		ByteBuilder bb=new ByteBuilder(256);
		boolean success=false;
		for(int i=0; i<10 && !success; i++){
			try {
				URL url=new URL(address);
				InputStream is=url.openStream();

				byte[] buffer=new byte[4096];
				for(int len=is.read(buffer); len>0; len=is.read(buffer)){
					bb.append(buffer, 0, len);
				}
				is.close();
				success=true;
			} catch (MalformedURLException e) {
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			} catch (IOException e) {
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			}
		}
        return bb;
    }
	
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceive(byte[] message, String address){
    	address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		HttpURLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=(HttpURLConnection) url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(40000); //For testing
			os=connection.getOutputStream();
		} catch (IOException e1) {
//			System.err.println("A:\t"+address+" -> "+url+" -> "+connection+" -> "+os);
			// TODO Auto-generated catch block
			if(!suppressErrors){e1.printStackTrace();}
		}
		
		if(os!=null){
			try {
				//TODO: It may be useful to set a timeout here.
				if(message!=null){os.write(message);}
			} catch (IOException e) {
//				System.err.println("B:\t"+connection+" -> "+os);
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
			try {
				os.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		String result=null;
		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try {
				is=connection.getInputStream();
			} catch (IOException e) {
				noInputStream=e;
			}
			if(is==null){is=connection.getErrorStream();}
			
			try {
				responseCode=connection.getResponseCode();
				if(!suppressErrors && (responseCode<200 || responseCode>299)){
					System.err.println("Error: Server returned response code "+responseCode);
				}
			} catch (IOException e) {
				if(!suppressErrors) {
					e.printStackTrace();
				}
			}
			
			if(is!=null){
				result=readStream(is);
				try {
//					System.err.println("C:\t"+connection+" -> "+os+" -> "+is+" -> "+(result!=null));
					is.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		
		return new StringNum(result, responseCode);
	}

	/** Read the body of an incoming HTTP session */
	public static String receive(HttpExchange t){
		InputStream is=t.getRequestBody();
		String s=readStream(is);
		return s;
	}
	
	/** Completely read an InputStream into a String */
	public static String readStream(InputStream is){
		if(is==null){return null;}
		try {
			byte[] buffer=new byte[256];
			int count=is.read(buffer);
			int next=0;
			while(count>-1){
				next+=count;
				if(next>=buffer.length){
					buffer=Arrays.copyOf(buffer, buffer.length*2);
				}
				count=is.read(buffer, next, buffer.length-next);
			}
			is.close();
			return new String(buffer, 0, next);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean reply(String response, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}

		return reply(response.getBytes(), type, t, false, code, close);
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean replyWithFile(String path, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+path);}
		
		byte[] response=null;
		try {
			response=ReadWrite.readRaw(path);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			response=new byte[0];
			code=400; //Node sure about this
		}
		
		return reply(response, type, t, false, code, close);
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean reply(byte[] response, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}
		
		{
			Headers h = t.getResponseHeaders();
//			String type="text/plain";
			h.add("Content-Type", type);
		}
		try {
			t.sendResponseHeaders(code, response.length);
			OutputStream os = t.getResponseBody();
			os.write(response);
			os.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			if(close){t.close();}
			return false;
		}
		if(close){t.close();}
		return true;
	}
	
	
	/**
	 * Wait for a set amount of time
	 * @param millis Time to wait
	 */
	public static void pause(long millis){
		Integer lock=new Integer(1);
		synchronized(lock){
			final long time=System.currentTimeMillis()+millis;
			while(System.currentTimeMillis()<time){
				try {
					lock.wait(Tools.max(100, time-System.currentTimeMillis()));
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	public static String getClientAddress(HttpExchange t) {
		
		InetSocketAddress client=t.getRemoteAddress();
//		InetSocketAddress server=t.getLocalAddress();
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
//		String ls=server.toString();
		
		if(clientAddress.contains("127.0.0.1")){
			Headers clientRequestHeaders=t.getRequestHeaders();
//			Headers resh=t.getResponseHeaders();
	
			String xff=clientRequestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){clientAddress=xff;}
			
//			System.err.println("\nRequest: ");
//			for(Entry<String, List<String>> entry : reqh.entrySet()){
//				System.err.println(entry.getKey()+" -> "+entry.getValue());
//			}
		}
		return clientAddress;
	}
	
	public static boolean isInternalQuery(HttpExchange t, String prefix, boolean allowLocalHost, boolean printIP, boolean printHeaders){
		
		InetSocketAddress client=t.getRemoteAddress();
		InetSocketAddress server=t.getLocalAddress();
		
		if(printIP){System.err.println(client+"\t"+server);}
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
		String serverAddress=server.toString();
		
		if(clientAddress.contains("127.0.0.1")){//TODO: contains versus startsWith?
			Headers requestHeaders=t.getRequestHeaders();
			
			if(printHeaders){
				Headers responseHeaders=t.getResponseHeaders();
				System.err.println("\nRequest: ");
				for(Entry<String, List<String>> entry : requestHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
				System.err.println("\nResponse: ");
				for(Entry<String, List<String>> entry : responseHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
			}
			
			final String xff=requestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){
				if(xff.startsWith(prefix)){return true;}
				clientAddress=xff;
			}else{
				return allowLocalHost;
			}
		}else{
			if(clientAddress.startsWith(prefix)){return true;}
		}
		
		//Makes sure they match up to the first delimiter
		//TODO: This needs to be reviewed
		for(int i=0, max=Tools.max(clientAddress.length(), serverAddress.length()); i<max; i++){
			char cc=clientAddress.charAt(i), sc=serverAddress.charAt(i);
			if(cc!=sc){break;}
			if(cc=='.'){//IPv4
				return true; 
			}else if(cc==':'){//IPv6; probably depends on how long the mask is
				return true;
			}
		}
	
		return false;
	}
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
}
