package fun;

import shared.Tools;

public class Palindrome {
	
	public static void main(String[] args){
		System.out.println(longestPalindrome(args[0]));
	}
	
	public static String longestPalindrome(String s){
		int longestLength=0;
		int longestStart=0;
		for(int i=0; i<s.length(); i++){
			int lenEven=palindromeLength(s, i, i+1);
			if(lenEven>longestLength){
				longestLength=lenEven;
				longestStart=i-lenEven/2+1;
			}
			int lenOdd=palindromeLength(s, i, i);
			if(lenOdd>longestLength){
				longestLength=lenOdd;
				longestStart=i-lenOdd/2;
			}
		}
		return s.substring(longestStart, longestStart+longestLength+1);
	}
	
	public static int palindromeLengthOdd(String s, int middle){
		int length=1;
		int a=middle-1, b=middle+1;
		while(a>=0 && b<s.length()){
			if(s.charAt(a)==s.charAt(b)){
				a--;
				b++;
				length+=2;
			}else{
				break;
			}
		}
		return length;
	}
	
	public static int palindromeLengthEven(String s, int middle){
		int length=0;
		int a=middle, b=middle+1;
		while(a>=0 && b<s.length()){
			if(s.charAt(a)==s.charAt(b)){
				a--;
				b++;
				length+=2;
			}else{
				break;
			}
		}
		return length;
	}
	
	public static int palindromeLength(String s, int a, int b){
		while(a>=0 && b<s.length()){
			if(s.charAt(a)!=s.charAt(b)){break;}
			a--;
			b++;
		}
		return Tools.max(0, b-a-2);
	}
	
	public static boolean isPalindrome(String s, int a, int b){
		while(a<b){
			if(s.charAt(a)!=s.charAt(b)){
				return false;
			}
		}
		return true;
	}
	
}
