/**
 * 
 */
package hangman;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * @author raag
 */

public class FindDup {
 
  private static BufferedReader in;

public static void countDupChars(String str){
 
    //HashMap 
    Map<Character, Integer> map = new HashMap<Character, Integer>(); 
 
    // String to char array
    char[] chars = str.toCharArray();

        
    for(Character ch:chars){
      if(map.containsKey(ch)){
         map.put(ch, map.get(ch)+1);
      } else {
         map.put(ch, 1);
        }
    }
    int len=str.length();
    
    //Retrieve set of keys
    Set<Character> keys = map.keySet();
    
    DecimalFormat df2 = new DecimalFormat(".##");
    df2.setRoundingMode(RoundingMode.UP);
    
    for(Character ch:keys){
      if(map.get(ch) > 1){
        System.out.println("Char "+ch+" "+ df2.format((double) map.get(ch)/len*100)+"%");
      }
    }
  }
 
 public static void main(String a[]) throws Exception{
	
	String fn ="";

	if (a.length>0) fn=a[0];
	fn = (new BufferedReader(new InputStreamReader(System.in))).readLine();
	
	in = new BufferedReader(new FileReader( fn ) );
    
    String line   = in.readLine();
    
    if( line.equals(null))
    	 {
        throw new IOException( fn + " is an empty file" );
         }
    else {
    	 if( line.length()>0)
    		 countDupChars(line);
    	}
    }
}
