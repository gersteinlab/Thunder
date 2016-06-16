package utils;

import java.text.SimpleDateFormat;
import java.util.Calendar;

public class IO_utils {


	public static String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}
	
	public static void printOut(String message){ System.out.print(getTime()+" "+message); }
	public static void printErr(String message){ System.err.print(getTime()+" "+message); }
	public static void printLineOut(String message){ System.out.println(getTime()+" "+message); }
	public static void printLineErr(String message){ System.err.println(getTime()+" "+message); }
	
	public static void printProgressBar(int percent){
		printProgressBar(percent+0.0);
	}
	public static void printProgressBar(double percent){
		percent = Math.round(percent*100.0)/100.0;
	    StringBuilder bar = new StringBuilder("[");
	    for(int i = 0; i < 50; i++){
	        if( i < (percent/2)){
	            bar.append("=");
	        }else if( i == (percent/2)){
	            bar.append(">");
	        }else{
	            bar.append(" ");
	        }
	    }
	    bar.append("]   " + percent + "%     ");
	    printErr(bar.toString()+"\r");
	    //if(percent == 100.0)
	    //	System.err.println();
	}
	
}
