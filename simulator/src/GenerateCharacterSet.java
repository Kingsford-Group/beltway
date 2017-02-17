import java.io.*;
import java.math.*;
import java.util.stream.*;
import java.util.concurrent.ThreadLocalRandom;

/* 
 * This generates n-sized characterset as well as associated weights for each entry in the set.
 *
 */
public class GenerateCharacterSet{
    
    /* 
     * @param n       : number of distinct strings to be built
     * @param minW    : minimum weight (inclusive) to be assigned to a character 
     * @param maxW    : maximum weight (exclusive) to be assigned to a character
     * @param intOnly : force the weights to be integers if true. 
     * @param outFile : output file name 
     */
    public static void generate(int n, double minW, double maxW, String outFile, String outDir, boolean intOnly){
	String expDir = (outDir.endsWith(File.separator) ? outDir : outDir + File.separator);
	if(intOnly){
	    int[] intWs = ThreadLocalRandom.current().ints((int) minW, (int) maxW).limit(n).toArray();
	    outToFile(intWs, expDir + outFile);
	}else{
	    double[] weights = ThreadLocalRandom.current().doubles(minW , maxW).limit(n).toArray();
	    outToFile(weights, expDir + outFile);
	}
    }
    
    /* DOUBLE VERSION
     * @param weights : array of n weights 
     * @param f       : output file name
     */
    public static void outToFile(double[] weights, String f){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(f));
	    for(int i=0; i < weights.length; i++)
		bw.write(charPrefix + i + "\t" + new BigDecimal(weights[i]).setScale(4, RoundingMode.DOWN) + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    System.exit(1);
	}  
    }
    
    /* INT VERSION
     * @param weights : array of n weights 
     * @param f       : output file name
     */
    public static void outToFile(int[] weights, String f){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(f));
	    for(int i=0; i < weights.length; i++)
		bw.write(charPrefix + i + "\t" + new BigDecimal(weights[i]).setScale(4, RoundingMode.DOWN) + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    System.exit(1);
	}
    }
    
    public static void main(String[] args){
	if(args.length == 6)
	    generate(Integer.parseInt(args[0]), Double.parseDouble(args[1]), Double.parseDouble(args[2]), args[3], args[4], (args[5].equals("Y") || args[5].equals("y") ? true : false ));
	else if(args.length == 5)
	    generate(Integer.parseInt(args[0]), Double.parseDouble(args[1]), Double.parseDouble(args[2]), args[3], args[4], false);
	else{
	    System.err.println("java GenerateCharacterSet <n> <minW> <maxW> <outFile> <oDir> [YyNn: INT?]");
	    System.err.println("Parameters:");
	    System.err.println("\tn       [INT]  \tSize of character set");
	    System.err.println("\tminW    [FLOAT]\tMinimum weight allowed");
	    System.err.println("\tmaxW    [FLOAT]\tMaximum weight allowed");
	    System.err.println("\toutFile [STR]  \tFilename to store a generated set of characters and their associated weights");
	    System.err.println("\toDir    [STR]  \tOutput directory (ex: data/int_simul1)");
	    System.err.println("\tInt?    [YyNn] \t Force weights to be integers. [default:N]");
	}
    }

    public static String charPrefix = "AA";
    
    public static String expDir;
}
