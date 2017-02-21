import java.util.*;
import java.io.*;
import java.math.*;
import java.util.stream.*;
import java.util.concurrent.ThreadLocalRandom;


public class Simulator{

    /* constructor */
    
    public Simulator(){
	this.character2weight = new Hashtable<String, BigDecimal>();
	this.theoSpectrum = new ArrayList<Peptide>();
    }

    public Simulator(int cLen, int numSim, double sampRatio, int numChars, String theoSpectrumOutFile, double minWeight, double maxWeight, String charsetOutFile, String outDir, String useIntegerWeights){
	this();
	this.expDir = (outDir.endsWith(File.separator) ? outDir : outDir + File.separator);
	this.clen = cLen;
	this.prefixW = new BigDecimal[cLen+1];
	if(this.cyclopeptide == null){
	    System.err.println("cyclopeptide is NULL");
	}
	boolean intOnly = ( (useIntegerWeights.equals("Y") || useIntegerWeights.equals("y")) ? true : false);
	GenerateCharacterSet.generate(numChars, minWeight, maxWeight, charsetOutFile, this.expDir, intOnly);
	
	this.loadChar2WeightHash(charsetOutFile);
	this.loadCyclopeptide(cLen, numSim);
	
	this.generateTheoreticalSpectrum(theoSpectrumOutFile);
	this.sampleSpectrum(numSim, sampRatio, "cycoLen" + cLen);
    }
    
    public Simulator(int cLen, int numSim, double sampRatio, String theoSpectrumOutFile, String charsetFile, String outDir, String cpepf){
	this();
	System.err.println(cpepf);
	this.expDir = (outDir.endsWith(File.separator) ? outDir : outDir + File.separator);
	this.clen = cLen;
	this.prefixW = new BigDecimal[cLen+1];
	//this.cyclopeptide = cpep;
	this.loadChar2WeightHash(charsetFile);
	this.loadCyclopeptide(cLen, numSim, cpepf);

	if(cpepf == null)
	    this.generateTheoreticalSpectrum(theoSpectrumOutFile);
	else
	    this.loadTheoreticalSpectrum(theoSpectrumOutFile);
	this.sampleSpectrum(numSim, sampRatio, "cycloLen" + cLen);
    }
    
    private void loadTheoreticalSpectrum(String inF){
	
	BufferedReader br = null;
	try{
	    if(new File(inF).isFile()){
		br = new BufferedReader(new FileReader(inF));
		String curline = "";
		while((curline=br.readLine())!=null){
		    String[] tokens = curline.split("\\t");
		    this.theoSpectrum.add(new Peptide(tokens[1].split(Simulator.delim), new BigDecimal(tokens[0]).setScale(4, RoundingMode.DOWN)));
		}
	    }else{
		this.generateTheoreticalSpectrum(inF);
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    /*
     * generates theoretical spectrum for cyclopeptide
     *
     */
    private void generateTheoreticalSpectrum(String outF){
	String[] aas = this.cyclopeptide.getAAs();
	
	String[] tmpSeq;
	for(int i=0; i<this.clen; i++){
	    for(int j=i+1; j<=this.clen;j++){
		tmpSeq = new String[j-i];
		for(int k=i; k<j; k++)
		    tmpSeq[k-i] = aas[k];

		this.theoSpectrum.add(new Peptide(tmpSeq, this.prefixW[j].subtract(this.prefixW[i])));

		if(i > 0 && j < this.clen){
		    tmpSeq = new String[ aas.length-j + i ];
		    int offset = 0;
		    for(int k=j;k<aas.length;k++){
			tmpSeq[k-j] = aas[k];
			offset = k-j;
		    }for(int k=0;k<i;k++)
			tmpSeq[offset+1+k] = aas[k];
		    
		    this.theoSpectrum.add(new Peptide(tmpSeq, this.prefixW[this.prefixW.length-1].subtract(this.prefixW[j].subtract(this.prefixW[i]))));
		    
		}
	    }
	}
	
	//add all cycles of length l
	for(int i=1; i<this.clen;i++){
	    tmpSeq = new String[this.clen];
	    for(int k=i; k<this.clen; k++)
		tmpSeq[k-i] = aas[k];
	    for(int k=0; k<i; k++)
		tmpSeq[this.clen-i+k] = aas[k];
	    //spectrum.add(new Double(this.prefixW[this.prefixW.length-1]));
	    this.theoSpectrum.add(new Peptide(tmpSeq, this.prefixW[this.prefixW.length-1]));
	}
	
	//for(Peptide p : this.theoSpectrum){
	//    System.out.println(p.toString());
	//}
	
	Collections.sort(this.theoSpectrum);
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.expDir + outF));
	    for(Peptide p : this.theoSpectrum)
		bw.write(p.toString() + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    

    /*
     * This simply loads character2weight hashtable from a input file
     * containing characterset and its associated weights.
     *
     * @param char2weightf : 2-column tab-delimited file. Column 1: aminoacid name, Column2: its associated weight
     *
     */
    private void loadChar2WeightHash(String char2weightf){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(this.expDir + char2weightf));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		this.character2weight.put(tokens[0], new BigDecimal(tokens[1]).setScale(4, RoundingMode.DOWN));
	    }
	    br.close();
	}catch(IOException ioe){
	    System.err.println("Error parsing char2weight file : " + char2weightf);
	    ioe.printStackTrace();
	    System.exit(1);
	}
    }
    
    private boolean loadCyclopeptide(int l, int n){
	return loadCyclopeptide(l, n, null);
    }

    /*
     * If cyclopeptide sequence is not given, it generates a random cyclopeptide of length l 
     * using characterset defined in character2weight. Then it calls checkCyclopeptideAndLoadPrefixWeights(l) 
     * to generate all prefix peptide weights
     *
     * @param l : lenght of cyclopeptide
     * @param n : size of chracterset
     * @param cpef : cyclopeptide file
     * 
     */
    private boolean loadCyclopeptide(int l, int n, String cpepf){
	String cpep;
	if(cpepf == null){
	    System.err.println("cyclo null");
	    Random rand = new Random();
	    int[] indices = new int[l];
	    for(int i=0;i<indices.length;i++){
		indices[i] = rand.nextInt(this.character2weight.size());
	    }

	    for(int i =0; i<indices.length; i++){
		System.err.println("index\t" + indices[i]);
	    }
	    String[] chars = new String[0];
	    System.err.println("hashSize:\t" + this.character2weight.size());
	    chars = this.character2weight.keySet().toArray(chars);
	    System.err.println("charsLen:\t" + chars.length);
	    if(chars.length > 0 && indices.length >0){
		StringBuffer sb = new StringBuffer(chars[indices[0]]);
		for(int i=1; i<indices.length; i++)
		    sb.append(Simulator.delim + chars[indices[i]]);
		System.err.println("Setting cyclopeptide as:\t" + sb.toString());
		cpep = sb.toString();
		
		BufferedWriter bw = null;
		try{
		    bw = new BufferedWriter(new FileWriter(expDir + "cycloLen" + l +".txt"));
		    bw.write(cpep + "\n");
		    bw.close();
		}catch(IOException ioe){
		    ioe.printStackTrace();
		}

	    }else{
		System.err.println("cyclopeptide length [l] must be GREATER than 0");
		return false;
	    }
	}else
	    cpep = getCpepFromInputFile(cpepf);
	
	return this.checkCyclopeptideAndLoadPrefixWeights(l, cpep);
    }
    
    /* grabs the first line containing the cyclopeptide sequence from a input file */
    private String getCpepFromInputFile(String cpepf){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(cpepf));
	    String curline = br.readLine();;
	    if(curline == null){
		br.close();
		throw new Exception("cyclopeptide file [" + cpepf + "] EMPTY." );
	    }else{
		br.close();
		return curline;
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}catch(Exception e){
	    e.printStackTrace();
	}
	return null;
    }

    /*
     * @param l : length of cyclopeptide
     * @return  True if cyclopeptide sequence is a valid peptides 
     *          (length must be l and only amino acids from characterset are used. False, otherwise.
     */
    private boolean checkCyclopeptideAndLoadPrefixWeights(int l, String cpep){

	String[] aas = cpep.split(Simulator.delim);
	//check length
	if(aas.length != l){
	    System.err.println("The length of the input cyclopeptide sequence is NOT " + l);
	    return false;
	}
	
	/* generate prefix weights */
	int i = 0;
	this.prefixW[i] = new BigDecimal("0").setScale(4, RoundingMode.DOWN);
	for(String aa : aas){
	    i++;
	    BigDecimal w = character2weight.get(aa);
	    if(w == null){ 
		System.err.println("UNKNOWN amino acid [" + aa + "] present in cyclopeptide sequence");
		return false;
	    }
	    prefixW[i] = prefixW[i-1].add(w);
	}
	
	this.cyclopeptide = new Peptide(aas, prefixW[prefixW.length-1]);
	return true;
    }


    /*
     * randomly samples a fraction [samplingRatio] from theorectical spectrum in r repetitions.
     */
    private void sampleSpectrum(int r, double samplingRatio, String expName){
	//int[] indices = new int[this.theoSpectrum.size()];
	ArrayList<Integer> indices = new ArrayList<>();
	for(int i=0;i<this.theoSpectrum.size();i++)
	    indices.add(new Integer(i));
	if(samplingRatio > 0 && samplingRatio <= 1){
	    //Random rand = new Random();
	    for(int i=0; i<r; i++)
		this.sampleSpectrumSingle(samplingRatio, indices, expName + "_" + i);
	}
    }
    
    /*
     * randomly samples a fraction [samplingRatio] from theorectical spectrum.
     */
    private void sampleSpectrumSingle(double samplingRatio, ArrayList<Integer> indices, String filePrefix){
	StringBuffer bf = new StringBuffer("#n = " + this.clen + "\n");
	StringBuffer answerbf = new StringBuffer("#n = " + this.clen + "\n");
	int ns = (int) Math.ceil(this.theoSpectrum.size() * samplingRatio);
	Collections.shuffle(indices);
	//int[] indices = rand.ints(ns, 0, this.theoSpectrum.size()).toArray();
	
	for(int i=0; i<ns; i++){
	    bf.append(this.theoSpectrum.get(indices.get(i).intValue()).getWeight() + "\n");
	    answerbf.append(this.theoSpectrum.get(indices.get(i).intValue()).toString() + "\n");
	}

	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.expDir + filePrefix + "_" + ns + "_" + this.theoSpectrum.size() + ".txt"));
	    bw.write(bf.toString());
	    bw.close();
	    bw = new BufferedWriter(new FileWriter(this.expDir + filePrefix + "_" + ns + "_" + this.theoSpectrum.size() + ".wPeptide.txt"));
	    bw.write(answerbf.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}

    }


    public static void main(String[] args){
	
	//generate characterSet + cyclopeptide and use float weights and then simulate
	if(args.length == 9){
	    new Simulator(Integer.parseInt(args[0]) , Integer.parseInt(args[1]) , Double.parseDouble(args[2]) , Integer.parseInt(args[3]) , args[4], Double.parseDouble(args[5]) , Double.parseDouble(args[6]) , args[7] , args[8], "N");
	}
	//generate characterSet + cyclopeptide and use integer weights and then simulate
	else if(args.length == 10){
	    new Simulator(Integer.parseInt(args[0]) , Integer.parseInt(args[1]) , Double.parseDouble(args[2]) , Integer.parseInt(args[3]) , args[4], Double.parseDouble(args[5]) , Double.parseDouble(args[6]) , args[7] , args[8], args[9]);
	}
	
	//Generate cyclopeptide sequence using the input characterset and then simulate
	else if(args.length == 6){
	    new Simulator(Integer.parseInt(args[0]) , Integer.parseInt(args[1]) , Double.parseDouble(args[2]) , args[3], args[4], args[5], null);
	}
	//Using the input characterset as well as cyclopeptide, simulate
	else if(args.length == 7){
	    new Simulator(Integer.parseInt(args[0]) , Integer.parseInt(args[1]) , Double.parseDouble(args[2]) , args[3], args[4], args[5], args[6]);
	}else{
	    System.err.println("USAGE:\n" 
			       + "To generate characterset and simulate data:\n"
			       + "        java Simulator <l> <r> <f> <n> <sf> <minW> <maxW> <cf> <oDir> [YyNn: INT?]\n" 
			       + "\n" 
			       + "To use input chracterset and simulate data based on the input:\n"
			       + "        java Simulator <l> <r> <f> <sf> <chrf> <oDir> [cpep]\n");
	    System.err.println("Parameters:");
	    System.err.println("\tl        [INT]  \tLength of cyclopeptide");	
	    System.err.println("\tr        [INT]  \tNumber of simulations");
	    System.err.println("\tf        [FLOAT]\tSampling ratio. How much of theoretical spectrum you "
			       +"\n\t                \twant to sample. (0-1.0]");
	    System.err.println("\tsf       [STR]  \tFilename for theoretical spectrum. If <cpepf> is given," 
			       +"\n\t                \tthis is an input file containing theoretical spectrum."
			       +"\n\t                \tIf <cpepf> is given and <sf> is NOT found, it will overwrite.");
	    System.err.println("\tchrf     [STR]  \tCharacterset file");
	    System.err.println("\toDir     [STR]  \tOutput directory (ex: data/int_simul1)");
	    System.err.println("\tcpepf    [STR]  \tFile containing a single line : '-' delimited cyclopeptide"
			       +"\n\t                \tof length l (ex: AA15-AA3-AA5-AA8-AA18-AA13");
	    System.err.println("\tn        [INT]  \tSize of character set");
	    System.err.println("\tminW     [FLOAT]\tMinimum weight allowed");
	    System.err.println("\tmaxW     [FLOAT]\tMaximum weight allowed");
	    System.err.println("\tcf       [STR]  \tFilename to output a generated set of characters and their" 
			       +"\n\t                \tassociated weights");
	    System.err.println("\tInt?     [YyNn] \t Force weights to be integers. [default:N]");
	}
		
    }
    
    
    
    /* class fields */
    public Hashtable<String, BigDecimal> character2weight; // KEY: character(amino acid), VALUE: its associated weight
    public ArrayList<Peptide> theoSpectrum; //theoretical spectrum. it includes 0 as well as all l full-length peptide
    public Peptide cyclopeptide; //holds the sequence of cyclopeptide
    public int clen;
    public BigDecimal[] prefixW; //prefixWeights
    public String expDir;
    public static String delim = "-";
    
}
