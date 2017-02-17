import java.math.*;

public class Peptide implements Comparable<Peptide>{

    public Peptide(String[] aas, BigDecimal w){
	this.aas = aas;
	this.weight = w;
    }
    
    /*
    public Peptide(String[] prev, String nextChar, double w){
	this.seqs = new String[prev.length+1];
	for(int i=0; i<prev.length; i++)
	    this.seqs[i] = prev[i];
	this.seqs[prev.length] = nextChar;
	this.weight = w;
	}*/
    
    public String toString(){
	if(aas.length < 1)
	    return null;
	StringBuffer bf = new StringBuffer(weight + "\t" + aas[0]);
	for(int i=1; i<aas.length; i++)
	    bf.append(Simulator.delim + aas[i]);
	
	return bf.toString();
    }

    public int compareTo(Peptide p){
	BigDecimal diff = this.weight.subtract(p.getWeight());
	if(diff.doubleValue() > 0)
	    return 1;
	else if(diff.doubleValue() < 0)
	    return -1;
	return 0;
    }

    public String[] getAAs(){
	return aas;
    }

    public BigDecimal getWeight(){
	return weight;
    }

    private String[] aas;
    private BigDecimal weight;
    
}
