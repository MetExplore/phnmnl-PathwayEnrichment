package fr.inra.toulouse.metexplore;

public class MappingElement implements Comparable <MappingElement>{

    protected Boolean isMapped;
    protected String inFileName, sbmlName, ID, inFileVal, sbmlVal;

    public MappingElement(Boolean isMapped, String inFileName, String sbmlName, String ID,
                          String inFileVal, String sbmlVal) {
        this.isMapped = isMapped;
        this.inFileName = inFileName;
        this.sbmlName = sbmlName;
        this.ID = ID;
        this.inFileVal = inFileVal;
        this.sbmlVal = sbmlVal;
    }
    public int compareTo(MappingElement m) {
        if (this.isMapped == m.isMapped) return this.inFileName.compareToIgnoreCase(m.inFileName);
        else if (this.isMapped && !m.isMapped) return -1;
        return 1;
    }

    public String toString(){
        return this.isMapped + "\t" + this.inFileName + "\t" + this.sbmlName + "\t" + this.ID + "\t" + this.inFileVal + "\t" + this.sbmlVal + "\n";
    }
}