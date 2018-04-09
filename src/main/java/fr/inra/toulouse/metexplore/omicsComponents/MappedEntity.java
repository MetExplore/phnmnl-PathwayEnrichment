package fr.inra.toulouse.metexplore.omicsComponents;

public class MappedEntity implements Comparable <MappedEntity>{

    protected Boolean isMapped;
    protected String inFileName, sbmlName, ID, inFileVal, sbmlVal;

    public MappedEntity(Boolean isMapped, String inFileName, String sbmlName, String ID,
                        String inFileVal, String sbmlVal) {
        this.isMapped = isMapped;
        this.inFileName = inFileName;
        this.sbmlName = sbmlName;
        this.ID = ID;
        this.inFileVal = inFileVal;
        this.sbmlVal = sbmlVal;
    }
    public int compareTo(MappedEntity m) {
        if (this.isMapped == m.isMapped) return this.inFileName.compareToIgnoreCase(m.inFileName);
        else if (this.isMapped && !m.isMapped) return -1;
        return 1;
    }

    public String toString(){
        return this.isMapped + "\t" + this.inFileName + "\t" + this.sbmlName + "\t" + this.ID + "\t" + this.inFileVal + "\t" + this.sbmlVal + "\n";
    }
}
