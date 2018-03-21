package fr.inra.toulouse.metexplore;

import java.util.List;

public class PathwayEnrichmentElement implements Comparable <PathwayEnrichmentElement> {

    protected String pathName, mappedMetabolitesSBML, mappedMetabolitesFingerprint, mappedMetabolitesID, coverage,
            p_value, q_value_Bonf, q_value_BenHoc;
    protected int nb_mapped;
    protected int nb_unmappedInPathway;
    protected int nb_unmappedInFingerprint;

    public void setNb_unmappedInPathway(int nb_unmappedInPathway) {
        this.nb_unmappedInPathway = nb_unmappedInPathway;
    }

    public void setNb_unmappedInFingerprint(int nb_unmappedInFingerprint) {
        this.nb_unmappedInFingerprint = nb_unmappedInFingerprint;
    }

    public void setNb_remainingInNetwork(int nb_remainingInNetwork) {
        this.nb_remainingInNetwork = nb_remainingInNetwork;
    }

    protected int nb_remainingInNetwork;
    protected WritingComportment write;

    public PathwayEnrichmentElement(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                                    List<String> mappedMetabolitesSBML, List<String> mappedMetabolitesFingerprint,
                                    List<String> mappedMetabolitesID, int nb_mapped, String coverage) {

        this.write = new WritingComportment("");
        this.pathName = pathName;
        this.p_value = write.removeSciNot(p_value);
        this.q_value_Bonf = write.removeSciNot(q_value_Bonf);
        this.q_value_BenHoc = write.removeSciNot(q_value_BenHoc);
        this.mappedMetabolitesFingerprint = String.join(";", mappedMetabolitesFingerprint);
        this.mappedMetabolitesSBML = String.join(";", mappedMetabolitesSBML);
        this.mappedMetabolitesID = String.join(";", mappedMetabolitesID);
        this.nb_mapped = nb_mapped;
        this.coverage = coverage;
    }

    //TODO: hashmap pour list metab et id et sort functionality

    public int compareTo(PathwayEnrichmentElement p){
        return (this.pathName).compareToIgnoreCase(p.pathName);
    }

    public String toString() {
        String line = this.pathName + "\t" + this.p_value + "\t" + this.q_value_Bonf + "\t" + this.q_value_BenHoc
                + "\t" + this.mappedMetabolitesSBML + "\t" + this.mappedMetabolitesFingerprint + "\t" + this.mappedMetabolitesID + "\t" + this.nb_mapped + "\t" + this.coverage;
        if (!write.galaxy.equals(""))
            line += "\t" + this.nb_unmappedInPathway + "\t" + this.nb_unmappedInFingerprint + "\t" + this.nb_remainingInNetwork;
        return line + "\n";
    }
}
