package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioPathway;

import java.util.List;

public class PathwayEnrichmentElement implements Comparable <PathwayEnrichmentElement> {

    protected String pathName, mappedMetabolites, mappedMetabolitesID, coverage,
            p_value, q_value_Bonf, q_value_BenHoc;
    protected int nb_mapped, nb_unmappedInPathway, nb_unmappedInFingerprint, nb_remainingInNetwork;
    protected Boolean ifGalaxy;

    public PathwayEnrichmentElement(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                                    List<String> mappedMetabolites, List<String> mappedMetabolitesID,
                                    int nb_mapped, String coverage, Boolean ifGalaxy) {
        this.pathName = pathName;
        this.p_value = (new WritingComportment()).removeSciNot(p_value);
        this.q_value_Bonf = (new WritingComportment()).removeSciNot(q_value_Bonf);
        this.q_value_BenHoc = (new WritingComportment()).removeSciNot(q_value_BenHoc);
        this.mappedMetabolites = String.join(";", mappedMetabolites);
        this.mappedMetabolitesID = String.join(";", mappedMetabolitesID);
        this.nb_mapped = nb_mapped;
        this.coverage = coverage;
        this.ifGalaxy = ifGalaxy;
    }

    //TODO: hashmap pour list metab et id et sort functionality

    public int compareTo(PathwayEnrichmentElement p){
        return (this.pathName).compareToIgnoreCase(p.pathName);
    }

    public String toString() {
        String line = this.pathName + "\t" + this.p_value + "\t" + this.q_value_Bonf + "\t" + this.q_value_BenHoc
                + "\t" + this.mappedMetabolites.toString() + "\t" + this.mappedMetabolitesID.toString() + "\t" + this.nb_mapped + "\t" + this.coverage;
        if (this.ifGalaxy)
            line += "\t" + this.nb_unmappedInPathway + "\t" + this.nb_unmappedInFingerprint + "\t" + this.nb_remainingInNetwork;
        return line + "\n";
    }
}
