package fr.inra.toulouse.metexplore.omicsComponents;

import fr.inra.toulouse.metexplore.io.WritingComportment;

import java.util.List;

public class EnrichedEntity implements Comparable <EnrichedEntity>, WritingComportment {

    protected String pathName, mappedMetabolitesSBML, mappedMetabolitesFingerprint, mappedMetabolitesID, coverage,
            p_value, q_value_Bonf, q_value_BenHoc;
    protected int nb_mapped, nb_unmappedInPathway, nb_unmappedInFingerprint, nb_remainingInNetwork;

    public int getNb_mapped() {
        return nb_mapped;
    }

    public void setNb_unmappedInPathway(int nb_unmappedInPathway) {
        this.nb_unmappedInPathway = nb_unmappedInPathway;
    }

    public void setNb_unmappedInFingerprint(int nb_unmappedInFingerprint) {
        this.nb_unmappedInFingerprint = nb_unmappedInFingerprint;
    }

    public void setNb_remainingInNetwork(int nb_remainingInNetwork) {
        this.nb_remainingInNetwork = nb_remainingInNetwork;
    }

    public EnrichedEntity(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                          List<String> mappedMetabolitesSBML, List<String> mappedMetabolitesFingerprint,
                          List<String> mappedMetabolitesID, int nb_mapped, String coverage) {

        this.pathName = pathName;
        this.p_value = removeSciNot(p_value);
        this.q_value_Bonf = removeSciNot(q_value_Bonf);
        this.q_value_BenHoc = removeSciNot(q_value_BenHoc);
        this.mappedMetabolitesFingerprint = String.join(";", mappedMetabolitesFingerprint);
        this.mappedMetabolitesSBML = String.join(";", mappedMetabolitesSBML);
        this.mappedMetabolitesID = String.join(";", mappedMetabolitesID);
        this.nb_mapped = nb_mapped;
        this.coverage = coverage;
    }

    //TODO: hashmap pour list metab et id et sort functionality

    public int compareTo(EnrichedEntity p){
        return (this.pathName).compareToIgnoreCase(p.pathName);
    }

    public String toString(Boolean galaxyCompliance) {
        String line = this.pathName + "\t" + this.coverage + "\t" + this.nb_mapped + "\t" + this.p_value +
                "\t" + this.q_value_Bonf + "\t" + this.q_value_BenHoc + "\t" + this.mappedMetabolitesSBML
                + "\t" + this.mappedMetabolitesFingerprint + "\t" + this.mappedMetabolitesID;
        if (galaxyCompliance)
            line += "\t" + this.nb_unmappedInPathway + "\t" + this.nb_unmappedInFingerprint + "\t" + this.nb_remainingInNetwork;
        return line + "\n";
    }
}