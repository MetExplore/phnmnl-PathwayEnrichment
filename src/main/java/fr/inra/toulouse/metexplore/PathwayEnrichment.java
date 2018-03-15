package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PathwayEnrichment extends Omics{

    protected List<HashMap<BioPathway, Double>> list_pathwayEnr = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    protected String outFilePathEnr;

    public PathwayEnrichment (BioNetwork network, ArrayList<String[]> list_fingerprint, HashMap <BioEntity, String> list_mappedEntities,
                              String outFilePathEnr, String galaxy , int bioEntityType) throws IOException {
        super(galaxy, list_fingerprint, list_mappedEntities, network, bioEntityType);
        this.outFilePathEnr=outFilePathEnr;
        this.computeEnrichmentWithCorrections();
    }

    public void computeEnrichmentWithCorrections() throws IOException {
        System.out.println("Pathway enrichment in progress...");

        if (bioEntityType == 5) this.parseProtList4Genes();
        System.out.println("MapSize: " + list_mappedEntities.size());
        PathwayEnrichmentCalculation pathEnr = new PathwayEnrichmentCalculation(this.network, this.list_mappedEntities,this.bioEntityType);
        HashMap<BioPathway, Double> pathEnrWhithPval = pathEnr.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioPathway, Double> pathEnrBenHoc = pathEnr.benjaminiHochbergCorrection(pathEnrWhithPval);

        this.list_pathwayEnr.add(sortPathByBenHoc(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval,
        // need to do the same here to join with it
        this.list_pathwayEnr.add(sortPathByBenHoc(pathEnr.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        this.list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg
        write.writeLog(this.list_pathwayEnr.get(0).size() + " pathways are concerned among the network (on " + this.network.getPathwayList().size() + " in the network).");
        writeOutputPathEnr();
    }

    public HashMap<BioPathway, Double> sortPathByBenHoc(HashMap<BioPathway, Double> disorderedPathEnr, HashMap<BioPathway, Double> pathEnrBenHoc) {
        ArrayList<BioPathway> pathBenHoc = new ArrayList(pathEnrBenHoc.keySet());
        HashMap<BioPathway, Double> orderedPathEnr = new HashMap();
        for (int i = 0; i < pathBenHoc.size(); ++i) {
            BioPathway path = (BioPathway) pathBenHoc.get(i);
            double pval = (Double) disorderedPathEnr.get(path);
            orderedPathEnr.put(path, pval);

        }
        return orderedPathEnr;
    }

    public void writeOutputPathEnr() throws IOException{

        File fo2 = new File(this.outFilePathEnr);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <PathwayEnrichmentElement> list_pathwayEnrElement = new ArrayList<PathwayEnrichmentElement>(); //list of pathway enrichment instantiation for sorting
        HashMap <String, String> hm_pathwayMetabolitesFingerprint;
        List <String> list_pathwayMetabolitesID;
        List <String> list_pathwayMetabolitesSBML;
        List <String> list_pathwayMetabolitesFingerprint;

        f.write("Pathway_name\tp-value\tBonferroni_corrected_p_value\tBH_corrected_p_value\tMapped_entities_SBML\tMapped_entities_Fingerprint\t" +
                "Mapped_entities_ID\tNb. of mapped\tCoverage (%)");
        if (this.write.galaxy !="")  f.write("\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network\n");
        else f.write("\n");

        HashMap<BioPathway, Double> result = this.list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = this.list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = this.list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            hm_pathwayMetabolitesFingerprint = new HashMap<String, String>();
            list_pathwayMetabolitesSBML = new ArrayList<String>();
            list_pathwayMetabolitesFingerprint = new ArrayList<String>();
            BioPathway path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (Map.Entry<BioEntity, String> entry : this.list_mappedEntities.entrySet()) {
                BioEntity bpe = entry.getKey();
                if (omics.getEntitySetInPathway(path).contains(bpe)) {
                    hm_pathwayMetabolitesFingerprint.put(bpe.getId(), entry.getValue());
                    j++;
                }
            }
            list_pathwayMetabolitesID = new ArrayList<String>(hm_pathwayMetabolitesFingerprint.keySet());
            Collections.sort(list_pathwayMetabolitesID);
            for(String id_bpe : list_pathwayMetabolitesID){
                list_pathwayMetabolitesSBML.add(((BioEntity)omics.getEntitySetInNetwork().get(id_bpe)).getName());
                list_pathwayMetabolitesFingerprint.add(hm_pathwayMetabolitesFingerprint.get(id_bpe));
            }
            String coverage = this.write.round((double) j / (double) omics.getEntitySetInPathway(path).size() * (double) 100);
            PathwayEnrichmentElement pathEnrElement = new PathwayEnrichmentElement(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),
                    (double)itBonCorr.next(),(double)itBenHocCorr.next(),list_pathwayMetabolitesSBML,list_pathwayMetabolitesFingerprint, list_pathwayMetabolitesID,j,coverage, this.write.galaxy);
            if (this.write.galaxy != "") this.settings4Galaxy(path, pathEnrElement);
            list_pathwayEnrElement.add(pathEnrElement);
        }

        Collections.sort(list_pathwayEnrElement);
        for (int i=0;i< list_pathwayEnrElement.size();i++){
            f.write(list_pathwayEnrElement.get(i).toString());
        }
        if (f != null) {
            f.close();
        }
    }

    public void settings4Galaxy(BioPathway pathway, PathwayEnrichmentElement pathEl) {
        int fisherTestParameters[] = omics.getFisherTestParameters(pathway);
        pathEl.nb_unmappedInFingerprint = fisherTestParameters[1];
        pathEl.nb_unmappedInPathway = fisherTestParameters[2];
        pathEl.nb_remainingInNetwork = fisherTestParameters[3];
    }
}