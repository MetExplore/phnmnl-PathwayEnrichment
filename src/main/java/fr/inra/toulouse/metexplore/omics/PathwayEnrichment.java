package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.PathwayEnrichmentCalculation;
import fr.inra.toulouse.metexplore.EnrichedEntity;
import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PathwayEnrichment extends Omics{

    protected List<HashMap<BioEntity, Double>> list_pathwayEnr = new ArrayList<HashMap<BioEntity, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    protected String outFilePathEnr;
    protected PathwayEnrichmentCalculation pathEnr;
    protected int entityType2Enrich;
    protected String typeOfEnrichedEntity;

    public PathwayEnrichment (String logContent, BioNetwork network, ArrayList<String[]> list_fingerprint, HashMap <BioEntity, String> list_mappedEntities,
                              String outFilePathEnr, String galaxy , int entityType2Map, int entityType2Enrich) throws IOException {
        super(logContent, galaxy, list_fingerprint, list_mappedEntities, network, entityType2Map);
        this.outFilePathEnr = outFilePathEnr;
        this.entityType2Enrich = entityType2Enrich;
        this.typeOfEnrichedEntity = getTypeOfEntity(entityType2Enrich);
        this.computeEnrichmentWithCorrections();
    }

    public void computeEnrichmentWithCorrections() throws IOException {
        System.out.println(this.typeOfEnrichedEntity + " enrichment in progress...");

        this.pathEnr = new PathwayEnrichmentCalculation(this.network, this.list_mappedEntities,this.bioEntityType, this.entityType2Enrich);
        HashMap<BioEntity, Double> pathEnrWhithPval = this.pathEnr.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioEntity, Double> pathEnrBenHoc = this.pathEnr.benjaminiHochbergCorrection(pathEnrWhithPval);

        this.list_pathwayEnr.add(sortPathByBenHoc(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval,
        // need to do the same here to join with it
        this.list_pathwayEnr.add(sortPathByBenHoc(this.pathEnr.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        this.list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg

        int nbEnriched = this.list_pathwayEnr.get(0).size();
        int nbEnrichedTypeInNetwork = getEntitySetInNetwork(network,this.entityType2Enrich).size();
        String plural = (nbEnriched > 1) ? "s are": " is";

        this.logContent = writeLog(logContent,nbEnriched + " " + typeOfEnrichedEntity.toLowerCase() + plural + " concerned among the network (on " + nbEnrichedTypeInNetwork + " in the network; " + calculPercent(nbEnriched,nbEnrichedTypeInNetwork) + "%).");
        writeOutputPathEnr();
    }

    public HashMap<BioEntity, Double> sortPathByBenHoc(HashMap<BioEntity, Double> disorderedPathEnr, HashMap<BioEntity, Double> pathEnrBenHoc) {
        ArrayList<BioEntity> pathBenHoc = new ArrayList(pathEnrBenHoc.keySet());
        HashMap<BioEntity, Double> orderedPathEnr = new HashMap();
        for (BioEntity path : pathBenHoc) {
            double pval = disorderedPathEnr.get(path);
            orderedPathEnr.put(path, pval);
        }
        return orderedPathEnr;
    }

    public void writeOutputPathEnr() throws IOException{

        File fo2 = new File(this.outFilePathEnr);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <EnrichedEntity> list_enrichedEntity = new ArrayList<EnrichedEntity>(); //list of pathway enrichment instantiation for sorting
        HashMap <String, String> hm_pathwayMetabolitesFingerprint;
        List <String> list_pathwayMetabolitesID;
        List <String> list_pathwayMetabolitesSBML;
        List <String> list_pathwayMetabolitesFingerprint;

        f.write(this.typeOfEnrichedEntity + "_name\tp-value\tBonferroni_corrected_p_value\tBH_corrected_p_value\tMapped_" + typeOfMappedEntity + "s_(SBML)\tMapped_" + typeOfMappedEntity + "s_(fingerprint)\t" +
                "Mapped_" + typeOfMappedEntity + "s_ID\tNb. of mapped\tCoverage (%)");
        if (!this.galaxy.equals(""))  f.write("\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network\n");
        else f.write("\n");

        HashMap<BioEntity, Double> result = this.list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = this.list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = this.list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioEntity, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            hm_pathwayMetabolitesFingerprint = new HashMap<>();
            list_pathwayMetabolitesSBML = new ArrayList<>();
            list_pathwayMetabolitesFingerprint = new ArrayList<>();
            BioEntity path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (Map.Entry<BioEntity, String> entry : this.list_mappedEntities.entrySet()) {
                BioEntity bpe = entry.getKey();
                if (this.pathEnr.getEntityListInEnrichedEntity(path).contains(bpe)) {
                    hm_pathwayMetabolitesFingerprint.put(bpe.getId(), entry.getValue());
                    j++;
                }
            }
            list_pathwayMetabolitesID = new ArrayList<>(hm_pathwayMetabolitesFingerprint.keySet());
            Collections.sort(list_pathwayMetabolitesID);
            for(String id_bpe : list_pathwayMetabolitesID){
                list_pathwayMetabolitesSBML.add(((BioEntity)getEntitySetInNetwork(network,this.bioEntityType).get(id_bpe)).getName());
                list_pathwayMetabolitesFingerprint.add(hm_pathwayMetabolitesFingerprint.get(id_bpe));
            }
            String coverage = round((double) j / (double) this.pathEnr.getEntityListInEnrichedEntity(path).size() * (double) 100);
            EnrichedEntity pathEnrElement = new EnrichedEntity(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),
                    (double)itBonCorr.next(),(double)itBenHocCorr.next(),list_pathwayMetabolitesSBML,list_pathwayMetabolitesFingerprint, list_pathwayMetabolitesID,j,coverage);
            if (!this.galaxy.equals("")) this.settings4Galaxy(path, pathEnrElement);
            list_enrichedEntity.add(pathEnrElement);
        }

        Collections.sort(list_enrichedEntity);
        for (EnrichedEntity enrichedEntity : list_enrichedEntity) {
            f.write(enrichedEntity.toString(!galaxy.equals("")));
        }
        if (f != null) {
            f.close();
        }
    }

    public void settings4Galaxy(BioEntity pathway, EnrichedEntity pathEl) {
        int fisherTestParameters[] = this.pathEnr.getFisherTestParameters(pathway);
        pathEl.setNb_unmappedInFingerprint(fisherTestParameters[1]);
        pathEl.setNb_unmappedInPathway(fisherTestParameters[2]);
        pathEl.setNb_remainingInNetwork(fisherTestParameters[3]);
    }
}