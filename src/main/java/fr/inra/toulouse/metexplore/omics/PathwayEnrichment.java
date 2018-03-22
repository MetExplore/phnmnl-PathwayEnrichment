package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.omicsComponents.EnrichedEntity;
import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;

import static java.util.Collections.sort;

public class PathwayEnrichment extends Omics{

    protected List<HashMap<BioEntity, Double>> list_pathwayEnr = new ArrayList<HashMap<BioEntity, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    protected String outFilePathEnr;
    protected int entityType2Enrich;
    protected String typeOfEnrichedEntity;
    protected  Set<BioChemicalReaction> reactionSet;

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

        HashMap<BioEntity, Double> pathEnrWhithPval = this.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioEntity, Double> pathEnrBenHoc = this.benjaminiHochbergCorrection(pathEnrWhithPval);

        this.list_pathwayEnr.add(sortPathByBenHoc(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval,
        // need to do the same here to join with it
        this.list_pathwayEnr.add(sortPathByBenHoc(this.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        this.list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg

        int nbEnriched = this.list_pathwayEnr.get(0).size();
        int nbEnrichedTypeInNetwork = getEntitySetInNetwork(this.entityType2Enrich).size();
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

        f.write(this.typeOfEnrichedEntity + " name\tCoverage (%)\tNb. of mapped\tP-value\tBonferroni corrected p-value" +
                "\tBH corrected p-value\tMapped " + typeOfMappedEntity + "s (SBML)\tMapped " + typeOfMappedEntity + "s (fingerprint)\t" +
                "Mapped " + typeOfMappedEntity + "s ID");
        if (!this.galaxy.equals(""))  f.write("\tNb. of unmapped (pathway)\tNb. of unmapped (fingerprint)\tNb. of remaining (network)\n");
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
                if (this.getEntityListInEnrichedEntity(path).contains(bpe)) {
                    hm_pathwayMetabolitesFingerprint.put(bpe.getId(), entry.getValue());
                    j++;
                }
            }
            list_pathwayMetabolitesID = new ArrayList<>(hm_pathwayMetabolitesFingerprint.keySet());
            Collections.sort(list_pathwayMetabolitesID);
            for(String id_bpe : list_pathwayMetabolitesID){
                list_pathwayMetabolitesSBML.add(((BioEntity)getEntitySetInNetwork().get(id_bpe)).getName());
                list_pathwayMetabolitesFingerprint.add(hm_pathwayMetabolitesFingerprint.get(id_bpe));
            }
            String coverage = round((double) j / (double) this.getEntityListInEnrichedEntity(path).size() * (double) 100);
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
        int fisherTestParameters[] = this.getFisherTestParameters(pathway);
        pathEl.setNb_unmappedInFingerprint(fisherTestParameters[1]);
        pathEl.setNb_unmappedInPathway(fisherTestParameters[2]);
        pathEl.setNb_remainingInNetwork(fisherTestParameters[3]);
    }

    public void setReactionSet(Set<? extends BioEntity> BioEntitySet) {
        this.reactionSet = new HashSet<>();
        for(BioEntity e : BioEntitySet) {
            if (e instanceof BioChemicalReaction) {
                this.reactionSet.add((BioChemicalReaction) e);
            /*}else if(this.methods.bioEntityTYpe == 4) {
                for (BioChemicalReaction r : network.getBiochemicalReactionList().values()) {
                    if (r.getEnzList().values().contains(e)) {
                        System.out.println("\t"  + ". " + r.getId());
                        for (BioPathway p : network.getPathwayList().values()) {
                            if(p.getReactions().values().contains(r)){
                                System.out.println("\t\t" + p.getId());
                            }
                        }
                        reactionSet.add(r);
                    }
                }*/
                //BUG: TODO: see in JSBML2BioNetwork why only one enzyme is associated to a reaction
                // (instead of multiple for protein and genes)
            }else if(this.entityType2Map == 4 || e instanceof BioProtein) {
                BioProtein p = network.getProteinList().get(e.getId());
                HashMap<String, BioGene> list_genes = p.getGeneList();
                for (BioGene g : list_genes.values()){
                    addReactions(network.getReactionsFromGene(g.getId()));
                }
            } else if (e instanceof BioGene) {
                BioGene g = (BioGene) e;
                addReactions(network.getReactionsFromGene(g.getId()));
            } else if (e instanceof BioPhysicalEntity) {
                if (this.entityType2Map == 1) {
                    BioPhysicalEntity m = (BioPhysicalEntity) e;
                    this.reactionSet.addAll(m.getReactionsAsProduct().values());
                    this.reactionSet.addAll(m.getReactionsAsSubstrate().values());
                }
            }
        }
    }

    public void addReactions(Set <String> list_reactions_ID){

        for (String reac_ID : list_reactions_ID) {
            BioChemicalReaction r = network.getBiochemicalReactionList().get(reac_ID);
            this.reactionSet.add(r);
        }
    }

    public HashMap<BioEntity, Double> computeEnrichment() {
        HashSet<BioEntity> entityType2Enrich = new HashSet<>();
        this.parseProtList4Genes();
        setReactionSet(this.list_mappedEntities.keySet());
        for ( BioChemicalReaction r : this.reactionSet){
            //System.out.println("Reac: " + r.getName() + ": " + r.getListOfGeneNames().size());

            switch (this.entityType2Enrich){
                case 2:
                    entityType2Enrich.add(r);
                    break;
                case 3:
                    entityType2Enrich.addAll(r.getPathwayList().values());
                    break;
                case 4: case 5:
                    for (BioGene g : r.getListOfGenes().values()){
                        //Rare case of some proteins which do not have any mapped genes
                        // (sometimes it works in gene to protein list but not in protein to gene list
                        // for getEntityListInEnrichedEntity step : problem in JSBML2Bionetwork parsing?
                        // because of alternative splicing not taken in account?)
                        for(BioProtein p : g.proteinList.values()){
                            entityType2Enrich.add(p);
                        }
                    }
                    break;
                case 6:
                    entityType2Enrich.addAll(r.getListOfGenes().values());
                    break;
            }
        }

        HashMap<BioEntity, Double> res = new HashMap<>();
        for (BioEntity p : entityType2Enrich){
            res.put(p, this.getPvalue(p));
        }
        return res;
    }

    public Collection getEntityListInEnrichedEntity(BioEntity entityType2Enrich) {
        switch (this.entityType2Enrich){
            case 1:
                return getEntityListInMetabolite(entityType2Enrich);
            case 2:
                return getEntityListInReaction(entityType2Enrich);
            case 3:
                return getEntityListInPathway(entityType2Enrich);
            case 4: case 5:
                return getEntityListInProtein(entityType2Enrich);
            case 6:
                return getEntityListInGene(entityType2Enrich);
        }
        return null;
    }

    public Collection getEntityListInReaction(BioEntity enrichedEntity) {
        BioChemicalReaction reaction = (BioChemicalReaction) enrichedEntity;

        switch (entityType2Map) {
            case 1:
                Set <BioPhysicalEntity> mets = new HashSet<>(reaction.getListOfSubstrates().values());
                mets.addAll(reaction.getListOfProducts().values());
                return mets;
            case 4: case 5:
                Set <BioProtein> prots = new HashSet<>();
                for (BioGene g : reaction.getListOfGenes().values()){
                    prots.addAll(g.getProteinList().values());
                }
                return prots;
            case 6:
                return reaction.getListOfGenes().values();
        }
        return null;
    }

    public Collection getEntityListInPathway(BioEntity enrichedEntity) {
        BioPathway pathway = (BioPathway) enrichedEntity;
        switch (entityType2Map){
            case 1:
                return pathway.getListOfInvolvedMetabolite().values();
            case 2:
                return pathway.getReactions().values();
          /* case 4:
               Set <BioEntity> enzymes = new HashSet<BioEntity>();
               //System.out.println("ReacSize2: " + pathway.getReactions().size());
               for (BioChemicalReaction r : pathway.getReactions().values()){
                   //System.out.println("EnzSize2: " + r.getEnzList().size());
                   enzymes.addAll(r.getEnzList().values());
               }
               return enzymes;*/
            //BUG: TODO: see in JSBML2BioNetwork why only one enzyme is associated to a reaction
            // (instead of multiple for protein and genes)
            case 4 : case 5:
                Set <BioProtein> proteins = new HashSet<>();
                //System.out.println("GenesSize2: " + pathway.getGenes().size());
                for (BioGene g : pathway.getGenes()){
                    //System.out.println("ProtSize2: " + g.getProteinList().size());
                    proteins.addAll(g.getProteinList().values());
                }
                return proteins;
            case 6:
                return pathway.getGenes();
        }
        return null;
    }

    public Collection getEntityListInGene(BioEntity enrichedEntity) {
        BioGene g = (BioGene) enrichedEntity;
        Set <String> list_reactions_ID = network.getReactionsFromGene(g.getId());
        HashSet list_entities = new HashSet();

        switch (entityType2Map) {
            case 1:
                for (String reac_ID : list_reactions_ID) {
                    BioChemicalReaction reaction = network.getBiochemicalReactionList().get(reac_ID);
                    list_entities.addAll(reaction.getListOfSubstrates().values());
                    list_entities.addAll(reaction.getListOfProducts().values());
                }
            case 2:
                for (String reac_ID : list_reactions_ID) {
                    list_entities.add(network.getBiochemicalReactionList().get(reac_ID));
                }
            case 4: case 5:
                list_entities.addAll(g.getProteinList().values());
        }
        return list_entities;
    }

    public Collection getEntityListInProtein(BioEntity enrichedEntity) {
        BioProtein p = (BioProtein) enrichedEntity;
        HashMap<String, BioGene> list_genes = p.getGeneList();
        Set <BioEntity> list_entities = new HashSet<>();
        for (BioGene g : list_genes.values()) {
            Set <String> list_reactions_ID = network.getReactionsFromGene(g.getId());

            switch (entityType2Map) {
                case 1:
                    for (String reac_ID : list_reactions_ID) {
                        BioChemicalReaction reaction = network.getBiochemicalReactionList().get(reac_ID);
                        list_entities.addAll(reaction.getListOfSubstrates().values());
                        list_entities.addAll(reaction.getListOfProducts().values());
                    }
                case 2:
                    for (String reac_ID : list_reactions_ID) {
                        list_entities.add(network.getBiochemicalReactionList().get(reac_ID));
                    }
                case 6:
                    list_entities.add(g);
            }
        }
        return list_entities;
    }

    public Collection getEntityListInMetabolite(BioEntity enrichedEntity) {
        BioPhysicalEntity m = (BioPhysicalEntity) enrichedEntity;
        Set <BioChemicalReaction> reactionContainingMetabolites = new HashSet<>();
        reactionContainingMetabolites.addAll(m.getReactionsAsProduct().values());
        reactionContainingMetabolites.addAll(m.getReactionsAsSubstrate().values());
        Set <BioEntity> list_entities = new HashSet<>();
        for (BioChemicalReaction r : reactionContainingMetabolites) {

            switch (entityType2Map) {
                case 2:
                    list_entities.add(r);
                case 4: case 5:
                    for (BioGene g : r.getListOfGenes().values()) {
                        list_entities.addAll(g.proteinList.values());
                    }
                case 6:
                    list_entities.addAll(r.getListOfGenes().values());
            }
        }
        return list_entities;
    }

    public double getPvalue(BioEntity enrichedEntity) throws IllegalArgumentException {
        int fisherTestParameters[] = this.getFisherTestParameters(enrichedEntity);
        return exactFisherOneTailed(fisherTestParameters[0], fisherTestParameters[1],
                fisherTestParameters[2], fisherTestParameters[3]);
    }

    public int[] getFisherTestParameters(BioEntity enrichedEntity) {
        Collection entityListInEnrichedEntity = this.getEntityListInEnrichedEntity(enrichedEntity);
        //nb of mapped in the pathway
        int a = intersect(entityListInEnrichedEntity).size();
        //unmapped metabolites in the fingerprint
        int b = this.list_mappedEntities.size() - a;
        //unmapped metabolites in the pathway
        int c = entityListInEnrichedEntity.size() - a;
        //remaining metabolites in the network
        int d = getEntitySetInNetwork().size() - (a + b + c);

        //System.out.println(pathway.getName() + ": " + Arrays.toString(fisherTestParameters));
        return new int[]{a, b, c, d};
    }

    public HashMap<BioEntity, Double> benjaminiHochbergCorrection(HashMap<BioEntity, Double> pvalues) {
        ArrayList<BioEntity> orderedEnrichedEntities = this.sortPval(pvalues);
        HashMap<BioEntity, Double> adjPvalues = new HashMap<>();
        BioEntity e, e_next;
        double pval, adjPval;
        Boolean wasAlreadyEqual = false;
        HashSet <BioEntity> entityWithEqualsPval = new HashSet<>();
        int k_adj = 0;

        for(int k = 0; k < orderedEnrichedEntities.size(); ++k) {
            e = orderedEnrichedEntities.get(k);
            entityWithEqualsPval.add(e);
            pval = pvalues.get(e);
            //System.out.println("Pval : " + e.getName() + ": " + pval);

            //case of probability equality
            if (k + 1 < orderedEnrichedEntities.size()) {
                e_next = orderedEnrichedEntities.get(k + 1);
                double pval_next = pvalues.get(e_next);
                if (pval_next == pval) {
                    entityWithEqualsPval.add(e_next);
                    if (!wasAlreadyEqual) {
                        k_adj = k;
                    }
                    wasAlreadyEqual = true;
                    continue;
                } else if (wasAlreadyEqual) {
                    adjPval = pval * (double) pvalues.size() / (double) (k_adj + 1);
                    for (BioEntity e2 : entityWithEqualsPval) {
                        adjPvalues.put(e2, adjPval);
                    }
                    entityWithEqualsPval = new HashSet<>();
                    wasAlreadyEqual = false;
                    continue;
                }

            }
            //default case
            adjPval = pval * (double) pvalues.size() / (double) (k + 1);
            adjPvalues.put(e, adjPval);
        }
        return adjPvalues;
    }

    public static BigDecimal fact(int n) {
        BigDecimal fact = new BigDecimal("1");

        for(int i = 1; i <= n; ++i) {
            fact = fact.multiply(new BigDecimal(i + ""));
        }

        return fact;
    }

    public HashMap<BioEntity, Double> bonferroniCorrection(HashMap<BioEntity, Double> pvalues) {
        HashMap<BioEntity, Double> adjPvalues = new HashMap<>();
        for (BioEntity e : pvalues.keySet()){
            double pval = pvalues.get(e);
            double adjPval = pval * (double)pvalues.size();
            adjPvalues.put(e, adjPval);
        }
        return adjPvalues;
    }

    public ArrayList<BioEntity> sortPval(HashMap<BioEntity, Double> map) {
        ArrayList<BioEntity> orderedEnrichedEntities = new ArrayList<>(map.keySet());
        sort(orderedEnrichedEntities, new significanceComparator(map));
        return orderedEnrichedEntities;
    }

    static class significanceComparator implements Comparator<BioEntity> {
        HashMap<BioEntity, Double> pvalMap;

        public significanceComparator(HashMap<BioEntity, Double> pvalMap) {
            this.pvalMap = pvalMap;
        }

        public int compare(BioEntity o1, BioEntity o2) {
            return Double.compare(this.pvalMap.get(o1), this.pvalMap.get(o2));
        }
    }

    public double exactFisherOneTailed(int a, int b, int c, int d) {
        double res = 0.0D;
        int lim = Math.min(a + c, a + b);

        for(int i = 0; a + i <= lim; ++i) {
            res += getHypergeometricProba(a + i, b - i, c - i, d + i);
        }

        return res;
    }

    public static double getHypergeometricProba(int a, int b, int c, int d) {
        BigDecimal numerator = fact(a + b).multiply(fact(c + d)).multiply(fact(a + c)).multiply(fact(b + d));
        BigDecimal denominator = fact(a).multiply(fact(b)).multiply(fact(c)).multiply(fact(d)).multiply(fact(a + b + c + d));
        BigDecimal res = numerator.divide(denominator, MathContext.DECIMAL64);
        return res.doubleValue();
    }

    public void parseProtList4Genes(){
        HashMap<String, BioProtein> protList;
        for (BioGene g : this.network.getGeneList().values()) {
            protList = new HashMap<>();
            for (BioProtein p : this.network.getProteinList().values()) {
                if(p.getGeneList().values().contains(g)) protList.put(p.getId(),p);
            }
            g.setProteinList(protList);
        }
    }
}