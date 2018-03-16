package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;
import parsebionet.statistics.PathwayEnrichment;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;

public class PathwayEnrichmentCalculation {

    protected static final int BONFERRONI = 0;
    protected static final int BENJAMINIHOCHBERG = 1;
    protected static final int HOLMBONFERRONI = 2;
    protected  BioNetwork network;
    protected  HashMap <BioEntity, String> list_mappedEntities;
    protected  OmicsMethods methods;
    protected  Set<BioChemicalReaction> reactionSet;
    protected int entityType2Enrich;
    
    public PathwayEnrichmentCalculation(BioNetwork network, HashMap <BioEntity, String> list_mappedEntities, int bioEntityType, int entityType2Enrich){
        this.network = network;
        this.list_mappedEntities = list_mappedEntities;
        this.methods = new OmicsMethods(list_mappedEntities,network, bioEntityType);
        this.setReactionSet(list_mappedEntities.keySet());
        this.entityType2Enrich = entityType2Enrich;
    }

    public void setReactionSet(Set<? extends BioEntity> BioEntitySet) {
        this.reactionSet = new HashSet();
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
            }else if(this.methods.bioEntityTYpe == 4 || e instanceof BioProtein) {
                BioProtein p = network.getProteinList().get(e.getId());
                HashMap<String, BioGene> list_genes = p.getGeneList();
                for (BioGene g : list_genes.values()){
                    addReactions(network.getReactionsFromGene(g.getId()));
                }
            } else if (e instanceof BioGene) {
                BioGene g = (BioGene) e;
                addReactions(network.getReactionsFromGene(g.getId()));
            } else if (e instanceof BioPhysicalEntity) {
                if (this.methods.bioEntityTYpe == 1) {
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
            for (BioPathway p : network.getPathwayList().values()) {
                if(p.getReactions().values().contains(r)){
                }
            }
            this.reactionSet.add(r);
        }
    }

    public HashMap<BioEntity, Double> computeEnrichment() {
        HashSet<BioEntity> entityType2Enrich = new HashSet();

       for ( BioChemicalReaction r : this.reactionSet){
            //System.out.println("Reac: " + r.getName() + ": " + r.getListOfGeneNames().size());
           switch (this.entityType2Enrich){
               case 2:
                   entityType2Enrich.add(r);
                   break;
               case 3:
                   entityType2Enrich.addAll(r.getPathwayList().values());
                   break;
               /*case 6:
                   entityType2Enrich.addAll(r.getListOfGenes().values());
                   break;*/
           };
        }

       HashMap<BioEntity, Double> res = new HashMap();
       for (BioEntity p : entityType2Enrich){
            //System.out.println("PathName: " + p.getName() + ": " + p.getGenes().size());
            res.put(p, this.getPvalue(p));
        }
        return res;
    }

    public Collection getMappedEntityInEnrichedEntity(BioEntity entityType2Enrich) {
        switch (this.entityType2Enrich){
            /*case 1:
                return ;*/
            case 2:
                return getMappedEntityInReaction(entityType2Enrich);
            case 3:
                return getMappedEntityInPathway(entityType2Enrich);
            /*case 4: case 5:
                return getMappedEntityWithProtein(entityType2Enrich);*/
            case 6:
                return getMappedEntityWithGene(entityType2Enrich);
        };
        return null;
    }

    public Collection getMappedEntityInReaction(BioEntity enrichedEntity) {
        BioChemicalReaction reaction = (BioChemicalReaction) enrichedEntity;

        switch (this.methods.bioEntityTYpe) {
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
        };
        return null;
    }

    public Collection getMappedEntityInPathway(BioEntity enrichedEntity) {
        BioPathway pathway = (BioPathway) enrichedEntity;
        switch (this.methods.bioEntityTYpe){
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
                Set <BioProtein> proteins = new HashSet<BioProtein>();
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

    public Collection getMappedEntityWithGene(BioEntity enrichedEntity) {
        BioGene g = (BioGene) enrichedEntity;
        Set <String> list_reactions_ID = network.getReactionsFromGene(g.getId());
        Set <BioEntity> list_mapped = new HashSet();

            switch (this.methods.bioEntityTYpe) {
                case 1:
                    for (String reac_ID : list_reactions_ID) {
                        BioChemicalReaction reaction = network.getBiochemicalReactionList().get(reac_ID);
                        list_mapped.addAll(reaction.getListOfSubstrates().values());
                        list_mapped.addAll(reaction.getListOfProducts().values());
                    }
                case 2:
                    for (String reac_ID : list_reactions_ID) {
                        list_mapped.add(network.getBiochemicalReactionList().get(reac_ID));
                    }
                case 4: case 5:
                   list_mapped.addAll(g.getProteinList().values());
            };
        return null;
    }

    public Collection getMappedEntityWithProtein(BioEntity enrichedEntity) {
        BioProtein p = (BioProtein) enrichedEntity;
        Set <String> list_reactions_ID = network.getReactionsFromGene(g.getId());
        Set <BioEntity> list_mapped = new HashSet();

        switch (this.methods.bioEntityTYpe) {
            case 1:
                for (String reac_ID : list_reactions_ID) {
                    BioChemicalReaction reaction = network.getBiochemicalReactionList().get(reac_ID);
                    list_mapped.addAll(reaction.getListOfSubstrates().values());
                    list_mapped.addAll(reaction.getListOfProducts().values());
                }
            case 2:
                for (String reac_ID : list_reactions_ID) {
                    list_mapped.add(network.getBiochemicalReactionList().get(reac_ID));
                }
            case 4: case 5:
                list_mapped.addAll(g.getProteinList().values());
        };
        return null;
    }

    public double getPvalue(BioEntity enrichedEntity) throws IllegalArgumentException {
        int fisherTestParameters[] = this.getFisherTestParameters(enrichedEntity);
            return this.exactFisherOneTailed(fisherTestParameters[0], fisherTestParameters[1],
                    fisherTestParameters[2], fisherTestParameters[3]);
    }

    public int[] getFisherTestParameters(BioEntity enrichedEntity) {
        Collection mappedEntityInEnrichedEntity = this.getMappedEntityInEnrichedEntity(enrichedEntity);
        //nb of mapped in the pathway
        int a = methods.intersect(mappedEntityInEnrichedEntity).size();
        //unmapped metabolites in the fingerprint
        int b = this.list_mappedEntities.size() - a;
        //unmapped metabolites in the pathway
        int c = mappedEntityInEnrichedEntity.size() - a;
        //remaining metabolites in the network
        int d = methods.getEntitySetInNetwork().size() - (a + b + c);

        int fisherTestParameters[] = {a, b, c, d};
        //System.out.println(pathway.getName() + ": " + Arrays.toString(fisherTestParameters));
        return fisherTestParameters;
    }

        public HashMap<BioEntity, Double> benjaminiHochbergCorrection(HashMap<BioEntity, Double> pvalues) {
        ArrayList<BioEntity> orderedEnrichedEntities = this.sortPval(pvalues);
        HashMap<BioEntity, Double> adjPvalues = new HashMap();
        BioEntity p, p_next;
        double pval, adjPval;

        for(int k = 0; k < orderedEnrichedEntities.size(); ++k) {
            p = (BioEntity)orderedEnrichedEntities.get(k);
            pval = ((Double)pvalues.get(p)).doubleValue();
           // System.out.println("Pval : " + p.getName() + ": " + pval);

            //case of probability equality
            if (k+1 < orderedEnrichedEntities.size()){
                p_next = (BioEntity)orderedEnrichedEntities.get(k+1);
                double pval_next = ((Double)pvalues.get(p_next)).doubleValue();
                if (pval_next == pval) {
                    //System.out.println(true);
                    adjPval = pval * (double) pvalues.size() / ((double) k + (double) 1.5);
                    adjPvalues.put(p, new Double(adjPval));
                    adjPvalues.put(p_next, new Double(adjPval));
                    continue;
                }
            }

            //default case
            adjPval = pval * (double) pvalues.size() / (double) (k + 1);
            adjPvalues.put(p, new Double(adjPval));
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
        HashMap<BioEntity, Double> adjPvalues = new HashMap();

        for (BioEntity p : pvalues.keySet()){
            double pval = (Double)pvalues.get(p);
            double adjPval = pval * (double)pvalues.size();
            adjPvalues.put(p, new Double(adjPval));
        }

        return adjPvalues;
    }

    public ArrayList<BioEntity> sortPval(HashMap<BioEntity, Double> map) {
        ArrayList<BioEntity> orderedEnrichedEntities = new ArrayList(map.keySet());
        Collections.sort(orderedEnrichedEntities, new PathwayEnrichmentCalculation.significanceComparator(map));
        return orderedEnrichedEntities;
    }

    static class significanceComparator implements Comparator<BioEntity> {
        HashMap<BioEntity, Double> pvalMap;

        public significanceComparator(HashMap<BioEntity, Double> pvalMap) {
            this.pvalMap = pvalMap;
        }

        public int compare(BioEntity o1, BioEntity o2) {
            return Double.compare((Double)this.pvalMap.get(o1), (Double)this.pvalMap.get(o2));
        }
    }
    
    public static double exactFisherOneTailed(int a, int b, int c, int d) {
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
}