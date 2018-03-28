package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.io.WritingComportment;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public abstract class Omics implements WritingComportment{
    protected ArrayList<String[]> list_fingerprint;//input file after formatting and filtering

    public HashMap<BioEntity, String> getList_mappedEntities() {
        return list_mappedEntities;
    }

    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;

    public WritingComportment getWrite() {
        return write;
    }

    protected WritingComportment write;
    protected int entityType2Map;
    protected String typeOfMappedEntity;
    protected String galaxy;

    public String getLogContent() {
        return logContent;
    }

    protected String logContent;

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int entityType2Map){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.entityType2Map=entityType2Map;
        this.typeOfMappedEntity = getTypeOfEntity(entityType2Map).toLowerCase();
    }

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int entityType2Map){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<>();
        this.network = network;
        this.entityType2Map=entityType2Map;
        this.typeOfMappedEntity = getTypeOfEntity(entityType2Map).toLowerCase();
    }

    public HashMap getEntitySetInNetwork(int bioEntityType) {
        switch (bioEntityType){
            case 1:
                return this.network.getPhysicalEntityList();
            case 2:
                return this.network.getBiochemicalReactionList();
            case 3:
                return this.network.getPathwayList();
            case 4:
                return this.network.getEnzymeList();
            case 5:
                return this.network.getProteinList();
            case 6:
                return this.network.getGeneList();
        }
        return null;
    }

    public HashMap getEntitySetInNetwork() {
       return getEntitySetInNetwork(this.entityType2Map);
    }

    public HashSet<BioEntity> intersect(Collection<BioEntity> set2) {
        HashSet<BioEntity> inter = new HashSet<>();
        for (BioEntity bpe : set2){
            if (this.list_mappedEntities.keySet().contains(bpe)) inter.add(bpe);
        }
        return inter;
    }

    public String getTypeOfEntity(int entityType){
        switch (entityType) {
            case 1:
                return "Metabolite";
            case 2:
                return "Reaction";
            case 3:
                return "Pathway";
            case 4:
                return "Enzyme";
            case 5:
                return "Protein";
            case 6:
                return "Gene";
        }
        return null;
    }

    public void writeLog(String warning){
        this.logContent = writeLog(this.logContent,warning);
    }
}