package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.OmicsMethods;
import fr.inra.toulouse.metexplore.WritingComportment;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.util.ArrayList;
import java.util.HashMap;

public abstract class Omics implements WritingComportment, OmicsMethods{
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
    protected OmicsMethods omics;
    protected int bioEntityType;
    protected String typeOfMappedEntity;
    protected String galaxy;

    public String getLogContent() {
        return logContent;
    }

    protected String logContent;

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.bioEntityType=bioEntityType;
        this.typeOfMappedEntity = getTypeOfEntity(bioEntityType).toLowerCase();
    }

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int bioEntityType){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<>();
        this.network = network;
        this.bioEntityType=bioEntityType;
        this.typeOfMappedEntity = getTypeOfEntity(bioEntityType).toLowerCase();
    }

}