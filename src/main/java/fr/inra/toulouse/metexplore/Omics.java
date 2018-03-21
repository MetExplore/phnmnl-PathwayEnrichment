package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.util.ArrayList;
import java.util.HashMap;

public abstract class Omics {
    protected ArrayList<String[]> list_fingerprint;//input file after formatting and filtering
    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;
    protected WritingComportment write;
    protected OmicsMethods omics;
    protected int bioEntityType;
    protected String typeOfMappedEntity;


    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
        this.typeOfMappedEntity = omics.getTypeOfEntity(bioEntityType).toLowerCase();
    }

    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int bioEntityType){
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<>();
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
        this.typeOfMappedEntity = omics.getTypeOfEntity(bioEntityType).toLowerCase();
    }

}