package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public abstract class Omics {
    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;
    protected WritingComportment write;
    protected OmicsMethods omics;
    protected int bioEntityType;


    public Omics (String galaxy, HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }

    public Omics (String galaxy, BioNetwork network, int bioEntityType){
        this.list_mappedEntities = new HashMap<BioEntity, String>();
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }
}