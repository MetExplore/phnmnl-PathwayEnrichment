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
    protected String galaxy, text4outputFileInfo;
    protected ArrayList<String[]> list_fingerprint;//input file after formatting and filtering
    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;
    protected WritingComportment writingComportment;
    protected OmicsMethods methods;
    protected int bioEntityType;
    protected File log;
    static int nbInstance=0;


    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        nbInstance++;
        this.galaxy = galaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.writingComportment = new WritingComportment();
        this.bioEntityType=bioEntityType;
        this.methods = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }

    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int bioEntityType){
        nbInstance++;
        this.galaxy = galaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<BioEntity, String>();
        this.network = network;
        this.writingComportment = new WritingComportment();
        this.bioEntityType=bioEntityType;
        this.methods = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }

    public void writeOutputInfo() throws IOException {
        if (this.galaxy != "") {//if "writing console log in a file" functionality is activated
            log = new File(this.galaxy);
            if(nbInstance==1){
                //avoid to erase the file if Omics class is called more than one time
                if(log.isFile()) {
                    //write a new file if already exists
                    log.delete();
                }
                log.createNewFile();
            }
            BufferedWriter b = new BufferedWriter(new FileWriter(log, true));
            b.write(text4outputFileInfo);
            b.close();
        }
    }

    public void writeLog(String message) {
        System.out.println(message.replaceAll("\n", ""));
        text4outputFileInfo += message;
    }
}