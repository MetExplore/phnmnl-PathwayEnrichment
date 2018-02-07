package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public abstract class Omics {

    protected Boolean ifGalaxy;
    protected String text4outputFileInfo;
    protected HashMap<String, String[]> list_fingerprint;//input file after formatting and filtering
    protected Set<BioEntity> list_mappedEntities; //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    protected BioNetwork network;
    protected WritingComportment writingComportment;

    public Omics (Boolean ifGalaxy, HashMap<String, String[]> list_fingerprint,
                       Set<BioEntity> list_mappedEntities, BioNetwork network){
        this.ifGalaxy =ifGalaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.writingComportment = new WritingComportment();
    }

    public Omics (Boolean ifGalaxy, HashMap<String, String[]> list_fingerprint,
                  BioNetwork network){
        this.ifGalaxy =ifGalaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashSet<BioEntity>();
        this.network = network;
        this.writingComportment = new WritingComportment();
    }

    public void writeOutputInfo() throws IOException {
        if (this.ifGalaxy) {//if "writing console output in a file" functionality is activated
            File f = new File("information.tsv");
            f.createNewFile();
            BufferedWriter b = new BufferedWriter(new FileWriter(f));
            b.write(text4outputFileInfo);
            b.close();
        }
    }

    public void writeLog(String message) {
        System.out.println(message.replaceAll("\n", ""));
        text4outputFileInfo += message;
    }
}
