package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.Fingerprint;
import fr.inra.toulouse.metexplore.JSBML2Bionetwork4Galaxy;
import fr.inra.toulouse.metexplore.omics.Mapping;
import fr.inra.toulouse.metexplore.omics.Omics;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import parsebionet.biodata.BioNetwork;
import java.io.IOException;

public class Launcher_mapping extends Launcher_Fingerprint {

    @Option(name = "-o2", aliases = "--outMap", usage = "Output file name for mapping result (by default: mapping.tsv).")
    protected String outFileMapping = "mapping.tsv";

    @Option(name = "-s", aliases = "--sbml", usage = "SBML file name (by default: Recon v2.02).")
    protected String sbml = "data/recon2.02_without_compartment.xml";

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", aliases = "--type", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: metabolites).")
    protected int entityType2Map = 1;

    @Option(name = "-name", usage = "Activate this option for a name mapping .")
    protected Boolean nameMapping = false;

    protected BioNetwork network;
    protected Fingerprint fingerprint;

    @SuppressWarnings("deprecation")
    public void printInfo(CmdLineParser parser, String[] args) throws CmdLineException {

        super.printInfo(parser, args);

        if (this.entityType2Map < 1 || this.entityType2Map > 6) {
            throw new CmdLineException("Type of biological object must be between 1 and 6.");
        }
    }

    public Omics exec(CmdLineParser parser, String[] args) throws IOException {

        this.fingerprint = (Fingerprint) super.exec(parser, args);
        this.network = (new JSBML2Bionetwork4Galaxy(this.sbml)).getBioNetwork();

        return new Mapping(network, fingerprint.getList_entities(), this.tab_inchiLayers,
                this.nameMapping, this.outFileMapping, this.galaxy, this.entityType2Map);
    }

    public static void main(String[] args) {
        Launcher_mapping launch = new Launcher_mapping();
        try {
            launch.exec(new CmdLineParser(launch), args);
            launch.write.writeOutputInfo();
        } catch (IOException e) {
            e.printStackTrace();
        }
        timeCalculation(System.nanoTime() - startTime);
    }
}
