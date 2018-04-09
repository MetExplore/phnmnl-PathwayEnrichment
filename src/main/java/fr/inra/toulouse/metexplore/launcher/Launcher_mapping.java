package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.io.Fingerprint;
import fr.inra.toulouse.metexplore.io.JSBML2Bionetwork;
import fr.inra.toulouse.metexplore.omics.Mapping;
import fr.inra.toulouse.metexplore.omics.Omics;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import parsebionet.biodata.BioNetwork;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

public class Launcher_mapping extends Launcher_Fingerprint {

    @Option(name = "-o2", aliases = "-outMap", usage = "Output file name for mapping result (by default: mapping.tsv).")
    protected String outFileMapping = "mapping.tsv";

    @Option(name = "-s", aliases = "-sbml", usage = "SBML file name (by default: Recon v2.02).")
    protected String sbml = "data/recon2.02_without_compartment.xml";

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", aliases = "-type", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: metabolites).")
    protected int entityType2Map = 1;

    @Option(name = "-name", usage = "Activate this option for a name mapping .")
    protected int nameMapping = -1;

    @Option(name = "-prec", aliases = "-precision", usage = "Indicate the allowed error in ppm (used in mass mapping).")
    protected int weightPrecision = 2;

    protected BioNetwork network;
    protected Fingerprint fingerprint;
    protected String[] warn_type = {"Type of ", " entity must be comprise between 1 and 6."};

    @SuppressWarnings("deprecation")

    public void testParameters(String[] args) throws CmdLineException2{
        if (this.entityType2Map < 1 || this.entityType2Map > 6) {
            throw new CmdLineException2(warn_type[0] + "mapped" + warn_type[1]);
        }

        if (this.weightPrecision > 100 || this.weightPrecision < 1) {
            throw new CmdLineException2("Weight precision must be comprise between 1 and 100.");
        }

        if (this.nameColumn != this.nameMapping && this.nameMapping > 0) {
            this.nameColumn = this.nameMapping;
            if (Arrays.asList(args).contains("-nameCol")) {
                writeLog("[WARNING] You have set both name column"
                        + " and name mapping parameters and with different parameters.\n" +
                        "[WARNING] By default, the name mapping is activated with the column number of this parameter.\n");
            }
        }

        for (String arg : args) {
            if (Pattern.matches("-+prec.*", arg)) {
                if (this.weightColumn < 0) {
                    this.weightColumn = 2;
                    writeLog("[WARNING] Weight precision has been set without specify isotopic mass column in the fingerprint.\n" +
                            "[WARNING] By default, it has been set to the 2nd column.");
                }
                break;
            }
        }

        super.testParameters(args);
    }

    public Omics analyse(CmdLineParser parser, String[] args) throws IOException {

        this.fingerprint = (Fingerprint) super.analyse(parser, args);
        this.network = (new JSBML2Bionetwork(this.sbml)).getBioNetwork();
        Mapping map = new Mapping(logContent, network, fingerprint.getList_entities(), this.tab_inchiLayers,
                this.nameMapping, this.weightPrecision, this.outFileMapping, this.galaxyFile, this.entityType2Map);
        this.logContent = map.getLogContent();
        return map;
    }

    public static void main(String[] args) {
        exec(new Launcher_mapping(), args);
    }
}
