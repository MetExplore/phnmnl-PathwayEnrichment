package fr.inra.toulouse.metexplore;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import parsebionet.biodata.*;

import java.util.HashSet;
import java.util.Set;

import static java.lang.System.exit;

public class FilterbyBioData {
    @Option(name = "-h", usage = "Prints this help.")
    public boolean phelp = false;

    @Option(name = "-s", usage = "SBML file name.")
    public String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name = "-t", usage = "Type of biological object selected : 1 for metabolites, 2 for reactions, 3 for pathway (by default: metabolites).")
    public int typeOfBioEntity = 1;

    @Option(name = "-id", usage = "[REQUIRED] List containing the id of the object(s) in the network: separated by comma (e.g., M_taur,M_mal_L,M_pnto_R).")
    public String selectedID;

    public BioNetwork network;
    public Set<String> list_selectedID = new HashSet<String>();


    public void setParameters(String sbml, String selectedID) {
        this.network = (new JSBML2Bionetwork4Galaxy(sbml)).getBioNetwork();

        for (String id : selectedID.split(",")) {
            this.list_selectedID.add(id);
        }
    }

    public void filterBioEntityByID() {

        System.out.println(this.list_selectedID);

        if (this.typeOfBioEntity == 1) {
            Set<String> list_BioEntity = this.network.getPhysicalEntityList().keySet();
            int initialSizeInNetwork = list_BioEntity.size();
            this.network.filterByMetabolites(this.list_selectedID);
            System.out.println(list_BioEntity);
            System.out.println(list_BioEntity.size());
            this.generateErrorIfExists(list_BioEntity, initialSizeInNetwork);

        } else if (this.typeOfBioEntity == 2) {
            Set<String> list_BioEntity = this.network.getBiochemicalReactionList().keySet();
            int initialSizeInNetwork = list_BioEntity.size();
            this.network.filterByReactions(this.list_selectedID);
            this.generateErrorIfExists(list_BioEntity, initialSizeInNetwork);
        } else {
            System.out.println("Here !");
            Set<String> list_BioEntity = this.network.getPathwayList().keySet();
            //this.list_selectedID=new HashSet<String>();
            //this.list_selectedID.add("Alkaloid synthesis");
            int initialSizeInNetwork = list_BioEntity.size();
            this.network.filterByPathways(this.list_selectedID);
            System.out.println(list_BioEntity.size());
            this.generateErrorIfExists(list_BioEntity, initialSizeInNetwork);
        }
    }

    public void generateErrorIfExists(Set <String> obtainedList, int initialSizeInNetwork){
        int nbExpected = this.list_selectedID.size();
        if (nbExpected != obtainedList.size()) {
            if (obtainedList.size()==initialSizeInNetwork)
                System.err.println("None of the specified ID have been found.\nPlease, check both the type of your biological object and your ID list");
            else {
                this.list_selectedID.removeAll(obtainedList);
                String missingID=this.list_selectedID.toString();
                //TODO: remove []
                System.err.println(missingID + " have not been found.\nPlease, check the corresponding ID.");
            }
            exit(1);
        }
    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        FilterbyBioData launch = new FilterbyBioData();
        CmdLineParser parser = new CmdLineParser(launch);

        try {
            parser.parseArgument(args);

            //Print help
            if (launch.phelp) {
                System.out.println("Options:");
                parser.printUsage(System.out);
                exit(0);
            }

            if (launch.typeOfBioEntity < 1 || launch.typeOfBioEntity > 3) {
                throw new CmdLineException("Type of biological object must be 1, 2 or 3.");
            }

        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("Options:");
            parser.printUsage(System.err);
            exit(1);
        }

        launch.setParameters(launch.sbml, launch.selectedID);
        launch.filterBioEntityByID();
    }

    public void test(String[] args){
       FilterbyBioData launch = new FilterbyBioData();
       launch.main(args);
    }

    public void test() {
        String[][] list_args = {{"-t 4"},{"-id M_taur,M_mal_L,M_pnto_R,M_malt"},{"-t 2 -id R_ILETA"},{"-t 3 -id Unassigned"},{"-t 3 -id Alkaloid synthesis"}};
        for (String[] args : list_args) {
            test(args);
        }
        //TODO: M_mal_L unrocognize
        //TODO: Unknown
        //TODO: space allowed for pathway
    }
}
