package fr.inra.toulouse.metexplore;

import fr.inra.toulouse.metexplore.io.JSBML2Bionetwork;
import fr.inra.toulouse.metexplore.io.WritingComportment;
import fr.inra.toulouse.metexplore.omics.PathwayEnrichment;
import junit.framework.TestCase;
import org.apache.commons.io.FileUtils;

import java.io.*;

import java.security.Permission;

public class Test extends TestCase implements WritingComportment {
    protected String galaxy, logContent, logFile = "temp/information.txt";
    protected File file;
    protected BufferedReader buffer;
    protected PathwayEnrichment pathEnr;


    public void setUp() throws Exception {
        //Initialization of the parameters before each tests
        super.setUp();
        //this.setSecurityManager();
        (new File("temp")).mkdir();
        this.galaxy = "";
        this.logContent = "";
    }

    protected void tearDown() throws Exception {
        super.tearDown();
        File folder = new File("temp");
        FileUtils.forceDelete(folder);
        folder.delete();
    }

    /*******************************
     *             SETTINGS
     *****************************/

    public static void setSecurityManager() {
        //Allow to catch exit(1) when object is called in a main function (or by other object)
        //but not when running this class as testunit
        System.setSecurityManager(new SecurityManager() {

            @Override
            public void checkPermission(Permission perm) {
            }

            @Override
            public void checkExit(int status) {
                throw new ThreadDeath();
            }

        });
    }

    /************setWriteOutput************/
    public void setBufferReader(String outFile) {
        this.file = new File(outFile);
        try {
            this.buffer = new BufferedReader(new FileReader(this.file));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setBufferTest(String fileName, String header, String line) {
        setBufferReader(fileName);
        try {
            if (fileName.equals(this.logFile)) {
                this.logContent = this.pathEnr.getLogContent();
                writeOutput(this.logContent,this.file);
            }
            assertEquals(buffer.readLine(), header);
            assertEquals(buffer.readLine(), line);
            assertEquals(buffer.readLine(), null);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



}
