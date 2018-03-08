package fr.inra.toulouse.metexplore;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class WritingComportment {

    protected String galaxy;
    static String text4outputFileInfo = "";
    protected File log;
    static int nbInstance = 0;

    public WritingComportment (String galaxyOut) {
        nbInstance++;
        this.galaxy = galaxyOut;
        if(this.galaxy != "") {
            this.log = new File(this.galaxy);
            if (nbInstance == 1) {
                //avoid to erase the file if this class is called more than one time
                if (this.log.isFile()) {
                    //write a new file if already exists
                    this.log.delete();
                }
                try {
                    this.log.createNewFile();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public void writeOutputInfo() throws IOException {
        if (this.galaxy != "") {//if "writing console log in a file" functionality is activated
            this.log = new File(this.galaxy);
            BufferedWriter b = new BufferedWriter(new FileWriter(this.log, true));
            b.write(text4outputFileInfo);
            b.close();
        }
    }

    public void writeLog(String message) {
        //Avoid to skip a line
        System.out.println(message.replaceAll("\n$", ""));
        //TODO: regex, enlever celui de la fin
        text4outputFileInfo += message;
    }

    public String removeSciNot(double value) {
        String tmp = (String.valueOf(value));
        if (value < 1e-3) {
            String[] splitting = tmp.split("E");
            String power = splitting[1];
            if (power.charAt(1) == '0') {//remove the O after the E if exists
                power = power.substring(2, power.length());
            } else power = power.substring(1, power.length());
            String[] number = splitting[0].split("\\.");//obtain the integer and the decimal parts
            return "0." + (new String(new char[Integer.parseInt(power) - 1]).replace("\0", "0")) + number[0] + number[1];
        }
        return tmp;
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }
}
