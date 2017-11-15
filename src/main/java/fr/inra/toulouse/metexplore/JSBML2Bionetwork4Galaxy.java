package fr.inra.toulouse.metexplore;


import org.apache.log4j.Level;
import org.w3c.dom.CharacterData;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import parsebionet.io.PreProcessSBML;
import parsebionet.utils.JSBMLUtils;
import parsebionet.utils.StringUtils;

import javax.xml.stream.XMLStreamException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.sbml.jsbml.AbstractNamedSBase;
import org.sbml.jsbml.Compartment;
import org.sbml.jsbml.Constraint;
import org.sbml.jsbml.Event;
import org.sbml.jsbml.FunctionDefinition;
import org.sbml.jsbml.InitialAssignment;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.ListOf;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.ModifierSpeciesReference;
import org.sbml.jsbml.Parameter;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.Rule;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLException;
import org.sbml.jsbml.SBMLReader;
import org.sbml.jsbml.SBMLWriter;
import org.sbml.jsbml.SBase;
import org.sbml.jsbml.Species;
import org.sbml.jsbml.SpeciesReference;
import org.sbml.jsbml.SpeciesType;
import org.sbml.jsbml.Unit;
import org.sbml.jsbml.UnitDefinition;
import parsebionet.biodata.BioAnnotation;
import parsebionet.biodata.BioChemicalReaction;
import parsebionet.biodata.BioCompartment;
import parsebionet.biodata.BioCompartmentType;
import parsebionet.biodata.BioComplex;
import parsebionet.biodata.BioEntity;
import parsebionet.biodata.BioGene;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.biodata.BioPhysicalEntityParticipant;
import parsebionet.biodata.BioProtein;
import parsebionet.biodata.BioRef;
import parsebionet.biodata.BioUnitDefinition;
import parsebionet.biodata.Comment;
import parsebionet.biodata.Flux;
import parsebionet.biodata.Notes;
import parsebionet.biodata.UnitSbml;

public class JSBML2Bionetwork4Galaxy{
    private BioNetwork bioNetwork;
    private Set<String> Warnings = new HashSet();
    private String errorThread;
    private SBMLDocument SBMLdoc;
    private Model jSBMLmodel;
    public ArrayList<String> errorSet = new ArrayList();
    private boolean processing = true;
    private HashMap<String, String> IdCorrespondance = new HashMap();

    public static void main(String[] args) {
        String inputFile = "/home/bmerlet/Documents/PathwayEnrichment/recon2.v03_ext_noCompartment_noTransport_v2.xml";
        JSBML2Bionetwork4Galaxy jtb = new JSBML2Bionetwork4Galaxy(inputFile);
        BioNetwork network = jtb.getBioNetwork();
        network.printBioNetworkSizeToErr();
    }

    public JSBML2Bionetwork4Galaxy(String inputFile) {
        System.out.println("Pre-processing File...");
        try {
            this.setLogMode();
        } catch (IOException var6) {
            var6.printStackTrace();
        }

        File PreProcesFile = null;

        try {
            PreProcesFile = PreProcessSBML.process(inputFile, this.IdCorrespondance);
        } catch (IOException var5) {
            var5.printStackTrace();
        }

        this.bioNetwork = new BioNetwork();

        try {
            this.setSBMLdoc(SBMLReader.read(PreProcesFile));
            this.errorSet = this.getErrorsFromLogs();
        } catch (XMLStreamException | IOException var4) {
            var4.printStackTrace();
            this.bioNetwork = null;
        }

        if (this.bioNetwork != null && this.errorSet.isEmpty()) {
            this.convertToBioNetwork();
        } else {
            System.err.println("Error(s) while parsing the the SBML file: ");
        }

    }

    public void convertToBioNetwork() {
        long startTime = System.nanoTime();
        this.jSBMLmodel = this.SBMLdoc.getModel();
        this.bioNetwork.setId(this.getOriginalSiD(this.jSBMLmodel));
        this.bioNetwork.setName(this.jSBMLmodel.getName());
        this.bioNetwork.setType("sbml" + this.jSBMLmodel.getLevel() + "." + this.jSBMLmodel.getVersion());
        if (!StringUtils.isVoid(this.jSBMLmodel.getMetaId())) {
            try {
                this.bioNetwork.setModelAnnot(new BioAnnotation(this.getOriginalMetaiD(this.jSBMLmodel), this.jSBMLmodel.getAnnotationString()));
                this.bioNetwork.setModelNotes(new Notes(this.jSBMLmodel.getNotesString()));
            } catch (XMLStreamException var6) {
                var6.printStackTrace();
            }
        }

        try {
            this.parseSbmlListOfUnitDefinitions(this.jSBMLmodel);
            this.parseSbmlListOfCompartments(this.jSBMLmodel);
            this.parseSbmlListOfReactions(this.jSBMLmodel);
            this.parseForUnusedSpecies(this.jSBMLmodel);
        } catch (Exception var5) {
            this.errorThread = var5.toString() + "\n" + var5.getStackTrace()[0].toString();
            System.err.println(this.errorThread);
            System.exit(1);
        }

    }

    public void parseSbmlListOfUnitDefinitions(Model jSBMLmodel) {
        long numUnitDef = (long)jSBMLmodel.getNumUnitDefinitions();
        if (numUnitDef != 0L) {
            for(int i = 0; (long)i < numUnitDef; ++i) {
                UnitDefinition jSBMLUD = jSBMLmodel.getUnitDefinition(i);
                BioUnitDefinition bionetUD = new BioUnitDefinition(this.getOriginalSiD(jSBMLUD), jSBMLUD.getName());
                if (jSBMLUD.getName().isEmpty()) {
                    bionetUD.setName(this.getOriginalSiD(jSBMLUD));
                }

                ListOf<Unit> listofunits = jSBMLUD.getListOfUnits();
                if (listofunits.size() != 0) {
                    for(int n = 0; n < listofunits.size(); ++n) {
                        Unit jSBMLUnit = (Unit)listofunits.get(n);
                        Unit.Kind kind = jSBMLUnit.getKind();
                        String Exp = "" + (int)jSBMLUnit.getExponent();
                        String Scale = "" + jSBMLUnit.getScale();
                        String Multiplier = "" + jSBMLUnit.getMultiplier();
                        UnitSbml bionetUnit = new UnitSbml(kind.getName().toUpperCase(), Exp, Scale, Multiplier);
                        bionetUD.addUnit(bionetUnit);
                    }
                }

                this.getBioNetwork().addUnitDefinition(bionetUD);
            }
        } else {
            this.Warnings.add("Warnings, no list of UnitDefinition, please check your SBML file");
        }

    }

    public void parseSbmlListOfCompartments(Model jSBMLmodel) {
        long numCompart = (long)jSBMLmodel.getNumCompartments();
        if (numCompart != 0L) {
            for(int i = 0; (long)i < numCompart; ++i) {
                Compartment jSBMLCompart = jSBMLmodel.getCompartment(i);
                String compartId = this.getOriginalSiD(jSBMLCompart);
                String compartName = jSBMLCompart.getName();
                if (compartName.isEmpty()) {
                    compartName = compartId;
                }

                BioCompartment bionetCompart = this.getBioNetwork().findbioCompartmentInList(compartId);
                if (bionetCompart == null) {
                    bionetCompart = new BioCompartment(compartName, compartId);
                    this.getBioNetwork().addCompartment(bionetCompart);
                }

                if (jSBMLCompart.isSetOutside()) {
                    Compartment outsideJSBMLComp = jSBMLmodel.getCompartment(jSBMLCompart.getOutside());
                    BioCompartment outsideCompart = this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(outsideJSBMLComp));
                    if (outsideCompart == null) {
                        outsideCompart = new BioCompartment(outsideJSBMLComp.getName(), this.getOriginalSiD(outsideJSBMLComp));
                        this.getBioNetwork().addCompartment(outsideCompart);
                    }

                    bionetCompart.setOutsideCompartment(outsideCompart);
                }

                if (jSBMLCompart.isSetCompartmentType()) {
                    BioCompartmentType bionetCompartType = this.getBioNetwork().findbioCompartmentTypeInList(jSBMLCompart.getCompartmentType());
                    if (bionetCompartType == null) {
                        bionetCompartType = new BioCompartmentType(jSBMLCompart.getCompartmentType(), jSBMLmodel.getCompartmentType(jSBMLCompart.getCompartmentType()).getName());
                        this.getBioNetwork().addCompartmentType(bionetCompartType);
                    }

                    bionetCompart.setCompartmentType(bionetCompartType);
                }

                if (jSBMLCompart.isSetUnits()) {
                    UnitDefinition JsbmlUnitDef = jSBMLmodel.getUnitDefinition(jSBMLCompart.getUnits());
                    BioUnitDefinition bionetUnitDef = this.getBioNetwork().findUnitInUnitDefinituin(this.getOriginalSiD(JsbmlUnitDef));
                    if (bionetUnitDef == null) {
                        System.err.println("SBML Error: The compartment " + jSBMLCompart.getName() + " reference units that are not defined in the ListOfUnitDefinition.");
                    } else {
                        bionetCompart.setUnit(bionetUnitDef);
                    }
                }

                bionetCompart.setSboterm(jSBMLCompart.getSBOTermID());
                bionetCompart.setConstant(jSBMLCompart.getConstant());
                if (jSBMLCompart.isSetSize()) {
                    bionetCompart.setSize(jSBMLCompart.getSize());
                } else {
                    bionetCompart.setSize(1.0D);
                }

                bionetCompart.setSpatialDimensions((int)jSBMLCompart.getSpatialDimensions());
                if (jSBMLCompart.isSetAnnotation()) {
                    BioAnnotation bionetCompartAnnot = null;

                    try {
                        bionetCompartAnnot = new BioAnnotation(this.getOriginalMetaiD(jSBMLCompart), jSBMLCompart.getAnnotationString());
                    } catch (XMLStreamException var12) {
                        var12.printStackTrace();
                    }

                    bionetCompart.setCompartAnnot(bionetCompartAnnot);
                }

                if (jSBMLCompart.isSetNotes()) {
                    Notes bionetCompartNotes = null;

                    try {
                        bionetCompartNotes = new Notes(jSBMLCompart.getNotesString());
                    } catch (XMLStreamException var11) {
                        var11.printStackTrace();
                    }

                    bionetCompart.setCompartNotes(bionetCompartNotes);
                }
            }
        } else {
            this.Warnings.add("Errors while parsing the list of compartments, please check your SBML file");
        }

    }

    public void parseSbmlListOfReactions(Model jSBMLmodel) throws XMLStreamException {
        long numReaction = (long)jSBMLmodel.getNumReactions();
        if (numReaction != 0L) {
            for(int i = 0; (long)i < numReaction; ++i) {
                Reaction jSBMLReaction = jSBMLmodel.getReaction(i);
                String reactionId = this.getOriginalSiD(jSBMLReaction);
                String reactionName = jSBMLReaction.getName();
                if (reactionName.isEmpty()) {
                    reactionName = reactionId;
                }

                BioChemicalReaction bionetReaction = new BioChemicalReaction(reactionId, reactionName);
                bionetReaction.setSboterm(jSBMLReaction.getSBOTermID());
                if (jSBMLReaction.isSetAnnotation()) {
                    bionetReaction.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLReaction), jSBMLReaction.getAnnotationString()));
                }

                if (jSBMLReaction.isSetFast()) {
                    bionetReaction.setSpontaneous((new Boolean(jSBMLReaction.getFast())).toString());
                } else {
                    bionetReaction.setSpontaneous("false");
                }

                if (jSBMLReaction.isSetReversible() && !jSBMLReaction.getReversible()) {
                    bionetReaction.setReversibility(false);
                } else {
                    bionetReaction.setReversibility(true);
                }

                int n;
                SpeciesReference jSBMLProduct;
                Species jSBMLModifierSpecies;
                BioPhysicalEntity modifier;
                String bionetMetabId;
                String Stoechio;
                String tempParticipantId;
                BioPhysicalEntityParticipant bionetProdPart;
                for(n = 0; n < jSBMLReaction.getNumReactants(); ++n) {
                    jSBMLProduct = jSBMLReaction.getReactant(n);
                    jSBMLModifierSpecies = jSBMLmodel.getSpecies(jSBMLProduct.getSpecies());
                    modifier = this.parseParticipantSpecies(jSBMLModifierSpecies);
                    bionetMetabId = modifier.getId();
                    Stoechio = jSBMLReaction.getReactant(n).getStoichiometry() + "";
                    tempParticipantId = bionetReaction.getId() + "__With__" + bionetMetabId;
                    bionetProdPart = new BioPhysicalEntityParticipant(tempParticipantId, modifier, Stoechio, modifier.getCompartment());
                    if (jSBMLReaction.getReactant(n).isSetConstant()) {
                        bionetProdPart.setIsConstant(jSBMLReaction.getReactant(n).getConstant());
                    } else {
                        bionetProdPart.setIsConstant(false);
                    }

                    if (jSBMLReaction.getReactant(n).isSetId()) {
                        bionetProdPart.setId(this.getOriginalSiD(jSBMLProduct));
                    }

                    bionetReaction.addLeftParticipant(bionetProdPart);
                }

                for(n = 0; n < jSBMLReaction.getNumProducts(); ++n) {
                    jSBMLProduct = jSBMLReaction.getProduct(n);
                    jSBMLModifierSpecies = jSBMLmodel.getSpecies(jSBMLProduct.getSpecies());
                    modifier = this.parseParticipantSpecies(jSBMLModifierSpecies);
                    bionetMetabId = modifier.getId();
                    Stoechio = jSBMLReaction.getProduct(n).getStoichiometry() + "";
                    tempParticipantId = bionetReaction.getId() + "__With__" + bionetMetabId;
                    bionetProdPart = new BioPhysicalEntityParticipant(tempParticipantId, modifier, Stoechio, modifier.getCompartment());
                    if (jSBMLReaction.getProduct(n).isSetConstant()) {
                        bionetProdPart.setIsConstant(jSBMLReaction.getProduct(n).getConstant());
                    } else {
                        bionetProdPart.setIsConstant(false);
                    }

                    if (jSBMLReaction.getProduct(n).isSetId()) {
                        bionetProdPart.setId(this.getOriginalSiD(jSBMLProduct));
                    }

                    bionetReaction.addRightParticipant(bionetProdPart);
                }

                bionetReaction.setListOfSubstrates();
                bionetReaction.setListOfProducts();

                for(n = 0; n < jSBMLReaction.getNumModifiers(); ++n) {
                    ModifierSpeciesReference jSBMLModifier = jSBMLReaction.getModifier(n);
                    jSBMLModifierSpecies = jSBMLmodel.getSpecies(jSBMLModifier.getSpecies());
                    modifier = this.parseModifierSpecies(jSBMLModifierSpecies);
                    bionetReaction.addEnz(modifier);
                    this.getBioNetwork().getEnzList().put(modifier.getId(), modifier);
                }

                KineticLaw kine = jSBMLReaction.getKineticLaw();
                if (kine != null) {
                    if (kine.getFormula() != null) {
                        bionetReaction.setKineticFormula(kine.getFormula());
                    }

                    boolean hasBounds = false;
                    n = new Integer (null);
                    UnitDefinition jsbmlUnit;
                    BioUnitDefinition UD;
                    Flux newflux;
                    if (jSBMLmodel.getLevel() < 3) {
                        for(n = 0; n < kine.getNumParameters(); ++n) {
                            jsbmlUnit = jSBMLmodel.getUnitDefinition(kine.getParameter(n).getUnits());
                            if (jsbmlUnit == null) {
                                UD = new BioUnitDefinition("dimensionless", "dimensionless");
                                newflux = new Flux("" + kine.getParameter(n).getValue(), UD);
                                bionetReaction.addFluxParam(kine.getParameter(n).getId(), newflux);
                            } else {
                                UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get(this.getOriginalSiD(jsbmlUnit));
                                if (!kine.getParameter(n).getId().equalsIgnoreCase("UPPER_BOUND") && !kine.getParameter(n).getName().equalsIgnoreCase("UPPER_BOUND")) {
                                    if (!kine.getParameter(n).getId().equalsIgnoreCase("LOWER_BOUND") && !kine.getParameter(n).getName().equalsIgnoreCase("LOWER_BOUND")) {
                                        if (UD != null) {
                                            newflux = new Flux("" + kine.getParameter(n).getValue(), UD);
                                            bionetReaction.addFluxParam(kine.getParameter(n).getId(), newflux);
                                        } else if (UD == null && kine.getParameter(n).getUnits().equalsIgnoreCase("dimensionless")) {
                                            UD = new BioUnitDefinition("dimensionless", "dimensionless");
                                            newflux = new Flux("" + kine.getParameter(n).getValue(), UD);
                                            bionetReaction.addFluxParam(kine.getParameter(n).getId(), newflux);
                                        }
                                    } else {
                                        if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("mmol_per_gDW_per_hr")) {
                                            UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("mmol_per_gDW_per_hr");
                                        } else if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("FLUX_UNIT")) {
                                            UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("FLUX_UNIT");
                                        }

                                        if (UD == null) {
                                            UD = new BioUnitDefinition();
                                            UD.setDefault();
                                            this.getBioNetwork().getUnitDefinitions().put(UD.getId(), UD);
                                        }

                                        newflux = new Flux("" + kine.getParameter(n).getValue(), UD);
                                        bionetReaction.setLowerBound(newflux);
                                        hasBounds = true;
                                    }
                                } else {
                                    if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("mmol_per_gDW_per_hr")) {
                                        UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("mmol_per_gDW_per_hr");
                                    } else if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("FLUX_UNIT")) {
                                        UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("FLUX_UNIT");
                                    }

                                    if (UD == null) {
                                        UD = new BioUnitDefinition();
                                        UD.setDefault();
                                        this.getBioNetwork().getUnitDefinitions().put(UD.getId(), UD);
                                    }

                                    newflux = new Flux("" + kine.getParameter(n).getValue(), UD);
                                    bionetReaction.setUpperBound(newflux);
                                    hasBounds = true;
                                }
                            }
                        }
                    } else {
                        for(n = 0; n < kine.getLocalParameterCount(); ++n) {
                            jsbmlUnit = jSBMLmodel.getUnitDefinition(kine.getLocalParameter(n).getUnits());
                            if (jsbmlUnit == null) {
                                UD = new BioUnitDefinition("dimensionless", "dimensionless");
                                newflux = new Flux("" + kine.getLocalParameter(n).getValue(), UD);
                                bionetReaction.addFluxParam(kine.getLocalParameter(n).getId(), newflux);
                            } else {
                                UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get(this.getOriginalSiD(jsbmlUnit));
                                if (!kine.getLocalParameter(n).getId().equalsIgnoreCase("UPPER_BOUND") && !kine.getLocalParameter(n).getName().equalsIgnoreCase("UPPER_BOUND")) {
                                    if (!kine.getLocalParameter(n).getId().equalsIgnoreCase("LOWER_BOUND") && !kine.getLocalParameter(n).getName().equalsIgnoreCase("LOWER_BOUND")) {
                                        if (UD != null) {
                                            newflux = new Flux("" + kine.getLocalParameter(n).getValue(), UD);
                                            bionetReaction.addFluxParam(kine.getLocalParameter(n).getId(), newflux);
                                        } else if (UD == null && kine.getLocalParameter(n).getUnits().equalsIgnoreCase("dimensionless")) {
                                            UD = new BioUnitDefinition("dimensionless", "dimensionless");
                                            newflux = new Flux("" + kine.getLocalParameter(n).getValue(), UD);
                                            bionetReaction.addFluxParam(kine.getLocalParameter(n).getId(), newflux);
                                        }
                                    } else {
                                        if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("mmol_per_gDW_per_hr")) {
                                            UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("mmol_per_gDW_per_hr");
                                        } else if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("FLUX_UNIT")) {
                                            UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("FLUX_UNIT");
                                        }

                                        if (UD == null) {
                                            UD = new BioUnitDefinition();
                                            UD.setDefault();
                                            this.getBioNetwork().getUnitDefinitions().put(UD.getId(), UD);
                                        }

                                        newflux = new Flux("" + kine.getLocalParameter(n).getValue(), UD);
                                        bionetReaction.setLowerBound(newflux);
                                        hasBounds = true;
                                    }
                                } else {
                                    if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("mmol_per_gDW_per_hr")) {
                                        UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("mmol_per_gDW_per_hr");
                                    } else if (UD == null && this.getBioNetwork().getUnitDefinitions().containsKey("FLUX_UNIT")) {
                                        UD = (BioUnitDefinition)this.getBioNetwork().getUnitDefinitions().get("FLUX_UNIT");
                                    }

                                    if (UD == null) {
                                        UD = new BioUnitDefinition();
                                        UD.setDefault();
                                        this.getBioNetwork().getUnitDefinitions().put(UD.getId(), UD);
                                    }

                                    newflux = new Flux("" + kine.getLocalParameter(n).getValue(), UD);
                                    bionetReaction.setUpperBound(newflux);
                                    hasBounds = true;
                                }
                            }
                        }
                    }

                    if (!hasBounds) {
                        bionetReaction.setUpperBound((Flux)null);
                        bionetReaction.setLowerBound((Flux)null);
                    }
                }

                if (jSBMLReaction.isSetNotes()) {
                    bionetReaction.setEntityNotes(new Notes(jSBMLReaction.getNotesString()));
                    this.parseReactionNotes(bionetReaction);
                }

                bionetReaction.setCompartment((BioCompartment)null);
                this.getBioNetwork().addBiochemicalReaction(bionetReaction);
            }
        } else {
            this.Warnings.add("Error while parsing the list of reaction. No reaction in the file, please check your SBML file");
        }

    }

    private void parseForUnusedSpecies(Model jSBMLmodel) throws XMLStreamException {
        long numspecies = (long)jSBMLmodel.getNumSpecies();

        for(int i = 0; (long)i < numspecies; ++i) {
            Species jSBMLSpecie = jSBMLmodel.getSpecies(i);
            String specieId = this.getOriginalSiD(jSBMLSpecie);
            if (!this.getBioNetwork().getPhysicalEntityList().containsKey(specieId) && !this.getBioNetwork().getProteinList().containsKey(specieId) && !this.getBioNetwork().getComplexList().containsKey(specieId)) {
                String specieName = jSBMLSpecie.getName();
                if (specieName.isEmpty()) {
                    specieName = specieId;
                }

                if (jSBMLSpecie.getSBOTerm() != 2 && jSBMLSpecie.getSBOTerm() != 247) {
                    if (jSBMLSpecie.getSBOTerm() == 252 || jSBMLSpecie.getSBOTerm() == 297) {
                        HashSet<BioProtein> Compo = this.parseComplexComponent(specieId);
                        Compartment jsbmlComp;
                        if (Compo.size() == 1) {
                            BioProtein unusedProtein = (BioProtein)Compo.toArray()[0];
                            unusedProtein.setId(specieId);
                            unusedProtein.setName(specieName);
                            unusedProtein.setBoundaryCondition(jSBMLSpecie.getBoundaryCondition());
                            unusedProtein.setConstant(jSBMLSpecie.getConstant());
                            unusedProtein.setSubstanceUnits(jSBMLSpecie.getSubstanceUnits());
                            unusedProtein.setSboterm(jSBMLSpecie.getSBOTermID());
                            jsbmlComp = jSBMLmodel.getCompartment(jSBMLSpecie.getCompartment());
                            unusedProtein.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
                            if (jSBMLSpecie.isSetInitialAmount()) {
                                unusedProtein.addInitialQuantity("amount", jSBMLSpecie.getInitialAmount());
                            } else if (jSBMLSpecie.isSetInitialConcentration()) {
                                unusedProtein.addInitialQuantity("concentration", jSBMLSpecie.getInitialConcentration());
                            }

                            if (jSBMLSpecie.isSetAnnotation()) {
                                unusedProtein.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLSpecie), jSBMLSpecie.getAnnotationString()));
                                this.parseMetaboliteAnnot(unusedProtein);
                            }

                            this.parseIDForGene(unusedProtein);
                            if (jSBMLSpecie.isSetNotes()) {
                                unusedProtein.setEntityNotes(new Notes(jSBMLSpecie.getNotesString()));
                                this.parseModifierNote((BioPhysicalEntity)unusedProtein);
                            }

                            this.AddNote(unusedProtein, "Unused");
                        } else {
                            BioComplex unusedComplex = new BioComplex(specieId, specieName);
                            unusedComplex.setBoundaryCondition(jSBMLSpecie.getBoundaryCondition());
                            unusedComplex.setConstant(jSBMLSpecie.getConstant());
                            unusedComplex.setSubstanceUnits(jSBMLSpecie.getSubstanceUnits());
                            unusedComplex.setSboterm(jSBMLSpecie.getSBOTermID());
                            jsbmlComp = jSBMLmodel.getCompartment(jSBMLSpecie.getCompartment());
                            unusedComplex.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
                            if (jSBMLSpecie.isSetInitialAmount()) {
                                unusedComplex.addInitialQuantity("amount", jSBMLSpecie.getInitialAmount());
                            } else if (jSBMLSpecie.isSetInitialConcentration()) {
                                unusedComplex.addInitialQuantity("Concentration", jSBMLSpecie.getInitialConcentration());
                            }

                            Iterator var11 = Compo.iterator();

                            while(var11.hasNext()) {
                                BioProtein protcompo = (BioProtein)var11.next();
                                protcompo.setCompartment(unusedComplex.getCompartment());
                                BioPhysicalEntityParticipant bionetParticipant = new BioPhysicalEntityParticipant(specieId + "__With__" + protcompo.getId(), protcompo);
                                unusedComplex.addComponent(bionetParticipant);
                            }

                            if (jSBMLSpecie.isSetAnnotation()) {
                                unusedComplex.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLSpecie), jSBMLSpecie.getAnnotationString()));
                                this.parseMetaboliteAnnot(unusedComplex);
                            }

                            if (jSBMLSpecie.isSetNotes()) {
                                unusedComplex.setEntityNotes(new Notes(jSBMLSpecie.getNotesString()));
                                this.parseModifierNote(unusedComplex);
                            }

                            this.AddNote(unusedComplex, "Unused");
                            this.getBioNetwork().addComplex(unusedComplex);
                        }
                    }
                } else {
                    this.parseParticipantSpecies(jSBMLSpecie);
                }
            }
        }

    }

    private BioPhysicalEntity parseParticipantSpecies(Species jSBMLSpecie) throws XMLStreamException {
        String specieId = this.getOriginalSiD(jSBMLSpecie);
        BioPhysicalEntity bionetSpecies = this.getBioNetwork().getBioPhysicalEntityById(specieId);
        if (bionetSpecies != null) {
            return bionetSpecies;
        } else {
            String specieName = jSBMLSpecie.getName();
            if (specieName.isEmpty()) {
                specieName = specieId;
            }

            bionetSpecies = new BioPhysicalEntity(specieId, specieName);
            bionetSpecies.setBoundaryCondition(jSBMLSpecie.getBoundaryCondition());
            bionetSpecies.setConstant(jSBMLSpecie.getConstant());
            bionetSpecies.setSubstanceUnits(jSBMLSpecie.getSubstanceUnits());
            bionetSpecies.setSboterm(jSBMLSpecie.getSBOTermID());
            Compartment jsbmlComp = this.jSBMLmodel.getCompartment(jSBMLSpecie.getCompartment());
            bionetSpecies.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
            if (jSBMLSpecie.isSetInitialAmount()) {
                bionetSpecies.addInitialQuantity("amount", jSBMLSpecie.getInitialAmount());
            } else if (jSBMLSpecie.isSetInitialConcentration()) {
                bionetSpecies.addInitialQuantity("Concentration", jSBMLSpecie.getInitialConcentration());
            }

            if (jSBMLSpecie.isSetAnnotation()) {
                bionetSpecies.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLSpecie), jSBMLSpecie.getAnnotationString()));
                this.parseMetaboliteAnnot(bionetSpecies);
            }

            if (jSBMLSpecie.isSetCharge()) {
                bionetSpecies.setCharge(jSBMLSpecie.getCharge() + "");
            }

            if (jSBMLSpecie.isSetNotes()) {
                bionetSpecies.setEntityNotes(new Notes(jSBMLSpecie.getNotesString()));
                this.parseMetaboliteNote(bionetSpecies);
            }

            if (jSBMLSpecie.isSetSpeciesType()) {
                SpeciesType stype = jSBMLSpecie.getSpeciesTypeInstance();
                if (stype.isSetSBOTerm()) {
                    bionetSpecies.setSboterm(stype.getSBOTermID());
                }

                if (stype.isSetAnnotation()) {
                    bionetSpecies.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(stype), stype.getAnnotationString()));
                    this.parseMetaboliteAnnot(bionetSpecies);
                }

                if (stype.isSetNotes()) {
                    bionetSpecies.setEntityNotes(new Notes(stype.getNotesString()));
                    this.parseMetaboliteNote(bionetSpecies);
                }
            }

            this.getBioNetwork().addPhysicalEntity(bionetSpecies);
            return bionetSpecies;
        }
    }

    private BioPhysicalEntity parseModifierSpecies(Species jSBMLModifierSpecies) throws XMLStreamException {
        String specieId = this.getOriginalSiD(jSBMLModifierSpecies);
        if (this.getBioNetwork().getProteinList().containsKey(specieId)) {
            return (BioPhysicalEntity)this.getBioNetwork().getProteinList().get(specieId);
        } else if (this.getBioNetwork().getComplexList().containsKey(specieId)) {
            return (BioPhysicalEntity)this.getBioNetwork().getComplexList().get(specieId);
        } else if (this.getBioNetwork().getPhysicalEntityList().containsKey(specieId)) {
            return (BioPhysicalEntity)this.getBioNetwork().getPhysicalEntityList().get(specieId);
        } else {
            String specieName = jSBMLModifierSpecies.getName();
            if (specieName.isEmpty()) {
                specieName = specieId;
            }

            HashSet<BioProtein> Compo = this.parseComplexComponent(specieId);
            Compartment jsbmlComp;
            if (Compo.size() == 1) {
                BioProtein proteinModifier = (BioProtein)Compo.toArray()[0];
                proteinModifier.setId(specieId);
                proteinModifier.setName(specieName);
                proteinModifier.setBoundaryCondition(jSBMLModifierSpecies.getBoundaryCondition());
                proteinModifier.setConstant(jSBMLModifierSpecies.getConstant());
                proteinModifier.setSubstanceUnits(jSBMLModifierSpecies.getSubstanceUnits());
                proteinModifier.setSboterm(jSBMLModifierSpecies.getSBOTermID());
                jsbmlComp = this.jSBMLmodel.getCompartment(jSBMLModifierSpecies.getCompartment());
                proteinModifier.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
                if (jSBMLModifierSpecies.isSetInitialAmount()) {
                    proteinModifier.addInitialQuantity("amount", jSBMLModifierSpecies.getInitialAmount());
                } else if (jSBMLModifierSpecies.isSetInitialConcentration()) {
                    proteinModifier.addInitialQuantity("concentration", jSBMLModifierSpecies.getInitialConcentration());
                }

                if (jSBMLModifierSpecies.isSetAnnotation()) {
                    proteinModifier.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLModifierSpecies), jSBMLModifierSpecies.getAnnotationString()));
                    this.parseMetaboliteAnnot(proteinModifier);
                }

                this.parseIDForGene(proteinModifier);
                if (jSBMLModifierSpecies.isSetNotes()) {
                    proteinModifier.setEntityNotes(new Notes(jSBMLModifierSpecies.getNotesString()));
                    this.parseModifierNote((BioPhysicalEntity)proteinModifier);
                }

                return proteinModifier;
            } else {
                BioComplex complexModifier = new BioComplex(specieId, specieName);
                complexModifier.setBoundaryCondition(jSBMLModifierSpecies.getBoundaryCondition());
                complexModifier.setConstant(jSBMLModifierSpecies.getConstant());
                complexModifier.setSubstanceUnits(jSBMLModifierSpecies.getSubstanceUnits());
                complexModifier.setSboterm(jSBMLModifierSpecies.getSBOTermID());
                jsbmlComp = this.jSBMLmodel.getCompartment(jSBMLModifierSpecies.getCompartment());
                complexModifier.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
                if (jSBMLModifierSpecies.isSetInitialAmount()) {
                    complexModifier.addInitialQuantity("amount", jSBMLModifierSpecies.getInitialAmount());
                } else if (jSBMLModifierSpecies.isSetInitialConcentration()) {
                    complexModifier.addInitialQuantity("Concentration", jSBMLModifierSpecies.getInitialConcentration());
                }

                Iterator var7 = Compo.iterator();

                while(var7.hasNext()) {
                    BioProtein protcompo = (BioProtein)var7.next();
                    protcompo.setCompartment(complexModifier.getCompartment());
                    BioPhysicalEntityParticipant bionetParticipant = new BioPhysicalEntityParticipant(specieId + "__With__" + protcompo.getId(), protcompo);
                    complexModifier.addComponent(bionetParticipant);
                }

                if (jSBMLModifierSpecies.isSetAnnotation()) {
                    complexModifier.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLModifierSpecies), jSBMLModifierSpecies.getAnnotationString()));
                    this.parseMetaboliteAnnot(complexModifier);
                }

                complexModifier.getIsGeneticallyPossible();
                if (jSBMLModifierSpecies.isSetNotes()) {
                    complexModifier.setEntityNotes(new Notes(jSBMLModifierSpecies.getNotesString()));
                    this.parseModifierNote(complexModifier);
                }

                this.getBioNetwork().addComplex(complexModifier);
                return complexModifier;
            }
        }
    }

    public void parseUnusedData(Model jSBMLmodel) throws SBMLException, XMLStreamException {
        SBMLWriter tmpWriter = new SBMLWriter();
        SBMLDocument tmpDoc = new SBMLDocument(jSBMLmodel.getLevel(), jSBMLmodel.getVersion());
        Model tmpModel;
        Iterator var5;
        if (jSBMLmodel.getListOfConstraints().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfConstraints().iterator();

            while(var5.hasNext()) {
                Constraint constraint = (Constraint)var5.next();
                tmpModel.addConstraint(constraint.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfConstraints", tmpWriter.writeSBMLToString(tmpDoc));
        }

        if (jSBMLmodel.getListOfEvents().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfEvents().iterator();

            while(var5.hasNext()) {
                Event event = (Event)var5.next();
                tmpModel.addEvent(event.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfEvents", tmpWriter.writeSBMLToString(tmpDoc));
        }

        if (jSBMLmodel.getListOfFunctionDefinitions().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfFunctionDefinitions().iterator();

            while(var5.hasNext()) {
                FunctionDefinition fctDef = (FunctionDefinition)var5.next();
                tmpModel.addFunctionDefinition(fctDef.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfFunctionDefinitions", tmpWriter.writeSBMLToString(tmpDoc));
        }

        if (jSBMLmodel.getListOfInitialAssignments().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfInitialAssignments().iterator();

            while(var5.hasNext()) {
                InitialAssignment initAss = (InitialAssignment)var5.next();
                tmpModel.addInitialAssignment(initAss.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfInitialAssignments", tmpWriter.writeSBMLToString(tmpDoc));
        }

        if (jSBMLmodel.getListOfParameters().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfParameters().iterator();

            while(var5.hasNext()) {
                Parameter param = (Parameter)var5.next();
                tmpModel.addParameter(param.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfParameters", tmpWriter.writeSBMLToString(tmpDoc));
        }

        if (jSBMLmodel.getListOfRules().size() != 0) {
            tmpModel = tmpDoc.createModel();
            var5 = jSBMLmodel.getListOfRules().iterator();

            while(var5.hasNext()) {
                Rule rule = (Rule)var5.next();
                tmpModel.addRule(rule.clone());
            }

            this.getBioNetwork().addUnusedSBMLdata("ListOfRules", tmpWriter.writeSBMLToString(tmpDoc));
        }

    }

    private String getOriginalSiD(AbstractNamedSBase JSBMLComponent) {
        String modifiedSiD = JSBMLComponent.getId();
        return (String)this.IdCorrespondance.get(modifiedSiD);
    }

    private String getOriginalMetaiD(AbstractNamedSBase JSBMLComponent) {
        String modifiedSiD = JSBMLComponent.getMetaId();
        return (String)this.IdCorrespondance.get(modifiedSiD);
    }

    public String ListOfToString(ListOf<? extends SBase> listOf) {
        String list = "";
        int i = 0;

        for(int n = listOf.getChildCount(); i < n; ++i) {
            list = list + listOf.getChildAt(i).toString();
        }

        return list;
    }

    public HashSet<BioProtein> parseComplexComponent(String complexID) throws XMLStreamException {
        HashSet<BioProtein> ComponentList = new HashSet();
        String regexCompartIds = "[";

        String compartId;
        for(Iterator var4 = this.getBioNetwork().getCompartments().keySet().iterator(); var4.hasNext(); regexCompartIds = regexCompartIds + "(" + compartId + ")") {
            compartId = (String)var4.next();
        }

        regexCompartIds = regexCompartIds + "]";
        String[] simpleIDs = complexID.split("(?<=_" + regexCompartIds + "{0,1})(?=_)()");
        if (simpleIDs.length == 1) {
            BioProtein bionetProt = (BioProtein)this.getBioNetwork().getProteinList().get(simpleIDs[0]);
            if (bionetProt != null) {
                ComponentList.add(bionetProt);
            } else {
                String geneId = simpleIDs[0];
                bionetProt = new BioProtein(geneId);
                this.parseIDForGene(bionetProt);
                ComponentList.add(bionetProt);
                this.getBioNetwork().addProtein(bionetProt);
            }
        } else {
            String[] var17 = simpleIDs;
            int var18 = simpleIDs.length;

            for(int var7 = 0; var7 < var18; ++var7) {
                String value = var17[var7];
                String encodedValue = StringUtils.sbmlEncode(value);
                Species jSBMLProt = this.jSBMLmodel.getSpecies(encodedValue);
                BioProtein bionetProt = (BioProtein)this.getBioNetwork().getProteinList().get(value);
                if (bionetProt != null) {
                    ComponentList.add(bionetProt);
                } else if (jSBMLProt != null) {
                    bionetProt = new BioProtein(value, jSBMLProt.getName());
                    bionetProt.setBoundaryCondition(jSBMLProt.getBoundaryCondition());
                    bionetProt.setConstant(jSBMLProt.getConstant());
                    bionetProt.setSubstanceUnits(jSBMLProt.getSubstanceUnits());
                    bionetProt.setSboterm(jSBMLProt.getSBOTermID());
                    Compartment jsbmlComp = this.jSBMLmodel.getCompartment(jSBMLProt.getCompartment());
                    bionetProt.setCompartment(this.getBioNetwork().findbioCompartmentInList(this.getOriginalSiD(jsbmlComp)));
                    if (jSBMLProt.isSetInitialAmount()) {
                        bionetProt.addInitialQuantity("amount", jSBMLProt.getInitialAmount());
                    } else if (jSBMLProt.isSetInitialConcentration()) {
                        bionetProt.addInitialQuantity("Concentration", jSBMLProt.getInitialConcentration());
                    }

                    if (jSBMLProt.isSetAnnotation()) {
                        bionetProt.setEntityAnnot(new BioAnnotation(this.getOriginalMetaiD(jSBMLProt), jSBMLProt.getAnnotationString()));
                        this.parseMetaboliteAnnot(bionetProt);
                    }

                    this.parseIDForGene(bionetProt);
                    if (jSBMLProt.isSetNotes()) {
                        bionetProt.setEntityNotes(new Notes(jSBMLProt.getNotesString()));
                        this.parseModifierNote((BioPhysicalEntity)bionetProt);
                    }

                    this.getBioNetwork().addProtein(bionetProt);
                    ComponentList.add(bionetProt);
                } else {
                    bionetProt = new BioProtein(value, value + "_TH");
                    bionetProt.setBoundaryCondition(true);
                    bionetProt.setConstant(false);
                    bionetProt.setSboterm("SBO:0000252");
                    bionetProt.setSubstanceUnits("false");
                    String protNotes = "<notes><body xmlns=\"http://www.w3.org/1999/xhtml\"><p> METEXPLORE_NOTE: This theoretical protein is a component of the protein complexe " + complexID + "</p></body></notes>";
                    bionetProt.setEntityNotes(new Notes(protNotes));
                    String geneId = value;
                    if (value.startsWith("_")) {
                        geneId = value.substring(1, value.length());
                    }

                    if (Pattern.compile("^.+_\\w$", 2).matcher(geneId).matches()) {
                        geneId = geneId.substring(0, geneId.length() - 2);
                    }

                    geneId = geneId.replaceAll("_", ".");
                    BioGene newGene;
                    if (this.getBioNetwork().getGeneList().containsKey(geneId)) {
                        newGene = (BioGene)this.getBioNetwork().getGeneList().get(geneId);
                    } else {
                        newGene = new BioGene(geneId, geneId);
                        this.getBioNetwork().addGene(newGene);
                    }

                    bionetProt.addGene(newGene);
                    this.getBioNetwork().addProtein(bionetProt);
                    ComponentList.add(bionetProt);
                }
            }
        }

        return ComponentList;
    }

    public void parseIDForGene(BioProtein bionetProt) {
        if (bionetProt.getId().equals("")) {
            System.out.println("Enpty gene value for :" + bionetProt.getName());
        }

        String geneId = bionetProt.getId();
        if (geneId.startsWith("_")) {
            geneId = geneId.substring(1, geneId.length());
        }

        if (Pattern.compile("^.+_\\w$", 2).matcher(geneId).matches()) {
            geneId = geneId.substring(0, geneId.length() - 2);
        }

        geneId = geneId.replaceAll("_", ".");
        String[] ids;
        int var16;
        if (bionetProt.getEntityAnnot() != null) {
            HashSet<String> geneNames = bionetProt.getEntityAnnot().getEncodedBy();
            BioGene newGene;
            if (!geneNames.isEmpty()) {
                for(Iterator var4 = geneNames.iterator(); var4.hasNext(); bionetProt.addGene(newGene)) {
                    String name = (String)var4.next();
                    newGene = (BioGene)this.getBioNetwork().getGeneList().get(geneId);
                    if (newGene == null) {
                        newGene = new BioGene(geneId, name);
                        newGene.setSboterm("SBO:0000335");
                        this.getBioNetwork().addGene(newGene);
                    } else {
                        newGene.setName(name);
                    }
                }
            } else if (bionetProt.getId().contains("_and_")) {
                ids = bionetProt.getId().split("_and_");
                String[] var14 = ids;
                var16 = ids.length;

                for(int var7 = 0; var7 < var16; ++var7) {
                    String id = var14[var7];
                    newGene = (BioGene)this.getBioNetwork().getGeneList().get(id);
                    if (newGene == null) {
                        newGene = new BioGene(id, id);
                        newGene.setSboterm("SBO:0000335");
                        this.getBioNetwork().addGene(newGene);
                    }

                    bionetProt.addGene(newGene);
                }
            } else {
                newGene = (BioGene)this.getBioNetwork().getGeneList().get(geneId);
                if (newGene == null) {
                    newGene = new BioGene(geneId, geneId);
                    newGene.setSboterm("SBO:0000335");
                    this.getBioNetwork().addGene(newGene);
                }

                bionetProt.addGene(newGene);
            }
        } else if (bionetProt.getId().contains("_and_")) {
            ids = bionetProt.getId().split("_and_");
            ids = ids;
            int var15 = ids.length;

            for(var16 = 0; var16 < var15; ++var16) {
                String id = ids[var16];
                BioGene newGene = (BioGene)this.getBioNetwork().getGeneList().get(id);
                if (newGene == null) {
                    newGene = new BioGene(id, id);
                    newGene.setSboterm("SBO:0000335");
                    this.getBioNetwork().addGene(newGene);
                }

                bionetProt.addGene(newGene);
            }
        } else {
            BioGene newGene = (BioGene)this.getBioNetwork().getGeneList().get(geneId);
            if (newGene == null) {
                newGene = new BioGene(geneId, geneId);
                newGene.setSboterm("SBO:0000335");
                this.getBioNetwork().addGene(newGene);
            }

            bionetProt.addGene(newGene);
        }

    }

    public void parseReactionNotes(BioChemicalReaction bionetReaction) {
        String reactionNotes = bionetReaction.getEntityNotes().getXHTMLasString();
        reactionNotes = reactionNotes.replaceAll(">\\s+<", "><").replaceAll("[\\n\\r]", "");
        Matcher m = Pattern.compile(".*[> ]+SUBSYSTEM:\\s+([^<]+)<.*").matcher(reactionNotes);
        String[] GeneAssoList;
        String[] var6;
        int var7;
        int var8;
        String value;
        if (m.matches()) {
            GeneAssoList = m.group(1).split(StringUtils.escapeSpecialRegexChars(" || "));
            var6 = GeneAssoList;
            var7 = GeneAssoList.length;

            for(var8 = 0; var8 < var7; ++var8) {
                value = var6[var8];
                value = value.replaceAll("[^\\p{ASCII}]", "");
                if (this.getBioNetwork().getPathwayList().containsKey(value)) {
                    bionetReaction.addPathway((BioPathway)this.getBioNetwork().getPathwayList().get(value));
                } else {
                    BioPathway bionetPath = new BioPathway(value, value);
                    this.getBioNetwork().addPathway(bionetPath);
                    bionetReaction.addPathway(bionetPath);
                }
            }
        }

        m = Pattern.compile(".*[> ]+PROTEIN.CLASS:\\s+([^<]+)<.*").matcher(reactionNotes);
        Matcher m2 = Pattern.compile(".*[> ]+EC.Number:\\s+([^<]+)<.*").matcher(reactionNotes);
        value=null;
        if (m.matches()) {
            value = m.group(1);
            if (value == "") {
                value = "NA";
            }

            bionetReaction.setEcNumber(value);
        } else if (m2.matches()) {
            value = m2.group(1);
            if (value == "") {
                value = "NA";
            }

            bionetReaction.setEcNumber(value);
        } else {
            bionetReaction.setEcNumber("NA");
        }

        m = Pattern.compile(".*[> ]+SCORE:\\s+([^<]+)<.*").matcher(reactionNotes);
        if (m.matches()) {
            value = m.group(1);
            bionetReaction.setScore(value);
        }

        m = Pattern.compile(".*[> ]+STATUS:\\s+([^<]+)<.*").matcher(reactionNotes);
        if (m.matches()) {
            value = m.group(1);
            bionetReaction.setStatus(value);
        }

        m = Pattern.compile(".*[> ]+AUTHORS:\\s+([^<]+)<.*").matcher(reactionNotes);
        if (m.matches()) {
            GeneAssoList = m.group(1).split("(,\\s?)(?=PMID:)");
            var6 = GeneAssoList;
            var7 = GeneAssoList.length;

            for(var8 = 0; var8 < var7; ++var8) {
                value = var6[var8];
                m = Pattern.compile("PMID:\\s?(\\d++).*").matcher(value);
                if (m.matches()) {
                    bionetReaction.addPmid(m.group(1));
                }
            }
        }

        m = Pattern.compile(".*[> ]+NOTES:\\s+([^<]+)<.*").matcher(reactionNotes);
        m2 = Pattern.compile(".*[> ]+COMMENTS:\\s+([^<]+)<.*").matcher(reactionNotes);
        if (m.matches()) {
            value = m.group(1);
            if (value == "") {
                value = "NA";
            }

            bionetReaction.addComment(new Comment(value, "NA"));
        }

        if (m2.matches()) {
            value = m2.group(1);
            if (value == "") {
                value = "NA";
            }

            bionetReaction.addComment(new Comment(value, "NA"));
        }

        m = Pattern.compile(".*[> ]+GENE.{0,1}ASSOCIATION:\\s+([^<]+)<.*").matcher(reactionNotes);
        if (m.matches()) {
            GeneAssoList = m.group(1).replaceAll("[\\(\\)]", "").split("(?i) or ");
            if (GeneAssoList.length > bionetReaction.getEnzList().size()) {
                String regexCompartIds = "[" + this.getReactionCompart(bionetReaction) + "]";
                BioCompartment DefaultCompart;
                if (this.getBioNetwork().findbioCompartmentInList("fake_compartment") == null) {
                    DefaultCompart = new BioCompartment();
                    DefaultCompart.setAsFakeCompartment();
                    this.getBioNetwork().addCompartment(DefaultCompart);
                } else {
                    DefaultCompart = this.getBioNetwork().findbioCompartmentInList("fake_compartment");
                }

                String[] var31 = GeneAssoList;
                int var32 = GeneAssoList.length;

                label202:
                for(int var33 = 0; var33 < var32; ++var33) {
                    String GeneAsso = var31[var33];
                    if (!GeneAsso.contains("...")) {
                        String[] geneList = GeneAsso.replaceAll(" ", "").split("(?i)and");
                        String enzID_Regex = "";
                        String newEnzId = "";
                        HashMap<String, String> protIDs_REGEX = new HashMap();
                        String[] foundEnzyme = geneList;
                        int var17 = geneList.length;

                        String geneN;
                        String prot_REGEX;
                        for(int var18 = 0; var18 < var17; ++var18) {
                            geneN = foundEnzyme[var18];
                            newEnzId = newEnzId + "_" + geneN.replaceAll("\\.", "_");
                            prot_REGEX = "_" + geneN.replaceAll("\\.", "_") + "_" + regexCompartIds + "{0,1}";
                            enzID_Regex = enzID_Regex + prot_REGEX;
                            protIDs_REGEX.put(geneN, prot_REGEX);
                        }

                        Iterator var37;
                        if (geneList.length == 1) {
                            String geneName = geneList[0];
                            var37 = bionetReaction.getEnzList().values().iterator();

                            BioPhysicalEntity networkEnz;
                            do {
                                if (!var37.hasNext()) {
                                    var37 = this.getBioNetwork().getEnzList().values().iterator();

                                    do {
                                        if (!var37.hasNext()) {
                                            var37 = this.getBioNetwork().getProteinList().values().iterator();

                                            BioProtein prot;
                                            do {
                                                if (!var37.hasNext()) {
                                                    BioGene gene;
                                                    if (this.getBioNetwork().getGeneList().containsKey(geneName)) {
                                                        gene = (BioGene)this.getBioNetwork().getGeneList().get(geneName);
                                                    } else {
                                                        gene = new BioGene(geneName);
                                                        this.getBioNetwork().addGene(gene);
                                                    }

                                                    prot = new BioProtein("_" + geneName.replaceAll("\\.", "_").toUpperCase(), geneName + " (TH)");
                                                    prot.setBoundaryCondition(true);
                                                    prot.setConstant(false);
                                                    prot.setSboterm("SBO:0000252");
                                                    prot.setSubstanceUnits("false");
                                                    prot.setCompartment(DefaultCompart);
                                                    geneN = "<notes><body xmlns=\"http://www.w3.org/1999/xhtml\"><p> METEXPLORE_NOTE: This is a theoretical protein created through the Gene Association</p></body></notes>";
                                                    prot.setEntityNotes(new Notes(geneN));
                                                    prot.addGene(gene);
                                                    this.getBioNetwork().addProtein(prot);
                                                    bionetReaction.addEnz(prot);
                                                    this.getBioNetwork().getEnzList().put(prot.getId(), prot);
                                                    continue label202;
                                                }

                                                prot = (BioProtein)var37.next();
                                            } while(!Pattern.compile(enzID_Regex, 2).matcher(prot.getId()).matches());

                                            bionetReaction.addEnz(prot);
                                            this.getBioNetwork().getEnzList().put(prot.getId(), prot);
                                            continue label202;
                                        }

                                        networkEnz = (BioPhysicalEntity)var37.next();
                                    } while(!Pattern.compile(enzID_Regex, 2).matcher(networkEnz.getId()).matches());

                                    bionetReaction.addEnz(networkEnz);
                                    break;
                                }

                                networkEnz = (BioPhysicalEntity)var37.next();
                            } while(!Pattern.compile(enzID_Regex, 2).matcher(networkEnz.getId()).matches());
                        } else if (geneList.length > 1) {
                            foundEnzyme = null;
                            var37 = this.getBioNetwork().getComplexList().values().iterator();

                            int needed;
                            int found;
                            BioComplex existingCplx;
                            label193:
                            do {
                                label191:
                                while(true) {
                                    Collection setRegex;
                                    Collection compolist;
                                    do {
                                        if (!var37.hasNext()) {
                                            if (foundEnzyme != null) {
                                                continue label202;
                                            }

                                            BioComplex foundEnzyme2 = new BioComplex(newEnzId, newEnzId);
                                            foundEnzyme2.setBoundaryCondition(true);
                                            foundEnzyme2.setConstant(false);
                                            foundEnzyme2.setSboterm("SBO:0000297");
                                            foundEnzyme2.setSubstanceUnits("false");
                                            foundEnzyme2.setCompartment(DefaultCompart);
                                            this.getBioNetwork().addComplex(foundEnzyme2);
                                            var37 = protIDs_REGEX.entrySet().iterator();

                                            while(true) {
                                                label167:
                                                while(var37.hasNext()) {
                                                    Map.Entry<String, String> geneProtRegex = (Map.Entry)var37.next();
                                                    geneN = (String)geneProtRegex.getKey();
                                                    prot_REGEX = (String)geneProtRegex.getValue();
                                                    Iterator var45 = this.getBioNetwork().getProteinList().values().iterator();

                                                    BioProtein prot;
                                                    while(var45.hasNext()) {
                                                        prot = (BioProtein)var45.next();
                                                        if (Pattern.compile(prot_REGEX, 2).matcher(prot.getId()).matches()) {
                                                            foundEnzyme2.addComponent(new BioPhysicalEntityParticipant(newEnzId + "__With__" + prot.getId(), prot));
                                                            continue label167;
                                                        }
                                                    }

                                                    BioGene gene;
                                                    if (this.getBioNetwork().getGeneList().containsKey(geneN)) {
                                                        gene = (BioGene)this.getBioNetwork().getGeneList().get(geneN);
                                                    } else {
                                                        gene = new BioGene(geneN);
                                                        this.getBioNetwork().getGeneList().put(geneN, gene);
                                                    }

                                                    prot = new BioProtein("_" + geneN.replaceAll("\\.", "_").toUpperCase(), geneN + " (TH)");
                                                    prot.setBoundaryCondition(true);
                                                    prot.setConstant(false);
                                                    prot.setSboterm("SBO:0000252");
                                                    prot.setSubstanceUnits("false");
                                                    prot.setCompartment(DefaultCompart);
                                                    String protNotes = "<notes><body xmlns=\"http://www.w3.org/1999/xhtml\"><p> METEXPLORE_NOTE: This theoretical protein was created through the Gene Association and is a component of the protein complexe " + newEnzId + "</p></body></notes>";
                                                    prot.setEntityNotes(new Notes(protNotes));
                                                    prot.addGene(gene);
                                                    this.getBioNetwork().addProtein(prot);
                                                    foundEnzyme2.addComponent(new BioPhysicalEntityParticipant(newEnzId + "__With__" + prot.getId(), prot));
                                                }

                                                bionetReaction.addEnz(foundEnzyme2);
                                                this.getBioNetwork().getEnzList().put(foundEnzyme2.getId(), foundEnzyme2);
                                                continue label202;
                                            }
                                        }

                                        existingCplx = (BioComplex)var37.next();
                                        setRegex = protIDs_REGEX.values();
                                        compolist = existingCplx.getAllComponentList().values();
                                        needed = setRegex.size();
                                        found = 0;
                                    } while(setRegex.size() != compolist.size());

                                    Iterator var23 = compolist.iterator();

                                    while(true) {
                                        label187:
                                        while(true) {
                                            if (!var23.hasNext()) {
                                                continue label193;
                                            }

                                            BioPhysicalEntity compo = (BioPhysicalEntity)var23.next();
                                            boolean notFound = true;
                                            Iterator var26 = setRegex.iterator();

                                            while(var26.hasNext()) {
                                                String regex = (String)var26.next();
                                                if (Pattern.compile(regex, 2).matcher(compo.getId()).matches()) {
                                                    ++found;
                                                    continue label187;
                                                }
                                            }

                                            if (notFound) {
                                                continue label191;
                                            }
                                        }
                                    }
                                }
                            } while(found != needed);

                            if (!this.getBioNetwork().getEnzList().containsValue(existingCplx)) {
                                this.getBioNetwork().getEnzList().put(existingCplx.getId(), existingCplx);
                            }

                            if (!bionetReaction.getEnzList().containsValue(existingCplx)) {
                                bionetReaction.addEnz(existingCplx);
                            }
                        }
                    }
                }
            }
        }

    }

    private String getReactionCompart(BioChemicalReaction bionetReaction) {
        String reactionComparts = "";
        HashSet<String> ids = new HashSet();
        Iterator var4 = bionetReaction.leftParticipantList.values().iterator();

        BioPhysicalEntityParticipant parti;
        while(var4.hasNext()) {
            parti = (BioPhysicalEntityParticipant)var4.next();
            ids.add("(" + parti.getLocation().getId() + ")");
        }

        var4 = bionetReaction.rightParticipantList.values().iterator();

        while(var4.hasNext()) {
            parti = (BioPhysicalEntityParticipant)var4.next();
            ids.add("(" + parti.getLocation().getId() + ")");
        }

        String id;
        for(var4 = ids.iterator(); var4.hasNext(); reactionComparts = reactionComparts + id) {
            id = (String)var4.next();
        }

        return reactionComparts;
    }

    public void parseMetaboliteNote(BioPhysicalEntity bionetSpecies) {
        String metaboNotes = bionetSpecies.getEntityNotes().getXHTMLasString();
        metaboNotes = metaboNotes.replaceAll(">\\s+<", "><");
        Matcher m = Pattern.compile(".*FORMULA:\\s*([^<]+)<.*").matcher(metaboNotes);
        String value;
        if (m.matches()) {
            value = m.group(1);
            bionetSpecies.setChemicalFormula(value);
        }

        m = Pattern.compile(".*INCHI:\\s*([^<\\s]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            value = m.group(1);
            bionetSpecies.setInchi(value);
            if (!bionetSpecies.hasRef("inchi", value)) {
                bionetSpecies.addRef(new BioRef("SBML File", "inchi", value, 1));
            }
        }

        m = Pattern.compile(".*SMILES:\\s*([^<\\s]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            value = m.group(1);
            bionetSpecies.setSmiles(value);
        }

        m = Pattern.compile(".*INCHIKEY:\\s*([^<\\s]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            value = m.group(1);
            if (!bionetSpecies.hasRef("inchikey", value)) {
                bionetSpecies.addRef(new BioRef("SBML File", "inchikey", value, 1));
            }
        }

        m = Pattern.compile(".*KEGG.COMPOUND:\\s*([^<]+)<.*").matcher(metaboNotes);
        String[] var5;
        int var6;
        int var7;
        String[] ids;
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("kegg.compound", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "kegg.compound", value, 1));
                }
            }
        }

        m = Pattern.compile(".*HMDB:\\s*([^<]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("hmdb", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "hmdb", value, 1));
                }
            }
        }

        m = Pattern.compile(".*CHEBI:\\s*([^<]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("chebi", "CHEBI:" + value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "chebi", value, 1));
                }
            }
        }

        m = Pattern.compile(".*PUBCHEM.COMPOUND:\\s*([^<]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("pubchem.compound", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "pubchem.compound", value, 1));
                }
            }
        }

        m = Pattern.compile(".*PUBCHEM.SUBSTANCE:\\s*([^<]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("pubchem.substance", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "pubchem.substance", value, 1));
                }
            }
        }

        m = Pattern.compile(".*KEGG.GENES:\\s*([^<\\s]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("kegg.genes", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "kegg.genes", value, 1));
                }
            }
        }

        m = Pattern.compile(".*UNIPROT:\\s*([^<\\s]+)<.*").matcher(metaboNotes);
        if (m.matches()) {
            ids = m.group(1).split(" \\|\\| ");
            var5 = ids;
            var6 = ids.length;

            for(var7 = 0; var7 < var6; ++var7) {
                value = var5[var7];
                if (!bionetSpecies.hasRef("uniprot", value)) {
                    bionetSpecies.addRef(new BioRef("SBML File", "uniprot", value, 1));
                }
            }
        }

        if (bionetSpecies.getCharge() != null) {
            m = Pattern.compile(".*CHARGE:\\s*([^<]+)<.*").matcher(metaboNotes);
            if (m.matches()) {
                value = m.group(1);
                bionetSpecies.setCharge(value);
            }
        }

    }

    public void parseMetaboliteAnnot(BioPhysicalEntity bionetSpecies) {
        BioAnnotation annot = bionetSpecies.getEntityAnnot();
        if (annot.getXMLasString() != null) {
            String ref = annot.getXMLasString().replaceAll(">\\s*<", "><");
            String regex1 = ".*?<bqbiol:([^>]+)><rdf:Bag>(.*?)</rdf:Bag>.*";
            String regex2 = ".*?<rdf:li rdf:resource=\"http://identifiers.org/([^/]+)/([^\"]+)\"/>.*";

            Matcher m;
            String relation;
            for(m = Pattern.compile(regex1).matcher(ref); m.matches(); m = Pattern.compile(regex1).matcher(ref)) {
                relation = m.group(1);
                String links = m.group(2);

                for(Matcher m2 = Pattern.compile(regex2).matcher(links); m2.matches(); m2 = Pattern.compile(regex2).matcher(links)) {
                    String dbName = m2.group(1);
                    String dbId = m2.group(2);
                    if (!bionetSpecies.hasRef(dbName, dbId)) {
                        bionetSpecies.addRef(dbName, dbId, 1, relation, "SBML File");
                    }

                    if (dbName.equals("inchi")) {
                        bionetSpecies.setInchi(dbId);
                    }

                    links = links.replace(dbName + "/" + dbId, "");
                }

                ref = ref.replaceFirst("<bqbiol:([^>]+)><rdf:Bag>", "");
            }

            relation = ".*?<in:inchi xmlns:in=[^>]+>([^<]+)</in:inchi>.*";
            m = Pattern.compile(relation).matcher(ref);
            if (m.matches() && m.group(1).startsWith("InChI=")) {
                bionetSpecies.addRef("inchi", m.group(1), 1, "is", "SBML File");
                bionetSpecies.setInchi(m.group(1));
            }
        }

    }

    public void AddNote(BioEntity unusedSpecie, String note) {
        String ComplementaryNote = "<p>METEXPLORE_NOTE: " + note + "</p>\n</body>\n</notes>";
        String newNotes;
        if (unusedSpecie.getEntityNotes() != null && !unusedSpecie.getEntityNotes().getXHTMLasString().equals("")) {
            newNotes = unusedSpecie.getEntityNotes().getXHTMLasString();
            newNotes = newNotes.replaceAll("</body>\\s*</notes>", ComplementaryNote);
            unusedSpecie.getEntityNotes().setXHTMLasString(newNotes);
        } else {
            newNotes = "<notes>\n<body xmlns=\"http://www.w3.org/1999/xhtml\">" + ComplementaryNote;
            unusedSpecie.setEntityNotes(new Notes(newNotes));
        }

    }

    public void parseModifierNote(BioPhysicalEntity bionetSpecies) {
    }

    public void parseModifierNote(BioComplex bionetSpecies) {
    }

    public static String getCharacterDataFromElement(Element e) {
        Node child = e.getLastChild();
        if (child instanceof CharacterData) {
            CharacterData cd = (CharacterData)child;
            return cd.getData();
        } else {
            return "";
        }
    }

    private void setLogMode() throws FileNotFoundException, IOException {
        JSBMLUtils.setDefaultLog4jConfiguration();
        JSBMLUtils.setJSBMLlogToLevel(Level.ERROR);
        JSBMLUtils.addFileLog("/tmp/importLog.log");
    }

    private ArrayList<String> getErrorsFromLogs() throws IOException {
        ArrayList<String> errors = new ArrayList();
        File logFile = new File("/tmp/importLog.log");
        BufferedReader br = new BufferedReader(new FileReader(logFile));

        String sCurrentLine;
        while((sCurrentLine = br.readLine()) != null) {
            Matcher m = Pattern.compile(".* (ERROR .*)").matcher(sCurrentLine);
            if (m.matches()) {
                errors.add(m.group(1));
            }
        }

        br.close();
        logFile.delete();
        return errors;
    }

    public BioNetwork getBioNetwork() {
        return this.bioNetwork;
    }

    public void setBioNetwork(BioNetwork bioNetwork) {
        this.bioNetwork = bioNetwork;
    }

    public SBMLDocument getSBMLdoc() {
        return this.SBMLdoc;
    }

    public void setSBMLdoc(SBMLDocument sBMLdoc) {
        this.SBMLdoc = sBMLdoc;
    }

    public Set<String> getWarnings() {
        return this.Warnings;
    }

    public String getErrorThread() {
        return this.errorThread;
    }

    public void setWarnings(Set<String> warnings) {
        this.Warnings = warnings;
    }

    public void setErrorThread(String errorThread) {
        this.errorThread = errorThread;
    }

    public Model getjSBMLmodel() {
        return this.jSBMLmodel;
    }

    public void setjSBMLmodel(Model jSBMLmodel) {
        this.jSBMLmodel = jSBMLmodel;
    }

    public boolean isProccessing() {
        return this.processing;
    }

    public void setProccessing(boolean proccessing) {
        this.processing = proccessing;
    }
}

