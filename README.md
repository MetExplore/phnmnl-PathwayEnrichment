![LOGO](Logo_Metexplore.png)
# PathwayEnrichement

Version: 1.1

## Short description
Predicts pathway enrichment into a (human) metabolic network (Recon v2.02)

## Description
Metabolites (or other bio-entities, e.g., reaction, genes) belonging to the fingerprint (input file) are mapped on Recon v2.02 network (Thiele et al., 2013). For that purpose, ID information are used among: name of the entities, SBML ID, ChEBI, InChI, InChIKey, SMILES, CSID, PubChem and HMDB). Pathway enrichment is calculated with an exact Fisher one-tailed test corrected by Bonferroni and Benjamini Hochberg methods. This tool is part of the MetExplore's project consisting in a web server dedicated to the analysis of omics data in the context of genome scale metabolic networks (Cottret et al., 2010).

### Input files
- a ```fingerprint``` (tsv or tabular format): required, composed by at least a column containing ID values to map on the network file. Multi-mapping (i.e., mapping on different ID) could be performed if these three kind of values are included each in a separate column. Optionally, this program could filter empty values from a designated column (e.g., non-significant bio-entities).
- a ```metabolic network``` (SBML) : optional, by default Recon v2.02 SBML file (without compartments) (Thiele et al., 2013).

### Output files
- ```mapping.tsv```: each line corresponds to the bio-entities from the dataset file: the success or the failure of the mapping, their names in the dataset, those of one or several elements of the network in case of matching.
- ```pathwayEnrichment.tsv```: contains for each pathway associated with the mapped bio-entity: their names, the Fisher's p value of enrichment, the Bonferroni and the Benjamini-Hochberg corrections, the list of the mapped bio-entities and their corresponding ID in the SBML, the number of mapping and its coverage on the total of bio-entities contained in the studied pathway.
- ```information.txt``` (if -gal or --galaxy parameter is activated): contains general information about mapping and pathway enrichment results. This file contains the total number of mapped bio-entities with the coverage in the dataset and in the network, the total number of enriched pathways and the coverage in the network. Eventually, warnings alert about doublets in mapping which are discarded from the pathway analysis. In this case, the user must choose the corresponded bio-entity's ID in the network in order to add them in a new column of the fingerprint dataset. Then, the program must be relaunched by using the SBML ID mapping only (-id <columnNumber>).

## Key features
- Metabolic network
- Modeling
- Pathway analysis
- Prediction

## Functionality
- Post-processing
- Statistical Analysis

## Tool Authors
- MetExplore Group contact-metexplore@inra.fr

## Git Repository
- https://github.com/phnmnl/container-pathwayEnrichment.git

## Installation
For local installation of the container:
```
VERSION=1.1
git clone --depth 1 --single-branch --branch $VERSION https://github.com/MetExplore/phnmnl-PathwayEnrichment.git pathwayEnrichment
cd pathwayEnrichment/
git checkout $VERSION
mvn install:install-file install:install-file -Dfile=parseBioNet.jar -DgroupId=fr.inra.toulouse.metexplore -DartifactId=parseBioNet -Dversion=0.0.1 -Dpackaging=jar && \
mvn install
mv target/PathwayEnrichment-$VERSION-jar-with-dependencies.jar pathwayEnrichment.jar
rm -rf target
```

## Usage Instructions
Test the program:

```
java -jar pathwayEnrichment.jar -i  data/recon2.02_without_compartment.xml
```

For direct usage:

```
java -jar pathwayenrichment.jar -i  <input_file> [-s <sbml_file>] [-o1 <mapping_output_file>] [-o2 <pathway_enrichment_output_file>] [-f <filtered_column>] [-sep <separator>] [--header] [-t <bio_entity>] [-name <name_column>] [-idSBML <SBML_ID_column>] [-chebi <ChEBI_ID_column>] [-kegg <KEGG_ID_column>], [-hmdb <HMDB_ID_column>] [-csid <ChemSpider_ID_column>] [-inchikey <InChIKey_column>] [-inchi <InChI_column>] [-l <layer_selection>] [-h]
```

- ```-h (--help)``` for printing the help.

##### Files parameters
- ```-i (--inFile)```, ```-s (--sbml)``` (STRING) are used to specify the inputs files. Only ```-i``` - corresponding to the dataset of fingerprint - is required. ```-s``` - the sbml file where the bionetwork is extracted - used Recon2.02 network by default (Thiele et al., 2013).
- ```-o1 (--outMap)``` and ```-o2 (--outPath)``` (STRING) could be used to specify the output file. They contains respectively results from mapping (default: "mapping.tsv") and pathway enrichment (default: "pathwayEnrichment.tsv") 

##### Parsing parameters
- ```-f (--filter)``` (NUMERICAL) the specified column is used to discard from the analysis lines containing empty values (for example : after a statistical pre-selection among the fingerprint). 
- ```--header``` could be used to discard first line of the file to the analysis.
- ```-sep (--separator)``` (STRING) specify the character used to separate the column in the fingerprint dataset (e.g., ";"). By default, the program uses "\\t" for tabulation separators.

##### Mapping parameters
- ```-name```, ```-idSBML```, ```-inchi```, ```-inchikey```, ```-chebi```,```-kegg```, ```-hmdb``` and ```-csid (--chemspider)``` (NUMERICAL) options point out specific column numbers in the fingerprint file. ```-name``` indicates the name of the bio-entities. ```-idSBML``` corresponds to the ID of the bio-entities in the SBML. At least one mapping parameter should be activated; multiple parameters could be activated together. By default, ```-name``` and ```-id``` are set respectively to 1 and 2 ("1" is considered as the first column); the other parameters are disabled.
- ```-l``` (STRING) selects the InChI's layers used to map metabolites. Nine layers could be selected: connections (c), hydrogens atoms (h), charge (q), protons (p), double-bond stereo (b), tetrahedral sp3 stereo (t), isotopic atoms (i), fixed hydrogens (f) and reconnected layers (r). These layers must be set as a list of characters containing the letter of each selected layer separated by commas. By default, this option is set on ```c,h``` for a mapping only on connections and hydrogens atom layers. A whole layers selected must be ```c,h,q,p,b,t,i,f,r```. For mapping on formula only, select ```-l``` option with no parameters after.
- ```-t``` (NUMERICAL) specify the bio-entities used for pathway analysis: 1 for metabolites (by default), 2 for reactions. 


## References
- Thiele I., Swainston N., Fleming R.M.T., et al. A community-driven global reconstruction of human metabolism (2013). Nature biotechnology 31(5):10. doi:10.1038/nbt.2488.
- Cottret L. , Wildridge D., Vinson F., Barrett M.P., Charles H., Sagot M.-F. and Jourdan F. MetExplore: a web server to link metabolomic experiments and genome-scale metabolic networks. (2010) Nucleic Acids Research 1:38 Suppl:W132-7. doi:10.1093/nar/gkq312.