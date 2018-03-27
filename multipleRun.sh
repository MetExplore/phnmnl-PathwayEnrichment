#!/bin/bash

VERSION=1.1
#TESTS=('')
TESTS=('-t 2' '-t 4' '-t 5' ' -t 6 -idSBML 3' '-t 5 -tEnr 6' '-t 4 -tEnr 6' '-t 6 -idSBML 3 -tEnr 5' '-t 1' '' '-tEnr 2' '-f 13' '-idSBML 2' '-chebi 3' '-inchi 4' '-inchi 4 -l c,h,t' '-inchi 4 -l t,h,c' '-l c,h,t' '-inchi 4 -l' '-l' '-inchikey 5' '-kegg 6' '-pubchem 7' '-hmdb 8' '-csid 9' '-mass 12 -prec 2' '-prec 2' '-name 1 -nameCol 2' '-nameCol -1 -name 1' '-s data/recon2.02_without_compartment.xml' '-header' '-sep \t' '-sepID ;' '-inchi 4 -l c,h,t -lWarn -o1 temp/checking_format.tsv' '-inchi 4 -o1 temp/checking_format.tsv' '-noCheck -o1 temp/checking_format2.tsv')
echo '' > resultRuns.log
NB_RUN=0;

#TODO: verifier les output similaires ex: -t 1 avec -s without comp

#'without inchi layers'
#TESTS=('-o1 mappingX.tsv' '-o2 mappingY.tsv')

#smiles

#cp target/PathwayEnrichment-$VERSION-jar-with-dependencies.jar pathwayEnrichment.jar
cd temp
mkdir map enr info
cd ..

run(){
	printf "\n\n$NB_RUN. ${TESTS[NB_RUN]}\n" >> resultRuns.log 2>&1
	java -jar pathwayEnrichment.jar -o2 temp/map/mapping.tsv$NB_RUN -o3 temp/enr/pathwayEnrichment.tsv$NB_RUN  -gal temp/info/information.txt$NB_RUN $@ >> resultRuns.log 2>&1
	let NB_RUN+=1
}

for i in `seq 0 $((${#TESTS[@]} -1))`; do
	if [ $i -eq 0 ]; then
		run -i /home/user/bin/phnmnl-PathwayEnrichment/data/reactions_recon2.02.tsv ${TESTS[i]}
	elif [[ $i -ge 1 ]] && [[ $i -le 6 ]]; then
		run -i /home/user/bin/phnmnl-PathwayEnrichment/data/gpr_recon2.02.tsv -s data/recon2.02.xml ${TESTS[i]}
	else
		run -i /home/user/bin/phnmnl-PathwayEnrichment/data/sacurine_workflow_output.tsv ${TESTS[i]}
	fi
done
