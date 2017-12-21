----------

## High Resolution and Hierachical Chromosome Modelling 

----------

#### Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory, 
#### University of Missouri, Columbia MO 65211

----------

#### Developer: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tuan Trieu <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: tuantrieu@mail.missouri.edu <br/>

#### Contact: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Jianlin Cheng, PhD <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: chengji@missouri.edu <br/>





## 1. Content of folders:
- src: source code in java
- output: all experimental data
- executable: executable jar file

## 2. Usage ##

To run the tool, type: `java -jar HierarchicalModeller.jar chr_id resolution observed_contact_data normalized_contact_data domain_file output_folder`


- Parameters:
	+ **chr_id**: eg. 1, 2, ..
	+ **resolution**: e.g 5000
	+ **observed_contact_data**: observed hi-C contact file, each line contains 3 numbers (separated by a space) of a contact, position_1 position_2 interaction_frequencies (input/chr10_5kb.RAWobserved)(can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)	
	+ **normalized_contact_data**: normalized hi-C contact file, each line contains 3 numbers (separated by a space) of a contact, position_1 position_2 interaction_frequencies (chr10_5kb_gm12878_list.txt) (can be downloaded and normalized from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)
	+ **domain_file**: file contains domains identified by Juicer (input/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)	
	+ **output_folder**: output folder

- Typically, the input is several GBs in size and therefore, the program requires a lot of RAM memory to run. We ran our experiment in a server with 120 GB RAM and 80 cores. 

## 4. Disclaimer ##

The executable software and the source code of is distributed free of 
charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.


