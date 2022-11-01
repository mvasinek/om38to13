# om38to13
Tool for annotation and filtering of optical mapping structural variants. The tool detects structural variants detected due to misassembly of GRCh38 to T2T-CHM13

# Python scripts
* om38to13.py - simple console based application.
* To run: 
  - Annotation - provides summary over all variants in input SMAP file.
    `python om38to13.py annotate variants.smap`
  - Filter - the operation evaluates whether there are induced structural variants and removes these variants from the SMAP file.
    `python3 om38to13.py filter variants.smap -o variants.filtered.smap`
  - View - operation searches for the input interval and prints all the information to the console.
    `python3 view chr1:1000000-2000000`

# Web app
* Offers simple user interface accessible online.
* http://app.olgen.cz/om38to13/