## Code to preprocess and load data

Code used to build the AD Atlas. To refill the database BUILD SCRIPTS #1 and #2 need to be executed. For more infromation please see the [documentation](https://sysmet.pages.hzdr.de/adatlas-backend/.).

    .
    ├── preprocessing                   # Code to preprocess and load data
    │   ├── BiGG                        # Parse and format BiGG models downloaded from www.vmh.life
    │   ├── GTEx                        # Download and format GTEx data
    │   ├── data-format-scripts         # Scripts used to format raw data
    │   ├── gene-hgnc-ensembl           # Download and format gene information
    │   ├── misc-scripts                # Misc scripts - need to be sorted
    ├── setup                           # Scripts to load data into neo4j      [BUILD SCRIPTS #1]
    ├── setup_ADatlas                   # Scripts to build/summarize AD Atlas  [BUILD SCRIPTS #2]
    ├── quality_control.R               # <= needs to be adapted        
    └── README.md

> Data files and temporary files used to build the atlas are NOT included in this repository as they may include sentsitive data. 
> For access to the data please contact the owner of this repository. 
