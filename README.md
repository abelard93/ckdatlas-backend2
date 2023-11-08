## Git repository
All the code present in this repository comes from Project ADatals (https://adatlas.org). There are just a few small modifications.

This repository contains code for a two-stage process in building the atlas. The first stage involves loading all the data into the Neo4j database. In the second stage, data is abstracted and summarized using a combination of Neo4j and R scripts to generate the final CKD Atlas data model.


### 

    .
    ├── code                   # Code to preprocess and load data
    ├── docs                   # Documentation files 
    ├── misc                   # Miscellaneous files - need to be sorted
    ├── schema                 # Html representation of current data model 
    ├── .gitignore
    ├── .gitlab-ci.yml         
    └── README.md

> Data files and temporary files used to build the atlas are NOT included in this repository as they may include sentsitive data. 
> For access to the data please contact the owner of this repository. 

## Note
For more information on the data model, please visit the ADatlas repository:
[https://github.com/compneurobio/adatlas-backend/](https://sysmet.pages.hzdr.de/adatlas-backend/).

<!--
## Documentation

For a more comprehensive documentation of the AD Atlas backend including a guide for developers and file format specification please visit [https://sysmet.pages.hzdr.de/adatlas-backend/](https://sysmet.pages.hzdr.de/adatlas-backend/).


## Data model 

General           |  AD Atlas 
:-------------------------:|:-------------------------:
![Data_model](misc/data_model/data_model_current.png?raw=true) | ![Data_model](misc/data_model/data_model_adatlas_current.png?raw=true)

> Click image to view details. 
-->
