# Supporting data and code for the paper:

"Siderophores as an iron source for picocyanobacteria in deep chlorophyll maximum layers of the oligotrophic ocean"

[Published version available Open Access from the ISME Journal](https://doi.org/10.1038/s41396-022-01215-w)

[Preprint available from bioRxiv](https://doi.org/10.1101/2021.11.13.468467)

Data and code here is provided under the MIT License. Feel free to use or remix as you see fit.

## NOTE ON REPRODUCIBILITY
The [GEOTRACES](https://www.bodc.ac.uk/geotraces/data/idp2017/) data license is not GPL3 and you must request access and be approved before being allowed to download the data. Therefore, I can not include the GEOTRACES data here since it would be in violation of the license agreement. Unfortunately, this means that the workflows in the code directory are not truly reproducible since they do not include the GEOTRACES data. You can fill in the missing data for yourself by downloading the necessary metadata and placing it in the `01_env_data` directory. You will then simply need to replace the completed data files with the onces provided in this repository here which contain NA values for all GEOTRACES chemical measurements. However, all measurements for Prochlorococcus and SAR11 ecotype abundance are included.

The MARMICRODB reference database version used in this analysis and the instructions for its use are available from https://dx.doi.org/10.5281/zenodo.3520509.

## ATTRIBUTION
If you use any of the data mentioned or described in this repository it is ESSENTIAL!! that you properly cite the data sources. These include (but are not limited to):

Biller, S. J., P. M. Berube, K. Dooley, and others. 2018. Marine microbial metagenomes sampled across space and time. [Scientific Data 5: 180176.](https://www.nature.com/articles/sdata2018176)

Schlitzer, R., R. F. Anderson, E. M. Dodas, and others. 2018. The GEOTRACES Intermediate Data Product 2017. [Chem. Geol. 493: 210–223.](https://doi.org/10.1016/j.chemgeo.2018.05.040)

The license for GEOTRACES data does not allow it to be distributed outside of the original. You are encouraged to also cite specific papers describing individual measurement types and/or reach out to those folks to include them as collaborators.

> GEOTRACES Intermediate Data Product (IDP) Download Agreement
> The GEOTRACES programme is keen to ensure that the very significant effort and expertise involved in making trace-element and isotope measurements is acknowledged as fully as possibly in subsequent publications.
>
>Users of the GEOTRACES Intermediate Data Product are expected to abide to the following rules regarding citation
>
>To the greatest extent possible, please cite all relevant publications from researchers that made the measurements you use in your work. Details of publications that should be cited are provided point-by-point in the IDP dataset (in the ODV and ASCII versions) and will be updated on the online database as new papers are published. Where your research deals particularly with data measured by a single group of data originators, you are invited to please contact that group to discuss your work prior to publication and jointly consider the synergy and mutual benefit of co-authorship where appropriate.
>
>Where other constraints prevent citation of all relevant publications, for instance where there is a journal limitation on the maximum number of publications that can be cited, or if the dataset is only used in a minor supportive way, please cite the data compilation itself (as below). In such cases, also please cite any original individual papers that you rely on particularly heavily for your research and interpretations.
>
>Where using data from the IDP2017 product in your publication, please cite the data compilation as: Schlitzer, R., Anderson, R. F., Masferrer Dodas, E, et al., The GEOTRACES Intermediate Data Product 2017, Chem. Geol. (2018), https://doi.org/10.1016/j.chemgeo.2018.05.040 .
>
>Where using data from the IDP2014 product in your publication, please cite the data compilation as: Mawji, E., et al., The GEOTRACES Intermediate Data Product 2014, Mar. Chem. (2015), http://dx.doi.org/10.1016/j.marchem.2015.04.005 .
>
>Users of the GEOTRACES Intermediate Data Product shall not distribute downloaded data to any third party.

## WORKFLOW DESCRIPTION
- `01_env_data` contains the workflow for matching GEOTRACES and Tara Oceans chemical, biological, and hydrographic data to metagenome sample identifiers. 
  - Also calculates the Deep Chloropyll Maximum depth at each geographic location and determines whether one of the metagenome samples overlaps with the DCM.
  - To match internal sample identifiers (in form `S0031`) to actual NCBI sequence read archive identifiers you need the files found [here](https://static-content.springer.com/esm/art%3A10.1038%2Fsdata.2018.176/MediaObjects/41597_2018_BFsdata2018176_MOESM325_ESM.zip)

- `01_tbdt_sequence_similarity` contains code for reproducing Supplementary Fig. S2 and for looking at sequence similarity networks of the tonB dependent transporter in marine cyanobacteria.

- `03_marker_tbdt_counts` contains code for formatting counts of Prochlorococcus marker genes and also tonB dependent transporter counts and for calculating relative abundance of TBDTs.

- `04_sag_genome_biogeography` contains code for reproducing Fig. 1 from the main text

- `05_hot_bats` contains code for reproducing Fig. 2, Supplementary Fig. S3, and Supplmentary Table S1.

- `06_global_metagenome_biogeography` contains code for reproducing Fig. 3.

- `07_global_metagenome_models` contains code for reproducing Fig. 4, Fig, 5, Fig. 6, Supplementary Fig. S4, and Supplementary Table S2.

