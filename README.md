This serves as a detailed description of all scripts, packages, algorithms, extractions, and sanity checks associated with Historical Lightcurves of the R Coronae Borealis Stars (Gutierrez & Clayton 2026) [Will be linked] conducted during the Maria Mitchell Association Research Experience for Undergraduates (REU). The descriptions below are meant to help the reader understand each directory and/or package found in this repository for its future use.

[literature_And_Previous_Work] : This includes citations and papers of relevant pre-existing work relating to convolutional neural networks, photographic plates, R Coronae Borealis variables, and astrophysical time series analysis practices. This also includes an older version of DeepDisc titled [Astro RCNN](https://github.com/burke86/astro_rcnn).

[project_Notes] : My notes taken throughout the project, including meeting notes, example plots, etc.

[test_CNN] : The test run of the algorithm survey, performed on R CrB.
    - Within [test_CNN/notebooks], there exist [.ipynb] files (Python Notebooks). Those with a number at the start of their names are associated with the main pipeline of this project. Other notebooks are draft/concept notebooks used to develop the main pipeline. :
        [01_Load_Cutout.ipynb] : This notebook downloads the photographic plate cutouts for R CrB and includes a function that concatenates surrounding cutouts when more stars are necessary. This
        also fetches lim_mag_apass and lim_mag_atlas, the limiting magnitude of a plate, for each plate.
        [02_A_Source_Identification_And_Sorting_Algorithm.ipynb] : This
        [02_B_PyTorch_Algorithm] : 
        [02_C_DeepDisc_Algorithm] : 
    - Within [test_CNN/scripts], there exist [.py] script files that help the [.ipynb] files at [test_CNN/notebooks] run properly.
