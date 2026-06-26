# rcb_Star_Dust_Survey
Maria Mitchell Association REU Project, Dr. Geoffrey Clayton, Isaac Gutierrez

# Project Directions

- Started with R CrB photographic plates, based on work of /literature_And_Previous_Work/rnaaspdf.pdf.
- test_CNN, used the R CrB data as unit tests:
    - 01_Load_Cutout.ipynb - Loads available plates. To train Deepdisc (later step), I downloaded all available plates on the object due to many bad plates.
    - 02_Source_Identification_And_Sorting_By_Algorithm.ipynb - Uses three algorithms for plates downloaded successfully (with proper imaging solutions) to locate sources: 
        (1) find_peaks, 
        (2) DAOStarFinder
        (3) deepdisc (a deep learning algorithm, Merz et al. 2023). I was originally given astro_rcnn for (3) by mentor (found here: https://github.com/burke86/astro_rcnn/blob/master/README.md, Burke et al. 2019, MNRAS, 490 3952.) but found that it is built on the old Matterport Mask R-CNN framework which uses TensorFlow 1.x/Keras. It's quite dated and getting it running can be painful. The authors of astro_rcnn made an updated version called DeepDISC (burke86/deepdisc) which is built on the modern detectron2 framework and is much easier to install. (3) Requires a separate environment, which is setup and implemented through 02_DeepDisc_Algorithm.ipynb.
        (4) Pytorch. Another deep learning algorithm, as I am trying to compare this to (3) and see which is most efficient/practical to start with.

    - 03_Photometric_Survey.ipynb - 