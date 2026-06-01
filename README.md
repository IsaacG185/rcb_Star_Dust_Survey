# rcb_Star_Dust_Survey
Maria Mitchell Association REU Project, Dr. Geoffrey Clayton, Isaac Gutierrez

# Project Directions

- Started with V854 Cen photographic plates, based on work of /literature_And_Previous_Work/rnaaspdf.pdf.
- test_CNN, used the V854 Cen data as unit tests:
    - 01_Load_Cutout.ipynb - Loads available plates.
    - 02_Source_Identification_And_Sorting_By_Algorithm.ipynb - Uses three algorithms for plates downloaded successfully (with proper imaging solutions): (1) find_peaks, (2) DAOStarFinder, (3) deepdisc [Deep Learning Algorithm]. I was originally given astro_rcnn for (3) by mentor (found here: https://github.com/burke86/astro_rcnn/blob/master/README.md) but found that it is built on the old Matterport Mask R-CNN framework which uses TensorFlow 1.x/Keras. It's quite dated and getting it running can be painful. The authors of astro_rcnn made an updated version called DeepDISC (burke86/deepdisc) which is built on the modern detectron2 framework and is much easier to install. (3) Requires a separate environment, which is setup and implemented through 02_DeepDisc_Algorithm.ipynb.