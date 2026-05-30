# rcb_Star_Dust_Survey
Maria Mitchell Association REU Project, Dr. Geoffrey Clayton, Isaac Gutierrez

# Project Directions

- Started with V854 Cen photographic plates, based on work of /literature_And_Previous_Work/rnaaspdf.pdf.
- test_CNN, used the V854 Cen data as unit tests:
    - 01_Load_Cutout.ipynb - Loads available plates.
    - 02_Source_Identification_And_Sorting_By_Algorithm.ipynb - Uses three algorithms for plates downloaded successfully (with proper imaging solutions): (1) find_peaks, (2) DAOStarFinder, (3) astro_rcnn [Deep Learning Algorithm].