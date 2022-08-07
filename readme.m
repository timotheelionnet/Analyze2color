This code was used to match single molecule trajectories in two colors in the following publication:
"Mapping translation 'hot-spots' in live cells by tracking single molecules of mRNA and ribosomes"
Zachary B Katz, Brian P English, Timothée Lionnet, Young J Yoon, Nilah Monnier, Ben Ovryn, Mark Bathe, Robert H Singer
Elife 2016 https://doi.org/10.7554/eLife.10415

requirements:
Matlab with the following toolboxes installed:
Curve Fitting Toolbox
Statistics and Machine Learning Toolbox

usage:
1 Set up the 'code' subfolder as Matlab's working directory.
2 > enter Analyze2color_diatrack3(); in Matlab command line.
3 set parameters in the following GUI:
  - Pixel size in um is the physical size of a camera pixel (trajectory spatial coordinates are expected in pixel units)
  - Time between frames (s) is self explanatory (trajectory time coordinates are expected in frame units)
  - check the diffusion analysis radio button to compute various
    metrics of diffusion on each trajectory (increases computation time)
  - Matched trajcetories minimum number of co-detections:
    required minimum number of frames where the two tracjetories are within a threshold distance
  - Maximum distance allowed between codected spots (pix):
    threshold distance used to call colocalized spots in above matching test.
  - Maximum distance allowed between trajectories (pix):
    maximum distance between the centroid of trajectories required to call a match.
  - Auto shift (recommended) will iteratively correct any systematic x,y offset between trajectories
    by computing the average offset between trajectories and subtracting that value
4 use the file pop-ups to enter in the following order:
  4.1 the file holding the trajectories for channel 1
  4.2 the file holding the trajectories for channel 2
  4.3 the output folder where the results will be saved.
the analysis should then proceed automatically and load various GUIs.

the 'out' subfolde contains the results of the anaysis of the
trajectories in the 'sample_trajectories' subfolder

*****
Diatrack trajectory file input format:
trajectory j occupies 3 columns: (j-1)*3 +1 , (j-1)*3 + 2, (j-1)*3 + 3
row 1 column (j-1)*3 +1: Index t1 of the first frame where the trajectory is detected
row 1 of the next two columns are set to zero
row 2 of the three columns contains the coordinates: (y(t1) x(t1) z(t1))
...
row n+1 of the three columns contains the coordinates: (y(tn) x(tn) z(tn)),
Zeros pad the end of each column once the trajectory has ended.
Beware! the number of tracks is limited by the width of the widest text file on your machine.
Rowwise export preferred.

example for 4 trajectories, the first three starting at frame 1, the fourth on frame 16.
Because this is 2D tracking all z coordinates are set to 1.
1	0	0	1	0	0	1	0	0	16	0	0
16.13	14.05	1	36.9	100.03	1	44.15	40.23	1	37.15	100.82	1
15.04	13.79	1	37.11	101.07	1	42.99	40.4	1	36.21	100.23	1
15.24	14.83	1	36.93	98.87	1	44.11	40.21	1	36.91	98.94	1
14.96	15.16	1	37.35	99.88	1	45.07	40.25	1	36.44	99.12	1
14.97	14.37	1	37.99	101.13	1	44.16	40.23	1	37.37	97.17	1
14.95	14.13	1	37.85	99.26	1	43.11	40.25	1	36.03	99.82	1
15.05	14.25	1	37.33	99.93	1	44.74	41.25	1	37.64	99.96	1
14.4	14.84	1	37.93	100.18	1	43.82	40.9	1	38.63	99.17	1
15.25	14.23	1	37.99	99.66	1	44.25	41.02	1	37.15	99.73	1
14.95	13.82	1	36.94	99.76	1	44.1	40.86	1	36.32	98.28	1
14.94	14.13	1	37.17	100.02	1	44.11	40.28	1	36.91	99.13	1
15.07	14.8	1	37.61	98.83	1	45.06	41.61	1	36.98	97.93	1
15.4	14.85	1	37.24	99.71	1	44.18	41.64	1	37.81	99.31	1
15.2	14.21	1	36.94	99.69	1	43.99	40.66	1	36.27	99.93	1
15.9	14.25	1	36.26	99.11	1	44.14	40.17	1	37.26	98.79	1
15.29	15.01	1	36.63	98	1	43.22	41.12	1	37.05	99.19	1
15.24	13.92	1	0	0	0	44.84	40.91	1	36.95	99.12	1
15.3	15.21	1	0	0	0	43.76	40.13	1	37	98.12	1
15.21	15.14	1	0	0	0	43.31	40.83	1	37.79	97.16	1
15.02	13.84	1	0	0	0	42.94	40.31	1	37.72	98.18	1
15.83	14.34	1	0	0	0	44.11	40.35	1	37.38	98.99	1
15.7	14.14	1	0	0	0	0	0	0	36.91	98.2	1
15.15	15.17	1	0	0	0	0	0	0	36.21	97.8	1
15.83	13.83	1	0	0	0	0	0	0	37.41	98.07	1
15.77	14.28	1	0	0	0	0	0	0	37.27	98.12	1
15.23	14.19	1	0	0	0	0	0	0	37	97.86	1
15.09	14.73	1	0	0	0	0	0	0	37.64	98.09	1
15.7	14.15	1	0	0	0	0	0	0	37.7	96.94	1
16.12	14.85	1	0	0	0	0	0	0	37.26	97.75	1
15.1	14.26	1	0	0	0	0	0	0	37.11	97.19	1
15.62	13.97	1	0	0	0	0	0	0	36.92	98.11	1
15.76	13.32	1	0	0	0	0	0	0	36.97	97.28	1

*****
Code written by Timothee Lionnet, HHMI Janelia Research Campus.

Example tracks supplied by Brian P. English, HHMI Janelia Research Campus.

Otsu code from:
-- Damien Garcia -- 2007/08, revised 2010/03 http://www.biomecardio.com/matlab/otsu.html
