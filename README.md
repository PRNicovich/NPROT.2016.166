#NPROT.2016.166

Code associated with publication 

"Turning single-molecule localization microscopy into a quantitative bioanalytical tool" <br>
Philip R Nicovich, Dylan M Owen & Katharina Gaus<br>
*Nature Protocols* volume 12, pages 453–460 (2017)<br>

* Figure 1 panels created with NatProtTimelineFigure.m.
* Figure 2 panels created with NatProtFigure.m.
* Figure 3 panels created with DrawNatProtKineticsFigure.m.

Each can generate new SMLM data for illustration.  Data generated with GenerateNatProtData.m.  Option in code to share datasets by defining previously-saved file for this purpose. 

NatProtFigure.m calls examples of nearest-neighbor, Ripley's K, pair correlation, DBSCAN, Delaunay, and Voronoi clustering analysis all on the same dataset.  See blocks of this code for examples of implementing each of these analyses.

DrawNatProtKineticsFigure.m generates single spatial cluster of points which is then simulated as pentamer of underlying points.  Includes Hidden Markov model simulation of fluorophore kinetics.