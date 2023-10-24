# fesom-c
FESOM-C model

The FESOM-C model is a coastal branch of the global Finite volumE sea ice–Ocean Model (FESOM2) (Danilov et al., 2017). In addition to the partially common interface of the models, FESOM-C has speciﬁc features that are important for our work. The model was originally developed for applications with a high horizontal resolution as ﬁne as several meters (Neder et al., 2022; Kuznetsov et al., 2020; Fofonova et al., 2019). This model uses the discretization of cells and vertices of a ﬁnite volume, which allows the use of unstructured computational grids. We use this function to move the boundary of the computational domain away from the region of interest, the "core" of the model grid, without creating a system of nested grids. At the same time, the horizontal resolution outside the core is quite coarse, which allows us to reduce the inﬂuence of the boundary on our solution inside the core. The most important distinguishing feature from the global FESOM2 model is the possibility of using hybrid grids consisting of triangles and  squares. The effectiveness of this approach in enhancing stability and using larger time steps is shown by Danilov and Androsov (2015) and Androsov et al. (2019).
