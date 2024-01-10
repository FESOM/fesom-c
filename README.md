# FESOM-C

The FESOM-C, the coastal branch of the Finite-volumE Sea ice – Ocean Model (FESOM2). This model retains many numerical characteristics of FESOM2, notably its finite-volume cell-vertex discretization. However, FESOM-C's dynamical core is distinct, particularly in its approach to time stepping, adoption of a terrain-following vertical coordinate, and the development of hybrid meshes that integrate triangles and quads. These first two features were essential in establishing FESOM-C as a separate branch. The implementation of hybrid mesh technology enhances computational efficiency. This is because quadrilateral cells, which have fewer edges than their triangular counterparts, are less prone to the spurious inertial modes associated with triangular cell-vertex discretization and require reduced dissipation. The hybrid mesh design enables the use of quasi-quadrilateral unstructured meshes, where triangular cells are employed solely to connect quadrilateral sections of varying resolutions or in places where quadrilateral cells are excessively distorted.

# References

Androsov, A., Fofonova, V., Kuznetsov, I., Danilov, S., Rakowsky, N., Harig, S., Brix, H., and Helen Wiltshire, K.: FESOM-C v.2: Coastal dynamics on hybrid unstructured meshes, Geoscientific Model Development, 12, 1009–1028, https://doi.org/10.5194/gmd-12-1009-2019, 2019.

Danilov, S. and Androsov, A.: Cell-vertex discretization of shallow water equations on mixed unstructured meshes, Ocean Dynamics, 65, 33–47, https://doi.org/10.1007/s10236-014-0790-x, 2015.

Fofonova, V., Androsov, A., Sander, L., Kuznetsov, I., Amorim, F., Hass, H. C., and Wiltshire, K. H.: Non-linear aspects of the tidal dynamics in the Sylt-Rømø Bight, south-eastern North Sea, Ocean Science, 15, 1761–1782, https://doi.org/10.5194/os-15-1761-2019, 2019.

Kuznetsov, I., Androsov, A., Fofonova, V., Danilov, S., Rakowsky, N., Harig, S., and Wiltshire, K. H.: Evaluation and Application of Newly Designed Finite Volume Coastal Model FESOM-C, Effect of Variable Resolution in the Southeastern North Sea, Water, 12, https://doi.org/10.3390/w12051412, 2020.

Nerger, L., Tang, Q., and Mu, L.: Efficient ensemble data assimilation for coupled models with the Parallel Data Assimilation Framework: example of AWI-CM (AWI-CM-PDAF 1.0), Geoscientific Model Development, 13, 4305–4321, https://doi.org/10.5194/gmd-13-4305-2020, 2020.

Rubinetti, S.; Kuznetsov, I.; Fofonova, V.; Androsov, A.; Gnesotto, M.; Rubino, A.; Zanchettin, D. Water Mass Transport Changes through the Venice Lagoon Inlets from Projected Sea-Level Changes under a Climate Warming Scenario. Water, 15, 3221. https://doi.org/10.3390/w15183221 , 2023.

Pham, H.V. Dal Barco, M.K. Cadau, M. Harris, R. Furlan, E. Torresan, S. Rubinetti, S. Zanchettin, D. Rubino, A. Kuznetsov, I., Barbariol, F., Benetazzo, A., Sclavo, M., Critto, A. Multi-Model Chain for Climate Change Scenario Analysis to Support Coastal Erosion and Water Quality Risk Management for the Metropolitan City of Venice. Sci. Total Environ., 904, 166310. https://doi.org/10.1016/j.scitotenv.2023.166310 , 2023.

Kuznetsov, I., Rabe, B., Androsov, A., Fang, Y.-C., Hoppmann, M., Quintanilla-Zurita, A., Harig, S., Tippenhauer, S., Schulz, K., Mohrholz, V., Fer, I., Fofonova, V., and Janout, M.: Dynamical reconstruction of the upper-ocean state in the Central Arctic during the winter period of the MOSAiC Expedition, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1353, 2023. 
