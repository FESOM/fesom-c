# FESOM-C

The FESOM-C, the coastal branch of the Finite-volumE Sea ice – Ocean Model (FESOM2). FESOM-C solves 3D primitive equations in the Boussinesq, hydrostatic, and traditional approximations for the momentum, continuity, and density constituents. This model retains many numerical characteristics of FESOM2, notably its finite-volume cell-vertex discretization. However, FESOM-C's dynamical core is distinct, particularly in its approach to time stepping, adoption of a terrain-following vertical coordinate, the inclusion of a wetting/drying option, and the usage of hybrid meshes that integrate triangles and quads and allow high computational efficiency and geometrical flexibility. The numerical scheme of FESOM-C splits the fast and slow motions into barotropic and baroclinic modes. It uses an asynchronous time stepping, assuming that the integration of temperature and salinity is half-step-shifted with respect to momentum. For modeling wetting and drying we use the method proposed by Stelling and Duinmeijer, the algorithm to account for wetting and drying (if an appropriate flag is on) is performed on each time step of the external mode. Three schemes have been implemented in FESOM-C for horizontal advection: upwind, Miura, and MUSCL-type.
## 
 **Website:** [FESOM-C on fesom.de](https://fesom.de/models/fesom-c/)
## Contacts

- Vera Sidorenko: `Vera.Sidorenko@awi.de`
- Alexey Androsov: `Alexey.Androsov@awi.de`
- Ivan Kuznetsov: `Ivan.Kuznetsov@awi.de`


# References

Androsov, A., Fofonova, V., Kuznetsov, I., Danilov, S., Rakowsky, N., Harig, S., Brix, H., and Helen Wiltshire, K.: FESOM-C v.2: Coastal dynamics on hybrid unstructured meshes, Geoscientific Model Development, 12, 1009–1028, https://doi.org/10.5194/gmd-12-1009-2019, 2019.

Danilov, S. and Androsov, A.: Cell-vertex discretization of shallow water equations on mixed unstructured meshes, Ocean Dynamics, 65, 33–47, https://doi.org/10.1007/s10236-014-0790-x, 2015.

Fofonova, V., Androsov, A., Sander, L., Kuznetsov, I., Amorim, F., Hass, H. C., and Wiltshire, K. H.: Non-linear aspects of the tidal dynamics in the Sylt-Rømø Bight, south-eastern North Sea, Ocean Science, 15, 1761–1782, https://doi.org/10.5194/os-15-1761-2019, 2019.

Kuznetsov, I., Androsov, A., Fofonova, V., Danilov, S., Rakowsky, N., Harig, S., and Wiltshire, K. H.: Evaluation and Application of Newly Designed Finite Volume Coastal Model FESOM-C, Effect of Variable Resolution in the Southeastern North Sea, Water, 12, https://doi.org/10.3390/w12051412, 2020.

Nerger, L., Tang, Q., and Mu, L.: Efficient ensemble data assimilation for coupled models with the Parallel Data Assimilation Framework: example of AWI-CM (AWI-CM-PDAF 1.0), Geoscientific Model Development, 13, 4305–4321, https://doi.org/10.5194/gmd-13-4305-2020, 2020.

Fofonova​​​​​​​, V., Kärnä, T., Klingbeil, K., Androsov, A., Kuznetsov, I., Sidorenko, D., Danilov, S., Burchard, H., and Wiltshire, K. H.: Plume spreading test case for coastal ocean models, Geosci. Model Dev., 14, 6945–6975, https://doi.org/10.5194/gmd-14-6945-2021, 2021.

Galvez, D.S.; Papenmeier, S.; Sander, L.; Hass, H.C.; Fofonova, V.; Bartholomä, A.; Wiltshire, K.H. Ensemble Mapping and Change Analysis of the Seafloor Sediment Distribution in the Sylt Outer Reef, German North Sea from 2016 to 2018. Water 2021, 13, 2254. https://doi.org/10.3390/w13162254

Neder, C., Fofonova, V., Androsov, A., Kuznetsov, I., Abele, D., Falk, U., et al.. Modelling suspended particulate matter dynamics at an Antarctic fjord impacted by glacier melt. Journal of Marine Systems 2022, 231, 103734. https://doi.org/10.1016/j.jmarsys.2022.103734.

Danilov, S., Mehlmann, C., Fofonova, V. On discretizing sea-ice dynamics on triangular meshes using vertex, cell or edge velocities, Ocean Modelling, 170, 2022, 101937, https://doi.org/10.1016/j.ocemod.2021.101937.

Rubinetti, S.; Kuznetsov, I.; Fofonova, V.; Androsov, A.; Gnesotto, M.; Rubino, A.; Zanchettin, D. Water Mass Transport Changes through the Venice Lagoon Inlets from Projected Sea-Level Changes under a Climate Warming Scenario. Water, 15, 3221. https://doi.org/10.3390/w15183221 , 2023.

Pham, H.V. Dal Barco, M.K. Cadau, M. Harris, R. Furlan, E. Torresan, S. Rubinetti, S. Zanchettin, D. Rubino, A. Kuznetsov, I., Barbariol, F., Benetazzo, A., Sclavo, M., Critto, A. Multi-Model Chain for Climate Change Scenario Analysis to Support Coastal Erosion and Water Quality Risk Management for the Metropolitan City of Venice. Sci. Total Environ., 904, 166310. https://doi.org/10.1016/j.scitotenv.2023.166310 , 2023.

Kuznetsov, I., Rabe, B., Androsov, A., Fang, Y.-C., Hoppmann, M., Quintanilla-Zurita, A., Harig, S., Tippenhauer, S., Schulz, K., Mohrholz, V., Fer, I., Fofonova, V., and Janout, M.: Dynamical reconstruction of the upper-ocean state in the Central Arctic during the winter period of the MOSAiC Expedition, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1353, 2023. 

Bonthond, G., Beermann, J., Gutow, L. et al. Benthic microbial biogeographic trends in the North Sea are shaped by an interplay of environmental drivers and bottom trawling effort. ISME COMMUN. 2023, 3, 132. https://doi.org/10.1038/s43705-023-00336-3.
