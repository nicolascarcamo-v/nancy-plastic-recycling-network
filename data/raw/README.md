# Raw database

`model/build_instances.py` reads a single source file, `Base_de_datos_FINAL_v2.xlsx`
(~70 MB), and serializes the geographic core of each zone (buildings, hubs, candidate
facilities, coordinates, admissible arcs, monthly waste supply `w_imt` and parameters) to
`data/instances/zone_*.json`.

The Excel file is not versioned in git because of its size. The precomputed instances in
`data/instances/` are all that is needed to run the model and reproduce every result.
To regenerate them, place `Base_de_datos_FINAL_v2.xlsx` in this folder and run, from the
repository root:

    python model/build_instances.py

Main sheets: `Parameters`, `Polymer`, `Buildings` (with monthly `W_1..W_12`), `w_imt`,
`Hubs`, `Facilities`, `D_Build_Hub` (building→hub road distance, m) and `D_Cand_Hub`
(facility→hub road distance, m).
