# Pi0_Reconstruction

This Repo contains any protodune analysis work geared towards pi0 reconstuction. Main python file to analyse root files is EventSeection.py, to run:

```
python EventSelection.py --file <root file> --outName <file name to append> --outDir <output directory> --parameter <0 - 8> --conditional < l, g, ne, e > --cut <cut value> --beam <1/0> --plot <single/both>
```
---
options for what selection quantities and cuts can be made, see EventSelection.py.
