# FlexibleStructuralAlignment
Flexible Structural Alignment : Flexible protein alignment by a method of
hierarchical segmentation

A reccurent issue in structural alignment is to compare protein in the same folding family with different conformations. In fact, proteins are not rigid body.

In this method, a protein is cut into Protein Units (PUs) and these PUs are aligned on the second protein and conversely. These PUs are produced by the _Protein Peeling_ algorithm. Then, for each cutting, these PUs are aligned and a TMscore is calculated by TMalign. 


## Install and execute
To install it, first install all the requiered module with correct version with this command:
```
pip install -r requirements.txt
```

Then, it's possible to run the programm with this command:
```
./main.py -i <file> -p <file2>
```


## Running the test



## Authors

Sophie LEMATRE
