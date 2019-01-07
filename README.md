# FlexibleStructuralAlignment
Flexible Structural Alignment : Flexible protein alignment by a method of
hierarchical segmentation

A reccurent issue in structural alignment is to compare protein in the same folding family with different conformations. In fact, proteins are not rigid body.

In this method, a protein is cut into Protein Units (PUs) and these PUs are aligned on the second protein and conversely. These PUs are produced by the _Protein Peeling_ algorithm. Then, for each cutting, these PUs are aligned and a TMscore is calculated by TMalign.

To compare this method to others, TMalign, a classic rigid algorithm, is implemented and parMATT, a flexible algorithm, too.


## Install and execute
To install it, first install all the requiered module with correct version with this command:
```
pip install -r requirements.txt
```

Then, it's possible to run the programm with this command:
```
./main.py -i <file> -c <chain1> -f <file2> -d <chain2>
```


## Running the test
For this program, a benchmark is given and saved in data folder. If you want to test it you can use this command :
```
./bench.py
```

## Example output
![exemple_output](https://github.com/sophielem/FlexibleStructuralAlignment/blob/dev/results/example_output.png "Example output")
```
		****************
		* 1l5b vs 1l5e *
		****************
Number of PUs :  2                 Score : 0.9767
Number of PUs :  3                 Score : 0.9771
Number of PUs :  4                 Score : 0.9782

		****************
		* 1l5e vs 1l5b *
		****************
Number of PUs :  2                 Score : 0.9664
Number of PUs :  3                 Score : 0.9668
Number of PUs :  4                 Score : 0.9688

		  ***********
		  * TMalign *
		  ***********
Score : 0.55845


		  ***********
		  * parMATT *
		  ***********
Score : 0.5625
```

## Authors

Sophie LEMATRE
