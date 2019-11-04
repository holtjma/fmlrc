# Examples
## Run _E. coli_ Example
### Commands
1. Go to a directory to store fmlrc and clone this repository:
```bash
git clone https://github.com/holtjma/fmlrc.git
```
2. Change into the fmlrc example directory:
```bash
cd fmlrc/example
```
3. Run the examples script.  This will automatically install ropebwt2 to a local path, download short- and long-read _E. coli_ data, build the multi-string BWT using ropebwt2 and fmlrc-convert, and lastly run use fmlrc to correct the first 400 reads in the long-read data:
```bash
./run_example.sh
```

### Outputs
1. `./ecoli_comp_msbwt.npy` - this contains the run-length encoded multi-string BWT using the same encoding as the [msbwt](https://github.com/holtjma/msbwt) python package.  Note: if you wish to use this with msbwt, it will need to be placed in it's own directory and renamed to `comp_msbwt.npy` in that directory.  For more information, please refer to the msbwt wiki pages.
2. `./corrected_final.fa` - this contains the first 400 long reads after correction using fmlrc.