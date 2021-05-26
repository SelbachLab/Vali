# Vali
Vali analyzes PRM raw files acquired by Thermo-Fisher Orbitrap instruments. It is optimized to work with [Picky](https://picky.mdc-berlin.de) generated spectrum libraries and provides semi-automated peak determination. All peak assignments can be manually adjusted or corrected. Subsequently, quantified fragment peaks can be exported as tab delimited text files.
## Requirements
* R > 4.0.3
* [RawDiag Requirements](https://github.com/fgcz/rawDiag):
    * [Mono](https://www.mono-project.com)
    * [MS .Net Framework](https://dotnet.microsoft.com/download)
* Packages:
    * list follows
 * You must set '.' as decimal seperator in your system. 

## Getting Started
1. Start Vali
    * Vali.bat (Windows)
    * Vali.sh (Linux or MacOs)
2. Go to the RawScanner Tab.
3. Set the Path to the Main Analysis folder, which contains your raw files (PRM, tDIA or DIA).
4. Set the Picky Export folder path. This folder must be the exact same one you exported from Picky. 
5. Choose the correct analysis mode (PRM or DIA) and adjust the mass accuracy settings. Leave the Recompile Database unset.
6. Press GO and wait for the program to finalize the analysis.
## Manual adjustments
The Validation Tool Panel allows you to go through all peptides and readjust the peak selection.

1. Use the mouse over to define the RT range of the peak and double click. 
2. Press the Set button. The corresponding area will be highlighted in grey. You can remove a peak definition by pressing the Remove button.
3. (Optional) You can rate the quality of the peak with the Rating buttons.

## Export
Press the Export button in the upper right Corner. After all sequences have been analysed Vali will provide you with the following list of tables:
1. TransitionList.txt (Quantitative information of Transitions, use this table for downstream analysis).
2. Peaks.txt (for replotting the peaks, or downstream analysis such as lm based ratio analysis.)
4. Ratios.txt (for SILAC pairs)

## Known issues
* The MS1 read out function is at the moment not working due to substantial changes in the internal data management.
* The fragment selection and exclusion does not work properly. The current workarround is to document bad transitions in order to remove them manually in a downstream analysis step. 

## Citation
Picky:
Zauber et al. (2018). Nature Methods. http://doi.org/10.1038/nmeth.4607


