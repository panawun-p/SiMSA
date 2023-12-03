# SiMSA
*Code unmaintained after 2016, might not be working properly*  

Similarity and Identity Calculation Multiple Sequence Alignment  
(Developed in Ubuntu Linux with Python ver.2)

## Set Up
### Requirements
- Python 2.x.x  
  For Python Interpreter please refer to https://wiki.python.org/moin/BeginnersGuide/Download
- Tkinter python package (For SiMSA GUI version)  
  For Tkinter package please refer to http://tkinter.unpythonic.net/wiki/How_to_install_Tkinter for the instructions.

## Running Script Version
### Running the script
Script version can be run in the command line using:
```
python SiMSA.py (parameters)
```
User needs to specify parameters in the command. Please refer to the next section for details.
The output will not be displayed on the screen but will be automatically saved to the file. 
### Parameters
- ***Input Parameter***  
  *This is a required parameter.* This parameter can be specified by adding `-i (input file name)` to the command excluding the bracket:
  ```
  python SiMSA.py -i (input file name)
  ```
- ***Output Parameter***  
  Let user specify the path for the output file. If this option is not used, the program will output file to same folder containing 'SiMSA.py'. User should only specify the path and the prefix of desired file name, the script will automatically add the rest of file name. *Output folder should already exist before running the program.* Use `-o (output path with file name prefix)` should be added to the command.
  ```
  python SiMSA.py -i (input file name) -o (output path)
  ```
- ***Report Parameter***  
  Let user specify which type of output file they want. If this option is not specified the program will save both report file and matrix file as the output.
  Use: `-r rpt` or `-r report` to output only the report file with statistical summary.   
  ```
  python SiMSA.py -i (input file name) -r report
  ```
  Use: `-r mat` or `-r matrix` to output only the matrix file containing the similarity and identity value between sequences.  
  ```
  python SiMSA.py -i (input file name) -r matrix
  ```
