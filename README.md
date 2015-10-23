# ERV_Simulations
version 1.0 to simulate the evolution of ERVs

<<<<<<< HEAD
ERV_Simulations comprise 4 main Python codes named them as "Master_Mortal.py", "Master_Imortal.py", "Transposon_Mortal.py", and "Transposon_Immortal.py".
These 4 Python codes make use of three Python classes (Random_functions.py, Tree_utilities.py, and writeFiles.py) with functions necessary for the main codes.
=======
I have upload 4 main Python codes named as "Master_Mortal.py", "Master_Imortal.py", "Transposon_Mortal.py", and "Transposon_Immortal.py".
I have also uploaded three Python classes (Random_functions.py, Tree_utilities.py, and writeFiles.py) which contains functions necessary for the main codes.
>>>>>>> origin/master

The main codes simulates phylogenetic trees following variations of the Master and Transposon models of ERV/Transposable element evolution in host genomes.
A paper has been recently submitted for publication, explaining in details these simulations. Once it is accepted, citation will also be included in this README file. 


####################################################################################################################
These codes requires the following Python version and packages
Python version 2.7. These codes were not tested with Python version 3.0.

Package ETE2: It can be downloaded at http://etetoolkit.org/download/
Obs. after installing ETE2 (depending on the way you decide to do it). Sometimes when running the code, it will appear some messages saying that some packages are missing from the installation of ETE2 (for example package pyqt) . You can ignore these messages and the code should still be able to run.

####################################################################################################################
How to run the code?
After installing Python and ETE2 package, go to the terminal and move to the directory where all the python codes were downloaded from GitHub.
If file is a zip file, you should unzip the files first.

In the directory in which the codes are, you can type in the Terminal:
python *code_name_here*

Sometimes instead of python you should type python2 (it depends on how the installation was made in your computer).
Please substitute *code_name_here* for the name of any of the four codes you are interested in reproducing the simulations

####################################################################################################################

When running the code you will get some messages printed in the terminal window. 
These messages start with "Elapsed time of first simulation is: ..."
This is only to show that the code is running properly. You can uncomment the line responsible to print such messages in the terminal.

These codes will generate lots of individual files that may starts with the following words: "Info", "subs_per_site" and "ultrametric".
The "Info" files contain some information about the simulations

The "subs_per_site" files contain phylogenetic trees in which branch lengths are represented in substitutions per sites

And the "ultrametric" files contain phylogenetic trees in which branch lengths are in host generations.

If you prefer that files are organised in a different way, modifications of these codes are allowed as far as the paper and these original codes are acknowledged.  

####################################################################################################################
For any questions or comments you can e-mail me at thednainus@yahoo.com


___________________

Reproduction or any form of use of these codes should be acknowledged. 
