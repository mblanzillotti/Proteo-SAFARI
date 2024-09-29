PROTEO-SAFARI

Hello, and thank you for your interest in Proteo-SAFARI, an application built in R Shiny designed to enable fragment ion identification in averaged top-down MS/MS spectra of intact proteins based on their isotopic distributions directly in the m/z domain. This application can process a single-scan Thermo Orbitrap ".raw" file (typically averaged and exported using FreeStyle or QualBrowser) and a corresponding protein sequence, analogous to those approaches used in programs like ProSight Lite. Proteo-SAFARI includes support for highly specific post-translational modification assignment (either on the peptide backbone nitrogen, side chain, or carbonyl), and automated detection of hydrogen shifted fragment ions found commonly in ultraviolet photodissociation spectra. The following sections in this document will detail the installation process, and brief descriptions of each module. The "Install Prerequisites" file includes the required packages for Proteo-SAFARI, and should be utilized whenever using Proteo-SAFARI for the first time or in a new version of R.


INSTALLATION

Proteo-SAFARI is a Shiny application utilizing the module format - each "module" represents a self-contained piece of the application, and typically corresponds to a particular tab in the app. For the app to run properly, please include each of the "...Module.R" files in the same directory as the "app.R" file so each of the dependencies is able to load properly.

To install the requisite packages utilized within Proteo-SAFARI, please run the "Installation Prerequisites.R" file, which contains install.packages commands for each of the required packages for Proteo-SAFARI. Most notably, the "rawrr" package installed via Bioconductor requires agreement to to Thermo RawFileReader end user license, which is necessary for direct read-in of Thermo .raw files into the R environment. Additionally, ensuring the rawrr .dll files and .exe file runs properly is critical, and these steps are included in the Installation Prerequisites.

With all packages installed, you may now proceed with running the application! To do this, only open the "app.R" file in your preferred IDE, and run the file or click "Run App" if using RStudio. This should open a new window containing the application.


PREREQUISITES AND HELPERS

This file contains necessary reference tables for proper sequence parsing and mass generation. This includes a table of monoisotopic masses for the elements considered in Proteo-SAFARI and a table of symbols and their corresponding chemical formulae. Any modifications with additional heteroatoms must have the monoisotopic mass of that heteroatom added to the element_masses table, and be included as a new column in the symbol_compositions table for proper isotopic modelling.


SEQUENCE PARSING

This module contains the sequence parsing strategy that breaks up a protein into its nitrogen (n), residue (R), and carbonyl (c) components, which are used for generation of fragment ions masses later in the application. The main regular expression requires definition of a residue (single-letter codes of the 20 canonical amino acids), and optinally allows for assignment of a modification on the nitrogen, residue, or carbonyl by writing out the net chemical formula (e.g., C6H12O6) in parentheses after the appropriate symbol (e.g., a methylalanine would look like: A(Methyl) or A(C1H2). A methylaion of the backbone amide would appear like: n(Methyl)A). The resulting data.table of the parsed sequence is passed to the chemical compositions module.


CHEMICAL COMPOSITIONS

The parsed sequence is sent to this module for determination of the chemical composition of each row in the parsed sequence table. These chemical compositions are either joined from the symbol_compsitions table, or interpreted if they are explicitly written out (e.g. C6H12O6). Chemical compositions associated with each row are also organized by which part of the protein at each residue position that they belong to, either nitrogen (n), residue (R), or carbonyl (c), before being passed to the mass generation module. Additionally, the net precursor composition and the composition of each fragment ion (generated as a cumulative sum of each element) are generated in this step. This module also considers generation of neutral-loss fragments (currently b-H2O, y-H2O, b-NH3, and y-NH3) based on a toggle selection from the User Inputs module.


MASS GENERATION

In this module, the chemical compositions for the precursor and fragment ions are converted into neutral masses, first monoisotopic and the isotopic distributions using the IsoSpecR package. The IsoSpecify function is applied to each compsotion, generating a list of fine isotopic abundances, which are then condensed based on their nominal mass to generate an isotopic distribution compatible with the observed data in typical top-down mass spectrometry experiments. These neutral mass data are then sent to the m/z generation module. Neutral fragment ion mass values can be downloaded from within the application.


MZ GENERATION

The m/z generation module takes the neutral masses of each ion based on the precursor charge and number of charge states arguments from the User Inputs module. This approach constrains the number of fragment ion charge states considered by Proteo-SAFARI to only those closest in charge density (thus m/z) to the precursor ion, which for a denatured protein is a reasonable assumption. Charge-reduced precursor ions or native-like precursor ions may need a larger number of fragment ion charge states considered. m/z values are generated based on these user-defined parameters, and are by default considered with a proton adduct and from m/z 200-2000. These parameters can be changed in this module. Fragment ion m/z values can be downloaded from within the application.


DATA READIN

This module handles the uploaded filepath from the User Inputs module and reads in a subset of data from the rawfile, most importantly the centroid data stream and profile data stream, if available.


FRAGMENT IDENTIFICATION

Fragment ions are identified based on the isotopic distributions and their corresponding m/z values calculated previously within the application and the centroids loaded in from the .raw file. These theoretical m/z values are first matched within a user-defined m/z tolerance (10 ppm by default), where consecutive peaks amounting to over 70% of the predicted isotopic abundance must be identified before proceeding with fitting. Additionally, in this step potential hydrogen shift fragments are identified and if found are also moved to fitting with the parent fragment. Using an nnls approach, the fragment ion abundances are predicted, reporting a predicted total abundance for each fragment ion considered. The average preducted abundance error must be within 25% for the fragment ion to be considered a match.


VISUALIZATION

This module takes outputs from previous modules and utilizes a set of predefined visualization functions to show fragment ion identification maps, an interactive annotated spectrum, and summaries of the identified fragment ions, all of which can be saved from within the application.


USER INPUTS

This module coordinates user input data from the app's side panel with the modules that require these data.



Please direct any questions about Proteo-SAFARI to m.b.lanzillotti@gmail.com.


