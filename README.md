# Bayesian_Ground_Motion_Selection
These are a set of codes to perform a Bayesian analysis of the Conditional Spectrum approach for ground motion selection. A set of instructions to run the codes are provided below.

# Corresponding paper
"A Bayesian Treatment of the Conditional Spectrum Approach for Ground Motion Selection". Somayajulu Dhulipala (Virginia Tech), Jack Baker (Stanford), and Madeleine Flint (Virginia Tech).

Created by Somayajulu Dhulipala, 20 May 2017

Modified by Somayajulu Dhulipala, 23 December 2017

PLEASE REPORT ANY ISSUES OR PROBLEMS WITH THIS CODE TO
lakshd5@vt.edu

# Instructions
- Dowload all the codes to a single directory in your PC which is also the current directory in MATLAB.
- You can run the codes by directly executing Master.m and following the instructions displayed in the MATLAB command window.
- Two options after executing Master.m: 
   
   (1) Define the conditioning IMs and other options by your self => enter 0 for the question 'Do you have a ready to use, Bayesian amenable workspace?'
   
    - Either a ground motion database of your choice can be used after executing Master.m or the provided W_Curtailed_NGA_West2.mat database can be utilized.
    - If you'd like to modify the priors, simulate some ground motions, make the simulated ground motions Bayesian amenable by running the first part of Master.m where the workspace is set up, and select the right option for the question: 'Do you want to use simulated priors? (1-YES, 0-NO)'. And, don't forget to change the workspace name in 'Prior_Extraction.m' (line 41) to the one defined by you.
   
   (2) Stick to some pre-defined conditioning IMs => enter 1 for the question 'Do you have a ready to use, Bayesian amenable workspace?'
   
    - Then load anyone of the pre-defined workspaces provided here. 
    - For a description of these workspaces take a look at 'Workspaces description.txt'
    - If you'd like to modify the priors, select the right option for the question: 'Do you want to use simulated priors? (1-YES, 0-NO)'. The code automatically loads W_EXSIM_Priors1.mat for generating the priors.
    
   
- Of course, option (2) is much simpler to follow if you are running these codes for the first time.
