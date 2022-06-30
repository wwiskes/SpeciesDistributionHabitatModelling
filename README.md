# Species Distribution Habitat Modelling (SDHM)
 
In this repository are two scripts designed to create species distribution models from known presence locations. SDHMmain is the primary script where you input your data, set your parameters, select your predictors, and receive your output. SDHMfunctions is a script containing functions written to make the main script run, you should never need to alter the functions script. Also provided is a CSV of random points throughout the state of Utah, provided as a placeholder for your own data.

The package contained in this repository has everything needed to complete a test run of the functions.

This code is designed to be ran in a cloud environment (although it can be run locally just fine).

https://cloud.google.com/vertex-ai

To run this code in Google Cloud's Vertex AI first deploy a virtual machine running R. This script was built using 16cpu @ 32ram. On these specifications this script typically takes around an hour to run.
Once the notebook is running click "Git Clone" in the control bar and enter: 

https://github.com/wwiskes/SpeciesDistributionHabitatModelling.git

This will clone this repository into your virtual machine, anytime you want to sync to check for new code revisions simply open the "Git" tab and select "Pull from Remote (force)". This will ensure you have the most recent version of this script.


# Tutorial to run in the cloud

This tutorial assumes you have access to the Google Cloud Platform, are proficient in the R programming language, and are capable of installing the required R libraries in a linux environment.

First navigate to the GCP project you wish to use:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture.PNG)

Next search for Vertex AI:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture2.PNG)

Once in Vertex AI select the 'Workbench' tab and click 'Create Notebook':

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture3.PNG)

Give the notebook a creative name and select the region closest to you:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture4.PNG)

Select Debian as the operating system and R as the environment:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture5.PNG)

Select the machine type, I have found personally through testing that the bottle neck with species distribution modelling is not RAM (like most people often say), but CPU, so I always choose a high CPU machine. I highly recommend using a compute optimized machine (C2), over an E2 or N2, as C2 instances are tuned for the High-performance computing (HPC) that SDHMs require. By keeping the RAM low we can boost up the CPU while keeping costs reasonable. Next time you're running a SDHM open a linux terminal and type 'top' - this will display your resource consumption. 

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture6.PNG)

After you click 'Create' you should see your notebook in your workbench. Now click 'Open JupyterLab':

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture7.PNG?)

Now a page will open with your Jupyter environment. Click the Git symbol to connect to a repository. Paste in: https://github.com/wwiskes/SpeciesDistributionHabitatModelling.git

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture8.PNG)

Navigate to the script titled 'SDHMmain' and open it with a notebook:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture9.PNG)

Use 'Pull from Remote (force)' to pull in new revisions of this code:

![alt text](https://wwiskes.github.io/datadump/SDHMtut/Capture10.PNG)


Thanks and happy modelling!
