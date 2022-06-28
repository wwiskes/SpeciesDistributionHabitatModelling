# SpeciesDistributionHabitatModelling
 
In this repository are two scripts designed to create species distribution models from known presence locations. SDHMmain is the primary script where you input your data, set your parameters, select your predictors, and recieve your output. SDHMfunctions is a script containing functions written to make the main script run, you should never need to alter the functions script. Also provided is a CSV of random points throughout the state of Utah, provided as a placeholder for your own data.

The package contained in this repository has everything needed to complete a test run of the functions.

This code is designed to be ran in a cloud environment (although it can be run locally just fine).

https://cloud.google.com/vertex-ai

To run this code in Google Cloud's Vertex AI first deploy a virtual machine running R. This script was built using 16cpu @ 32ram. On these specifications this script typically takes around an hour to run.
Once the notebook is running click "Git Clone" in the control bar and enter: 

https://github.com/wwiskes/SpeciesDistributionHabitatModelling.git

This will clone this repository into your virtual machine, anytime you want to sync to check for new code revisions simply open the "Git" tab and select "Pull from Remote (force)". This will ensure you have the most recent version of this script.

Thanks and happy modelling. 
