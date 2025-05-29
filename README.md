# Runningwheel_Setup_with_automatic_reward
bachelor thesis project


This project was developed as part of my Bachelor's thesis. It was my first coding project. Sorry in Advance that the coad is not very appealing and easy to understand. Throughout the development process, OpenAI's ChatGPT-4o was used as a supportive tool for programming assistance. This included help with generating new lines of code, iterating on existing implementations, and refining overall functionality. All final decisions and code integration were made by the author

Starting with the code of the arduino setup with the name: Aduino_Code. There the miscalculation of the speed threshold is explained in the comments. In there all parameter changes regarding the setup can be made.

Contiuning with only_needed_singletrialanalysis: This file loads all the data where first the path /storage/share/matlab/labbox is needed. The LFP data is loaded with LoadBinary. When loading, manual renaming of the open-ephys session file is needed. In this file the filtering and interpolation is performed and at the bottom excel sheets are created of every session. 
plot_whaticutoutschaukeln is then used to ensure that the filtering worked. It is important to run only_needed_singletrialanalysis before plot_whaticutoutschaukeln. 

Cellloading loads the data from the excel files into matlab to compute everything needed for the plots. DeprxTime is used to create dictionaries for sessions for some of the stats out of the excel sheets. detect_rotary_encoder_activity and extractAnimalAndTrialFromFileName are needed functions. 

Afterwards all the plots were made who are in the directory storage3/manus/.

A metal running wheel (from Lafayette Instruments with the model name 80860W), which had a diameter of 36 centimeters from. Along with a rotary encoder (with the model LPD3806-600BM-G5-24C) were used
