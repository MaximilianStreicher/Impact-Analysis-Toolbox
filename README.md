# Impact-Analysis-Toolbox
Matlab scripts to do filtering, integration/summation and event analysis of load cell and pressure sensor signals 

The values in the beginning of each Matlab script (section "MANUAL INPUT"), have to be set manually by the user. The settings depend on the measured signals, filter desicions amongst others. Also, the settings as given match the example included in the IAT. 

"Input" directory:        here the pressure sensor and load cell measurement files are stored

"Output" directory:       here any outputs produced with the IAT will be stored

"Force_Filter.m":         Matlab script to do the filtering of each load cell signal

"Force_Total.m":          Matlab script to do the summation of individual load cells and convert the force into a line force [kN/m]

"Force_Event.m":          Matlab script to select the highest force per impact event

"Pressure_Filter.m":      Matlab script to do the filtering of each pressure sensor signal

"Pressure_Integration.m": Matlab script to do the integration of the pressure sensor array and result given in [kN/m]

"Pressure_Event.m":       Matlab script to select the highest integrated pressure (=force) per impact event

The Impact-Analysis-Toolbox (IAT) is beeing developed by Maximilian Streicher, PhD student at Ghent University. For questions and suggestions you may contact via email: maximilian.streicher@ugent.be (alternative email: str.max1@gmx.de).
