# Eavesdropping micropredators as dynamic limiters of sexual signal elaboration and intrasexual competition
## Authors: Leavell BC, Beaty LE, McNickle GG, Bernal XE

Dates of data collection: Summers of 2010 and 2012

Location of data collection: Gamboa, Panamá

Keywords: _communication, foraging, male-male competition, predator-prey, structural equation model_

## Files:
HighstatLibV6.R
*   Code from Zuur AF, Ieno EN, Walkner N, et al (2009) Mixed effects models and extensions in ecology with R. Springer New York, New York, NY
    -   Available at https://highstat.com/Books/BGS/GAMM/RCodeP2/HighstatLibV6.R
    -   Necessary to run part of analysis in "MicropredatorsLimitSignalElaboration_code.R"
    -   Downloaded 2019-07-20

MicropredatorsLimitSignalElaboration_code.R
-   Script for analyses
*   Variables
    -   date = date (day-month-year) of observation
    -   time = duration (seconds) from beginning of 1st call to beginning of 50th sequential call
    -   call_rate = # the total number of calls (50), minus one, divided by the time from the beginning of the first call to the beginning of the last call 
    -   chucks = total # of chucks over the 50 sequential calls
    -   midges = total # of frog-biting midges observed landing on focal frog over 50 sequential calls  
    -   swatcount = total # of swats observed over the 50 sequential call duration
    -   males_lessthan1m = # of neighbor male competitors present within 1 meter of focal frog
    -   males_morethan1m = level of perceived abundance of calling conspecifics beyond 1 meter. (0 = only focal frog heard calling, 1 = individual calling frogs could be counted, 2 = calls of frogs overlapping but individuals distinguishable, 3 = full chorus, cannot distinguish individuals).
-   Modified 2021-11-24 (Added author names)

MicropredatorsLimitSignalElaboration_data.csv
-   Dataset used in "MicropredatorsLimitSignalElaboration_code.R"
-   Created 2020-09-14

Movies (https://doi.org/10.5281/zenodo.5759167)
-   1A10M3.mp4, 1A10M4.mp4, 3A12M1.mp4, 3A12M3.mp4, 7A12M3.mp4, 8A12M2.mp4, 11J10M5.mp4, 12J10M1.mp4, 12J10M2.mp4, 29J12M3.mp4
-   Recordings from 2010 and 2012 of calling male túngara frogs (*Engystomops pustulosus*). Most movies feature males being attacked by frog-biting midges (Diptera: Corethrellidae) and their anti-midge defensive swats. The size of the inflated vocal sac relative to the length of the frog's arm prevents the male from swatting while calling.
-   File names include the following information, in order: Day of month, Month ("J" = July, "A" = August), Year ("10" = 2010, "12" = 2012), Sex of frog ("M" = Male), Identifier for individual male on that date (i.e., "1", "2", "3", ...)
-   Created on date denoted in file name

## Methods: 
In Gamboa, Panamá, during the summers of 2010 and 2012, we observed calling male túngara frogs (*Engystomops pustulosus*) in situ. For each observation, which took place between 1900 and 2300 h, we illuminated a calling male (hereafter, 'focal frog') with infrared LED arrays (Sima SL-201R) and recorded audio and video using a Sony DCR-SR220D 120GB Handycam Camcorder. Following an observation, the focal frog was toe-clipped to prevent pseudoreplication and released at point of capture in adherence to protocols established by the American Society of Ichthyologists and Herpetologists (<https://asih.org/animal-care-guidelines>). We conducted all research in compliance with Panamanian legal and ethical regulations (MiAmbiente collection permits: SE/A-67-10, SC/A-20-12) and following IACUC protocols (Texas Tech University: 11056-08; Smithsonian Tropical Research Institute: 2011-0616-2014-11). See "Methods" in manuscript for more details.

To analyze calling behaviors, we selected a sequence of the first 50 clean, consecutive calls, per focal frog. We used the duration of the sequence, which varied per frog, to derive the call rate--i.e. the total number of calls, minus one, divided by the time from the beginning of the first call to the beginning of the last call. The total number of chucks produced was counted over the same time. To assess the threat of frog-biting midge (Diptera: Corethrellidae) attacks, the total number of instances in which midges landed on the focal frog were counted from video playback during the 50-call sequence. We also counted the total number of times the focal frog swatted during the same sequence. Out of a total of 100 focal frogs that were recorded, we omitted data from fifteen individuals due to high background noise levels or poor video quality that precluded behavioral analysis. The observed dataset thus includes data from a total of 85 males.
