#!/bin/bash
##source /home/davidm/Desktop/Cygno/python395Environment/venv3.9.5/bin/activate                     
for i in {11301,11400,11450,11500,11960,11961,11962,11963,11964,11965,11966,11967,11968,11969,11970};
do QT_QPA_PLATFORM=wayland python pmt_reconstruction.py configFile_LNGS.txt --pdir ./images -j10 -r ${i} --tmp /run/media/davidm/DAVID_2Tb/LIME_data; done





##new background runs with fitted L
##{11301,11400,11450,11500,11960,11961,11962,11963,11964,11965,11966,11967,11968,11969,11970}
##folder: /run/media/davidm/DAVID_2Tb/LIME_data

##background runs
##{12000,12001,12002,12003,12004,12005,12006,12007,12008,12009,12100,12101,12102,12103,12104,12105,12106,12107,12108,12109};

##iron runs diff positions
## folder: dataset_55Fe_Different_Z
##{9137,9138,9139,9140,9141,9758,9759,9760,9761,9762,11283,11284,11285,11286,11287,11954,11955,11956,11957,11958,12170,12171,12172,12173,12174};
