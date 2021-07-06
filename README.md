# BeamSpot

Analysis of the beam position.

Author: F. Bossù (fbossu at jlab.org)<br>
Based on S. Stepanynan's work: CLAS12 Note 2020-003<br>
Special thanks to N. Baltzell

## How to use it

You can use the `run.sh` script

> ./run.sh

The help message already will give you some info:

>   Usage : BeamSpot  [input1] [input2] ....

>   Options :

>        -B : Batch mode, no graphics (default = 0)
>        -H : Interpret inputs as HIPO histogram files (instead of DSTs) and add them together (default = 0)
>        -O : String prefix for output file names (default = BeamSpot)
>        -T : Interpret input as a TXT histogram file (instead of HIPO DSTs) (default = 0)
>        -X : Run with no output files (default = 0)

Example: running over many DST files in bach mode and then merging the output for the final analysis<br>

> ./run.sh -O BeamSpot[number] /path/to/rec_clas_003219.evio.[number]*hipo

Changing the `[number]` at your covenience. Then it is possible to merge all the output histograms either reading the output txt or hipo files

> ./run.sh -T 1 *txt

or 

> ./run.sh -H 1 *hipo

The output CCDB table will be then produced.

## What it does

This class analyses electrons and it stores the z position and the $\phi$ at the vertex for various $\theta$ bins.<br>
The z position of the target foil is extracted for each $\phi$ bin and each $\theta$ bin and its modulation versus $\phi$ is fitted with the function  

<img src="https://render.githubusercontent.com/render/math?math=f\(\phi\\)=A-B\cos\(\phi-\phi_0\\)/\tan\(\theta\\)">

The (x,y) position of the beam is then computed using B and $\phi_0$. The average values over the $\theta$ bins are then saved as results.


## Tweaks

Any refined cut on the electron can be added in the functions `checkTrack` and `checkParticle`. <br>
If the target foil position fits seem to not work properly, you can change the initial paramters of the fit function inside the `analyse`a function.



