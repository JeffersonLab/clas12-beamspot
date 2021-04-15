# BeamSpot

Analysis of the beam position.

Author: F. Bossù (fbossu at jlab.org)<br>
Based on S. Stepanynan's work: CLAS12 Note 2020-003

## How to use it

You can use the `run.sh` script or<br>
`java -Xms1024m -cp ".:$CLAS12DIR/lib/clas/*:$CLAS12DIR/lib/plugins/*" BeamSpot [arguments]`

If you run it on some hipo files, it runs the beam spot analysis. It saves some canvases, the results and a txt file containing the histograms.

If you have run it already (let's say on the farm) and you want to re-analyse the output without the hipo files, then you can run it without arguments as long as the `h2_z_phi.txt` file is in the current directory.

*Examples*<br>
> ./run.sh *hipo

or

> ./run.sh

## What it does

This class analyses electrons and it stores the z position and the $\phi$ at the vertex for various $\theta$ bins.<br>
The z position of the target foil is extracted for each $\phi$ bin and each $\theta$ bin and its modulation versus $\phi$ is fitted with the function  

<img src="https://render.githubusercontent.com/render/math?math=f\(\phi\\)=A-B\cos\(\phi-\phi_0\\)/\tan\(\theta\\)">

The (x,y) position of the beam is then computed using B and $\phi_0$. The average values over the $\theta$ bins are then saved as results.

## What it outputs

The final results are saved in `beamspot_results.txt`<br>
All the canvases are saved in png files.<br>
The `h2_z_phi.txt` file contains all the 2D histograms in text format: to be saved for "offline" analysis.

## Tweaks

Any refined cut on the electron can be added in the functions `checkTrack` and `checkParticle`. <br>
If the target foil position fits seem to not work properly, you can change the initial paramters of the fit function inside the `analyse`a function.

## TODO

A better CLI

