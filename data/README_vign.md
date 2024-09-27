# Explaination on vignetting files
The vignetting maps in this folder are taken with different cameras and conditions. The vignetting maps are stored in the root files in this folder.< \br>
To understand which vignetting map to use, consider that you need the correct camera. The file named `Vign_cam_matcher.txt` contains a txt map with the vignette root filename and the camera used for it. The reco code will read this file and prevent you from choosing a vignette map with a camera different than the one used in the data file you want to analyse. < \br>
Once the right camera is chosen, the vignetting map may change depending on the camera setup used: lens type, lens aperture, gold spacers between lens and sensor, focus of the sensor, etc.
The best would be to have one vignette map per configuration. **When a different lens or a different aperture is used you need to have a dedicated vignette map.** < \br>
Gold spacers, focus of the sensor should have minor impact, so if a specific map is not available use one in similar conditions (but be careful).< \br>
The vignette maps are not stored automatically in the data folder of the reconstruction repository to avoid saving unnecessarily heavy files. You can find all the described vignette maps in cygno-analysis/Useful_for_reco/Vignette/ .For example one can dowload with

`wget https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/cygno-analysis/Useful_for_reco/Vignette/vignette_QEHD_085.root`

In the following, there is a brief description of the vignette maps available.

## vignette_run03806.root
Camera: Fusion
Lens: Xenon Schneider
Aperture: 0.95
Distance to focus, golden rings and focus of lens as of LIME during data taking
Performed with images with cosmics on LIME turned on, so includes disuniformitis of electric field. Suggested not to be used.

## vignette_run04117.root
Camera: Fusion
Lens: Xenon Schneider
Aperture: 0.95
Distance to focus, golden rings and focus of lens as of LIME during data taking
Performed with images with cosmics on LIME turned on, so includes disuniformitis of electric field. Suggested not to be used

## vignette_runs03930to03932.root
Camera: Fusion
Lens: Xenon Schneider
Aperture: 0.95
Distance to focus, golden rings and focus of lens as of LIME during data taking
Performed with images facing a illuminated wall. Used for LIME data taking and for most of MANGO and GIN data as lenses and camera are the same and there are no dedicated vignette maps.

## vignette_QX_095.root
Camera: Quest
Lens: Xenon Schneider
Aperture: 0.95
Distance to focus: 657 mm (lens to focal plane)
Golden rings: none 
Focus of lens: 0.6 
Performed with images facing an illuminated wall. Original 0.95 image was saturated.  Slightly lower aperture used. Checked that there was no significant change. Also this was not taken in darkness with light shone on it.

## vignette_QEHD_095.root
Camera: Quest
Lens: EHD 25085
Aperture: 0.95
Distance to focus: 666 mm (lens to focal plane)
Golden rings: none 
Focus of lens: 0.6 
Performed with images facing an illuminated wall (or white paper sheet). This was not taken in darkness with light shone on it.

## vignette_QEHD_085.root
Camera: Quest
Lens: EHD 25085
Aperture: 0.85
Distance to focus: 666 mm (lens to focal plane)
Golden rings: none 
Focus of lens: 0.6 
Performed with images facing an illuminated wall (or white paper sheet). This was not taken in darkness with light shone on it.