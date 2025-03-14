# Climate-CWatM (C-CWatM)

## A simplified version of the Community Water Model (CWatM) to enable direct coupling to climate models.

Climate-CWatM (C-CWatM) is a flexible modelling tool that can be coupled to (regional) climate modelling systems. It is based on the socio-hydrological model [CWatM](https://cwatm.iiasa.ac.at) and enables the simulation of various hydrological quantities, such as discharge, sectoral water use, groundwater levels, and water reservoirs, as part of climate model simulations. C-CWatM is designed as a standalone Python model that can be easily integrated into various climate modelling systems as it operates on standard climate and land surface variables.
C-CWatM is developed and maintained by the [UWaRes](https://ms.hereon.de/uwares/) team at the Helmholtz-Zentrum [Hereon](https://www.hereon.de/).

> [!WARNING]
> C-CWatM is still under development and not fully functional in its current form.

## Recent development

### DONE

- Include flag in settingsfile to turn off paddy irrigation. If turned off, paddy fields will be treated as nonpaddy fields.
- Check which modules are not necessary anymore, remove those modules

### TODO

- Remove unnecessary code from existing modules
- Include OASIS coupling
- Read in soil parameters (field capacity, etc.) as input data and remove all unnecessary input data from settingsfile.
- many other things

## Input Data

Please check the CWatM repository [CWatM-Earth-30min](https://github.com/iiasa/CWatM-Earth-30min), which contains input data for CWatM and C-CWatM at 30 arcminutes and further links to the 5 arcminutes resolution input data. A dedicated repository for C-CWatM will follow.
