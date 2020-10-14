# ProSpect (R package)

## Synopsis

Core package containing all the tools for simple and advanced spectral energy distribution analysis (generating and fitting).

## Getting Started

New to Synthetic Spectra? Then take a look here: [sedfitting.org](http://www.sedfitting.org/Models.html). There's a bunch of useful stuff there, and a lot of dead links :-(. But it is worth exploring in some detail. In particular for the stellar population models have a look at BASTI, BPASS, Galaxev (which is BC03 to most astronomers), MILES, Pegase, SLUG and Starburst99. For dust stuff have a look at Dale+Helou and Draine+Li, for UV-FIR Da Cunha (which is MagPhys), CIGALE and Grasil. Those are the most popular variants of what they do in their respective fields. Caveats abound about which is better, but these days they are all pretty sophisticated in their own way.

## A New ProSpect

**ProSpect** is a package that aims to help users explore star formation histories (SFH) and spectral energy distributions (SED). Beneath it all it makes use of the BC03 (?BC03) and EMILES (?EMILES) sythetic stellar population libraries. On top of this it uses the Charlot and Fall model for birth cloud and screen dust attenuation (?CF) and re-emits MIR to FIR flux via the incorporation of the most recent Dale dust templates that model the heating of dust by a radiation field and AGN (?Dale). **ProSpect** can handle general SFHs via arbitrary functional forms (**SFHfunc**).

Conceptually this probably all sounds a bit similar to MagPhys and Cegale, bar the different dust model (Dale rather than greybody, though take a look at ?greybody) and some additional SSP libraries (?EMILES). To a degree this is true, but the main reason for putting it all together is to allow proper generative creation of SEDs via the alteration of user accessible parameters and arbitrary SFHs. MagPhys and similar codes do not allow easy access to the under-the-hood generative functionality which might allow this. Doing this is with a longer term aim of incorporating **ProSpect** into **ProFit** for multiband morphological decomposition of galaxies, revealing their component-wise SFHs. Before the completion of that component of the project **ProSpect** offers a readily accessible interface to multiband fitting of SED in order to measure (e.g.) stellar and dust masses, and it can also applied to arbitary SFHs computed by other codes, e.g. create a realistic SED for a semi-analytic (SAM) model SFH.

For most uses getting up and running, having a play with the simplest 5 phase burst star formation (**massfunc_p5**) is probably the best way to start. Have a read of the documentation and Examples in ?SFHfunc. Fitting data is a more complex affair, but users can get a feel for how this works in practice by looking at ?SFHfunc. The basic idea of these likelihood functions is that they allow the user to provide observed Data and can marginalise over various types of SFH and dust models in order to constrain a set of target parameters. It is also a good idea to get acquainted with the main vignettes covering the BC03 and EMILES data compiled into the package (SSPchecks and Comparisons). These give a good idea how to access the data at a low-level, and allows for easy manipulation and potentially user defined functions for processing them.

The 5 phase model (**massfunc_b5**) covers 4 key phases of star formation and models them as having constant star formation in each of these periods. By default they cover 0 - 100 Myr (burst, generally considered the shortest phase that broad band photometry can be sensitive to); 100 Myr - 1 Gyr (young; dominated by hot young stars and violent phases of stellar evolution); 1 Gyr - 5 Gyr (mid), 5 Gyr - 9 Gyr (old), and 9 Gyr - 13 Gyr (ancient). It is not really possible to break the SFH into any more independent phases than this using broad band photometry alone, but a physically motivated functional model (e.g. an exponentially decling SFR, or constant SFR) can be used using the functional interface in **SFHfunc**.

For a typical GAMA galaxy using the 21 band photometry from LAMBDAR a full SFH plus dust model fit takes about one minute using the **optim** function with the Nelder-Mead minimisation algorithm. A full MCMC will typically take about ten times longer, but in cases of a reasonable fit this is probably over-kill since the Laplace Approximation can be used.

## Summary of Things Used

Property | ProSpect Type  |  Relevant Help
----------------- | -------------------- | ---------------
SSPs | BC03 or EMILES  |  ?BC03, ?EMILES, ?BPASS
IMF | Chabrier (2003) |  ?BC03, ?EMILES, ?BPASS
Isochrones | Padova (1994) for BC03, EMILES   |   ?BC03, ?EMILES
Isochrones | STARS for BPASS |   ?BPASS
Stellar Atmospheres  |  See relevant text in BC03, EMILES, BPASS papers |   ?BC03, ?EMILES, ?BPASS
Dust Attenuation Law | Charlot and Fall (2000) for birth cloud and screen   |  ?CF
Dust Model  |  Dale et al (2014)  |  ?Dale

## Installation

### Getting R

First things first, you will probably want to install a recent version of **R** that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-built binaries on CRAN you might be waiting much longer.

#### Mac

For Mac just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/macosx/>

#### Windows

For Windows just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/windows/>

#### Linux

Debian:	`sudo apt-get install r-base r-base-dev`

Fedora:	`sudo yum install R`

Suse:	More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu:	`sudo apt-get install r-base-dev`

All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

If you have a poorly supported version of Linux (e.g. CentOS) you will need to install **R** from source with the development flags (this bit is important). You can read more here: <https://cloud.r-project.org/sources.html>

Now you have the development version of **R** installed (hopefully) I would also suggest you get yourself **R-Studio**. It is a very popular and well maintained **R** IDE that gives you a lot of helpful shortcuts to scripting and analysing with **R**. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

If you wish to use the command line version of **R** on Mac (why?!) then you might need to separately install **XQuartz** and set the DISPLAY system variable via something like export DISPLAY=:0 (this is not an issue for most people however).

### Getting ProSpect

Source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/ProSpect")
library(ProSpect)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools (which I do not use as the default since it has a low more dependencies, and is tricky to install on HPCs):

```R
install.packages('devtools')
devtools::install_github("asgr/ProSpect")
library(ProSpect)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true”)
remotes::install_github("asgr/ProSpect")
```

I also have these options set by default in my .Rprofile, which seems to help with some of the remote install issues some people face:

```R
options(download.file.method = "libcurl")
options(repos="http://cran.rstudio.com/")
options(rpubs.upload.method = "internal")
```

If all of these do not work than the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install ProSpect_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

If none of the above works then you should consider burning your computer in sacrifice to the IO Gods. Then buy a newer better computer, and try all the above steps again.

Failing all of the above, please email me for help.

Assuming you have installed all of the packages that you need/want, you should now be able to load **ProSpect** within **R** with the usual:

```R
library(ProSpect)
```
