# Production

In order to submit grid jobs, you will first have to build the `sbndcode` you wish to use and create a tarball of it. To do this, follow [this guide](https://sbnsoftware.github.io/sbndcode_wiki/commissioning/SBND_Commissioning_Get_Started.html), which basically amounts to setting up with 

```bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v09_88_00_04 -q e26:prof
```

and then using `mrb` to get the local build. The `setup.sh` script does all the setup, so then you can do 

```bash
make_tar_sbnd.sh epelaez_Aug_28_2024_v09_88_00_04.tar
```

to create the tarball (with the appropiate date, username, and version). Note that `setup.sh` only has to be run the first time, and any subsequent times you will only have to run `activate.sh` to activate the `sbndcode` version you built. Once the tarball is build, move it under

```
/pnfs/sbnd/resilient/users/${USER}/
```

so `config.xml` can find it. Before running jobs on the grid, we want to make sure that they work. With your `sbndcode` activated, you can test the job locally using

```bash
lar -c cafmakerjob_sbnd_systtools_and_fluxwgt.fcl -s /pnfs/sbnd/persistent/users/twester/sbnd/v09_88_00_04/detsys_debug_scaling/cv/reco2-cv-cdea-0ce8-68ed-25af.root -n -1
```

And if this works correctly, you can submit jobs to the grid by doing

```bash
project.py --xml config_cv.xml --stage flatcaf-cv --submit
```

where you can replace `config_cv.xml` with the appropiate `xml` file and `flatcaf-cv` with the corresponding stage for each file.