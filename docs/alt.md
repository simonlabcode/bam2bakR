## Alternative Strategies for Running bam2bakR
The old way of running bam2bakR required cloning the entire bam2bakR repo locally. This makes updating bam2bakR a bit of a hassle, as it is non-trivial to pull any updates from the bam2bakR repo if you edited the config file. That being said, clonging the full repo can be useful if you need to make changes to the workflow due to idiosyncracies in your particular data. Therefore, this section discusses three different ways bam2bakR can be run via cloning the entire repo. These instructions are also here for historical reasons, as versions 1.0.0 and earlier of bam2bakR were only compatible with these strategies.

The common denominator of all of these usage routes is, of course, cloning the repo. Clone the bam2bakR repository to wherever you would like on your system. You will eventually be navigating to this repo directory in the terminal and running Snakemake from inside the directory, so make sure your chosen location is conducive to this. Navigate to the directory in the terminal and run:

``` bash
$ git clone https://github.com/simonlabcode/bam2bakR.git
$ cd bam2bakR
```
You should be in the bam2bakR repo directory now!

Step 1 of the Snakedeploy route (i.e., installing mamba/conda) is also the same for the Activating --use-conda and Running all rules inside one environment routes.

### Activating --use-conda

Version 1.0.0 and later of bam2bakR is compatible with Snakemake's [--use-conda option](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). This will cause Snakemake to automatically create and activate conda environments for each step of the workflow to run inside. If you want to use this functionality, you can start by creating a simple conda environment that contains snakemake, as such:

``` bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake
```

You would then run the pipeline with the `snakemake` environment activated with:

``` bash
snakemake cores all --use-conda
```

where `cores all` is a convenient way to tell Snakemake to make use of all available cpus (all can be replaced with an explicit number as was shown in the installation/pipeline running instructions above).

### Running all rules inside one environment

All versions of bam2bakR are compatible with creating a single conda environment that contains all of bam2bakR's dependencies and running bam2bakR from inside of this environment. 

Inside of the cloned bam2bakR repo, you will find a file named `pipeline_env.yaml`. This is a YAML file with the list of exact dependencies that need to be installed to run the Snakemake workflow. Luckily, installation of everything can be completed automatically thanks to Mamba/Conda! If you followed the nstructions in [Running bam2bakR](deploy.md) for installing Mamba, just run:

``` bash
$ conda activate base
$ mamba env create --file pipeline_env.yaml
```
The first line activates the so-called `base` conda environment, which is what allows you to call `mamba` in the next line. This will create a new environment called complete_pipeline. If you don't like that environment name, you can change it by editing the first line of `pipeline_env.yaml`, which looks like:
``` yaml
name: complete_pipeline
```
by replacing `complete_pipeline` with the name of your choice! 

To make sure we are all on the same page, an environment is a collection of installed software tucked away in its own directory, separate from everything else installed on your computer. When installing software with Mamba/Conda, I highly suggest installing all of the software you need for a particular task into an isolated environment like this. What's great about it is that you can have multiple versions of any given software installed at the same time! 

Now, whenever you want to run the bam2bakR pipeline, you can just run:

``` bash
$ conda activate complete_pipeline
```
replacing `complete_pipeline` with whatever name you ended up giving the environment, and voila, all of the dependencies are ready to be called upon.

Once you edited the config as described in [Running bam2bakR](deploy.md), all you need to do to run bam2bakR is activate the pipeline environment and call snakemake from the top of the bam2bakR directory as follows (replacing `complete_pipeline` with whatever you named the environment and `4` replaced with however many cpus you would like to use):
``` bash
$ conda activate complete_pipeline
$ snakemake --cores 4
```

### Running bam2bakR in a Docker Container

If you are having trouble installing mamba/conda or creating the required conda environment, you may like to try running bam2bakR inside of a [Docker](https://docs.docker.com/get-started/) container instead. This will automatically install mambaforge, create bam2bakR's conda environment, and even reproduce the exact operating system that I have predominantly tested bam2bakR in. The steps to get up and running with this alternative installation route are (after cloning the bam2bakR repo):

1. [Install Docker](https://docs.docker.com/get-docker/). If you are running bam2bakR on a system where you are not the admin (e.g., a shared computing cluster), then you will want to use [Singularity](https://apptainer.org/) (recently changed name to Apptainer) instead. If using Apptainer, just replace `docker` in any code snippet that follows with `apptainer` and you should be good to go!
1. Navigate to the directory where you want to run bam2bakR (i.e., the directory that you will call `snakemake --cores all` or something similar to start the pipeline)
1. Pull the bam2bakr Docker image from Docker Hub by running: `docker pull isaacvock/bam2bakr:v1.0.1`
1. Run the bam2bakr image with a [bind mount](https://docs.docker.com/get-started/06_bind_mounts/). This allows you to have access to all of the files inside your working directory while inside the image's software environment (i.e., Ubuntu + bam2bakR conda environment): 
   ``` bash
   docker run -it --mount type=bind,src="$(pwd)", target=/pipeline isaacvock/bam2bakr:v1.0.1
   ```
1. Activate the conda environment: `conda activate docker_pipeline`
1. Move into the working directory: `cd pipeline`
1. Run bam2bakR: `snakemake --cores all`

When you want to exit the Docker container, hit "Ctrl" + "d". Now everytime you want to run bam2bakR again, you just have to repeat steps 4 through 7.