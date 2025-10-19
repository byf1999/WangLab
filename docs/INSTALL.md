# INSTALL Guide For WangLab

Please check the following instructions to complete your installation.

## Install through PyPI

The easiest way to install WangLab in through PYPI. 

```bash
pip install WangLab
```

**Note**: Bowtie/samtools/MACS do not provide releases for **Windows** systems by the time this document is written, so users on __Windows__ OS can not use full function in these modules.

Linux OS users need to manually install these three packages.

```bash
# install bowtie
# conda:
conda install bowtie -c bioconda

# apt
sudo apt install bowtie

```

``` bash
# install samtools
conda install -c bioconda samtools
```

```bash
# install MACS through pip
pip install macs3
```

## Install from source

### For Linux users

To install WangLab from source, users first download and extract .zip file from Github. Then, place the whole file folder `WangLab` in `WangLab-main` into `[ANACONDA_PATH/lib/python3.x/site-packages`. In the meanwhile, place the script file `wanglab` in `WangLab-main/bin` into `[ANACONDA_PATH]/bin`. If you want to use WangLab in a specific python environment, substitute `ANACONDA_PATH` to `ANACONDA_PATH/envs/ENV_NAME`.

Finally, install required packages through `pip install -r WangLab-main/requirements.txt`

Now users can call commands in WangLab in any path under specific python environment.

__Note__: To use TIS/ChIP_seq modules, you need to mannually install __Bowtie__, __samtools__ and  __MACS__ through:

```bash
# install bowtie
# conda:
conda install bowtie -c bioconda

# apt
sudo apt install bowtie

```

``` bash
# install samtools
conda install -c bioconda samtools
```

```bash
# install MACS through pip
pip install macs3
```

### For Windows users

To install __WangLab__ from source, users first download and extract .zip file from Github. Then, place the whole file folder `WangLab` in `WangLab-main`into `[ANACONDA_PATH/lib/python3.x/site-packages`. In the meanwhile, place the script files `wanglab` and `wanglab.bat` in `WangLab-main/bin` into `[ANACONDA_PATH]/scripts`. If you want to use WangLab in a specific python environment, substitute `ANACONDA_PATH` to `ANACONDA_PATH/envs/ENV_NAME`.

Finally, install required packages through `pip install -r WangLab-main/requirements.txt`

Now users can call commands in WangLab in any path under specific python environment.

__Note__: On Windows OS, users are not able to use functions in TIS/ChIP_seq modules.