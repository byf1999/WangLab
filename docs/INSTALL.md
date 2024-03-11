# INSTALL Guide For WangLab

Please check the following instructions to complete your installation.

## Install through PyPI

The easiest way to install WangLab in through PYPI. 

```bash
pip install WangLab
```

**Note**: Bowtie/MACS do not provide releases for **Windows** systems by the time this document is written, so users on __Windows__ OS can not use full function in these modules.Install from source

To install WangLab from source, users first download and extract .zip file from Github. Then, place the whole file folder __WangLab__ in __WangLab-main__ into __[ANACONDA_PATH/lib/python3.x/site-packages__. In the meanwhile, place the script file __wanglab__ in __WangLab-main/bin__ into __[ANACONDA_PATH]/bin__. If you want to use WangLab in a specific python environment, substitute __ANACONDA_PATH__ to __ANACONDA_PATH/envs/ENV_NAME__.

Now users can call commands in WangLab in any path under specific python environment.

__Note__: To use TIS/ChIP_seq modules, you need to mannually install __Bowtie__ and __MACS__ through:

```bash
# install bowtie
# conda:
conda install bowtie -c bioconda

# apt
sudo apt install bowtie
```

```bash
# install MACS through pip
pip install macs3
```

### 