
## Post-processing tools for the NEMO-OASIS-WRF configurations


### Installation

To install the tools, you need to have all the required python libraries. If you do not have [Anaconda](https://conda.io/projects/conda/en/latest/) or [Miniconda](https://conda.io/en/latest/miniconda.html) installed, you should first download one of the two and install it in your working directory. We recommend to use Miniconda on HPCs as it is lighter and install only required packages.

To install Miniconda3 just run the following commands:
```bash 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Once Miniconda is installed, you can create the conda environment required for running **now-posprocess** using:
```bash 
conda env create --file requirements.yml
```

Then, just run:
```bash
python setup.py install
```