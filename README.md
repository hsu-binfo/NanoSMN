# NanoSMN



## Usage

### Step 0: Install conda and Snakemake

Please follow the instruction from [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install conda and Snakemake.



### Step 1: Clone workflow

To use NanoSMN, you need to clone the packages to your local machine.

```bash
git clone https://github.com/hsu-binfo/NanoSMN.git
```



### Step 2: Configure workflow

The configuration file is config.yml in the root of the project. You can configure the workflow according to your needs by editing the file. Please place all the samples in path `data/samples` relative to project root. 



### Step 3: Execute workflow

Test your configuration by performing a dry-run via:

```bash
snakemake -n
```



Execute the workflow locally via:

```bash
snakemake --use-conda --cores N
```

`N` is stand for the number of CPU cores you want to use.

