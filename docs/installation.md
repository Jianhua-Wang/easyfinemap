# Installation

easyfinemap is a Python package that can be installed directly using pip. However, since easyfinemap requires additional fine-mapping software to compute posterior probabilities, such as FINEMAP and PAINTOR, you will need to install these software as well. Users have the option to install them manually or clone the easyfinemap conda environment. Here are the specific steps:

1. Download the environment configuration file:
   ```bash
   wget https://github.com/Jianhua-Wang/easyfinemap/blob/main/environment.yml
   ```
2. Create a conda environment:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the conda environment:
   ```bash
   conda activate easyfinemap
   ```
Once the conda environment is activated, you will have easyfinemap and the required fine-mapping software installed and ready to use.