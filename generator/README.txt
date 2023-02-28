How to setup a conda environment for synthetic data generation
==============================================================

We strongly recommend the use of conda environments for python
development. Conda is able to manage package dependencies much
better than pip. Conda also has a replacement for pip.

Here are instructions based on Miniconda, a "lite" version of
the full Anaconda installation.

-------------------------------------------------------------

Download the latest Miniconda package:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Compute the SHA256 hash of the downloaded file:

    sha256sum Miniconda3-latest-Linux-x86_64.sh

Check the hash against those on this page:

    https://conda.io/en/latest/miniconda_hashes.html

Install Miniconda:

    bash Miniconda3-latest-Linux-x86_64.sh
    accept the default install location and answer yes to all the questions

Activate the modifications to your .bashrc file by closing the
terminal window and starting a new terminal window.

Test the installation by listing the installed conda packages:

    conda list

If your system cannot find the 'conda' executable, then something went
wrong with the .bashrc modifications to your PATH environment variable.

Update the installation:

    conda update conda

Now setup a conda environment for the synthetic data generation code.
The environment we create will be called 'dev'.

    conda create --name dev
    <accept defaults when prompted>

Activate the 'dev' environment so that the packages required by the synthetic
data generation code can be installed into it.

    conda activate dev

Install the required packages:

    conda install -c anaconda jupyter
    conda install -c anaconda scipy
    conda install -c anaconda statsmodels
    conda install -c conda-forge matplotlib

That should be all you need to run the synthetic data generation software.

When finished, deactivate the 'dev' environment:

    conda deactivate

Delete the Miniconda installer if you want:

    rm Miniconda3-latest-Linux-x86_64.sh

To run the notebook, open a terminal window and change directories to the
location of "IntegratedModel.ipynb".

Start Jupyter with this command:

    jupyter notebook

It may take a few seconds for Jupyter to start, especially the first time
you run it. A web browser window will appear once it starts.

Click on "IntegratedModel_NETSS.ipynb" from the file list.
The notebook should start.

Scroll down to section 5, "Optinal Inputs".

       set the NETSS_DIR path
       set the OUTPUT_DIR path (create that directory if it does not exist)

Then enter a jurisdiction, condition code, and set any other variables in
sections 1-5 at the top of the notebook.

After setting the values of the variables, run the notebook as follows:

       From the "Kernel" menu, select "Restart & Clear Output",
       then "Restart & Run All".




