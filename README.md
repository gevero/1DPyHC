# 1DPyHC

<img src="https://github.com/gevero/1DPyHC/blob/master/images/bragg.png" width="600">

A python code to calculate the optical properties of 1D Photonic Crystals. The code follows the guidelines outlined in:

Yeh P, Yariv A, Hong C-S (1977) **Electromagnetic propagation in periodic stratified media. I. General theory.** *J Opt Soc Am 67:423–438*
[doi: 10.1364/JOSA.67.000423](http://dx.doi.org/10.1364/JOSA.67.000423)

For the time being we are only dealing with TE modes (which are specially interesting for **Bloch Surface Waves**). TM mode calculations will be added along the way

## Installation

Installing **1DPyHC** should be fairly easy. First of all you need a python distribution on your system. The simplest thing to do if you are using **Windows** or **OSX** would be to install [Anaconda](https://store.continuum.io/cshop/anaconda/), a beautiful, free, and easy to use python distribution: the relevant installations instructions are found [here](http://docs.continuum.io/anaconda/install.html). If you are using **Linux** you probably already know how to install python on your system, nevertheless my advice is to also install [The IPython Notebook](http://ipython.org/notebook.html), which is already bundled in Anaconda.

The second step would be installing **1DPyHC** itself. If you are familiar with [git](http://git-scm.com/) and [github](https://github.com/) you can simply clone the repository, otherwise just [download](https://github.com/gevero/py-matrix/archive/master.zip) the zipped version of the repo and unpack it wherever you like.

## Usage

The best thing to do for the moment is to start from the **.ipynb** files in the [examples](https://github.com/gevero/1DPyHC/tree/master/examples) folder. You can load them in your local IPython Notebook instance: they should give you a fair idea about how to proceed for your calculations. Each function in the code features a detailed documentation easily accessible with the `Shift-Tab` tool-tip shortcut from the IPython Notebook interface.
