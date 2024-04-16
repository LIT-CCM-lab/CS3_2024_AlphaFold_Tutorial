# Strasbourg Summer School in Chemoinformatics 2024
## Multistate modelling using LIT-AlphaFold

In this tutorial we will present you how to predict protein structures using AlphaFold, and how to evalaute the results, as well as some of the state of the art techniques for multistate modelling.

All calculations are performed using the LIT-AlphaFold pipeline, which simplifies input generation and customization compared to the original AlphaFold implementation.

To participate to the tutorial click on this link to the Colab Notebook:

[Link to the tutorial](https://colab.research.google.com/github/LIT-CCM-lab/CS3_2024_AlphaFold_Tutorial/blob/main/CS3_AlphaFold_tutorial.ipynb)

## How to
Here we include a few general notes that might be usefull during the tutorial

### Google Colab
> Colab is a hosted Jupyter Notebook service that requires no setup to use and provides free access to computing resources, including GPUs and TPUs.

A Google account is required to access the free computing resources of Colab.

By using Colab you do not need to install anything on your machine, therefore there should be no issue regardless of the used machine.

### py3Dmol
In the tutorial the 3D structure of the generated proteins will be visualized using py3Dmol.

To navigate the images the following command might be usefull:
* Rotation: left-click and move the mouse
* Translation: center-click or ctrl+left-click, and move the mouse
* Zoom: Scrool wheel, or shift+left-click and move the mouse
* Show/hide atom labels: single left-click (this feature is not universal and is specific for this tutorial)

### LIT-AlphaFold
LIT-AlphaFold stores protein information (multiple sequences alignements and templates) as compressed pickle file (*.pkl.bz2*). These file are generally read by the LIT-AlphaFold prediction pipeline to perform calculations, or can be used in python scripts using the package *litaf*.

For more information please visit the [LIT-AlphaFold GitHub repository](https://github.com/LIT-CCM-lab/LIT-AlphaFold).