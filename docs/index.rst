.. LineageGRN documentation master file, created by
   sphinx-quickstart on Tue Oct 29 17:29:50 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LineageGRN documentation
================================================

LineageGRN
----------


LineageGRN aims to infer dynamic gene regulatory networks (GRNs) at a real-time scale using paired **single-cell lineage tracing** and **scRNA-seq** data obtained from the same cell, along with **scATAC-seq** datasets.
The inferred GRNs are dependency and linked by a cell fate map representing the developmental trajectory from a progenitor cell to their descendants (sampled cells), which is modeled by a Bayesian network.

.. image:: _static/algorithm.png
   :alt: Overview of LineageGRN
   :align: center
   :scale: 70 %

Features
---------
The currently available features include:

- Inferring time-scaled **cell lineage**

- Inferring **dynamic gene regulatory networks**

- Revealing how inferred GRNs **reconfigure** along the cell lineage

- Identifying **key regulatory genes** driving **cell differentiation**

- Identifying **key regulatory genes** steering **cell fate bias**

- Identifying specific and constitutive **regulatory interactions**.



Reference
---------

Deciphering dynamic gene regulatory network underlying cell fate decisions in cell lineage

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Tutorials/index
   API/index
