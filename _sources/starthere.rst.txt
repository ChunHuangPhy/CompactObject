Start Here!
==========

Welcome to **CompactObject**! This guide will help you get started after successfully installing the package and activating your virtual environment. Follow the steps below to achieve your scientific goals using the functionalities of the **CompactObject** package.

Prerequisites
-------------
- **Package Installation:** Ensure you have installed **CompactObject** successfully.
- **Virtual Environment:** Activate your virtual environment where **CompactObject** is installed.

Depends on Your Scientific Goal
-------------------------------

 1. **Understanding Basic Concepts**

If you are new to **CompactObject** and wish to familiarize yourself with its basic concepts and introduction:

- **Index Page:** Explore the `Index Page <https://chunhuangphy.github.io/CompactObject/index.html>`_ for an overview of the package.
- **README:** Review the `README <https://github.com/ChunHuangPhy/CompactObject/blob/main/README.md>`_ on GitHub to understand the package’s purpose and features.

 2. **Generating Equations of State (EOS)**

To use **CompactObject** for generating equations of state:

- **EOS Generators Tutorial:** Start with the `EOS Generators Tutorial <https://chunhuangphy.github.io/CompactObject/test_EOSgenerators.html>`_ which showcases all currently integrated EOS computations in the package. The supported EOS models include:
  - **Relativistic Mean Field (RMF) Model**
  - **Strangeon Star Model**
  - **Quark Star Model:** MIT Bag Model
  - **Polytropic Model**
  - **Speed of Sound Model**
  - **Density-Dependent EOS Model**
  - *...and more*

 3. **Using the Tolman–Oppenheimer–Volkoff (TOV) Solver**

To utilize the TOV solver integrated within **CompactObject**:

- **TOV Solver Tutorial:** Refer to the `TOV Solver Tutorial <https://chunhuangphy.github.io/CompactObject/test_TOVsolver.html>`_. This tutorial demonstrates how to perform TOV computations using table-based EOS.
- **Inputting an Array:** If you prefer to input an array instead of a table, refer back to the `EOS Generators Tutorial <https://chunhuangphy.github.io/CompactObject/test_EOSgenerators.html>`_. The quantities currently available are:
  - **Mass**
  - **Radius**
  - **Tidal Deformability**

 4. **Performing Bayesian Inference on the EOS**

**CompactObject** offers robust tools for Bayesian inference related to neutron star EOS. There are two primary scenarios:

 a. **Using Integrated EOS for Inference**

If you wish to perform Bayesian inference using the EOS models already integrated into **CompactObject**, follow these examples:

1. **RMF EOS Inference Pipeline:** Explore the `RMF EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Inference.html>`_.
2. **Strangeon EOS Inference Pipeline:** Check out the `Strangeon EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_Strangeon_EOS.html>`_.
3. *...and more examples*

 b. **Using Your Own EOS for Inference**

If you have developed your own EOS and wish to perform Bayesian inference:

1. **Contributor Instructions:** Follow the `Instructions for Contributors <https://chunhuangphy.github.io/CompactObject/Contributor.html>`_ to define your EOS without conflicts and understand the necessary tests to ensure computational accuracy.
2. **Constructing a New Pipeline:** Utilize the existing inference pipelines, such as the `RMF EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Inference.html>`_ and `Strangeon EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_Strangeon_EOS.html>`_, as reference points to build your own inference workflow based on **CompactObject**.
3. **Performance Considerations:** Ensure that the computation time for a single EOS is within the sub-second range to handle millions of samples required for Bayesian inference efficiently.

**Contribution Invitation:**  
Consider contributing your EOS to our community if you are using our package! This will boost the influence of your work and make you a collaborator on this project. We can also promote your results on our page. Please check the `Instructions for Contributors <https://chunhuangphy.github.io/CompactObject/Contributor.html>`_ if you wish to contribute.

 5. **Performing Bayesian Inference in Other Fields**

If you aim to use the inference modules of **CompactObject** for Bayesian inference outside of EOS studies, such as in high-energy physics or other astrophysical fields:

- **Example Notebook:** Refer to the `Sample Glitch Bayesian Inference Notebook <https://github.com/ChunHuangPhy/CompactObject/blob/main/Test_Case/Sample_glitchBayesian.ipynb>`_ for an example of how to apply **CompactObject** to different inference scenarios.

 6. **Advanced Usage: Modifying Solvers or APIs**

If the existing tutorials do not cover your specific needs, such as:

- **Modifying the TOV Solver:** For example, changing to a double-fluid TOV solver.
- **Investigating Additional EOS-Related Quantities:** Beyond pressure and density used in TOV solving.

You can:

- **API Documentation:** Consult the detailed API documentation in this documentation to identify and utilize the functions you need.
- **GitHub Repository:** Visit our `GitHub Page <https://github.com/ChunHuangPhy/CompactObject/tree/main>`_ to explore the source code and understand the internal workings of **CompactObject**.

**Note:**  
The key functions, including the **TOVsolver module** code, form the foundation for all inference and computations. Modifying them may require substantial code restructuring. Therefore, if you clone the repository locally and need to make changes to the **TOVsolver**, proceed with caution to avoid extensive code rebuilding. Only original contributors have permission to modify the **TOVsolver** module code, and such commits must undergo detailed reviews by the original members of the project.

Acknowledgements
----------------

We welcome feedback and contributions to expand the functionalities of **CompactObject**. Your support helps enhance the tool for the entire research community.

Contact
-------

For inquiries, contributions, or to be featured in our publications list, please contact us at `chun.h@wustl.edu <mailto:chun.h@wustl.edu>`_.
