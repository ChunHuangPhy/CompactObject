Start Here!
===========

Welcome to **CompactObject**! This guide will help you get started after successfully installing the package and activating your virtual environment. Follow the steps below to achieve your scientific goals using the functionalities of the **CompactObject** package.

Prerequisites
-------------

- **Package Installation:** Ensure you have installed **CompactObject** successfully, see installation instructions here: `Installation <https://chunhuangphy.github.io/CompactObject/installation.html>`_ .
- **Virtual Environment:** Activate your virtual environment where **CompactObject** is installed , also see installation instructions here: `Installation <https://chunhuangphy.github.io/CompactObject/installation.html>`_ .

Depending on Your Scientific Goal
---------------------------------

1. Understanding Basic Concepts
   ----------------------------

   If you are new to **CompactObject** and wish to familiarize yourself with its basic concepts and introduction:

   - **Index Page:** Explore the `Index Page <https://chunhuangphy.github.io/CompactObject/index.html>`_ for an overview of the package.
   - **README:** Review the `README <https://github.com/ChunHuangPhy/CompactObject/blob/main/README.md>`_ on GitHub to understand the package’s purpose and features.

2. Generating Equations of State (EOS)
   ------------------------------------

   To use **CompactObject** for generating equations of state:

   - **EOS Generators Tutorial:** Start with the `EOS Generators Tutorial <https://chunhuangphy.github.io/CompactObject/test_EOSgenerators.html>`_, which showcases all currently integrated EOS computations in the package. The supported EOS models include:

     - **Relativistic Mean Field (RMF) Model**: Here we adopted the RMF model from the `(Lalazissis et al. 1997) <https://ui.adsabs.harvard.edu/link_gateway/1997PhRvC..55..540L/doi:10.1103/PhysRevC.55.540>`_ and `(Tolos et al. 2016) <https://arxiv.org/abs/1610.00919>`_.
     - **Strangeon Star Model**: Original model proposed in `(Xu 2003) <https://ui.adsabs.harvard.edu/abs/2003ApJ...596L..59X/abstract>`_
     - **Quark Star Model:** MIT Bag Model in `(Chodos 1974) <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.9.3471>`_
     - **Polytropic Model** proposed in `(Hebeler et al. 2013) <https://iopscience.iop.org/article/10.1088/0004-637X/773/1/11>`_, here we adapted the original code definition of "eos_core_pp" from NEoST software `(Raaijmakers et al. 2025) <https://joss.theoj.org/papers/10.21105/joss.06003.pdf>`_, see the page in `Polytrope model <https://github.com/xpsi-group/neost/blob/main/neost/eos/polytropes.py>`_.
     - **Speed of Sound Model** proposed in `(Grief et al. 2019) <https://academic.oup.com/mnras/article/485/4/5363/5374528>`_.
     - **Density-Dependent EOS Model** proposed in `(Typel et al. 2010) <https://journals.aps.org/prc/abstract/10.1103/PhysRevC.81.015803>`_.
     - *...and more*

3. Using the Tolman–Oppenheimer–Volkoff (TOV) Solver
   -------------------------------------------------

   To utilize the TOV solver integrated within **CompactObject**:

   - **TOV Solver Tutorial:** Refer to the `TOV Solver Tutorial <https://chunhuangphy.github.io/CompactObject/test_TOVsolver.html>`_. This tutorial demonstrates how to perform TOV computations using table-based EOS.
   - **Inputting an Array:** If you prefer to input an array instead of a table, refer back to the `EOS Generators Tutorial <https://chunhuangphy.github.io/CompactObject/test_EOSgenerators.html>`_. The quantities currently available are:

     - **Mass**
     - **Radius**
     - **Tidal Deformability**

4. Performing Bayesian Inference on the EOS
   -----------------------------------------

   **CompactObject** offers robust tools for Bayesian inference related to neutron star EOS. There are two primary scenarios:

   a. **Using Integrated EOS for Inference**

      If you wish to perform Bayesian inference using the EOS models already integrated into **CompactObject**, follow these examples:

      1. **RMF EOS Inference Pipeline:** Explore the `RMF EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Inference.html>`_.
      2. **density-dependent Hartree EOS Inference Pipeline:** Explore the `DDH EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_DDH.html>`_.
      3. **Strangeon EOS Inference Pipeline:** Check out the `Strangeon EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_Strangeon_EOS.html>`_.
      4. **MIT Bag Model EOS Inference Pipeline:** Check out the `MIT Bag Model EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_MITbag_EOS.html>`_.
      5. **Speed of Sound EOS Inference Pipeline:** Check out the `Speed of Sound EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_SpeedOfSound_EOS.html>`_.
      6. **Polytrope EOS Inference Pipeline:** Check out the `Polytrope EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Inference_polytrope.html>`_.
      
      If you computer don't have OpenMPI and only want to parallelize your code in your labtop or a server, we also provide options for
      implementing the MCMC sampling method, to do the inference, which could easily scale up in a personal computer. Except for using 
      UltraNest, we use EMCEE in this version of pipeline, here below are some of the examples (not all the equation of state we provided
      the EMCEE version, since some of them could be time-consuming and we default recommend you to use UltraNest to scale up computation
      in a High performance computation facility)

      1. **RMF EOS Inference MCMC Pipeline:**, see `RMF EOS EMCEE Inference Pipeline <https://chunhuangphy.github.io/CompactObject/EMCEE_RMF.html>`_
      2. **Strangeon EOS MCMC Inference Pipeline:** Check out the `Strangeon EOS EMCEE Inference Pipeline <https://chunhuangphy.github.io/CompactObject/EMCEE_strangeon.html>`_.
      3. **MIT Bag Model EOS MCMC Inference Pipeline:** Check out the `MIT Bag Model EOS EMCEE Inference Pipeline <https://chunhuangphy.github.io/CompactObject/EMCEE_MIT_bag.html>`_.
      4. **Polytrope EOS MCMC Inference Pipeline:** Check out the `Polytrope EOS EMCEE Inference Pipeline <https://chunhuangphy.github.io/CompactObject/EMCEE_polytrope.html>`_.
   
   b. **Using Your Own EOS for Inference**

      If you have developed your own EOS and wish to perform Bayesian inference:

      1. **Contributor Instructions:** Follow the `Instructions for Contributors <https://chunhuangphy.github.io/CompactObject/Contributor.html>`_ to define your EOS without conflicts and understand the necessary tests to ensure computational accuracy.
      2. **Constructing a New Pipeline:** Utilize the existing inference pipelines, such as the `RMF EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Inference.html>`_ and the `Strangeon EOS Inference Pipeline <https://chunhuangphy.github.io/CompactObject/test_Bayesian_inference_Strangeon_EOS.html>`_, as reference points to build your own inference workflow based on **CompactObject**.
      3. **Performance Considerations:** Ensure that the computation time for a single EOS is within the sub-second range to handle millions of samples required for Bayesian inference efficiently.

   **Contribution Invitation:**  
   Consider contributing your EOS to our community if you are using our package! This will boost the influence of your work and make you a collaborator on this project. We can also promote your results on our page. Please check the `Instructions for Contributors <https://chunhuangphy.github.io/CompactObject/Contributor.html>`_ if you wish to contribute.

5. Performing Bayesian Inference in Other Fields
   ----------------------------------------------

   If you aim to use the inference modules of **CompactObject** for Bayesian inference outside of EOS studies, such as in high-energy physics or other astrophysical fields:

   - **Example Notebook:** Refer to the `Sample Glitch Bayesian Inference Notebook <https://github.com/ChunHuangPhy/CompactObject/blob/main/Test_Case/Sample_glitchBayesian.ipynb>`_ for an example of how to apply **CompactObject** to different inference scenarios.

6. Advanced Usage: Modifying Solvers or APIs
   ------------------------------------------

   If the existing tutorials do not cover your specific needs, such as:

   - **Modifying the TOV Solver:** For example, changing to a double-fluid TOV solver.
   - **Investigating Additional EOS-Related Quantities:** Beyond pressure and density used in TOV solving.

   You can:

   - **API Documentation:** Consult the detailed API documentation to identify and utilize the functions you need.
   - **GitHub Repository:** Visit our `GitHub Page <https://github.com/ChunHuangPhy/CompactObject/tree/main>`_ to explore the source code and understand the internal workings of **CompactObject**.

   **Note:**  
   The key functions, including the **TOVsolver module** code, form the foundation for all inference and computations. Modifying them may require substantial code restructuring. Therefore, if you clone the repository locally and need to make changes to the **TOVsolver**, proceed with caution to avoid extensive code rebuilding. Only original contributors have permission to modify the **TOVsolver** module code, and such commits must undergo detailed reviews by the original members of the project.

Computation requirements and citation
-------------------------------------

Bayesian inference could be quite computation expensive especially when the EOS parameters dimensions increase, So we provide a detailed `Computation guide <https://chunhuangphy.github.io/CompactObject/cluster.html>`_ on estimating how many computation resources are required if you would like to use the integrated EOS functions in our package and what is the computation speed requirement if you would like to implement a self-defined EOS functions. 

If you used our pipeline and code in your study, please consider to cite us, there is the list of publication you could cite in `Publications page <https://chunhuangphy.github.io/CompactObject/publications.html>`_


Acknowledgements
----------------

We welcome feedback and contributions to expand the functionalities of **CompactObject**. Your support helps enhance the tool for the entire research community.

Contact
-------

For inquiries, contributions, or to be featured in our publications list, please contact us at `chun.h@wustl.edu <mailto:chun.h@wustl.edu>`_.
