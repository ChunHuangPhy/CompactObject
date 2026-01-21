<!-- ABOUT THE PROJECT -->
# CompactObject
![Logo](logo.png)
**CompactObject** is an open-source package designed to perform Bayesian inference on neutron star equation of state (EOS) constraints. It offers a comprehensive workflow that integrates astrophysical observations and nuclear measurements to explore the interior composition of neutron stars. The package is built to be user-friendly, easily extendable, and thoroughly documented, making it an essential tool for researchers in nuclear astrophysics.

We have a detailed documentation please check here: [CompactObject Package Website](https://chunhuangphy.github.io/CompactObject/)

**Please star this repository if you find it is helpful!**

## Table of Contents

- [Core Functionality](#core-functionality)
- [Papers](#papers)
- [Sample Citation](#sample-citation)
- [Includes](#includes)
- [Installation](#installation)
- [Physics Notations](#physics-notations)
- [License](#license)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

## Core Functionality

CompactObject offers a range of functionalities essential for constraining the EOS of neutron stars:

1. **Equation of State (EOS) Generation**
    - Utilizes multiple physics/meta models, including Relativistic Mean Field (RMF), strange star, quark star, polytrope, and speed of sound models, to generate neutron star EOS.
    - [EOSgenerators](https://github.com/ChunHuangPhy/EoS_inference/blob/main/EOSgenerators) Package

2. **Tolman-Oppenheimer-Volkoff (TOV) Solver**
    - Solves the TOV equations for a spherically symmetric compact object based on a given neutron star EOS.
    - [TOVsolver](https://github.com/ChunHuangPhy/EoS_inference/blob/main/TOVsolver) Package

3. **Bayesian Inference Workflow**
    - Implements neutron star EOS inference using Nested Sampling.
    - Provide options for single machine users using MCMC sampling by emcee.
    - Integrates constraints from nuclear experiments, neutron star mass and/or radius observations (from X-ray timing and/or radio timing), and tidal measurements from gravitational wave detections.
    - [InferenceWorkflow](https://github.com/ChunHuangPhy/EoS_inference/blob/main/InferenceWorkflow) Package

## Papers

Please consider cite the following papers if you use CompactObject in your research:

Package paper:

1. **Huang, C., Malik, T., Cartaxo, J., et al.** (2024). *CompactObject: An open-source Python package for full-scope neutron star equation of state inference* [arXiv:2411.14615](https://arxiv.org/abs/2411.14615)

2. **Cartaxo, J., Huang, C., Malik, T., Sourav, S. , Yuan, W., Zhou, T., Liu, X., Providência, C.**, “Covariant Energy Density Functionals for Neutron Star Matter Equation of State Modeling: Cross-Comparison Analysis Using \texttt{CompactObject}”, *The Astrophysical Journal Supplement Series*, Volume 282, Number 2, [https://iopscience.iop.org/article/10.3847/1538-4365/ae2310](https://iopscience.iop.org/article/10.3847/1538-4365/ae2310)


Other papers:

1. **Huang, C., Raaijmakers, G., Watts, A. L., Tolos, L., & Providência, C.** (2024). *Constraining fundamental nuclear physics parameters using neutron star mass-radius measurements I: Nucleonic models*. *Monthly Notices of the Royal Astronomical Society*, 529. [DOI:10.1093/mnras/stae844](https://academic.oup.com/mnras/article/529/4/4650/7634362)

3. **Huang, C., & Zheng, X.-P.** (2024). *Bayesian Insights into post-Glitch Dynamics: Model comparison and parameter constraint from decades-long observation data of the Crab pulsar*, *Monthly Notices of the Royal Astronomical Society*, Volume 542, Issue 4, October 2025, Pages 3198–3205, [https://doi.org/10.1093/mnras/staf1415](https://doi.org/10.1093/mnras/staf1415).

5. **Huang, C., Sourav, S.**, *Constraining first-order phase transition inside neutron stars with application of Bayesian techniques on PSR J0437-4715 NICER data*, (2025), *The Astrophysical Journal*, Volume 983, Number 1, [https://iopscience.iop.org/article/10.3847/1538-4357/adbca6](https://iopscience.iop.org/article/10.3847/1538-4357/adbb67/meta) 
   
7.  **Huang, C.**, *Equation of State Independent Determination on the Radius of a 1.4 $M_{\odot}$ Neutron Star Using Mass-Radius Measurements* (2025), *The Astrophysical Journal Letters*, Volume 978, Number 1, [https://iopscience.iop.org/article/10.3847/2041-8213/ad9f3c](https://iopscience.iop.org/article/10.3847/2041-8213/ad9f3c).
   

## Sample Citation


"The inference conducted here relies on the framework in the *CompactObject* \cite{CompactObject} package\footnote{https://chunhuangphy.github.io/CompactObject/}. This is an open-source, comprehensive package designed to implement Bayesian constraints on the neutron star EOS. Other works based on this package include ..."


## Includes

CompactObject includes the following components to facilitate neutron star EOS inference analysis:

1. **EOS Output Validation Routine**
    - Compute various type of EOS by different models
    - Checks the validity of EOS inputs.

2. **Mass, Radius, and Tidal Deformability Calculator**
    - Returns mass, radius, tidal deformability, by solve TOV equation, and computes the corresponding speed of sound.

3. **Sample TOV Solver Notebook**
    - [Test_TOVsolver.ipynb](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_TOVsolver.ipynb)
    - Demonstrates how to solve the TOV equation with a given EOS.

4. **Sample EOS Generators Notebook**
    - [test_EOSgenerators.ipynb](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_EOSgenerators.ipynb)
    - Showcases all integrated EOS computations, including:
        - Polytrope
        - Speed of Sound Model
        - RMF Model
        - Strange Star Model
        - Quark Star Model

5. **Sample Analysis and Tutorial Notebook**
    - [test_Inference.ipynb](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_Inference.ipynb)
    - Demonstrates the entire pipeline of Bayesian inference using supported EOS models, constructing priors and likelihoods, and the types of likelihoods supported in this project. Also provide a MCMC based emcee example for people don't have access to High Performance Computer. This is specifically focus on the RMF EOS,
    However, **please check this notebook before all other notebook,**
    **since here we showcase all the likelihood** 
    - Other Inference pipline that using different EOS are
        - [MIT bag inference](https://github.com/ChunHuangPhy/CompactObject/blob/main/Test_Case/test_Bayesian_inference_MITbag_EOS.ipynb)
        - [Strangeon Star inference](https://github.com/ChunHuangPhy/CompactObject/blob/main/Test_Case/test_Bayesian_inference_Strangeon_EOS.ipynb)
        - [Polytrope inference](https://github.com/ChunHuangPhy/CompactObject/blob/main/Test_Case/test_Inference_polytrope.ipynb) 

> **Note:** Please review these notebooks before starting your own project to familiarize yourself with the coding routines.

## Installation

Below are the commands to install and update the CompactObject package, along with a link to PyPI.

### [PyPI - CompactObject-TOV](https://pypi.org/project/CompactObject-TOV/)

1. **Install the Package**
    ```sh
    pip install CompactObject-TOV
    ```

2. **Update the Package**
    ```sh
    pip install CompactObject-TOV --upgrade
    ```

### Importing the Package

- **For EOS Computation:**
    ```python
    import EOSgenerators
    ```

- **For TOV Solver:**
    ```python
    import TOVsolver
    ```

- **For Bayesian Inference:**
    ```python
    import InferenceWorkflow
    ```

## Physics Notations

1. **Units:**
    - **CGS Units** are used throughout the package.
    - **Pressure (P):** erg/cm³
    - **Energy Density (ρ):** g/cm³
    - **Default Unit system** any input in this package is follow the cgs unit

2. **Unit Conventions:**
    - Follow the unit system defined in the [Unit Conversion](https://chunhuangphy.github.io/CompactObject/UnitConventionForDeveloper.html) notebook.
    - **Basic Rules:**
        1. If using a quantity with a unit, multiply by the unit during computations.
        2. If plotting or demonstrating values in a specific unit, divide accordingly.
    - Detailed guidelines can be found in the [Unit Conversion](https://chunhuangphy.github.io/CompactObject/UnitConventionForDeveloper.html) documentation.

## License

Distributed under the GNU GPL-3.0 license. See `LICENSE.txt` for more information.

<p align="right">(<a href="#compactobject">back to top</a>)</p>

## Contact

* **Chun Huang** - [chun.h@wustl.edu](mailto:chun.h@wustl.edu)

**Documentation Link:** [CompactObject Documentation](https://chunhuangphy.github.io/CompactObject/)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.8145167.svg)](http://dx.doi.org/10.5281/zenodo.8145167)

<p align="right">(<a href="#compactobject">back to top</a>)</p>

## Acknowledgments

We would like to acknowledge the support of the Code/Astro workshop, which was instrumental in the development of this project. Special thanks to Professor Anna Watts, Dr. Geert Raaijmakers, and Jeannie Kuijper for their assistance with coding and providing foundational strategies to address the research challenges.

* [Code Astro](https://github.com/semaphoreP/codeastro)
* [Choose an Open Source License](https://choosealicense.com)
* [Read Me Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#compactobject">back to top</a>)</p>
```
