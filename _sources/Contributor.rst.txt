Call for Contributions
======================

Introduction
------------

If you are interested in contributing to the **CompactObject** project, please send an email to `chun.h@wustl.edu <mailto:chun.h@wustl.edu>`_ to request collaboration. You can also open an issue on our GitHub repository at `Issue page <https://github.com/ChunHuangPhy/CompactObject/issues>`_ so that the community can track the status of collaborations.

We are always eager to expand our scope to include the most equation of state (EOS) computations currently available in the community, as well as the most relevant constraints. If you have your own EOS and wish to contribute, please follow the steps outlined below.

Contributing an Equation of State (EOS)
---------------------------------------

1. **Open an Issue**  
   Visit our GitHub issues page at `Issue page <https://github.com/ChunHuangPhy/CompactObject/issues>`_ and open a new issue to suggest the EOS you would like to contribute.

2. **Request Collaboration**  
   Send an email to `chun.h@wustl.edu` to request collaborator access to the repository.

3. **Commit the EOS Script**  
   - Create a `.py` script that computes the EOS. The script should:
     - Input the EOS parameter values.
     - Output the pressure and density in CGS units.
   - For unit conventions, please refer to the `UnitConventionForDeveloper <https://chunhuangphy.github.io/CompactObject/UnitConventionForDeveloper.html>`_.
   - Name the main compute function as `compute_EOS`.
   - Name the file as `(YourEOSName)_EOS.py`.

4. **Create a Demonstration Notebook**  
   Develop a Jupyter notebook that demonstrates:
   - How to compute the EOS.
   - Defining your EOS parameters.
   - Integrating your EOS with our default TOV solver (`main.outputMR`).
   - Ensuring the TOV solver correctly solves the EOS and outputs a mass-radius curve plot.

Contributing Inference Components
---------------------------------

If you have developed your own likelihood, prior, sampler, or any component related to inference that is suitable for addition to our functionalities, we highly recommend contacting us. Follow these steps to contribute:

1. **Open an Issue**  
   Visit our GitHub issues page at `Issue page <https://github.com/ChunHuangPhy/CompactObject/issues>`_ and open a new issue to suggest changes to the inference components.

2. **Request Collaboration**  
   Send an email to `chun.h@wustl.edu` to request collaborator access to the repository.

3. **Commit the Inference Script**  
   - Create a `.py` script within the `InferenceWorkflow` directory.
   - Name the script after the function you are contributing.

4. **Create a Demonstration Notebook**  
   Develop a Jupyter notebook that demonstrates:
   - How to use your prior, likelihood, or sampler.
   - Include sampling outputs to ensure compatibility with the current structure of the package.

Guidelines for Key Functions and Permissions
--------------------------------------------

The **TOVsolver** module is fundamental to all inference and computation within **CompactObject**. Modifying this module requires careful consideration, as changes may necessitate significant code restructuring.

- **Permissions:**  
  Only original contributors have permission to modify the **TOVsolver** code. All commits related to the **TOVsolver** module require detailed reviews by the original members of the package.

- **Local Changes:**  
  If you clone the repository locally and need to make changes to the **TOVsolver**, proceed with caution to avoid extensive code rebuilding.

Acknowledgements
----------------

We welcome feedback and contributions to expand the functionalities of **CompactObject**. Your support helps enhance the tool for the entire research community.

Contact
-------

For inquiries, contributions, or to be featured in our publications list, please contact us at `chun.h@wustl.edu <mailto:chun.h@wustl.edu>`_.
