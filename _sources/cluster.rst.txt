Computation Guide
=================

Computation Time Scaling
------------------------

Fully converged Bayesian inference requires numerous iterations to explore the entire prior parameter space. This process demands substantial computation time, which can be very expensive depending on the Equation of State (EOS) you are using.

 **1. Using Integrated EOS Models**

If you are utilizing the EOS models integrated into this package, here is an overview of the computational resources required:

a. **Relativistic Mean Field (RMF) Model**

   - **Optimization:** Please use the Numba-accelerated version of EOS computation.
   - **Resources:** With one NICER observation constraint and utilizing **256 CPU cores**, the computation can be completed within approximately **1 hour**.
   - **Scalability:** The computation time will increase as you include more observations.

b. **Speed of Sound Model and Polytropic Model**

   - **Nature:** These are meta-models with explicit mathematical formulas, allowing for analytical computations.
   - **Resources:** Using **24 CPU cores**, the inference can be completed within **6 hours**.
   - **Scalability:** The computation time scales with the number of observations included.

c. **Strangeon Model and MIT Bag Model**

   - **Nature:** Although these are physics-based models, the computations are relatively straightforward and fast.
   - **Resources:** On a laptop with **24 cores**, the computation can be completed within **5-6 hours**.
   - **Scalability:** The computation time increases with the addition of more observations.

**2. Using Custom EOS Models**

If you are defining your own EOS and wish to perform inference using this package:

- **Performance Requirement:** Ensure that the computation time for your EOS is **less than 1 second** per evaluation. Otherwise, the overall computation will become prohibitively expensive due to the large number of samples required for Bayesian inference.

Parallelization
---------------

Given the computational intensity of Bayesian inference, parallelization strategies are essential to optimize performance. Below are the steps to parallelize computations using **UltraNest** without modifying the package's code.

**Requirements**

- **OpenMPI:** Ensure that OpenMPI is installed and available on your system.
- **Python Packages:** Install `h5py` and `mpi4py`.

**Steps to Parallelize**

1. **Set the Number of Threads**

   Control the number of threads per process by setting the `OMP_NUM_THREADS` environment variable. For example, to set it to 1:

   .. code-block:: bash

       export OMP_NUM_THREADS=1

2. **Run the Inference Script with MPI**

   Use `mpiexec` to execute your inference script across multiple cores. Replace `inference.py` with the name of your inference script.

   .. code-block:: bash

       mpiexec -np 4 python inference.py

   - **Explanation:** This command runs the inference using **4 CPU cores**. Adjust the `-np` value according to the number of available cores you wish to utilize.

**Recommendations**

- **Virtual Environment:** Regardless of where you are running the code, it is highly recommended to use a virtual environment to manage your installations. This practice helps maintain dependencies and avoids conflicts with other projects.

**Additional Resources**

For more detailed information on performance and parallelization strategies with **UltraNest**, refer to the `UltraNest Performance Page <https://johannesbuchner.github.io/UltraNest/performance.html>`_.

