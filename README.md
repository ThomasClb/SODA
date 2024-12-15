# SODA
A Stochastic Optimisation solver with Differential Algebra.

## Requirements
To use the SODA code, you need a C++11 compiler.

The code was mainly tested on Unix distributions but could be adapted with very little effort.

The only dependency on this code is with the DACE, which can be found at: https://github.com/dacelib/dace.
At installation on the DACE do not forget to use the `AlgebraicMatrix<T>` class in the `CMakeLists.txt`.
Moreover, we recommend commenting the 189th line of the file `dace/interfaces/cxx/DACEException.cpp` to avoid useless and frequent warnings.

## Setting up SODA
We recommend using CMake to use SODA and build the test cases.
To use the SODA library, simply clone this repository:
```
git clone https://github.com/ThomasClb/SODA.git
```
You can verify that SODA is correctly cloned by running the unit tests:
```
cd SODA
mkdir build_test
cd build_test
cmake ../test/
make
./unit_test
```
Then create a build directory run cmake, and then compile the test_cases:
```
cd SODA
mkdir build
cd build
cmake ../script/test_cases/
make
```
You might need to locate the DACE for cmake if it was not installed using:
```
sudo make install
```

## Running the test cases
Six test cases are given.

In the SODA folder run:
```
./build/test_case <test_case ID> <parameters> 
```
Leaving the parameters field empty will return an error message that will inform you on what are the available options.

These parameters consist of 9 arguments:
- Test case ID.
- The address of the SpacecraftParameter. (They consist of the dynamical system, followed by the thrust-to-mass ratio of the spacecraft)
- The DDP method to use. (recommended 1)
- The number of trajectory segments N?
- The time of flight in days.
- Perform robust optimisation. (1=True, 0=False)
- Target failure risk. (between 0.0 and 1.0)
- Perform fuel optimal optimisation. (1=True, 0=False)
- Perform solution polishing with Newton method. (1=True, 0=False)
- Save the data. (1=True, 0=False)
- Verbosity. (0=Full, 1=Medium, 2=None, 3=Benchmark)

For instance 
- The double integrator problem from [Lantoine and Russell 2012]:
	```
	./build/test_case 0 ./data/spacecraft_parameters/double_integrator_1e3.dat 1 11 0 1 0.05 0 0 1 0
	```
- The energy optimal Earth-Mars transfer from [Lantoine and Russell 2012]:
	```
	./build/test_case 1 ./data/spacecraft_parameters/tbp_SUN_5e-4.dat 1 40 348.79 1 0.05 0 1 1 0
	```
- The fuel optimal Earth-Mars transfer from [Lantoine and Russell 2012]:
	```
	./build/test_case 1 ./data/spacecraft_parameters/tbp_SUN_5e-4.dat 1 40 348.79 1 0.05 1 1 1 0
	```
- The fuel optimal Halo L2 to Halo L1 transfer from [Aziz et al. 2019]:
	```
	./build/test_case 2 ./data/spacecraft_parameters/cr3bp_EARTH_MOON_5e-4.dat 1 110 20 1 0.05 1 1 1 0
	```
The robust trajectory is saved in `./data/robust_trajectory`, and the test sample is (partly) saved in `./data/sample_trajectory`.

Reading the source code of the test cases is highly recommended to understand how they work and how to use SODA.


## Visualisation
The visualisation of the results is done using Python 3. Just run:
```
./source/visualization/main.py <test_case ID> <parameters>
```
Where the parameters are:
- Test case ID. (integer)
- Target failure risk. (float)
- The thrust-to-mass ratio in N/kg. (string)
- The time of flight in days. (integer)
- The DDP method used. (integer)
- Show the associated sample. (boolean)

	
For instance 
- The double integrator problem from [Lantoine and Russell 2012] with DDP method 1:
	```
	./source/visualisation/main.py 0 0.05 1e3 0 1 1
	```
- The fuel optimal Earth-Mars transfer from [Lantoine and Russell 2012] with DDP method 1:
	```
	./source/visualisation/main.py 1 0.05 5e-4 348 1 1
	```
- The fuel optimal Halo L2 to Halo L1 transfer from [Aziz et al. 2019] with DDP method 1:
	```
	./source/visualisation/main.py 2 0.05 5e-4 20 1 1
	```
The results are saved in `./data/plots`.
