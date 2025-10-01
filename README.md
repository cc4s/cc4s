Coupled Cluster for Solids
==========================

`CC4S` provides a cmake build system.

`CC4S` depends on two libraries:
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [Cyclops Tensor Framework](https://github.com/cyclops-community/ctf)

These dependencies are built automatically by the cmake program.

Super Quick start
-----------------

Load the relevant modules in your computer, you'll need
- BLAS
- lapack
- Scalapack

run the build script,

```sh
./build.sh
```

This will build by default the project in the directory `build`.

Building the project
--------------------

To make the project, create a folder and build it there as follows

```sh
mkdir build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=mpicxx
```

All the project will be built in the `build` directory, including its dependencies.
You may have several build directories if you intend to compile with different flags,
compilers or settings.

If you want to make the debug version of `CC4S`, you should pass to cmake the
`-DCMAKE_BUILD_TYPE=Debug` flag, for example

```sh
mkdir build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Debug
```

This will compile with the `-O0` flag all the source files.

Running
-------

-   a `cc4s` operation file, e.g. `mp2.yaml`, can be run with
    ```
    mpirun ./build/Cc4s -i mp2.yaml
    ```

Testing
-------

### Running tests

- Once you have built the project, you can run the tests by issuing the cmake test command
  ```
  ctest -V
  ```
- The test scripts are found in the test directory inside the build directory,
  you can edit them to have more fine control on which tests to run:

  - test/run -> A general script to run tests
  - test/run-integration-tests -> Run all the essential tests
  - test/run-integration-tests-run -> only run the tests, do not check the results
  - test/run-integration-tests-check -> only check the results of the test runs.
