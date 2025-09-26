Coupled Cluster for Solids
==========================

`CC4S` provides a cmake build system.

`CC4S` depends on two libraries:
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [Cyclops Tensor Framework](https://github.com/cyclops-community/ctf)

These dependencies are built automatically by the cmake program.

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

- Run the main Makefile with
  ```
  make test CONFIG=gcc CC4S=build/Cc4s
  ```
- this issues all tests of the given type for local build binary of the given
  build environment.

### Writing tests

TODO

Update dependencies
-------------------

If any of the dependencies are updated by their respective maintainer, you
can incorporate the latest version into cc4s. Note that this may lead to
incompatabilities and it must therefore be done with good care.

- if you intent to update the master, create a branch from the master.
  In case anything goes wrong the damage is controlled
- update the dependency at the top level of the cc4s directory structure:
  by editing the `CMakeLists.txt` file and replacing the commit hash
  by the new commit hash in the dependency's repository.
- run the `cc4s` `essental` tests for all configurations supported
- fix all bugs that emerged from advancing to the new version in `cc4s` or
  let them be fixed in the modules.

If all tests pass, `cc4s` may be advanced to the new dependency's  versions
by
```
git commit -m "Update dependency to version {...}"
```

- if you want to advance the master branch, merge your branch into it.

Note that you may commit changes to the your branch even if things do
not work. However, each commit will be visible in the history.
