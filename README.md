Coupled Cluster for Solids
==========================

Quick start for icc and Intel MPI library
-----------------------------------------

To build the required libraries and cc4s issue the following command
```
make CONFIG=icc_impi extern -j 4
make CONFIG=icc_impi cc4s -j 4
```

These commands use the configuration specified in the `./etc/config/icc_impi.mk` file.

Library dependencies
--------------------

- The dependencies that are built together with `cc4s` for a given configuration
  are built by issuing the command
  ```
  make CONFIG=<config> extern
  ```
  where `<config>` is the name of the configuration used.
- by default the configuration `gcc` is used.

Features
--------

Different features of the code can be controlled by the `config.mk`
file in the root directory.

For instance, if you want the library `atrip` to be built,
using the `icc_impi` configuration,
you can write in your `config.mk` file the following

```make
CONFIG = icc_impi
```

and you can simply build by doing

```sh
make extern -j4
make cc4s -j4
```

Note that you can use the silent mode of the Makefile by issuing
the make commands with the silent flag `-s`, i.e.,
```sh
make -s cc4s -j 4
```
or if you want to activate always the silent mode you can write in your
`config.mk` file

```make
.SILENT: cc4s
```

Building
--------

-   check `gcc` version. Need at least `gcc 4.7.2` (As tested, gcc 4.6 does
    not work).
-   write or edit a `config.<config>` file for your build environment.
-   the configuration files given in etc/config/ contain
    predifined retail environments for gnu- and intel-compiler, respectively,
    using full optimization and without debugging info.
-   for the gcc configuration the following additional libraries are
    required
    - OpenBLAS
    - scalapack
- make sure, the above library dependencies are built for your configuration
- run `make -j 8 [CONFIG=<config]` to build for the desired environment, by
  default for `gcc`. The `-j` option issues a parallel make on 8 processes.
- run `make install [CONFIG=<config>]` to copy the executable to the specified
  target directory. The default is `~/bin/cc4s/<config>`.
- the intermediate build files for each build environment can be found in the
  directory `build/<config>/`.

Running
-------

-   a `cc4s` operation file, e.g. `mp2.yaml`, can be run with
    ```
    ~/bin/cc4s/gcc/Cc4s -i mp2.yaml
    ```

Testing
-------

### Running tests

- Run the main Makefile with
  ```
  make test CONFIG=gcc
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
  by editing the `Cc4s.mk` file and replacing the commit hash
  by the new commit hash in the dependency's repository.
- build the dependency for all configurations supported as described above.
- build `cc4s` for all configurations supported
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
