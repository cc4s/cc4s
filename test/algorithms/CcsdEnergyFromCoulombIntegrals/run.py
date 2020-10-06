#!/usr/bin/env python3

from testis import call, read_yaml


call("{CC4S_PATH} -i ccsdReal.yaml")

out = read_yaml("cc4s.out")
assert out["dryRun"] == 0, "We should not be doing dryRuns now"
assert (out["steps"][0]["out"]["tensor"]["fermiEnergy"] == -4.4484944443468191,
       "Fermi energy is not correct")

call("{CC4S_PATH} -i ccsdComplex.yaml")
out = read_yaml("cc4s.out")