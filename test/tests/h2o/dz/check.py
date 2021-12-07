#!/usr/bin/env python3

from testis import read_yaml

out = read_yaml("cc4s.out")

assert out["dryRun"] == 0, "We should not be doing dryRuns now"
