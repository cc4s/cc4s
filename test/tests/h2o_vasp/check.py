#!/usr/bin/env python3
from testis import read_yaml, compare_energies

compare_energies("correct.out.yaml", "cc4s.out.yaml", accuracy=1e-7)
