import sys
import os
sys.path.append("d:/measuring")
sys.path.append("h:/My Documents/measuring")#only for local SvD
import numpy as np
import h5py

from matplotlib import pyplot as plt

from analysis.lib.tools import toolbox as tb
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro, mbi, sequence, pqsequence
from analysis.lib.nv import nvlevels
from analysis.lib.lde import tail_cts_per_shot_v4 as tail
from analysis.lib.pq import pq_tools, pq_plots
from analysis.lib.math import readout_correction as roc
from analysis.lib.math import error
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
reload(m2)
reload(tb)
reload(ssro)
reload(mbi)
reload(sequence)
reload(pqsequence)
reload(tail)
reload(pq_tools)
reload(pq_plots)