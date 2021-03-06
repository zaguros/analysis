# This file is part of QuTiP.
#
#    QuTiP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    QuTiP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QuTiP.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2013, Paul D. Nation & Robert J. Johansson
#
###########################################################################

import os
from numpy import amax
from numpy.testing import assert_, run_module_suite
import scipy

from qutip import *


class TestFileIO:
    """
    A test class for the QuTiP functions for writing and reading data to files.
    """

    def testRWRealDefault(self):
        "Read and write real valued default formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data)
        data2 = file_data_read("test.dat")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWRealDecimal(self):
        "Read and write real valued decimal formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data, "real", "decimal")
        data2 = file_data_read("test.dat", ",")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWRealExp(self):
        "Read and write real valued exp formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data, "real", "exp")
        data2 = file_data_read("test.dat", ",")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWComplexDefault(self):
        "Read and write complex valued default formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N)) + 1j * (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data)
        data2 = file_data_read("test.dat")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWComplexDecimal(self):
        "Read and write complex valued decimal formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N)) + 1j * (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data, "complex", "decimal")
        data2 = file_data_read("test.dat", ",")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWComplexExp(self):
        "Read and write complex valued exp formatted data"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N)) + 1j * (1 - 2 * scipy.rand(N, N))

        file_data_store("test.dat", data, "complex", "exp")
        data2 = file_data_read("test.dat", ",")
        # make sure the deviation is small:
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")

    def testRWSeparatorDetection(self):
        "Read and write with automatic separator detection"

        # create some random data
        N = 10
        data = (1 - 2 * scipy.rand(N, N)) + 1j * (1 - 2 * scipy.rand(N, N))

        # comma separated values
        file_data_store("test.dat", data, "complex", "exp", ",")
        data2 = file_data_read("test.dat")
        assert_(amax(abs((data - data2))) < 1e-8)

        # semicolon separated values
        file_data_store("test.dat", data, "complex", "exp", ";")
        data2 = file_data_read("test.dat")
        assert_(amax(abs((data - data2))) < 1e-8)

        # tab separated values
        file_data_store("test.dat", data, "complex", "exp", "\t")
        data2 = file_data_read("test.dat")
        assert_(amax(abs((data - data2))) < 1e-8)

        # space separated values
        file_data_store("test.dat", data, "complex", "exp", " ")
        data2 = file_data_read("test.dat")
        assert_(amax(abs((data - data2))) < 1e-8)

        # mixed-whitespace separated values
        file_data_store("test.dat", data, "complex", "exp", " \t ")
        data2 = file_data_read("test.dat")
        assert_(amax(abs((data - data2))) < 1e-8)
        os.remove("test.dat")


if __name__ == "__main__":
    run_module_suite()
