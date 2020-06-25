LARF
======================================
LARF is a tool for time-division-multiplexing (TDM) optimization in multi-FPGA systems. 
Given the timing graph and the routing resources between FPGAs, it will generate the TDM ratio for each inter-FPGA signal 
such that the system clock period (the maximum arrival time of the sinks of the timing graph).
More details may refer to the following papers,
* Chak-Wa Pui and Evangeline F. Y. Young, 
['Lagrangian Relaxation-Based Time-Division Multiplexing Optimization for Multi-FPGA Systems'](https://dl.acm.org/doi/10.1145/3377551), 
accepted by ACM Transactions on Design Automation of Electronic Systems (TODAES)
* Chak-Wa Pui and Evangeline F. Y. Young, 
['Lagrangian Relaxation-Based Time-Division Multiplexing Optimization for Multi-FPGA Systems'](https://ieeexplore.ieee.org/document/8942125), 
IEEE/ACM International Conference on Computer-Aided Design (ICCAD), Westminster, CO, USA, Nov. 4-7, 2019.

(LARF supports [ISPD'16](https://www.dropbox.com/sh/9c74a6f4o0rrd2t/AAA3V_fiP15pV20fV62apLoqa) benchmarks.
This version of code is consistent with the one accpected by TODAES.)

## 1. How to Build

**Step 1:** Download the source code. For example,
```bash
$ git clone https://github.com/cuhk-eda/larf
```

**Step 2:** Go to the project root and build by
```bash
$ cd larf/src
$ make
```

Note that this will generate a folder `bin` under the root, which contains binaries and auxiliary files.
More details are in [`Makefile`](src/Makefile).

### 1.1. Dependencies

* [GCC](https://gcc.gnu.org/) (version >= 4.8.0) or other working c++ compliers
* [Boost](https://www.boost.org/) (version >= 1.58)
* [Python](https://www.python.org/) (version 3, optional, for utility scripts)

## 2. How to Run

### 2.1. Toy Test

#### Run Binary Directly

Go to the `bin` directory and run the binary `larf` with a toy case `ispd16_f01`:
```bash
$ cd bin
$ mkdir f01
$ ../larf -aux ../../toys/ispd2016/FPGA01/design.aux -flow tdm_part -partition 5
$ ../larf -aux 0/design.aux -out FPGA01_0.pl -flow tdm_place
$ ../larf -aux 1/design.aux -out FPGA01_1.pl -flow tdm_place
$ ../larf -aux 2/design.aux -out FPGA01_2.pl -flow tdm_place
$ ../larf -aux 3/design.aux -out FPGA01_3.pl -flow tdm_place
$ ../larf -aux 4/design.aux -out FPGA01_4.pl -flow tdm_place
$ ../larf -aux ../../toys/ispd2016/FPGA01/design.aux -flow tdm_time -out f01.tdm -partition 5
```

#### Run with a Wrapping Script

Instead of running the binary directly, you may also use a wrapping script `run.sh` to save typing and do more:
```bash
$ export BENCHMARK_PATH=$(pwd)/toys
$ cd bin
$ ./run.sh tdm_part f01 5
$ ./run.sh tdm_place f01 5
$ ./run.sh tdm_time f01 5
```
More usage may refer to `scripts/run.sh` and `src/main.cpp`.

### 2.2. Batch Test

The benchmarks can be downloaded from [ISPD'16](https://www.dropbox.com/sh/9c74a6f4o0rrd2t/AAA3V_fiP15pV20fV62apLoqa).
You may let `run.sh` know the benchmark path by setting OS environmental variable `BENCHMARK_PATH`.
Then,
```bash
$ cd bin
$ ./run.sh tdm_part <benchmark_name...|all> <#partitions> [option...]
$ ./run.sh tdm_place <benchmark_name...|all> <#partitions> [option...]
$ ./run.sh tdm_time <benchmark_name...|all> <#partitions> [option...]
```

## 3. Modules

* `scripts`: utility python/bash scripts
* `src`: C++ source code
    * `alg`: external algorithm packages
    * `db`: database
    * `gp`: global placement
    * `tdm`: time-division-multiplexing optimization
    * `utils`: utilities
* `toys`: toy test cases


## 4. Results

Experiments are performed on a 64-bit Linux workstation with Intel Xeon Silver 4114 CPU (2.20GHz, 40 cores) and 256GB memory, where eight threads are used.

| Design | System Clock Period (ps) | Runtime (sec) |
|:------:|:------------------------:|:-------------:|
| FPGA01 |            174           |       31      |
| FPGA02 |            415           |       38      |
| FPGA03 |            863           |      151      |
| FPGA04 |           3153           |      190      |
| FPGA05 |           8325           |      192      |
| FPGA06 |           1151           |      324      |
| FPGA07 |           2960           |      358      |
| FPGA08 |           4214           |      162      |
| FPGA09 |           4695           |      358      |
| FPGA10 |           1066           |      257      |
| FPGA11 |           7268           |      262      |
| FPGA12 |           1288           |      520      |

## 5. License

READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY USING THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THE FOLLOWING AGREEMENT. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT LICENSEE OF THIS PRODUCT.



License Agreement for LARF



Copyright (c) 2020 by The Chinese University of Hong Kong



All rights reserved



CU-SD LICENSE (adapted from the original BSD license) Redistribution of the any code, with or without modification, are permitted provided that the conditions below are met.



Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.



Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.



Neither the name nor trademark of the copyright holder or the author may be used to endorse or promote products derived from this software without specific prior written permission.



Users are entirely responsible, to the exclusion of the author, for compliance with (a) regulations set by owners or administrators of employed equipment, (b) licensing terms of any other software, and (c) local, national, and international regulations regarding use, including those regarding import, export, and use of encryption software.



THIS FREE SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR ANY CONTRIBUTOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, EFFECTS OF UNAUTHORIZED OR MALICIOUS NETWORK ACCESS; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
