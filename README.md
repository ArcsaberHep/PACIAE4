<!--
// README.md is a part of the PACIAE event generator.
// Copyright (C) 2024 PACIAE Group.
// PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
// Open source: https://github.com/ArcsaberHep/PACIAE4
// Author: An-Ke Lei, January 2024 - January 2025.

// This is a README file for PACIAE.

//                                               By An-Ke at CCNU on 16/01/2024
//                                          Last updated by An-Ke on 17/01/2025
 -->

# The parton and hadron cascade model PACIAE 4

***PACIAE 4*** (***Parton And-hadron China Institute of Atomic Energy***) is a multipurpose Monte Carlo event generator developed to describe a wide range of collisions, including lepton-lepton, lepton-hadron, lepton-nucleus, hadron-hadron, hadron-nucleus, and nucleus-nucleus collisions. It is built based on **PYTHIA 6.4** & **PYTHIA 8.3**, and incorporates parton and hadron rescattering stages to address the nuclear medium effects. PACIAE 4 is the new generation of ***PACIAE*** model surpassing the version ***PACIAE 3***. In PACIAE 4, the old fixed-format FORTRAN 77 code has been refactored and rewritten by the free-format modern Fortran mixed with C++ languages. The C++-based **PYTHIA 8** is interfaced in. More physics and features are available.

## Installation

Just download the file and decompress it. Then you can get the PACIAE 4 source code directly.

## Prerequisites

 1. **`make`** tool, **`gfortran`** and **`g++`** compilers are required.

 2. **`PYTHIA 8`** library is indispensable. PYTHIA official website: [https://pythia.org/](https://pythia.org/) . Download then build PYTHIA 8 using the following command:
     ```
      ./configure
      make -j
     ```
    **`-j`** flag will compile PYTHIA 8 more quickly using multiple cores. 
     After the compilation and build completed, copy the **<font color=red> Makefile.inc </font>** file the PYTHIA 8 generated and paste it into the PACIAE 4 directory. If you have the installed PYTHIA 8 already, just copy the **<font color=red> Makefile.inc </font>** from your own installed PYTHIA 8 directory and paste it into the PACIAE 4 directory. **<font color=red> Makefile.inc </font>**  file is very important for the running of PACIAE 4. It specifies the PYTHIA 8 distribution path.

**NB**:
 1. The current Makefile is designed for linux system. The support for other OS will be updated in the next version.
 2. Released PACIAE 4.002 requires PYTHIA 8.312 and below because PYTHIA 8.313 modifies several code, especially the parts of Angantyr. From this version, PACIAE requires PYTHIA 8.313.

## Usage

We encourage users to run the program on LINUX.

Two ways to run the program are provided.
 1. Manual running. Use the **`make`** tool to compile the source code, link and build the PACIAE program. Then execute it by hand. **`GFortran`** and **`g++`** compilers have been pre-specified inside **Makefile**:
    - Compile code, link and build the programs by the command:
        ```
            make
        ```
        One can choose **`make -j`** to compile in parallel. The executable program **xPaciae.x** and the input files **usu.dat**, **pythia6_extra.cfg** and **pythia8_extra.cfg** will be generated in the folder **`sim`** (simulation). Enter **`sim`**.
    - Modify the input files according to your wish. The basic PACIAE-related and few PYTHIA-related parameters and switches are included in the **usu.dat** file. Moreover, user can add extra settings for PYTHIA 6/8 in **pythia6_extra.cfg** and **pythia6_extra.cfg** files, which requires the knowledge of PYTHIA parameters. See [PYTHIA 6.4 Physics and Manual, JHEP 05 (2006) 026](https://arxiv.org/abs/hep-ph/0603175) and [PYTHIA 8 online manual](https://pythia.org//latest-manual/Welcome.html) for more details. Thanks to these two files, PACIAE 4 can use almost all PYTHIA 6/8 features. 
    - Run the program by the command:
        ```
            ./xPaciae.x
        ```
      One can use the following command to run the program in the background and record the time and log information.
        ```
            nohup time ./xPaciae.x > paciae.log &
        ```
 2. Use the [PACIAE.sh](PACIAE.sh) shell-script to compile, link and build the program, generate input files and run the program automatically.
    - Modify the [PACIAE.sh](PACIAE.sh) file as needed. usu.dat, pythia6_extra.cfg and pythia6_extra.cfg have been integrated within [PACIAE.sh](PACIAE.sh) file. User can modify them in [PACIAE.sh](PACIAE.sh) directly. See **`USU_DAT_BLOCKTEXT`**, **`PY6_CFG_BLOCKTEXT`** and **`PY8_CFG_BLOCKTEXT`** parts inside it.
    - If executable permissions are missing, grant executable permissions to [PACIAE.sh](PACIAE.sh) by command (only once):
        ```
            chmod +x PACIAE.sh
        ```
    - Run the [PACIAE.sh](PACIAE.sh) script by the command:
        ```
            ./PACIAE.sh
        ```
    - **rms_analysis.f90** is a stand-alone program used in conjunction with the [PACIAE.sh](PACIAE.sh) script to average the internal analysis files **rms.out** generated by each PACIAE simulation.

    It is also worth mentioning that tasks can be submitted to computer clusters and supercomputers using the [PACIAE.sh](PACIAE.sh) script (It is currently available for SLURM, LSF and PBS scheduling systems only.). More detailed information and usage can be found in the **`SCHEDULE_SYSTEM_BLOCKTEXT`** part of the [PACIAE.sh](PACIAE.sh).file.

## Oscar format output

- The ASCII (text) human-readable output file **oscar.out** of the final particles following the OSCAR1997A standard can be generated via the switch **`nosc=1`**. A file snippet is as follows:
```
1  | OSC1997A
2  | final_id_p_x
3  | PACIAE   program_version   projectile+target   frame   energy   1
4  | Event_1 number_of_events number_of_particles b_parameter b_angle event_weight
5  |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1  (GeV & fm)
6  |   Particle_No_2  id2  px2  py2  pz2  E2  m2  x2  y2  z2  t2
7  |   ...
8  | Event_2 number_of_events number_of_particles b_parameter b_angle event_weight
9  |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
10 |   ...
```

- The full event history following the OSCAR199A standard including **(`STAGE_0`)** the initial state,**(`STAGE_1`)** the initial partonic state,**(`STAGE_2`)** the final partonic state,**(`STAGE_3`)** the initial hadronic state from the hadronization, and **(`STAGE_4`)** the final hadronic state, can be generated via **`nosc=2`**. A file snippet is as follows:
```
1  | # OSC1999A
2  | # full_event_history
3  | # PACIAE  program_version
4  | # projectile+target   frame   energy   1
5  | Event_1  N_events  N_particles  b_parameter  b_angle  STAGE_0  event_weight
6  |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1  (GeV & fm)
7  |   Particle_No_2  id2  px2  py2  pz2  E2  m2  x2  y2  z2  t2
8  |   ...
9  | Event_1  N_events  N_particles  b_parameter  b_angle  STAGE_1  event_weight
10 |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
11 |   ...
12 | Event_1  N_events  N_particles  b_parameter  b_angle  STAGE_2  event_weight
13 |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
14 |   ...
15 | Event_1  N_events  N_particles  b_parameter  b_angle  STAGE_3  event_weight
16 |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
17 |   ...
18 | Event_1  N_events  N_particles  b_parameter  b_angle  STAGE_4  event_weight
19 |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
20 |   ...
21 | Event_2  N_events  N_particles  b_parameter  b_angle  STAGE_9  event_weight
22 |   Particle_No_1  id1  px1  py1  pz1  E1  m1  x1  y1  z1  t1
23 |   ...
```

## Maintainers

[@ArcsaberHep](https://github.com/ArcsaberHep) (An-Ke Lei)
[@ArcsaberkxL](https://github.com/ArcsaberkxL)

## Contributing

Feel free to dive in! Any bug reports, comments and suggestions are welcome. Please do not hesitate to contact us.

## Contributors

 - [An-Ke Lei](https://inspirehep.net/authors/1965068), ankeleihep@gmail.com or ankelei@mails.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. -->
 - [Zhi-Lei She](https://inspirehep.net/authors/1903611), shezhilei@cug.edu.cn <!-- School of Mathematical and Physical Sciences, Wuhan Textile University, Wuhan 430200, China --> 
 - [Dai-Mei Zhou](https://inspirehep.net/authors/1030208), zhoudm@mail.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. -->
 - [Yu-Liang Yan](https://inspirehep.net/authors/1051028), yuliang86@yeah.net <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->
 - [Ben-Hao Sa](https://inspirehep.net/authors/990834), sabhliuym35@qq.com <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->

## License

[GPL v2.0](LICENSE) and any later version.

## Released programs

PACIAE 4.0 code are hosted on [https://github.com/ArcsaberHep/PACIAE4](https://github.com/ArcsaberHep/PACIAE4) and [https://gitee.com/arcsaberhep/PACIAE4](https://gitee.com/arcsaberhep/PACIAE4).

The released code are available on [https://github.com/ArcsaberHep/PACIAE4/releases](https://github.com/ArcsaberHep/PACIAE4/releases) and [https://gitee.com/arcsaberhep/PACIAE4/releases](https://gitee.com/arcsaberhep/PACIAE4/releases).

## Released papers

 - **Recent version:**
   - PACIAE 4.0: A brief introduction to PACIAE 4.0. (waiting...)
   - PACIAE 3.0: An introduction to the parton and hadron cascade model PACIAE 3.0, [Phys. Rev. C 108 (2023) 6, 064909](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.108.064909) or [2309.05110 [hep-ph]](https://arxiv.org/abs/2309.05110).
<br/>

 - **Past version:**
   - PACIAE 2.2.2: Revisiting the centrality definition and observable centrality dependence of relativistic heavy-ion collisions in PACIAE model, [Comput. Phys. Commun. 284 (2023) 108615](https://doi.org/10.1016/j.cpc.2022.108615) or [arXiv:2212.04087 [nucl-th]](https://doi.org/10.48550/arXiv.2212.04087).

   - PACIAE 2.2.1: An updated issue of the parton and hadron cascade model PACIAE 2.2, [Comput. Phys. Commun. 274 (2022) 108289](https://doi.org/10.1016/j.cpc.2022.108289). <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

   - PACIAE 2.2.0: Announcement for the replacement of the PACIAE 2.1 and PACIAE 2.2 series, [Comput. Phys. Commun. 224 (2018) 417-418](https://doi.org/10.1016/j.cpc.2017.10.006). <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

   - PACIAE 2.2.0: An upgraded issue of the parton and hadron cascade model, PACIAE 2.2, [Comput. Phys. Commun. 193 (2015) 89-94](https://doi.org/10.1016/j.cpc.2015.01.022) or [arXiv:1412.7579 [nucl-th]](https://doi.org/10.48550/arXiv.1412.7579).

   - PACIAE 2.1: An updated issue of the parton and hadron cascade model PACIAE 2.0, [Comput. Phys. Commun. 184 (2013) 1476-1479](https://doi.org/10.1016/j.cpc.2012.12.026) or [arXiv:1206.4795 [nucl-th]](https://doi.org/10.48550/arXiv.1206.4795).

   - PACIAE 2.0: An updated parton and hadron cascade model (program) for the relativistic nuclear collisions, [Comput. Phys. Commun. 183 (2012) 333-346](https://doi.org/10.1016/j.cpc.2011.08.021) or [arXiv:1104.1238 [nucl-th]](https://doi.org/10.48550/arXiv.1104.1238).

   - PACIAE 1.0: Charge particle universal rapidity scaling in e+e-, pbar + p and Au + Au collisions at relativistic energies and its partonic origin, [J. Phys. G 32 (2006) 243-250](https://doi.org/10.1088/0954-3899/32/3/001); Influence of the partonic Pauli blocking on the hadronic final state in relativistic nucleus-nucleus collisions, [Phys. Rev. C 70 (2004) 034904](https://doi.org/10.1103/PhysRevC.70.034904).
<br/>

 - **Ancient version:**
   - JPCIAE: J/psi dynamical suppression in a hadron and string cascade model (or Formation time effect on J/psi dynamical nuclear suppression), [Phys. Rev. C 59 (1999) 2728-2733](https://doi.org/10.1103/PhysRevC.59.2728) or [arXiv:nucl-th/9803033](https://arxiv.org/abs/nucl-th/9803033); J/psi normal and anomalous suppressions in a hadron and string cascade model, [J.Phys.G 25 (1999) 1123-1133](https://doi.org/10.1088/0954-3899/25/6/302) or [arXiv:nucl-th/9809020](https://arxiv.org/abs/nucl-th/9809020); Inclusive and direct photons in S + Au collisions at 200A GeV/c, [Phys. Rev. C 61 (2000) 064905](https://doi.org/10.1103/PhysRevC.61.064905) or [arXiv:nucl-th/9904035](https://arxiv.org/abs/nucl-th/9904035).
<br/>

 - **Origin version:**
   - LUCIAE 3.0: A New version of a computer program for firecracker model and rescattering in relativistic heavy ion collisions, [Comput. Phys. Commun. 116 (1999) 353](https://doi.org/10.1016/S0010-4655(98)00138-6) or [arXiv:nucl-th/9804001](https://doi.org/10.48550/arXiv.nucl-th/9804001)

   - LUCIAE 2.0: An Event generator for the firecracker model and the rescattering in high-energy pA and AA collisions: LUCIAE version 2.0, [Comput. Phys. Commun. 90 (1995) 121-140](https://doi.org/10.1016/0010-4655(95)00066-O); Final state interactions in the (nuclear) FRITIOF string interaction scenario, [Z. Phys. C 70 (1996) 499](https://doi.org/10.1007/s002880050127).
<br/>

 - **Special version:** 
   - HYDRO-PACIAE, a hydrodynamic and transport hybrid model for ultra-relativistic heavy ion collisions, [J.Phys.G 40 (2013) 025102](https://doi.org/10.1088/0954-3899/40/2/025102) or [arXiv:1110.6704 [nucl-th]](https://doi.org/10.48550/arXiv.1110.6704).

## Relevant papers

 - A comprehensive guide to the physics and usage of PYTHIA 8.3, [SciPost Phys.Codeb. 8 (2022)](https://scipost.org/10.21468/SciPostPhysCodeb.8) or [arXiv:2203.11601 [hep-ph]](https://arxiv.org/abs/2203.11601).

 - PYTHIA 6.4 Physics and Manual, [JHEP 05 (2006) 026](https://doi.org/10.1088/1126-6708/2006/05/026) or [arXiv:hep-ph/0603175 [hep-ph]](https://doi.org/10.48550/arXiv.hep-ph/0603175) (and its update notes: [https://pythia.org/download/pythia6/pythia6428.update](https://pythia.org/download/pythia6/pythia6428.update)).

## External links

- **[PYTHIA 8 (& 6)](https://pythia.org/)**

- **[EPS & EPPS PDFs](https://research.hip.fi/qcdtheory/nuclear-pdfs/)**, **[LHAPDF](https://lhapdf.hepforge.org/)**, **[CTEQ PDFs](https://cteq.gitlab.io/project/pdfs/)**

- **[ROOT](https://root.cern.ch/)**, **[CERNLIB](https://cernlib.web.cern.ch/cernlib/)**

- **[HEPForge](https://www.hepforge.org//)**

- **[MCnet](https://montecarlonet.org/)**

## Update notes:

<!----------------------------------------------------------------------------->
### 11/2024: In version PACIAE 4.0.02
- In "analy_40.f90" and "Rms_analysis.f90", fixed analysis code of mean pT \< pT \>. Analysis code of anisotropic flows were added.

<!----------------------------------------------------------------------------->
### <font color=red> 11/2024 PACIAE 4.0 is now released! </font>
- The released version if 4.0.02.
- The code have been refactored and rewrite by the free-format modern Fortran + C++ from the fixed-format FORTRAN77.
- PYTHIA 8 can be interfaced in now.
- The gluon splitting and quark deexcitation mechanisms, along with the coalescence hadronization model have been improved.
- The partonic and hadronic cascades have been improved.
- Some other aspects.

<!----------------------------------------------------------------------------->
### <font color=red> 01/2024  </font> PACIAE 4 project launched!

<!----------------------------------------------------------------------------->
### <font color=red> 12/2023 PACIAE 3.0 is now released! </font>