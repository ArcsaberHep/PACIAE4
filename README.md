<!-- This is a README file for usage of PACIAE.
     Written by Markdown language.
                 By Anke at UiO on 31/07/2024 
                                               
                    Last updated on 31/07/2024 
 -->

# The parton and hadron cascade model PACIAE 4 (waiting...)

 ***PACIAE 4*** (***Parton And-hadron China Institute of Atomic Energy 4***) is a new generation of ***PACIAE 3***. In PACIAE 4, the modern C++-based PYTHIA 8 is interfaced in. More physics and features are availble. Now PACIAE 4 is written in a mix of Fortran and C++ languages. PACIAE 4 is a multipurpose Monte Carlo event generator developed to describe a wide range of collisions, including lepton-lepton, lepton-hadron, lepton-nucleus, hadron-hadron, hadron-nucleus, and nucleus-nucleus collisions. It is built based on PYTHIA 6.428 & PYTHIA 8, and incorporates parton and hadron rescattering stages to take care of the nuclear medium effects. 

## Installation

Just download the file and decompress it. Then you can get the PACIAE 4 source code directly.

## Usage

We encourage users to run the program on LINUX.

First and formost, the PYTHIA 8 library is indispensable

...(waiting to update)

There are two ways to run the program.
 1. Manul running. Use the "make" tool to complie the soruce code, link and generate the PACIAE program. Then execute it. GFortran and g++ compilers have been pre-specified inside "Makefile"on LINUX as an example:
    - Compile code, link and generate the programs by the command:
        ```
            make
        ```
        One can choose `make -j` to compile in parallel. The executable program "xPaciae.x" and the input file "usu.dat" will be generated in folder "sim" (simulation). Enter the "sim".
    - Modify the input file of "usu.dat" according to your wish.
    - Run the program by the command:
        ```
            ./xPaciae.x
        ```
      One could use the following command to run the program in the background and record the time and log information.
        ```
            nohup time ./xPaciae.x > paciae.log &
        ```
 2. Use the "PACIAE.sh" shell-script to compile and link the code, generate usu.dat file and run program automatically.
    - Modify the "PACIAE.sh" file as needed.
    - If executable permissions are missing, grant executable permissions to "PACIAE.sh" by command (only once):
        ```
            chmod +x PACIAE.sh
        ```
    - Run the "PACIAE.sh" script by the command:
        ```
            ./PACIAE.sh
        ```
    -  "rms_analysis.f90" is a stand-alone program used in conjunction with the "PACIAE.sh" script to average the internal analysis files  "rms.out" generated by each PACIAE simulation.

    It is also worth mentioning that tasks can be submitted to computer clusters and supercomputers using the "PACIAE.sh" script (It is currently only available in SLRUM, LSF and PBS scheduling systems). More detailed information and usage please read the PACIAE.sh file.

## Maintainers

[@ArcsaberHep](https://github.com/ArcsaberHep)
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

The code now are hosted on https://github.com/ArcsaberHep/PACIAE4 and https://gitee.com/arcsaberhep/PACIAE4.
<br/>

The released code are available on https://github.com/ArcsaberHep/PACIAE4/releases and https://gitee.com/arcsaberhep/PACIAE4/releases.

## Released papers

 - **Recent version: A brief introduction to PACIAE 4.0. (waiting...)
 - ** An introduction to the parton and hadron cascade model PACIAE 3.0, [Phys. Rev. C 108 (2023) 6, 064909](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.108.064909) or [2309.05110 [hep-ph]](https://arxiv.org/abs/2309.05110)
<br/>

 - **Past version:** Revisiting the centrality definition and observable centrality dependence of relativistic heavy-ion collisions in PACIAE model, [Comput. Phys. Commun. 284 (2023) 108615](https://doi.org/10.1016/j.cpc.2022.108615) or [arXiv:2212.04087 [nucl-th]](https://doi.org/10.48550/arXiv.2212.04087)

 - PACIAE 2.2.1: An updated issue of the parton and hadron cascade model PACIAE 2.2, [Comput. Phys. Commun. 274 (2022) 108289](https://doi.org/10.1016/j.cpc.2022.108289) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - Announcement for the replacement of the PACIAE 2.1 and PACIAE 2.2 series, [Comput. Phys. Commun. 224 (2018) 417-418](https://doi.org/10.1016/j.cpc.2017.10.006) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - An upgraded issue of the parton and hadron cascade model, PACIAE 2.2, [Comput. Phys. Commun. 193 (2015) 89-94](https://doi.org/10.1016/j.cpc.2015.01.022) or [arXiv:1412.7579 [nucl-th]](https://doi.org/10.48550/arXiv.1412.7579)

 - PACIAE 2.1: An updated issue of the parton and hadron cascade model PACIAE 2.0, [Comput. Phys. Commun. 184 (2013) 1476-1479](https://doi.org/10.1016/j.cpc.2012.12.026) or [arXiv:1206.4795 [nucl-th]](https://doi.org/10.48550/arXiv.1206.4795)

 - PACIAE 2.0: An updated parton and hadron cascade model (program) for the relativistic nuclear collisions, [Comput. Phys. Commun. 183 (2012) 333-346](https://doi.org/10.1016/j.cpc.2011.08.021) or [arXiv:1104.1238 [nucl-th]](https://doi.org/10.48550/arXiv.1104.1238)

 - (PACIAE 1.0) Charge particle universal rapidity scaling in e+e-, pbar + p and Au + Au collisions at relativistic energies and its partonic origin, [J. Phys. G 32 (2006) 243-250](https://doi.org/10.1088/0954-3899/32/3/001); Influence of the partonic Pauli blocking on the hadronic final state in relativistic nucleus-nucleus collisions, [Phys. Rev. C 70 (2004) 034904](https://doi.org/10.1103/PhysRevC.70.034904)
<br/>

 - **Ancient version:** (JPCIAE) J/psi dynamical suppression in a hadron and string cascade model (or Formation time effect on J/psi dynamical nuclear suppression), [Phys. Rev. C 59 (1999) 2728-2733](https://doi.org/10.1103/PhysRevC.59.2728) or [arXiv:nucl-th/9803033](https://arxiv.org/abs/nucl-th/9803033); J/psi normal and anomalous suppressions in a hadron and string cascade model, [J.Phys.G 25 (1999) 1123-1133](https://doi.org/10.1088/0954-3899/25/6/302) or [arXiv:nucl-th/9809020](https://arxiv.org/abs/nucl-th/9809020); Inclusive and direct photons in S + Au collisions at 200A GeV/c, [Phys. Rev. C 61 (2000) 064905](https://doi.org/10.1103/PhysRevC.61.064905) or [arXiv:nucl-th/9904035](https://arxiv.org/abs/nucl-th/9904035).
<br/>

 - **Origin version:** LUCIAE 3.0: A New version of a computer program for firecracker model and rescattering in relativistic heavy ion collisions, [Comput. Phys. Commun. 116 (1999) 353](https://doi.org/10.1016/S0010-4655(98)00138-6) or [arXiv:nucl-th/9804001](https://doi.org/10.48550/arXiv.nucl-th/9804001)

 - An Event generator for the firecracker model and the rescattering in high-energy pA and AA collisions: LUCIAE version 2.0, [Comput. Phys. Commun. 90 (1995) 121-140](https://doi.org/10.1016/0010-4655(95)00066-O); Final state interactions in the (nuclear) FRITIOF string interaction scenario, [Z. Phys. C 70 (1996) 499](https://doi.org/10.1007/s002880050127).
<br/>

 - **Special version:** HYDRO-PACIAE, a hydrodynamic and transport hybrid model for ultra-relativistic heavy ion collisions, [J.Phys.G 40 (2013) 025102](https://doi.org/10.1088/0954-3899/40/2/025102) or [arXiv:1110.6704 [nucl-th]](https://doi.org/10.48550/arXiv.1110.6704)

## Relevant papers

 - A comprehensive guide to the physics and usage of PYTHIA 8.3, [SciPost Phys.Codeb. 2022 (2022) 8](https://scipost.org/10.21468/SciPostPhysCodeb.8) or [arXiv:2203.11601 [hep-ph]](https://arxiv.org/abs/2203.11601)

 - PYTHIA 6.4 Physics and Manual, [JHEP 05 (2006) 026](https://doi.org/10.1088/1126-6708/2006/05/026) or [arXiv:hep-ph/0603175 [hep-ph]](https://doi.org/10.48550/arXiv.hep-ph/0603175)

## Update Notes:

<!----------------------------------------------------------------------------->
### <font color=red> PACIAE 4 project started! </font>

<!----------------------------------------------------------------------------->
### <font color=red> PACIAE 3.0 is now released! </font>
