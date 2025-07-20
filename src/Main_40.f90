!! Main_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2025 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, ???? 2006 - July 2025.

!> This is the main program to administrate the Monte Carlo simulation for
!!  the relativistic hh, hA(Ah), AB, ll, lh(hl) and lA(Al) collisions.

!!                                             By Ben-Hao at CIAE on DD/MM/2006
!!                                  Last updated by An-Ke at CCNU on 20/07/2025



!ccccccccccccccccccccccccccc   Brief Introduction    ccccccccccccccccccccccccccc

!       PACIAE4.0 is developed from PACIAE3.0, which is from PACIAE2.0.
!        PACIAE2.0, in turn, is from PACIAE1.0, which is the rename of
!        JPCIAE, after JETSET integrated into PYTHIA. PACIAE1.0 is based on
!        JETSET7.4 and PYTHIA5.7 having implementation of partonic and hadronic
!        rescattering, respectively, before and after hadronization. PACIAE1.0
!        was next to LUCIAE, which is a extension of FRITIOF with performance
!        of hadronic rescattering after hadronization.
!       Ben-Hao Sa and An Tai began constructing the LUCIAE model at around
!        1993 with encouragement of Xu Cai. Subsequently, Zhong-Dao Lu,
!        Zhong-Qi Wang, Guang Song, Zhong-Han Feng, Xiao-Mei Li, Xiao-Long
!        Wang, Dai-Mei Zhou, Guo-Liang Ma, Zhi-Guang Tan, Bao-Guo Dong,
!        Yu-Liang Yan, Hai-Liang Ma, Sheng-Qin Feng, Gang Chen, Zhi-Lei She,
!        Yun Cheng, Du-Juan Wang, An-Ke Lei, and Liang Zheng, took part in
!        the PACIAE group one behind the other. Later on, Gao-Chan Yong,
!        Wen-Chao Zhang, Hua Zheng, and Li-Ling Zhu are joined with us.
!       An-Ke Lei is the main author of PACIAE4.0 and Dai-Mei Zhou is the
!        head of PACIAE group.
!                                                     Ben-Hao Sa on 05/04/2024

!ccccccccccccccccccccccccccc   Brief Introduction    ccccccccccccccccccccccccccc



        program main
!!      Administrates the MC simulation for relativistic hh, hA(Ah), AB,
!!       ll, lN(Nl) & lA(Al) collisions.
!
!       The program comprises Main_40.f90, Parini_40.f90, Parcas_40.f90,
!        Sfm_40.f90, Coales_40.f90, Hadcas_40.f90, Analy_40.f90, P_40.f,
!        Pythia8_fort_interface.f90, Pythia8CppInterface.cpp,
!        PaciaeUserHooks.cpp, Pythia8CppInterface.hpp, and
!        PaciaeUserHooks.hpp.
!
!       Main_40.f90: administrates the MC simulation.
!       Parini_40.f90: generates partonic initial state of colliding system.
!       Parcas_40.f90: performs parton rescattering, where 2->2 processes
!        are considered only, and the LO pQCD cross section or its regularized
!        approximation is used.
!       Sfm_40.f90: hadronization with Lund string fragmentation model
!       Coales_40.f90: hadronization with Monte Carlo coalescence model
!       Hadcas_40.f90: performs hadronic rescattering
!       Analy_40.f90: an example of event analysis subroutine, users are free
!        to replace it with your own one.
!       P_40.f (PYTHIA 6.4.28): for the generation of the partonic initial state
!        and/or string hadronization via PYTHIA 6.
!       Pythia8_fort_interface.f90: the Fortran interface of PACIAE
!        to C++-based PYTHIA 8.
!       Pythia8CppInterface.cpp: the C++ interface of PACIAE
!        to C++-based PYTHIA 8.
!       Pythia8CppInterface.hpp: the C++ interface header file.
!       PaciaeUserHooks.cpp: PACIAE-defined Userhooks and HIUserhooks for
!        the interaction with PYTHIA 8.
!       PaciaeUserHooks.hpp: PACIAE-defined Userhooks and HIUserhooks
!        header file.
!
!       Rms_analysis.f90 : standalone averaging program.
!
!       PYTHIA 8 has been interfaced in.
!       Both PYTHIA 6.4.28 and PYTHIA 8.3.1x are available now.
!
!       Note: for historical reasons, part of the Fortran programs
!        were written in case-insensitive mode. Type ": set ic + enter key"
!        before searching in Vi/Vim.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max

!       Final particle information will be stored in /PYJETS/ except for gammas,
!        which will be stored in /sgam/.

!       i_mode: = 1, low energy simulation A-framework;
!               = 2, PYTHIA-like simulation B-framework via PYTHIA 6;
!               = 3, PACIAE simulation C-framework via PYTHIA 6;
!               = 4, PACIAE simulation D-framework via PYTHIA 6;
!               = 5, PACIAE simulation B-framework via PYTHIA 8;
!               = 6, PACIAE simulation C-framework via PYTHIA 8;
!               = 7, dummy, reserved for future expansion;
!               = 8, PACIAE simulation B-framework via PYTHIA8/Angantyr;
!               = 9, PACIAE simulation C-framework via PYTHIA8/Angantyr.

!   See subroutine "read_input_parameter" and "usu.dat" file for more details.


!       Records the time the program starting running.
        call CPU_TIME( time_simulation_start )

        call read_input_parameter

        call identify_collision_system

        call open_output_file

        call set_simulation_parameter

        call initialize_global_event_variable

        call initialize_system

        call set_stable_particle

        call set_subprocess

        call calculate_optical_glauber

        call record_simulation_parameter


!*******************************************************************************
!----------------------------   Event Generating   -----------------------------
!       The current run count.
        iii = 0
!       Optional simulation stopping point.
        i_stop_point = INT( adj1(40) )

!       Loops over events.
100     iii = iii + 1

!       Generates an NN (hh, hA/Ah, AB, ll, lh/hl, lA/Al etc.) collision event.

!-------------------------------------------------------------------------------
        call initialize_single_event_variable
        call get_impact_parameter

!-------------------------------------------------------------------------------
        call do_initial_state_generation( i_error, i_pure_hadron )
        if( i_error == 1 ) goto 100

        if( i_stop_point == 1 .OR. i_stop_point == 3 )then
            call stop_event_at_this_point( i_stop_point )
            goto 666
        end if

!-------------------------------------------------------------------------------
        call do_partonic_rescattering( i_error, i_pure_hadron )
        if( i_error == 1 ) goto 100

        if( i_stop_point == 2 )then
            call stop_event_at_this_point( i_stop_point )
            goto 666
        end if

!-------------------------------------------------------------------------------
        call do_hadronization( i_error, i_pure_hadron )
        if( i_error == 1 ) goto 100

!-------------------------------------------------------------------------------
        call do_hadronic_rescattering( i_error )
        if( i_error == 1 ) goto 100

!-------------------------------------------------------------------------------
        call do_unstable_hadron_decay

!-------------------------------------------------------------------------------
666     continue
!       Moves partons from "PYJETS" to "sbe" if any, for the error check.
        if( i_stop_point > 2 .AND. i_mode /= 1 ) call remop

!-------------------------------------------------------------------------------
        call output_final_particle_information
        call analyze_event

!-------------------------------------------------------------------------------
!       Prompts the current number of events.
        if( MOD(iii,nout) == 0 )then
            write(*,*) "i-event =", iii
        end if
        open( 18, file = "nout.out", status = "unknown" )
        write(18,*) "iii=", iii
        close(18)

!       Simulates the next event.
        if( iii < neve ) goto 100
!----------------------------   Event Generating   -----------------------------
!*******************************************************************************


!-------------------------------------------------------------------------------
!       Statistics of subprocesses generated.
!       In Pythia8_fort_interface.f90.
        call PASTAT(1,0)

!-------------------------------------------------------------------------------
!       Releases memory from the heap for PYTHIA 8 objects.
        ! In Pythia8_fort_interface.f90.
        call PADELE

        call print_time_consumption( time_simulation_start )

        call close_output_file


        stop
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine read_input_parameter
!!      Opens the input file and reads parameters.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa7/ispmax,isdmax,iflmax,ispkf(100),n_bin_hist,asd(10), &
         afl(100,10,2),i_y_or_eta
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa29/i_color_reconnection,NCR,parp2
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sa38/ prob_ratio_q(6), am(6), amqq(6)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/coal1/bmrat,i_mm
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)
!       For the independent Optical Clauber calculation.
        common/Glauber1/ kjp23, kjp24
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        open( 11, file = "usu.dat", status = "unknown" )


!-------------------------------------------------------------------------------
!------------------------------   Input Reading   ------------------------------
!       Reads input variables for event generation.
        read(11,*) neve, nout, nosc
        read(11,*) KF_proj, KF_targ
        read(11,*) i_mode, ifram, nchan
        read(11,*) win, energy_B
        read(11,*) bmin, bmax, psno, nmax
            b_min = bmin
            b_max = bmax
        read(11,*) kjp21, x_ratio, decpro
!       Reads input variables for event analyses.
        read(11,*) ispmax, isdmax, iflmax, n_bin_hist, i_y_or_eta
            KF_woDecay = 0
            ispkf = 0
            afl = 0D0
        read(11,*) ( KF_woDecay(i), i=1,10,1 )
        read(11,*) ( KF_woDecay(i), i=11,ispmax,1 )
        read(11,*) ( ispkf(i), i=1,10,1 )
        read(11,*) ( ispkf(i), i=11,ispmax,1 )
        read(11,*) ( asd(i), i=1,isdmax,1 )
        if( iflmax > 0 )then
            do kk=1,ispmax
                do i=1,iflmax
                    read(11,*) ( afl(kk,i,j), j=1,2,1 )
                end do
            end do
        end if
!       Reads input variables for the event generation.
        read(11,*) ddt, i_mm
        read(11,*) iparres, i_inel_proc, i_time_shower
            para = 0D0
        read(11,*) para(7), ttaup, taujp, para(10)
        read(11,*) i_sigma_AQM, kjp20
        read(11,*) para1_1, para1_2, ( para(i), i=2,5,1 )
        read(11,*) ( para(i), i=13,16,1 )
        read(11,*) ( adj1(i), i=1,10,1 )
        read(11,*) ( adj1(i), i=11,20,1 )
        read(11,*) ( adj1(i), i=21,30,1 )
        read(11,*) ( adj1(i), i=31,40,1 )
        read(11,*) kjp22, kjp23, kjp24, i_color_reconnection, i_tune
        read(11,*) parecc, smadel, cp0, cr0, seco
        read(11,*) i_deex, n_deex_step, i_pT_coal, i_pT_endpoint
        read(11,*) a_FF, aPS_c, aPS_b, bmrat
        read(11,*) ( prob_ratio_q(i), i=1,6,1 )
!------------------------------   Input Reading   ------------------------------
!-------------------------------------------------------------------------------


        close(11)


!       neve: number of events to be generated
!       nout: output the event logs (main.out, rms.out) per nout events
!       nosc: = 0, no OSCAR output (oscar.out), see subroutine "oscar"
!             = 1, OSCAR1997A (final_id_p_x, just PACIAE final output)
!             = 2, OSCAR1999A (full_event_history)
!             = 3, OSCAR2013A (full_event_history), dummy now

!       KF_proj, KF_targ: standard PDG code (id/KF code) of the incoming
!                         projectile and target particles, l, h & A,
!                E.g. 11 = e-, -11 = e+, 2212 = p, 2112 = n,
!                     1000822080 = Pb(208,82), 10000791970 = Au(197,79), ...
!        Cf.
!        1. PDG RPP2024: Phys.Rev.D 110 (2024) 3, 030001
!         https://pdg.lbl.gov/2024/reviews/rpp2024-rev-monte-carlo-numbering.pdf
!        2. PYTHIA 6.4 munual:  JHEP 05 (2006) 026
!         https://arxiv.org/abs/hep-ph/0603175
!        3. PYTHIA 8 online manual:
!         https://pythia.org//latest-manual/ParticleData.html
!        4. The end part of "usu.dat".

!       i_mode: = 1, low energy simulation A-framework;
!               = 2, PYTHIA-like simulation B-framework via PYTHIA 6;
!               = 3, PACIAE simulation C-framework via PYTHIA 6;
!               = 4, PACIAE simulation D-framework via PYTHIA 6;
!               = 5, PACIAE simulation B-framework via PYTHIA 8;
!               = 6, PACIAE simulation C-framework via PYTHIA 8;
!               = 7, dummy, reserved for future expansion;
!               = 8, PACIAE simulation B-framework via PYTHIA8/Angantyr;
!               = 9, PACIAE simulation C-framework via PYTHIA8/Angantyr.

!       ifram: = 0, for fixed target (FXT)
!              = 1, for collider (CMS)
!              = 2, back-to-back asymmetric system (BKB)
!       win: = incident momentum +pz if ifram=0 (GeV)
!            = sqrt(s) if ifram=1 (GeV)
!            = energy_A if ifram=2, for back-to-back system (GeV)
!       energy_B: if ifram=2; energy_A, +Z direction; energy_B, -Z direction

!       nchan: = 0 , inelastic (INEL)
!              = 1 , non-single diffractive (NSD)
!              = 2 , Drell-Yan process
!              = 3 , J/psi production (color singlet)
!              = 4 , c-cbar and b-bbar production
!              = 5 , prompt photon
!              = 6 , soft QCD
!              = 7 , single W+/- production
!              = 8 , hard QCD
!              = 9 , single Z0 production
!              = 10, inelastic， non-diffractive (minimum-bias, MB)
!              = other, user-defined subprocessses in "pythia6_extra.cfg"
!                       for PYTHIA 6 or "pythia8_extra.cfg" for PYTHIA 8.

!       ifram: = 0, for fixed target (FXT)
!              = 1, for collider (CMS)
!              = 2, back-to-back asymmetric system (BKB)
!       win: = incident momentum +pz if ifram=0 (GeV)
!            = sqrt(s) if ifram=1 (GeV)
!            = energy_A if ifram=2, for back-to-back system (GeV)
!       energy_B: if ifram=2; energy_A, +Z direction; energy_B, -Z direction

!       bmin, bmax : minimum and maximum impact parameters, bmin=bmax means
!                    definite impact parameter
!       psno: b parameter sampling method for [bmin, bmax].
!             = 0, fixed b (unit-weight, b = bmin, phi = 0).
!               1, systematic sampling (unit-weight, phi = 0).
!               2, random sampling (unit-weight, phi = 0).
!               3, for Angantyr ion mode, original Gaussian sampling.
!                  This option will generate weighted events and set
!                  a large range of 0 < b < 20+ fm automatically
!                  (weighted, phi = [0, 2*pi] ).
!               4, for Angantyr ion mode, fixed b
!                  (unit-weight, b = bmin, phi = [0, 2*pi] ).
!               5, for Angantyr ion mode, systematic sampling
!                  (unit-weight, phi = [0, 2*pi] ).
!               6, for Angantyr ion mode, random sampling
!                  (unit-weight, phi = [0, 2*pi] ).
!       2*nmax: the number of intervals segmented in [bmin,bmax] for the
!               systematic sampling (psno=1)

!       kjp21: = 0, without hadron rescattering;
!              = 1, with hadron rescattering via PACIAE;
!       x_ratio: ratio of inela. cross section to total cross section,
!                i.e. PARAM(6) in hadcas_40.f90, with default value of
!                D=0.85 for the high energy and D=0.1 for the low energy
!       decpro is Delta particle decay probability in low energy A-framework

!       ispmax: the maximum number of particle KF code to be counted (max=100)
!       isdmax: the maximum number of statistical distributions (max=10)
!       iflmax: the maximum number of the kinematic cuts (max=2, pT & y/eta)
!       n_bin_hist: the maximum number of bins for histograms.
!       i_y_or_eta: select y or eta in partial phase-space statistics.
!                   = 0 , y
!                   = 1 , eta

!       KF_woDecay: KF code of particles specified not to decay (max=100).
!       ispkf(i): KF code of i-th particle to be counted (max=100)
!       asd(i=1,isdmax): bin width, interval of the i-th distribution
!           for pp, pbarp, pA(Ap), AB etc.
!               i=1: rapidity distribution (dN/dy v.s. y)
!                =2: invariant transverse momentum distribution
!                    (1/pT*dN/dpT v.s. pT)
!                =3: pseudorapidity distribution (dN/deta v.s. eta)
!                =4: transverse mass distribution (1/mT*dN/dmT v.s. mT)
!                =5: event-wise multiplicity distribution (dNev/dmult)
!                =6: transverse momentum distribution (dN/dpT v.s. pT)
!                =7: mean transverse momentum distribution (<pT> v.s. mult)
!           for ep, nu_ep, etc.
!               i=1: Q^2=-q^2 (fq2 in code) distribution
!                =2: W^2 (w21) distribution
!                =3: y (yyl) distribution
!                =4: p_h (pph) distribution
!                =5: z (zl) distribution
!       afl(j,i,1): lower-boundary cut of i-th window for j-th particle
!       afl(j,i,2): upper-boundary cut of i-th window for j-th particle
!           for pp, pbarp, pA(Ap), AB etc.
!               i=1, rapidity/pseudorapidity window (depends on i_y_or_eta)
!                =2, transverse momentum
!           for ep, nu_ep, etc.
!               i=1, Q^2=-q^2 window
!                =2, W^2
!                =3, y
!                =4, p_h (hadron momentum)
!                =5: z

!       ddt : time accuracy used in parini, i.e. the minimum distinguishable
!             collision time interval used in partonic initiation
!       i_mm: only the hadrons below numb(i_mm) (cf. 'filt' in parini.f90)
!        join in updating the hh collision time list (cf. parini.f90).
!        Originally i_mm=2, i.e. only p (proton) and n (neutron) join in
!        updating the collision time list in parton initialization;
!        if i_mm=6, p, n, pbar, nbar, pi+ and pi- will join in updating
!        the collision time list in parton initialization.

!       iparres: =0, considers elastic parton-parton collisions only
!                    in parcas_40.f90, with diqarks breaking-up.
!                =1, elastic + inelastic collisions, with diqarks breaking-up.
!                =2, elastic collisions only, without diqarks breaking-up.
!                =3, elastic + inelastic collisions, without diqarks breaking-up
!       i_inel_proc: =6, with inelastic processes 4, 6, and 7 (parcas_40.f90)
!                    =7, with inelastic process 7 only (parcas_40.f90)
!       i_time_shower: =0, w/o final state time-like parton shower if
!                          iparres=1/3
!                      =1, w/ final state time-like parton shower if iparres=1/3

!       para(7): proper formation time in rest-frame of particle
!       ttaup: extra factor of the formation time of particles
!              (protons in parini, all of hadrons in hadcas)
!       taujp: extra factor of the formation time of J/psi in parini
!       para(10): largest allowed size of partonic (hadronic) rescattering
!                 region which is product of para(10) and target radius

!       i_sigma_AQM: = 1, uses input total cross sections of piN, KN, pi+pi,
!                         Jpsi(psi')+N, Jpsi(psi')+N, Jpsi(psi')+pi/rho,
!                         and AQM cross sections for D+N/pi/rho in hadcas.
!                         Ignores the channels which are not well defined.
!                    = 2, uses input cross sections of piN, KN, ... ,
!                         AQM cross sections for D+N/pi/rho and the channels
!                         which are not well defined (treated as elastic).
!                    = 3, uses AQM cross sections of piN, KN,... and D+N/pi/rho.
!                         Ignores the channels which are not well defined.
!                    = 4, uses AQM cross sections of piN, KN, ..., D+N/pi/rho,
!                         and the channels which are not well defined (treated
!                         as elastic).

!       kjp20: = 1, constant cross sections for piN -> K+Sigma/Lambda/Delta,
!                   piN -> pi+Delta and NN -> N+Delta in hadcas.
!              = 0, energy dependent cross sections for them in hadcas.

!       para1_1: NN total cross section, used in parini_40.f90
!                <= 0, it will be evaluated in the program (PAXTOT)
!                > 0, accept the input one
!       para1_2: NN total cross section, used in hadcas_40.f90
!       para(2): total cross-section of pi + N
!       para(3): total cross-section of K + N
!       para(4): total cross-section of pi + pi
!       para(5): cross-section of pipi -> KKbar

!       para(13): total cross-section of J/psi + N
!       para(14): total cross-section of J/psi + pi/rho
!       para(15): total cross-section of psi' + N
!       para(16): total cross-section of psi' + pi/rho

!       adj1(i), i=
!       1: K factor used in parton rescattering; K=0: no parton rescattering
!       2: parameter alpha_s ('as' in program) in parton rescattering
!       3: parameter (mu)^2 ('tcut' in program): to avoid divergence in
!          calculating parton-parton differential cross section in parcas_40.f90
!       4: parameter idw, # of integration intervals in parcas_40.
!       5: choice of nuclear shadowing or nPDF for partons in nuclei,
!          available only for the ion-collisions.
!              =0, without nuclear shadowing or nPDF.
!          For PYTHIA 6:
!              =1, Wang's nuclear shadowing (PLB 527 (2002) 85).
!              =other, same as 1.
!          For PYTHIA 8:
!              =-1, only Isospin effect.
!              = 1, EPS09, LO nPDF.
!              = 2, EPS09, NLO nPDF.
!              = 3, EPPS16, NLO nPDF.
!              = other, preset 208Pb EPS09 LO nPDF.
!          EPS/EPPS nPDF requirs for grid files. There is only one default
!              EPS09 LO grid file of 208Pb in PYTHIA 8. One needs to add other
!              grid files manually for other nulei or NLO nPDF.
!       6: parameter 'a', i.e. PARJ(41), in Lund string fragmentation function
!       7: parameter 'b', i.e. PARJ(42), in Lund string fragmentation function
!       8: 'a' in the Lund deexcitation function of coal model. See adj1(29).
!       9: 'b' in the Lund deexcitation function of coal model. See adj1(29).
!       10: PARP(31), K factor in PYTHIA
!       11: time accuracy used in the hadronic rescattering.
!       12: model of hadronization:
!           =0, string fragmentation model (SFM).
!           =1, coalescence model (Coal);
!           =2, with gluon splitting & quark deexcitation before parcas and
!               hadronized by coalescence model. (same as adj1(12)=1 + kjp22=2)
!           =3, with gluon splitting & quark deexcitation before parcas and
!               hadronized by string fragmentation model  (not available for
!               PYTHIA8 mode).
!           The subclasses please see kjp22.
!       13: not used.
!       14: not used.
!       15: default string tension (D= 1 GeV/fm ~ 0.2 GeV^2)
!       16: allowable number of generations of q/qbar deexcitation
!           in coales_40.f90.
!           0 = without deexcitation.
!           1 = deexcites only the original 1-th q/qbars.
!           2 = same as 1, and deexcites the 2-th q/qbars excited from 1-th.
!           3/4/... = ...
!       17: threshold energy in deexcitation of energetic quarks in coales.f90
!       18: =0, rest partons hadronize by string fragmentation via PYTHIA 6
!           =1, rest partons hadronize by coalescence
!           <0, no treatment for the rest partons
!       19: time accuracy used in parcas_40.f90 ('dddt' in program)
!       20: =0 exact pQCD parton-parton cross section
!           =1, limited and regularized parton-parton cross section (B. Zhang,
!               Comput. Phys. Commun. 109 (1998) 193)
!           =2, the same as 0 but flat scattering angle distribution is assumed
!           =3, the same as 1 but flat scattering angle distribution is assumed
!       21: =0, without phase space requirement in coales_40.f90
!           =1, with complete phase space requirement
!           =2, with requirement in spatial phase space only
!           =3, with requirement in momentum phase space only
!       22: critical value of the product of radii both in coordinate and
!           momentum phase space (4 is assumed) used in coales_40.f90
!       23: switch for chiral magnetic effect (CME) induced charge seperation
!           =0: CME off
!           =1: CME on
!       24: = 'tl0' , the virtuality cut in time-like radiation in parcas_40.f90
!             4*tl0 is assumed
!       25: Lambda_QCD in parcas_40.f90
!       26: selection of random number seed
!           =0, default seed = 20240116, can be used for debug
!           =1, seed from the real-time colock
!           >1, sets seed as adj1(26)
!       27: largest momentum allowed for particle into rescatterings ('dpmax')
!       28: largest position allowed for particle in hadcas ( which is 1. in
!            usu.dat, but is recalculated in the running,
!            drmax=para(10)*max(rnt,rnp)*adj1(28) )
!       29: For coal, the function used to sample deexcited daughter qqbar-pair
!             energy fraction z taking from mother.
!            =0 : random z
!            =1 : Lund symmetric function, see PARJ(41) - PARJ(45)
!            =2 : Field–Feynman + Peterson/SLAC, see PARJ(51) PARJ(59)
!            =3 : Lund + Peterson/SLAC (light flavor + heavier)
!            =4 : Lund + Bowler
!            =5 : as = 4, but interpolate for c and b; see PARJ(46) and PARJ(47)
!           For coal, when using 11, 12 or 13, it is the simple Lund/FF/PS.
!       30: =1: distributes participant nucleons in overlapping areas forcely
!           =0: no more requirements
!       31: not used
!       32: not used
!       33: not used
!       34: PARJ(21) (Lund) width of px/py/pT sampling in PYPTDI.
!       35: width of px/py/pT sampling in the coalescence model. See i_pT_coal.
!       36: =0 without phenomenological parton energy loss in parcas_40.f90
!           =1 with (Do not use this option with the paronic rescattering
!                    openning altogether. If so, double-counting!)
!       37: the coefficient ('c') in phenomenological parton energy loss
!       38: pt cut in phenomenological parton energy loss
!       39: not used
!       40: =1 simulation ends after partonic initiation
!           =2 simulation ends after partonic rescattering
!           =3 hadronization by coalescence model follows partonic initiation
!           =4 simulation ends after the entire simulation

!       kjp22: hadronization subclass.
!         For Lund string fragmentation model (sfm):
!             =1, fragmentation string by string, with the single string
!                 structure induced variable string tension and PARJ(1) etc.
!             =2, fragmentation string by string, with the multiple string
!                 interaction induced variable string tension and PARJ(1) etc.
!             =3, fragmentation string by string, with the (single string
!                 structure + multiple string interaction) induced variable
!                 string tension and PARJ(1) etc.
!             =4, fragmentation of all strings at once with the
!                 default string tension and parj(1) etc.
!             =5, fragmentation NN pair by pair with the
!                 default string tension and parj(1) etc.
!           For PYTHIA 8 only:
!             =6, fragmentation by the Rope Hadronization with Flavour Ropes.
!             =7, fragmentation by the Rope Hadronization with String Shoving.
!             =8, fragmentation by Rope Hadronization with Flavour Ropes
!                 + String Shoving.
!             =9, Thermodynamical String Fragmentation.
!         For PACIAE coalescence model (coal):
!             =1, mode 1, gluon splitting and energetic quark
!                 deexcitation AFTER the partonic rescattering,
!                 via the traversal coalescence.
!             =2, mode 2, gluon splitting and energetic quark
!                 deexcitation BEFORE the partonic rescattering
!                 via the traversal coalescence (same as adj1(12)=2 ).
!             =3, mode 3, as = 1, but with only the quark deexcitation.
!             =4, mode 4, as = 2, but with only the quark deexcitation.
!             =5, as = 2, to keep the consistency of "kjp22" switch.
!             =6, mode 6, as = 1, but via the random coalescence.
!             =7, mode 7, as = 2, but via the random coalescence.
!             =8, mode 8, as = 3, but via the random coalescence.
!             =9, mode 9, as = 4, but via the random coalescence.

!       For the independent Optical Glauber model.
!       kjp23: = 1 Npart is calculated by geometric model
!       kjp23: = 2 Npart is calculated by Glauber model
!       kjp24: = 1 sharp sphere in Glauber model
!       kjp24: = 2 Woods-Saxon in Glauber model

!       i_color_reconnection: Selection of color-reconnection (CR) model.
!                   = 0, No color-reconnection.
!       In PYTHIA 6 == "MSTP(95)" :
!                   'S': Seattle model, 'P': Paquis model.
!                   = 1,  Old cut-and-paste style CR, handled in PYMIHK.
!                   = 2,  Annealing Type I(no gg loops); hadron-hadron only.
!                   = 3,  Annealing Type I(no gg loops); all beams.
!                   = 4,  Annealing Type II(gg loops)  ; hadron-hadron only.
!                   = 5,  Annealing Type II(gg loops)  ; all beams.
!                   = 6,  Annealing Type S             ; hadron-hadron only.
!                   = 7,  Annealing Type S             ; all beams.
!                   = 8,  Annealing Type P             ; hadron-hadron only.
!                   = 9,  Annealing Type P             ; all beams.
!                   = 11, Soft Color Interaction (SCI model).
!                   = 12, SCI without diffraction (no inter-remnant CR).
!                   = 13, Generalized Area Law (GAL model).
!                   = 14, GAL without diffraction (no inter-remnant CR).
!       In PYTHIA 8 and PYTHIA8/Angantyr :
!                   < 15, The default PYTHIA 8 MPI-based scheme. (MPI-CR)
!                   = 15, The new more QCD based scheme. (QCD-CR, QCD-CR-BLC)
!                   = 16, The new gluon-move model.
!                   = 17, The SK I e+e- CR model.
!                   = 18, The SK II e+e- CR model.

!       i_tune: basic PYTHIA parameter set, tune number.
!               =0, "Perugia 2011 Tune" as the basis for PYTHIA6;
!                   "Monash 2013 Tune" as the basis for PYTHIA8/Angantyr.
!               >0, correspinding tune "i_tune".

!       parecc: a parameter converting initial spatial space eccentricity
!               to final momentum space ellipticity
!       smadel: small perpurbation of ellipse relative to circle for pT
!       cp0,cr0: parameters in the parameterization of the multiple string
!                effect (kjp22= 2 / 3).
!       seco: parameter in the popcorn mechanism for the correction of PARJ(1)

!       i_deex: the deexcitation scheme used in the coalescence model
!               = 1, light-cone momentum scheme
!               = 2, energy scheme
!               = 3, logitudinal momentum pz scheme
!       n_deex_step: the number of deexcitation steps per q/qbar
!       i_pT_coal: the pT sampling method of the daughter qqbar pair in coal
!                  = 0, 0 pT (collinear)
!                  = 1, Gaussian pT (px and py) with width adj1(35)
!                  = 2, Exponential pT with width adj1(35)
!       i_pT_endpoint: generation of transverse momentum for endpoint/mother
!                      quarks in coal
!                      = 0, no transverse momentum for endpoint quarks
!                      = 1, endpoint quarks obtain transverse momenta like
!                           exited qqbars

!       a_FF: parameter for light u, d and s in Field-Feynman function
!       aPS_c: parameter for c in Petersono/SLAC
!       aPS_b: parameter for b in P/S function
!       bmrat: parameter related to the ratio of baryon to meson in Coal

!       prob_ratio_q: probability ratio of excited qqbar
!                     u-ubar : d-dbar : s-sbar : c-cbar : b-bbar : t-tbar
!                     used in the gluon splitting and the quark deexcitation.


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine identify_collision_system
!!      Identfies and sets the collision system and some
!!       mistake-proofing measures.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        LOGICAL IS_PYTHIA8,IS_NUCLEUS
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa30/vneump,vneumt,mstptj
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max

!       KF_proj, KF_targ: standard PDG code (id/KF code) of the incoming
!                         projectile and target particles, l, h & A,
!                E.g. 11 = e-, -11 = e+, 2212 = p, 2112 = n,
!                     1000822080 = Pb(208,82), 10000791970 = Au(197,79), ...
!
!                The standard PDG (hyper)nucleus codes are 10-digit:
!                        +-10LZZZAAAI.
!                +- denotes (hyper)nucleus and anti-one.
!                L = n_Lambda, total number of strange quarks (Lambda hyperons).
!                Z = n_p, total charges (protons).
!                A = n_p + n_n + n_Lambda, total baryon number.
!                    n_n, total number of neutrons.
!                I = isomer level, = 0: groud state,
!                                  > 0: excitations. m, n, p, q = 1, 2, 3, 4.
!                Cf. PDG RPP "Monte Carlo particle numbering scheme".
!                Here for normal ions, we use
!                    id = 100ZZZAAA0 = 1,000,000,000 + Z*10,000 + A*10 + 0.
!                For example, deuteron D/2H: 1000010020; 197Au: 1000791970;
!                                     208Pb: 1000822080, etc.

!       ipden: =0, if projectile is proton/anti-proton/neutron
!              =1, if projectile is nucleus
!              =2, for e+e- collisions  ! 180921 yan
!              =11, projectile is e- (e+)
!              =12, projectile is nu_e (nu_ebar)
!              =13, projectile is mu- (mu+)
!              =14, projectile is nu_mu (nu_mubar)
!              =15, projectile is tau- (tau+)
!              =16, projectile is nu_tau (nu_taubar)
!              =0, any other subatomic partices in principle
!       itden: =0, if target is proton/anti-proton/neutron
!              =1, if target is nucleus
!              =2, for e+e- collisions  ! 180921 yan
!              =...
!              =0, any other subatomic partices in principle

!       nap (nzp) : # of nucleons (protons) in projectile nucleus
!       nat (nzt) : # of nucleons (protons) in target nucleus
!       for e-A: formally sets nap=1,nzp=-1,ipden=11,itden=1, kf=11;
!       for e+A: formally sets nap=1,nzp=1,ipden=11,itden=1, kf=-11;
!       for nu_eA: formally sets nap=1,nzp=-1,ipden=12,itden=1, kf=12;
!       for nu_ebarA: formally sets nap=1,nzp=1,ipden=12,itden=1, kf=-12;
!       If projectile is lepton, set nap =1
!       Here formally setting nzp=-1 is presently, when one uses it later,
!        one should make nzp=iabs(nzp)

!       ipden, itden, nap, nzp, nat and nzt will be assigned automatically
!        according to the incoming beams.

!       mstptj: =0, sets MSTP(111) = 0 to PYTHIA,
!                for PACIAE simulation developed from partonic initial stage,
!                to partonic rescattering, hadronization, and to hadronic
!                rescttering stage
!       mstptj: =1, PYTHIA-like simulation without partonic
!                   rescatterings and low energy simulation


!-------------------------------------------------------------------------------
!------------------------------   Beams Setting   ------------------------------
        KF_proj_abs = ABS( KF_proj )
        KF_targ_abs = ABS( KF_targ )
        nap = 1
        nzp = 1
        ! Non-nucleus.
        if( .NOT.IS_NUCLEUS(KF_proj) ) nzp = SIGN(1, PYCHGE( KF_proj ))
        nat = 1
        nzt = 1
        ! Non-nucleus.
        if( .NOT.IS_NUCLEUS(KF_targ) ) nzt = SIGN(1, PYCHGE( KF_targ ))
        ! Flags of collision system the PACIAE defined.
        ipden = 0
        itden = 0

!       Projectile.
        select case( KF_proj_abs )
        ! e, nu_e, mu, nu_mu, tau and nu_tau
        case( 11:16 )
            ipden = KF_proj_abs
        ! e+ + e-
            if( KF_proj_abs == 11 .AND. KF_targ_abs == 11 &
            .AND. KF_targ*KF_targ < 0 )  ipden = 2
        ! A
        case( 1000000000: )
            ipden = 1
            nap = MOD( MOD( KF_proj_abs, 10000000 ), 10000 ) / 10
            nzp = MOD( KF_proj_abs, 10000000 ) / 10000
        end select

!       Target.
        select case( KF_targ_abs )
        ! e, nu_e, mu, nu_mu, tau, and nu_tau
        case( 11:16 )
            itden = KF_targ_abs
        ! e+ + e-
            if( KF_proj_abs == 11 .AND. KF_targ_abs == 11 &
            .AND. KF_targ*KF_targ < 0 )  itden = 2
        ! B
        case( 1000000000: )
            itden = 1
            nat = MOD( MOD( KF_targ_abs, 10000000 ), 10000 ) / 10
            nzt = MOD( KF_targ_abs, 10000000 ) / 10000
        end select
!------------------------------   Beams Setting   ------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Mistake-proofing   -----------------------------
!       Automatic assignment of values avoiding some manual input mistakes
!        in "usu.dat".

!       For subatomic-particle + subatomic-particle.
        if( .NOT.IS_NUCLEUS(KF_proj) .AND. .NOT.IS_NUCLEUS(KF_targ) )then
            bmin = 0D0
            bmax = 0D0
            psno = 0D0
            if( i_mode == 8 .OR. i_mode == 9 ) psno = 3D0
            ! Without nuclear shadowing.
            adj1(5) = 0D0
        end if

!       The hadronization model.
        i_had_model = INT( adj1(12) )

!       Lets Coal5 and special versions of SFM and Coal2 work well.
        if( i_had_model == 1 )then
            if( kjp22 < 0 .OR. kjp22 > 9 ) kjp22 = 2
            if( kjp22 == 5 ) kjp22 = 2
        else if( i_had_model == 2 )then
            adj1(12) = 1D0
            kjp22 = 2
        else if( i_had_model == 3 )then
            if( kjp22 < 0 .OR. kjp22 > 4 ) kjp22 = 4
            ! PYTHIA 6 only.
            i_mode = 3
        end if

!       Lets "kjp22" of SFM with PYTHIA 6 works well.
        if( i_had_model == 0 .AND. .NOT.IS_PYTHIA8( i_mode ) )then
            if( kjp22 < 0 .OR. kjp22 > 5 ) kjp22 = 5
!#TODO(Lei20241014): for SFM with PYTHIA 8 , kjp22 >= 5, temporarily.
        else if( i_had_model == 0 .AND. IS_PYTHIA8( i_mode ) )then
            if( kjp22 < 5 .OR. kjp22 > 9 ) kjp22 = 5
        end if

!       Lets D-framework with SFM work well. Special simulation modes.
        if( i_mode == 4 .AND. ( i_had_model == 0 .OR. i_had_model == 3 ) &
            .AND. ( kjp22 < 0 .OR. kjp22 > 4 ) ) kjp22 = 4

!       Fragmentes partons directly after partonic initialization, or not.
        mstptj = 0
        if( i_mode == 1 .OR. i_mode == 2 &
            .OR. i_mode == 5 .OR. i_mode == 8 )then
            mstptj = 1
            ! Without partonic rescattering.
            adj1(1) = 0D0
        end if

!       Simulation mode warning.
        if( i_mode < 1 .OR. i_mode == 7 .OR. i_mode > 9 )then
            write(*,*) "Error, unrecognized i_mode:", i_mode, ", STOP!"
            write(22,*) "Error, unrecognized i_mode:", i_mode, ", STOP!"
            stop
        end if
        if( i_mode == 4 .AND. ipden /= 1 .AND. itden /= 1 )then
            write(*,*) "D-framework is only available for the nucleus!"
            write(22,*) "D-framework is only available for the nucleus!"
            stop
        end if

!       Unneccesary in A- and D-framework. (?)
        if( i_mode == 1 .OR. i_mode == 4 ) adj1(30) = 0D0

!       Lets "psno=3" work well for non-Angantyr modes.
        if( INT(psno) == 3 .AND. i_mode /= 8 .AND. i_mode /= 9 ) psno = 2D0
        if( INT(psno) < 0 .OR. INT(psno) > 6 ) psno = 2D0

!       Elastic or el. + inel. collisions and w/ or w/o diquarks breaking-up.
        if( iparres < 0 .OR. iparres > 3 ) iparres = 0

!       For the debug.
        if( i_had_model == 1 .AND. i_deex > 99 )then
            if( iparres == 0 ) iparres = 2
            if( iparres == 1 ) iparres = 3
        end if
!----------------------------   Mistake-proofing   -----------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine open_output_file
!!      Opens the output files "main.out", "rms0.out", "oscar.out", etc.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc


!-------------------------------------------------------------------------------
!------------------------------   Files Openning   -----------------------------
!       Opens a PACIAE and PYTHIA 6 log file.
        open( 22, file = "main.out", status = "unknown" )
!       Sets the output file of PYTHIA 6 (original unit = 6).
        MSTU(11) = 22
!       Opens a PACIAE analysis output file.
        open( 19, file = "rms0.out", status = "unknown" )
!       Opens and prints the OSCAR file header if nosc > 0.
        if( nosc > 0 )then
            open( 34, file = "oscar.out", status = "unknown" )
            ! In main.f90.
            call oscar(-1)
        end if
!------------------------------   Files Openning   -----------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine close_output_file
!!      Closes the output files "main.out", "rms0.out", "oscar.out", etc.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc


!-------------------------------------------------------------------------------
!------------------------------   Files Closing   ------------------------------
        close(22)
        close(19)
        if(nosc > 0) close(34)
!------------------------------   Files Closing   ------------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine set_simulation_parameter
!!      Sets basic parameters of the simulation.
!!      Parameters are categorized to "PACIAE-" and "PYTHIA-" related.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        common/pycidat2/kfmaxt,nont2,PARAM(20),weigh(600)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa29/i_color_reconnection,NCR,parp2
        common/sa30/vneump,vneumt,mstptj
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        CHARACTER*100 PARAM_PYTHIA6, PACIAE_INPUT


!-------------------------------------------------------------------------------
!------------------------   PACIAE Parameters Setting   ------------------------
!       Calculates the total cross section of NN (pp) for pA/Ap/AB collisions.
        if( para1_1 <= 1D-5 )then
            ECM = win
            if( ifram == 0 )then
                pzA = win
                pzB = 0D0
                dm_proton = PYMASS(2212)
                eA = SQRT( pzA**2 + dm_proton**2 )
                eB = dm_proton
                ECM = SQRT( (eA + eB)**2 - ( pzA + pzB )**2 )
            else if( ifram == 2 )then
                eA = win
                eB = energy_B
                dm_proton = PYMASS(2212)
                if( eA < dm_proton ) eA = dm_proton
                if( eB < dm_proton ) eB = dm_proton
                pzA = 0D0
                pzB = 0D0
                pzA =  SQRT( eA**2 - dm_proton**2 )
                pzB = -SQRT( eB**2 - dm_proton**2 )
                ECM =  SQRT( (eA + eB)**2 - ( pzA + pzB )**2 )
            end if
            ! In Pythia_fort_interface.f90.
            call PAXTOT( 2212, 2212, ECM, sigma_NN_tot )
            para1_1 = sigma_NN_tot
        end if

        ! In hadcas_40.f90.
        call PYCIDATA
!       PARAM(1) = para(1)
        PARAM(2) = para(2)
        PARAM(4) = para(4)
        PARAM(6) = x_ratio
!       x_ratio: ratio of inela. cross section to total cross section,
!        i.e. PARAM(6) in hadcas_40.f90, with default value of
!        D=0.85 for high energy and D=0.1 for low energy
        if( i_mode == 1 )then   ! 310723 Lei
            if( win > 2.015 .AND. win < 3.0 )then
                PARAM(6) = 1.35 * (win - 2.015)**2 &
                         / ( 0.015 + (win - 2.015) ) * x_ratio * 10.0
            end if
!           x_ratio = 1.35*(win-2.015)**2 / (0.015 + (win-2.015)) revised by yan
!       default x_ratio=0.1 for A-frame, so there are x_ratio*10.0, 010423 yan
        end if

!       For the partonic initialization.
        PARAM(8)  = ddt
!       For the hadronic rescattering.
        PARAM(7)  = para(7)
        PARAM(10) = para(10)
        ! Up to D*0bar in hadcas, for time-saving.
        kfmaxt = 46
        ! Up to Sigma_bbar- in hadcas, almost full heavy-flavor ground states.
        if( i_sigma_AQM ==  2 .OR. i_sigma_AQM == 4 ) kfmaxt = 92
!       total cross-section of J/psi + N
        PARAM(13) = para(13)
!       total cross-section of J/psi + pi/rho
        PARAM(14) = para(14)
!       total cross-section of psi' + N
        PARAM(15) = para(15)
!       total cross-section of psi' + pi/rho
        PARAM(16) = para(16)
!       For coalesecence.
        imc = 30
        ibc = 45

        pio = 3.141592653589793D0
!------------------------   PACIAE Parameters Setting   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------   PYTHIA Parameters Setting   ------------------------
!       Does not print the title page of PYTHIA 6 when running other modes.
        if( i_mode == 1 ) call PYGIVE( "MSTU(12)=12345" )
!       SetS basic parameters for PYTHIA 6/8. (In Pythia_fort_interface.f90)
        call PATUNE
!       Intrinsic Settings:
!       Sets output file of PYTHIA 6 (original unit = 6).
        MSTU(11) = 6
        call PYGIVE( "MSTU(11)=22" )
!       Extends the number of lines available in the common block PYJETS.
        ! MSTU(4) = 4000 -> KSZJ
        WRITE( PACIAE_INPUT, "(I10)" ) KSZJ
        PARAM_PYTHIA6 = "MSTU(4)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
        ! MSTU(4) = 10000 -> KSZJ
        WRITE( PACIAE_INPUT, "(I10)" ) KSZJ
        PARAM_PYTHIA6 = "MSTU(5)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
!       Prints all lines (including KS <= 0) when calling PYLIST / PALIST.
        CALL PYGIVE( "MSTU(15)=2" )

        if( i_mode == 1 ) return

!       Checks on possible errors during execution but continues the generation.
        CALL PYGIVE( "MSTU(21)=1" )
!       parp2: lowest CM energy for calling 'PYTHIA' (D=10.)
        parp2 = 0D0
!       For hA/Ah/AB.
        if( IS_NUCLEUS(KF_proj) .OR. IS_NUCLEUS(KF_targ) )then
!       For PYTHIA 6.4.28, the following CM energies for a NN pair are required.
            if( .NOT.IS_PYTHIA8(i_mode) )then
                select case( nchan )
                case( 0, 1, 6, 8, 10 )
                    parp2 = 2.7D0
                case( 2 )
                    parp2 = 0D0
                case( 7, 9 )
                    parp2 = 6.8D0
                case( 3 )
                    parp2 = 9.5D0
                case( 4 )
                    parp2 = 15.6D0
                case( 5 )
                    parp2 = 6.8D0
                case default
                    parp2 = 3D0
                end select
!       For PYTHIA 8.3, the following CM energies for NN pair are required.
            else if( IS_PYTHIA8(i_mode) )then
                select case( nchan )
                case( 0, 1, 6, 10 )
                    parp2 = 3.9D0
                case( 8 )
                    parp2 = 6.5D0
                case( 2 )
                    parp2 = 0D0
                case( 7, 9 )
                    parp2 = 12.8D0
                case( 3 )
                    parp2 = 7.2D0
                case( 4 )
                    parp2 = 7.3D0
                case( 5 )
                    parp2 = 6.9D0
                case default
                    parp2 = 3D0
                end select
            end if
        end if
!       Lowest c.m. energy that PYTHIA could work.
        WRITE( PACIAE_INPUT, "(F11.5)" ) parp2
        PARAM_PYTHIA6 = "PARP(2)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
!       Whether or not to perform fragmentation in "partonic initialization".
        WRITE( PACIAE_INPUT, "(I5)" ) mstptj
        PARAM_PYTHIA6 = "MSTP(111)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )

!------------------------   String Tension Backing-up   ------------------------
!       Backing-up basic string tension parameters.
!#TODO(Lei20240930): Let's extend it to PYTHIA 8.
        parj1  = PARJ(1)
        parj2  = PARJ(2)
        parj3  = PARJ(3)
        parj4  = PARJ(4)
        parj21 = PARJ(21)
!       RecalculateS PARJ(1) with popcorn mechanism correction.
!       No changes were actually made to string tension and PARJ(1).
!       Just to prepare for the string fragmentation when calling
!        "sfm" with kjp22 = 1, 2 and 3, i.e. gets "parj1", "parj10", "ww0", etc.
        parj10   = parj1
        wx0      = parj3
        wy0      = parj4
        wrho0    = parj2
        wnumer   = 1D0 + 2D0*wx0*wrho0 + 9D0*wy0 &
                    + 6D0 * wx0 * wy0 * wrho0       &
                    + 3D0 * wx0*wx0 *wy0 * wrho0*wrho0
        wdenom   = 2D0 + wrho0
        ww0      = wnumer / wdenom
        wweff    = ww0
        ! I.e. PARJ(1) = parj1.
        PARJ(1)  = seco * wweff * ( parj10/seco/ww0 )
        akapa(1) = PARJ(1)
        akapa(2) = PARJ(2)
        akapa(3) = PARJ(3)
        akapa(4) = PARJ(4)
        akapa(5) = PARJ(21)
        ! Default effective string tension.
        ! akapa(6) = 1D0
        akapa(6) = adj1(15)
!------------------------   String Tension Backing-up   ------------------------

!------------------------   PYTHIA Parameters Setting   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialize_global_event_variable
!!      Initializes the global variables for the entire event generation.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa6_sum/sthroq,sthrob,sthroc,sthroe(4)
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
         pph,vnlep
        common/schuds/schun,schudn,schudsn,sfra
!       For PACIAE statistical analysis. ( f: full phase space; else partial )
        common/anly1_sum1/ san(100,10,100),  sbn(100,2), &
                          sanf(100,10,100), sbnf(100,2), &
                           sao(100,10,100),  sbo(100,2), &
                          saof(100,10,100), sbof(100,2)
        common/anly1_sum2/ dnmin, dnminf, dncha, dnchaf
        common/anly4/ distr(100,10,5), dMult(5,2), distr_f(100,10,5), &
         dMult_f(5,2), san_distr(100,10,5), sbn_mult(5,2), &
         san_distr_f(100,10,5), sbn_mult_f(5,2), sao_distr(100,10,5), &
         sbo_mult(5,2), sao_distr_f(100,10,5), sbo_mult_f(5,2)
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
        common/string_kapa2/ snnc,sadiv, sgtime,sitime,sgpmax, skapa(6)
!       For the independent Optical Clauber calculation.
        common/Glauber1/ kjp23, kjp24
        common/Glauber2/ bpp(20),acoll(20),acollp(20),acollt(20),sbp(20)
        common/Glauber3/ avb, avneu, astbp, astbt, aanbin
        common/Glauber4/ averb, psnon, psnop, psnot, anbin
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the parini related statistics.
        common/anly_parini1/ woun, swoun, snpctl0, snpctlm, snpar
!       For the parcas and hadcas related statistics.
        common/anly_parcas1/ sreac(20), rinel
        common/anly_hadcas1/ sel, sinel, dinel(600), einel(600)
!       For the Angantyr related statistics.
        common/anly_angantyr1/ weight_event, sum_weight_event
        common/anly_angantyr2/ sum_bParam_ANG(3), sum_sigma_ANG(8,2), &
                               sum_Ncoll_ANG(8), sum_Npart_ANG(4,2)
!       gamma flavor code 22: hardonic decay photon
!                         44: prompt direct photon ( <- PYTHIA )
!                         55: photon from parton-parton scattering
!                             qg -> q(gamma) and qqbar -> g(gamma)
!                         66: hardonic direct photon
!                             pi+pi -> rho+(gamma) and pi+rho -> pi+(gamma)


!       fraction for charge separation
        schun=0.
        schudn=0.
        schudsn=0.
        sfra=0.

!       initiation of event averaged variales
        dnmin=0.
        dnminf=0.
        dncha=0.
        dnchaf=0.
        sbn=0.
        sbnf=0.
        san=0.
        sanf=0.
        sbn_mult=0.
        sbn_mult_f=0.
        san_distr=0.
        san_distr_f=0.

        stime_ini=0.
        stime_par=0.
        stime_had=0.

        snnc=0.
        sgtime=0.
        sitime=0.
        sadiv=0.
        sgpmax=0.
        skapa=0.

        segam=0.
        segam1=0.
        segam2=0.
        segam3=0.

        spathn=0.
        sevbin=0.

        rinel=0.
        sinel=0.
        sel=0.

        swoun=0.
        snpctl0=0.
        snpctlm=0.

        snpar=0.
        dinel=0.
        einel=0.
        acoll=0.
        acollp=0.
        acollt=0.
        sbp=0.

        skpar=0.
        sknn=0.
        skpp=0.
        sknp=0.
        skep=0.

        sthroq=0.
        sthrob=0.
        sthroc=0.
        sthroe=0.
        nrel=0
        nreac=0
        sreac=0.

        averb=0.
!       N_bin in case of psno=2 or 6
        psnon=0.
!       parojectile N_part in case of psno=2 or 6
        psnop=0.
!       target N_part in case of psno=2 or 6
        psnot=0.
!       statistics of the number of studied leptons
        vnlep=0.d0

!       For Angantyr mode.
        weight_event     = 1D0
        sum_weight_event = 0D0
        sum_bParam_ANG   = 0D0
        sum_sigma_ANG    = 0D0
        sum_Ncoll_ANG    = 0D0
        sum_Npart_ANG    = 0D0

!       Initialies them avoiding bad display in rms.out.
        avb    = 0D0
        avneu  = 0D0
        astbp  = 0D0
        astbt  = 0D0
        aanbin = 0D0
        woun   = 2D0
        pirr   = 1D0
        tirr   = 1D0
        anbin  = 0D0

!       Initializes the global variables counting x-sections.
        ! In Pythia8_fort_interface.f90.
        call PASTAT(-2,0)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialize_system
!!      Initializes the simulation system.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!       In parini.f90
        call sysini(win)
        adj1(28) = para(10) * dmax1(rnt,rnp) * adj1(28)


!-------------------------------------------------------------------------------
!------------------------   PYTHIA 8 Initialization   --------------------------
!       Forward instantiation and initialization for PYTHIA 8.
        if( IS_PYTHIA8(i_mode) )then
            ! In Pythia8_fort_interface.f90.
            call PAINST
            i_frameType = ifram
            idA = KF_proj
            idB = KF_targ
            if( i_mode /= 8 .AND. i_mode /= 9 &
                .AND. ( IS_NUCLEUS(KF_proj) .OR. IS_NUCLEUS(KF_targ) ) )then
                    idA = 2212
                    idB = 2212
            end if
            eCM = win
            N = 2
            do i=1,5,1
                P(1,i) = 0D0
                P(2,i) = 0D0
            end do
            ! Fixed target (FXT) with the incident momentum +pz.
            if( ifram == 0 )then
                i_frameType = 3
                P(1,3) = win
                P(2,3) = 0D0
            ! CMS.
            else if( ifram == 1 )then
                eCM = win
            ! Back-to-back collisions with eA +Z and eB -Z.
            else if( ifram == 2 )then
                energy_A = win
                P(1,4) = energy_A
                P(2,4) = energy_B
            ! Free 3-momentum.
            ! else if( ifram == 3 )then
            !     P(1,1) = pxA
            !     P(1,2) = pyA
            !     P(1,3) = pzA
            !     P(2,1) = pxB
            !     P(2,2) = pyB
            !     P(2,3) = pzB
            end if
!       Note this "KF" version of "PAINIT".
            ! In Pythia8_fort_interface.f90.
            call PAINIT_KF( i_frameType, idA, idB, eCM )
        end if
!------------------------   PYTHIA 8 Initialization   --------------------------
!-------------------------------------------------------------------------------


!       Loads the hadron table for the coalecnece model.
        call tabhb


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine set_stable_particle
!!      Sets particle stable or not (decay).
!!      Note: user can use "pythia6_extra.cfg" and "pythia8_extra.cfg"
!!            to specify more stable/unstable particles
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        CHARACTER*100 PARAM_PYTHIA6, PACIAE_INPUT


!-------------------------------------------------------------------------------
!------------------------   Particle Decay Specifying   ------------------------
!       Forbidden decay of particle (sets MDCY(...)=0)
!       Primary long-lived particle definition in ALICE-PUBLIC-2017-005.
        CALL PYGIVE( "MDCY( C211,   1 ) = 0" )   ! pi+
        CALL PYGIVE( "MDCY( C-211,  1 ) = 0" )   ! pi-
        CALL PYGIVE( "MDCY( C130,   1 ) = 0" )   ! K0_L
        CALL PYGIVE( "MDCY( C310,   1 ) = 0" )   ! K0_S
        CALL PYGIVE( "MDCY( C311,   1 ) = 0" )   ! K0
        CALL PYGIVE( "MDCY( C-311,  1 ) = 0" )   ! Kbar0
        CALL PYGIVE( "MDCY( C321,   1 ) = 0" )   ! K+
        CALL PYGIVE( "MDCY( C-321,  1 ) = 0" )   ! K-
        CALL PYGIVE( "MDCY( C2212,  1 ) = 0" )   ! p+
        CALL PYGIVE( "MDCY( C-2212, 1 ) = 0" )   ! pbar-
        CALL PYGIVE( "MDCY( C2112,  1 ) = 0" )   ! n0
        CALL PYGIVE( "MDCY( C-2112, 1 ) = 0" )   ! nbar0
        CALL PYGIVE( "MDCY( C3122,  1 ) = 0" )   ! Lambda0
        CALL PYGIVE( "MDCY( C-3122, 1 ) = 0" )   ! Lambdabar0
        CALL PYGIVE( "MDCY( C3112,  1 ) = 0" )   ! Sigma-
        CALL PYGIVE( "MDCY( C-3112, 1 ) = 0" )   ! Sigmabar+
        CALL PYGIVE( "MDCY( C3222,  1 ) = 0" )   ! Sigma+
        CALL PYGIVE( "MDCY( C-3222, 1 ) = 0" )   ! Sigmabar-
        CALL PYGIVE( "MDCY( C3312,  1 ) = 0" )   ! Xi-
        CALL PYGIVE( "MDCY( C-3312, 1 ) = 0" )   ! Xibar+
        CALL PYGIVE( "MDCY( C3322,  1 ) = 0" )   ! Xi0
        CALL PYGIVE( "MDCY( C-3322, 1 ) = 0" )   ! Xibar0
        CALL PYGIVE( "MDCY( C3334,  1 ) = 0" )   ! Omega-
        CALL PYGIVE( "MDCY( C-3334, 1 ) = 0" )   ! Omegabar+
        CALL PYGIVE( "MDCY( C11,    1 ) = 0" )   ! e-
        CALL PYGIVE( "MDCY( C-11,   1 ) = 0" )   ! e+
        CALL PYGIVE( "MDCY( C13,    1 ) = 0" )   ! mu-
        CALL PYGIVE( "MDCY( C-13,   1 ) = 0" )   ! mu+
        CALL PYGIVE( "MDCY( C22,    1 ) = 0" )   ! gamma
!       MDCY( KC , 1 ) = 0 means the particle will not decay (stable).
        do i1=1,100,1
            KF = KF_woDecay(i1)
            if( KF /= 0 )then
                WRITE( PACIAE_INPUT, "(I20)" ) KF
                PARAM_PYTHIA6 = "MDCY( C"// TRIM(ADJUSTL( PACIAE_INPUT )) &
                            // ",1)=0"
                CALL PYGIVE( PARAM_PYTHIA6 )
            end if
        end do
!       MDCY( KC , 1 ) = 1 means the particle will decay (unstable).
        if(i_mode == 1)then   ! For A-framework, hardcode.
            CALL PYGIVE( "MDCY( C111,  1 ) = 1" )   ! pi0
            CALL PYGIVE( "MDCY( C1114, 1 ) = 1" )   ! Delta-
            CALL PYGIVE( "MDCY( C2114, 1 ) = 1" )   ! Delta0
            CALL PYGIVE( "MDCY( C2214, 1 ) = 1" )   ! Delta+
            CALL PYGIVE( "MDCY( C2224, 1 ) = 1" )   ! Delta++
        end if
!------------------------   Particle Decay Specifying   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine set_subprocess
!!      Sets subprocesses.
!!      Note: user can use "pythia6_extra.cfg" and "pythia8_extra.cfg"
!!            to select more subprocesses.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8, IS_LEPTON
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        CHARACTER*100 PARAM_PYTHIA6, PACIAE_INPUT, CHAR_I


        if( i_mode == 1 ) return


!-------------------------------------------------------------------------------
!-------------------------   Subprocesses Selecting   --------------------------
!       For PYTHIA 6.
!       Default PYTHIA, PACIAE-defined or user-defind subprocesses.
        CALL PYGIVE( "MSEL=0" )
!       Inelastic (INEL)
        if( nchan == 0 )then
            CALL PYGIVE( "MSUB(11)=1" )   ! Hard QCD, f_i + f_j -> f_i + f_j
            CALL PYGIVE( "MSUB(12)=1" )   ! Hard QCD, f_i+f_i^bar -> f_k+f_k^bar
            CALL PYGIVE( "MSUB(13)=1" )   ! Hard QCD, f_i + f_i^bar -> g + g
            CALL PYGIVE( "MSUB(28)=1" )   ! Hard QCD, f_i + g -> f_i + g
            CALL PYGIVE( "MSUB(53)=1" )   ! Hard QCD, g + g -> f_k + f_k^bar
            CALL PYGIVE( "MSUB(68)=1" )   ! Hard QCD, g + g -> g + g
          ! CALL PYGIVE( "MSUB(91)=1" )   ! Soft QCD, elastic scattering
            CALL PYGIVE( "MSUB(92)=1" )   ! Soft QCD, single diffraction(AB->XB)
            CALL PYGIVE( "MSUB(93)=1" )   ! Soft QCD, single diffraction(AB->AX)
            CALL PYGIVE( "MSUB(94)=1" )   ! Soft QCD, double diffraction
            CALL PYGIVE( "MSUB(95)=1" )   ! Soft QCD, low_pT production
!       Non-Single Difractive (NSD)
        else if( nchan == 1 )then
            CALL PYGIVE( "MSUB(11)=1" )   ! Hard QCD, f_i + f_j -> f_i + f_j
            CALL PYGIVE( "MSUB(12)=1" )   ! Hard QCD, f_i+f_i^bar -> f_k+f_k^bar
            CALL PYGIVE( "MSUB(13)=1" )   ! Hard QCD, f_i + f_i^bar -> g + g
            CALL PYGIVE( "MSUB(28)=1" )   ! Hard QCD, f_i + g -> f_i + g
            CALL PYGIVE( "MSUB(53)=1" )   ! Hard QCD, g + g -> f_k + f_k^bar
            CALL PYGIVE( "MSUB(68)=1" )   ! Hard QCD, g + g -> g + g
          ! CALL PYGIVE( "MSUB(91)=1" )   ! Soft QCD, elastic scattering
          ! CALL PYGIVE( "MSUB(92)=1" )   ! Soft QCD, single diffraction(AB->XB)
          ! CALL PYGIVE( "MSUB(93)=1" )   ! Soft QCD, single diffraction(AB->AX)
            CALL PYGIVE( "MSUB(94)=1" )   ! Soft QCD, double diffraction
            CALL PYGIVE( "MSUB(95)=1" )   ! Soft QCD, low_pT production
!       Used to generate Drell-Yan
        else if( nchan == 2 )then
            CALL PYGIVE( "MSUB(1)=1" )   ! q + qbar -> gamma* -> l+ + l-
            CALL PYGIVE( "MSTP(43)=1" )
!       J/psi production (color singlet)
        else if( nchan == 3 )then
            CALL PYGIVE( "MSUB(86)=1" )    ! g + g -> J/psi + g
            CALL PYGIVE( "MSUB(106)=1")    ! g + g -> J/psi + gamma
!       c-cbar and b-bbar production
        else if( nchan == 4 )then
            CALL PYGIVE( "MSUB(81)=1" )    ! Open HF, f_i+f_i^bar -> Q_k+Q_k^bar
            CALL PYGIVE( "MSUB(82)=1" )    ! Open HF, g + g -> Q_k + Q_k^bar
          ! CALL PYGIVE( "MSUB(83)=1" )    ! Open HF, q_i + f_i -> Q_k + f_l
!       Prompt photon
        else if( nchan == 5 )then
            CALL PYGIVE( "MSUB(14)=1" )   !  Prompt photon, f_i+f_i^bar->g+gamma
            CALL PYGIVE( "MSUB(18)=1" )   !  Prompt photon, f_i+f_i^bar->2*gamma
            CALL PYGIVE( "MSUB(29)=1" )   !  Prompt photon, f_i +g -> f_i +gamma
            CALL PYGIVE( "MSUB(114)=1")   !  Prompt photon, g+g -> gamma +gamma
            CALL PYGIVE( "MSUB(115)=1")   !  Prompt photon, g + g -> g + gamma
!       Soft QCD
        else if( nchan == 6 )then
            CALL PYGIVE( "MSUB(91)=1" )   ! Soft QCD, elastic scattering
            CALL PYGIVE( "MSUB(92)=1" )   ! Soft QCD, single diffraction(AB->XB)
            CALL PYGIVE( "MSUB(93)=1" )   ! Soft QCD, single diffraction(AB->AX)
            CALL PYGIVE( "MSUB(94)=1" )   ! Soft QCD, double diffraction
            CALL PYGIVE( "MSUB(95)=1" )   ! Soft QCD, low_pT production
!       Single W+/- production
        else if( nchan == 7 )then
            CALL PYGIVE( "MSUB(2)=1" )   ! Single W, f_i + f_j^bar -> W
!       Hard QCD
        else if( nchan == 8 )then
            CALL PYGIVE( "MSEL=1" )
        else if( nchan == 9 )then
!       Single Z0 production
            CALL PYGIVE( "MSUB(1)=1" )   ! Single Z, f_i + f_i^bar -> Z0
            CALL PYGIVE( "MSTP(43)=2" )   ! Only Z0 included.
!       Inelastic nondiffrative (Minimum Bias)
        else if( nchan == 10 )then
            CALL PYGIVE( "MSUB(11)=1" )   ! Hard QCD, f_i + f_j -> f_i + f_j
            CALL PYGIVE( "MSUB(12)=1" )   ! Hard QCD, f_i+f_i^bar -> f_k+f_k^bar
            CALL PYGIVE( "MSUB(13)=1" )   ! Hard QCD, f_i + f_i^bar -> g + g
            CALL PYGIVE( "MSUB(28)=1" )   ! Hard QCD, f_i + g -> f_i + g
            CALL PYGIVE( "MSUB(53)=1" )   ! Hard QCD, g + g -> f_k + f_k^bar
            CALL PYGIVE( "MSUB(68)=1" )   ! Hard QCD, g + g -> g + g
          ! CALL PYGIVE( "MSUB(91)=1" )   ! Soft QCD, elastic scattering
          ! CALL PYGIVE( "MSUB(92)=1" )   ! Soft QCD, single diffraction(AB->XB)
          ! CALL PYGIVE( "MSUB(93)=1" )   ! Soft QCD, single diffraction(AB->AX)
          ! CALL PYGIVE( "MSUB(94)=1" )   ! Soft QCD, double diffraction
            CALL PYGIVE( "MSUB(95)=1" )   ! Soft QCD, low_pT production
        else
!       Other: user-defined subprocessses in "pythia6_extra.cfg" for
!              PYTHIA 6 or "pythia8_extra.cfg" for PYTHIA 8.
            write(22,*)
            write(22,*) "PACIAE Warning: You have chosen an undefined nchan =",&
                        nchan, "."
            write(22,*) "Please specify subprocesses in " &
                     // """pythia6_extra.cfg"" or ""pythia8_extra.cfg""."
            write(22,*)
        end if

!       For l+l- annihilation, e.g. e+e-. l+l- -> gamma*/Z -> qqbar.
        if( IS_LEPTON( KF_proj ) .AND. IS_LEPTON( KF_targ ) &
            .AND. KF_proj*KF_targ < 0 .AND. nchan <= 10 )then
            CALL PYGIVE( "MSEL=0" )
!           Hard process, with hadronic Z decays.
            CALL PYGIVE( "MSUB(1)=1" )
!           Only allow Z0 decay to quarks (i.e. no leptonic final states).
            do i = MDCY(23,2), MDCY(23,2) + MDCY(23,3) - 1, 1
                if( IABS(KFDP(i,1)) >= 6 )then
                    MDME(i,1) = MIN( 0, MDME(i,1) )
                    ME=MIN( 0, MDME(i,1) )
                    WRITE( CHAR_I, "(I20)" ) i
                    WRITE( PACIAE_INPUT, "(I20)" ) ME
                    PARAM_PYTHIA6 = "MDME(" // TRIM(ADJUSTL( CHAR_I )) &
                                 // ",1)="  // TRIM(ADJUSTL( PACIAE_INPUT ))
                    CALL PYGIVE( PARAM_PYTHIA6 )
                end if
            end do
!       For lh(hl) and lA(Al) deeply inelastic scatterings (DIS).
        else if( IS_LEPTON( KF_proj ) .AND. IS_LEPTON( KF_targ ) &
                .AND. KF_proj*KF_targ < 0 .AND. nchan <= 10 )then
            
        end if

!       For PYTHIA 8, they will be done in Pythia8CppInterface.cpp.
!-------------------------   Subprocesses Selecting   --------------------------
!-------------------------------------------------------------------------------


!       nPDF warning.
        if( IS_PYTHIA8(i_mode) .AND. ( adj1(5) < -1 .OR. adj1(5) > 3 ) )then
            write(22,*)
            write(22,*) "PACIAE Warning: You are using the preset 208Pb " &
                    //  "EPS09 LO nPDF."
            write(22,*) "To set other nPDF, modify ""adj1(5)"" in " &
                    //  """usu.dat"" of PACIAE,"
            write(22,*) " and place the appropriate nPDF grid file into"
            write(22,*) " your-PYTHIA8-directory/share/Pythia8/pdfdata/"
            write(22,*)
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine calculate_optical_glauber
!!      Independent Optical Glauber calculations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa30/vneump,vneumt,mstptj
        common/sa31/rmax,bbb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the independent Optical Clauber calculation.
        common/Glauber1/ kjp23, kjp24
        common/Glauber2/ bpp(20),acoll(20),acollp(20),acollt(20),sbp(20)
        common/Glauber3/ avb, avneu, astbp, astbt, aanbin
        common/Glauber4/ averb, psnon, psnop, psnot, anbin
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max

!       acoll: a array, demension of which should be larger
!        than or equal to 'nmax' (interval # of impact parameter)
!       the dimension of 'bpp' must be small than or equal to 'nmax'


!       # of segments used in integration
        idw = INT( adj1(4) )


!-------------------------------------------------------------------------------
!------------   Independent Optical Initial Geometry Calculation   -------------
!       calculates the nuclear overlap function etc.
!       if((ipden == 0 .and. itden == 0) .or. ipden >= 11) goto 80001
        if( (ipden == 0 .and. itden == 0) .or. (ipden == 2 .and. &
         itden == 2) .or. ipden >= 11 ) goto 80001
        csnn1=csnn*10   ! csnn in fm^2 csnn1 in mb
        idw1=idw/50   ! *100
        if( ipden < 2 ) call xOVERLAP(nap,nat,rnp,rnt,csnn1,kjp23,kjp24, &
         rou0,idw1)

80001   continue

!       calculate the impact parameter etc.
!       for given b (impact parameter)

!---------------------------------   Fixed b   ---------------------------------
        if( INT(psno) == 0 .OR. INT(psno) == 4 )then
            bp=bmin
            r4=rnp
            if(rnt > rnp) r4=rnt
            rr4=bp/r4
            vneu=dexp(-rr4*rr4)
!       calculates the overlap region of two nuclei at given b by
!        interpolation method
            if( (ipden == 0 .and. itden == 0) .or. (ipden == 2 .and. &
             itden == 2) .or. ipden >= 11 ) goto 80002
            ibpp=int(bp/0.1+1.0)
            ibpp=min(ibpp,200)
!           if(ipden == 1 .and. itden == 1)then   ! A+B
!               overlap function of A+B (1/fm^2)
                anbin=ta1a2(ibpp)
!           elseif(ipden == 0 .and. itden == 1)then   ! p+A
!               overlap function of B nucleus (1/fm^2)
!               anbin=ta2(ibpp)
!           elseif(ipden == 1 .and. itden == 0)then   ! A+p
!               overlap function of A nucleus (1/fm2)
!               anbin=ta1(ibpp)
!           else
!           endif
            pir=part1(ibpp)
            tir=part2(ibpp)
!           Nbin in current event
            evbin=anbin*csnn
            pirr=pir
            tirr=tir
!       endif

80002       if( (ipden == 0 .and. itden == 0) .or. (ipden == 2 .and. &
             itden == 2) .or. ipden >= 11 )then
                pir=1.
                tir=1.
                evbin=1.
                pirr=1.
                tirr=1.
            endif
            vneump=pir
            vneumt=tir
            write(19,*) "#! Fixed b = bmin"
            write(19,*) "#! psno, b, N_part_p, N_part_t, N_bin " // &
                        "(optical Glauber)=", psno,bp,vneump,vneumt,evbin
            return
!---------------------------------   Fixed b   ---------------------------------

!--------------------------    Systematic Sampling   ---------------------------
        else if( INT(psno) == 1 .OR. INT(psno) == 5 )then
!       systematic sampling method for given interval of b according to b**2 law
            nmax2=2*nmax
            bmaxn=(bmax*bmax-bmin*bmin)/(2D0*nmax)
            bmin2=bmin*bmin
            i2=0
            do i1=1,nmax2,2
                i2=i2+1
                bpp(i2)=dsqrt(i1*bmaxn+bmin2)
            enddo
            write(19,*) "#! Systematic sampling of b from (bmin,bmax)"
            write(19,*) "nmax=", nmax
            write(19,*) "b=", (bpp(i1),i1=1,i2)

            stab=0.
            stb=0.
            stbp=0.
            stbt=0.

            do i1=1,i2
                bp=bpp(i1)
                r4=rnp
                if(rnt > rnp) r4=rnt
                rr4=bp/r4
                acoll(i1)=dexp(-rr4*rr4)
!       calculate the overlap region of two nuclei at given b
                ibpp=int(bp/0.1+1.0)
                ibpp=min(ibpp,200)
!               if(ipden == 1 .and. itden == 1)then   ! A+B
!                   overlap function of A+B (1/fm^2)
                    anbin=ta1a2(ibpp)
!               elseif(ipden == 0 .and. itden == 1)then   ! p+A
!                   overlap function of B nucleus (1/fm^2)
!                   anbin=ta2(ibpp)
!               elseif(ipden == 1 .and. itden == 0)then   ! A+p
!                   overlap function of A nucleus (1/fm^2) 020718
!                   anbin=ta1(ibpp)
!               else
!               endif
                pir=part1(ibpp)
                tir=part2(ibpp)
                stab=stab+bp
                acoll(i1)=anbin
                acollp(i1)=pir
                acollt(i1)=tir
                stbp=stbp+acollp(i1)
                stbt=stbt+acollt(i1)
                stb=stb+acoll(i1)
            enddo
            stab=stab/float(i2)
            aneump=stbp/dfloat(i2)
            aneumt=stbt/dfloat(i2)
            vneum=stb/dfloat(i2)
            write(19,*) "psno, ave. b=", psno,stab
            write(19,*) "N_bin=", (acoll(i1)*csnn,i1=1,i2)
            write(19,*) "(N_part)_p=", (acollp(i1),i1=1,i2)
            write(19,*) "(N_part)_t=", (acollt(i1),i1=1,i2)
            write(19,*) "ave. N_part_p, N_part_t, N_bin " // &
                        "(optical Glauber)=", &
             aneump,aneumt,vneum*csnn

!       average b in [bmin,bmax]
            avb=2./3.*(bmin+bmax)
!       above equation is correct when bmin=0 only
            r4=rnp
            if(rnt > rnp) r4=rnt
            rr4=avb/r4
            avneu=dexp(-rr4*rr4)
!       calculate the overlap region of two nuclei at given b
            ibpp=int(avb/0.1+1.0)
            ibpp=min(ibpp,200)
!           if(ipden == 1 .and. itden == 1)then   ! A+B
!               overlap function of A+B (1/fm^2)
                anbin=ta1a2(ibpp)
!           elseif(ipden == 0 .and. itden == 1)then   ! p+A
!               overlap function of B nucleus (1/fm^2)
!               anbin=ta2(ibpp)
!           elseif(ipden == 1 .and. itden == 0)then   ! A+p
!               overlap function of A nucleus (1/fm2)
!               anbin=ta1(ibpp)
!           else
!           endif
            pir=part1(ibpp)
            tir=part2(ibpp)
            aanbin=anbin
            astbp=pir
            astbt=tir
!--------------------------    Systematic Sampling   ---------------------------

!-----------------------   Angantyr Gaussian Sampling   ------------------------
        else if( INT(psno) == 3 )then
            write(19,*) "#! Angantyr weighted events via Gaussian sampling " &
                     // "of b from (bmin,bmax)"
!-----------------------   Angantyr Gaussian Sampling   ------------------------

!-----------------------------   Random Sampling   -----------------------------
        else
            write(19,*) "#! Random sampling of b from (bmin,bmax)"
        end if
!-----------------------------   Random Sampling   -----------------------------

!------------   Independent Optical Initial Geometry Calculation   -------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine record_simulation_parameter
!!      Records basic parameters of the simulation. Some of them might
!!       have been re-evaluated.
!!      Parameters are categorized to "PACIAE-" and "PYTHIA-" related.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDATR/MRPY(6),RRPY(100)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
        disbe(100,100)
        common/sa7/ispmax,isdmax,iflmax,ispkf(100),n_bin_hist,asd(10), &
         afl(100,10,2),i_y_or_eta
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa29/i_color_reconnection,NCR,parp2
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sa38/ prob_ratio_q(6), am(6), amqq(6)
        common/count/isinel(600)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/coal1/bmrat,i_mm
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For the independent Optical Clauber calculation.
        common/Glauber1/ kjp23, kjp24
!       For cross sections of hh collisions etc.
        common/para_h1/ para(20)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!-------------------------------------------------------------------------------
!---------------------------   Parameter Recording   ---------------------------
        write(19,*)
        write(19,*)
        write(19,*) "#!-------------------------- Input Recording " // &
                    "--------------------------"
!       Variables for the event generation.
        write(19,*) "#! neve, nout, nosc:"
        write(19,*)     neve, nout, nosc
        write(19,*) "#! KF_proj, KF_targ:"
        write(19,*)     KF_proj, KF_targ
        write(19,*) "#! nap, nzp, nat, nzt, ipden, itden:"
        write(19,*)     nap, nzp, nat, nzt, ipden, itden
        write(19,*) "#! i_mode, ifram, nchan:"
        write(19,*)     i_mode, ifram, nchan
        write(19,*) "#! win, energy_B:"
        write(19,*)     win, energy_B
        write(19,*) "#! bmin, bmax, psno, nmax:"
        write(19,*)     bmin, bmax, psno, nmax
        write(19,*) "#! kjp21, x_ratio, decpro:"
        write(19,*)     kjp21, x_ratio, decpro
!       Variables for the event analysis.
        write(19,*) "#! ispmax, isdmax, iflmax, n_bin_hist, i_y_or_eta:"
        write(19,*)     ispmax, isdmax, iflmax, n_bin_hist, i_y_or_eta
        write(19,*) "#! KF_woDecay:"
        write(19,*)   ( KF_woDecay(i), i=1,10,1 )
        write(19,*)   ( KF_woDecay(i), i=11,ispmax,1 )
        write(19,*) "#! ispkf:"
        write(19,*)   ( ispkf(i), i=1,10,1 )
        write(19,*)   ( ispkf(i), i=11,ispmax,1 )
        write(19,*) "#! asd:"
        write(19,*)   ( asd(i), i=1,isdmax,1 )
        if( iflmax > 0 )then
            write(19,*) "#! afl:"
            do kk=1,ispmax
                do i=1,iflmax
                    write(19,*) ( afl(kk,i,j), j=1,2,1 )
                end do
            end do
        end if
!       Variables for the event generation.
        write(19,*) "#! ddt, i_mm:"
        write(19,*)     ddt, i_mm
        write(19,*) "#! iparres, i_inel_proc, i_time_shower:"
        write(19,*)     iparres, i_inel_proc, i_time_shower
        write(19,*) "#! para(7), ttaup, taujp, para(10):"
        write(19,*)     para(7), ttaup, taujp, para(10)
        write(19,*) "#! i_sigma_AQM, kjp20:"
        write(19,*)     i_sigma_AQM, kjp20
        write(19,*) "#! para1_1, para1_2, para(2), para(4), para(5):"
        write(19,*)     para1_1, para1_2, ( para(i), i=2,5,1 )
        write(19,*) "#! para(13), para(14), para(15), para(16):"
        write(19,*)     ( para(i), i=13,16,1 )
        write(19,*) "#! adj1(i):"
        write(19,*)   ( adj1(i), i=1,10,1  )
        write(19,*)   ( adj1(i), i=11,20,1 )
        write(19,*)   ( adj1(i), i=21,30,1 )
        write(19,*)   ( adj1(i), i=31,40,1 )
        write(19,*) "#! kjp22-24, i_color_reconnection, i_tune:"
        write(19,*)     kjp22,kjp23,kjp24,i_color_reconnection, i_tune
        write(19,*) "#! parecc, smadel, cp0, cr0, seco:"
        write(19,*)     parecc, smadel, cp0, cr0, seco
        write(19,*) "#! i_deex, n_deex_step, i_pT_coal, i_pT_endpoint:"
        write(19,*)     i_deex, n_deex_step, i_pT_coal, i_pT_endpoint
        write(19,*) "#! a_FF, aPS_c, aPS_b, bmrat:"
        write(19,*)     a_FF, aPS_c, aPS_b, bmrat
        write(19,*) "#! prob_ratio_q(i):"
        write(19,*)   ( prob_ratio_q(i), i=1,6,1 )
        write(19,*) "#!-------------------------- Input Recording " // &
                   "--------------------------"
        write(19,*)
        write(19,*)
!       Records the seed of random number generator.
        write(19,*) "#! Seed (PACIAE default=20240116) =", MRPY(1)
        write(19,*)
        write(19,*)
        write(19,*) "cspipiKK, t0, ddt, dep =", cspipiKK, t0, ddt, dep
        write(19,*) "rou0, rao, rnp, rnt =", rou0, rao, rnp, rnt
        write(19,*) "csnn, cspin, cskn =", csnn, cspin, cskn
        write(19,*) "cspipi, cspsn, cspsm =", cspipi, cspsn, cspsm
        write(19,*) "ifram, rcsit, kfmax, ipden, itden =", &
                    ifram, rcsit, kfmax, ipden, itden
        n_kfmax = kfmax / 10
        write(19,*) "#! Particle KF considered in the hadronic" // &
                    " rescattering, kfaco ="
        do i_kfmax = 1, n_kfmax+1, 1
            i_low = i_kfmax*10 - 9
            i_upp = i_kfmax*10
            if( i_kfmax == (n_kfmax+1) ) i_upp = kfmax
            write(19,*) ( kfaco(i), i = i_low, i_upp, 1 )
        end do
        write(19,*) "#! Allowable minimum distance (not used now)," // &
                    " disbe ="
        do i_kfmax = 1, n_kfmax+1, 1
            i_low = i_kfmax*10 - 9
            i_upp = i_kfmax*10
            if( i_kfmax == (n_kfmax+1) ) i_upp = kfmax
            write(19,*) ( disbe(i,i), i=i_low,i_upp,1 )
        end do
        write(19,*) "#! isinel (open # inel. channel or not," // &
                    " list ref. subroutine coinel in hadcas.f90) ="
        write(19,600) isinel
600     format( 25(1x,i2)/ )
        write(19,*)
!---------------------------   Parameter Recording   ---------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initialize_single_event_variable
!!      Initializes the variables for one single event generation.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/delt/ndel,nodel,kdel(kszj,5),pdel(kszj,5),vdel(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/ctllist/npctl,npinel(600),npctl0,npctlm
        common/ctllist_h/nctl,noinel(600),nctl0,noel
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For the independent Optical Clauber calculation.
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the charge and 4-momentum check.
        common/cp_check/ p_init(4), p_final(4)
!       Arrays of particle information of one NN pair for PYTHIA 8.
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
!       Arrays of particle information of a AA collision (total NN pairs).
        COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
               K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
!       Arrays of freeze-out particle information.
        COMMON/PYJETS_FO/ N_FO, NPAD_FO, &
                          K_FO(KSZJ,8), P_FO(KSZJ,7), V_FO(KSZJ,5)


!-------------------------------------------------------------------------------
!-------------------   Single Event Variable Initializing   --------------------
!       Particle data.
        N    = 0
        do j=1,5,1
            do i=1,2,1
                K(i,j) = 0
                P(i,j) = 0D0
                V(i,j) = 0D0
            end do
        end do
        nsa  = 0
        nth  = 0
        nbe  = 0
        naf  = 0
        naff = 0
        nbh  = 0
        nn   = 0
        ndel = 0
        ngam = 0
        ntrs = 0
        N_PY8        = 0
        N_AA         = 0
        N_TOT_AA     = 0
        N_COLL_NN_AA = 0
        N_FO = 0
!       Diquark data.
        idi   = 0
        idio  = 0
        ndiq  = 0
        npt   = 0
        ifcom = 0
!       String list.
        nstr = 0
        nstr0 = 0
        nstr1 = 0
!       For parini.
        npctlm = 0
!       For parcas.
        ishp = 0
        tau  = 0D0
!       For sfm.
!       For coales.
!       For hadcas.
!       For the evolution time.
        time_ini = 0D0
        time_par = 0D0
        time_had = 0D0

        pathn   = 0D0
!       For number of elastic and inelasitc collisions in parcas.
        npinel  = 0
!       For number of elastic and inelasitc collisions in hadcas.
        noel    = 0
        noinel  = 0
!       For lost 4-momentum and charge.
        ithroq   = 0
        ithrob   = 0
        ithroc   = 0
        non6_p   = 0
        throe    = 0D0
        ithroq_p = 0
        ithrob_p = 0
        ich_p    = 0
        non6_c   = 123456
        throe_p  = 0D0
        ppsa     = 0D0
        nnc      = 0D0
        p_init   = 0D0
        p_final  = 0D0

!       Initializes the single-event variables couting cross sections.
        call PASTAT(-1,0)
!-------------------   Single Event Variable Initializing   --------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_impact_parameter
!!      Gets an impact parameter for a nucleus-induced collision.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_NUCLEUS
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa16/x_ratio,decpro,dni(10),dpi(10),edi(10),bmin,bmax, &
         bar(10),abar(10),barf(10),abarf(10),emin(10),eminf(10), &
         eplu(10),epluf(10),nmax
        common/sa30/vneump,vneumt,mstptj
        common/sa31/rmax,bbb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)
!       For the independent Optical Clauber calculation.
        common/Glauber2/ bpp(20),acoll(20),acollp(20),acollt(20),sbp(20)
        common/Glauber4/ averb, psnon, psnop, psnot, anbin
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        data iii_0 /0/
        data bp_0 /0D0/
        data anbin_0 /0D0/
        data vneump_0 /0D0/
        data vneumt_0 /0D0/
        data jjj / 0 /
!       bp: current impact parameter


!-------------------------------------------------------------------------------
!------------------------   Impact Parameter Sampling   ------------------------
!       Fixed b.
        if( INT(psno) == 0 .OR. INT(psno) == 4 )then
            bp = bmin

!       Systematic sampling method.
        else if( INT(psno) == 1 .OR. INT(psno) == 5 )then
            jjj = jjj + 1
            ! For the error event.
            if( iii == iii_0 ) jjj = jjj - 1
            sbp(jjj) = sbp(jjj) + 1
            bp       = bpp(jjj)
            vneump   = acollp(jjj)
            vneumt   = acollt(jjj)
            evbin    = acoll(jjj)*csnn
            pirr     = vneump
            tirr     = vneumt
!           10: total number of impact paremeters.
            if( jjj == 10 ) jjj = 0

!       Angantyr Gaussian sampling method.
        else if( INT(psno) == 3 )then

!       Random sampling method.
        else
            bmin2 = bmin*bmin   ! 190620
            bp    = sqrt( pyr(1)*( bmax*bmax - bmin2 ) + bmin2 )
!       calculate the overlap region of two nuclei at given bp
            ibpp  = int( bp/0.1 + 1.0 )
            ibpp  = min( ibpp, 200 )
!           if(ipden == 1 .and. itden == 1)then   ! A+B
!               overlap function of A+B (1/fm^2)
                anbin = ta1a2(ibpp)
!           elseif(ipden == 0 .and. itden == 1)then   ! p+A
!               overlap function of B nucleus (1/fm^2)
!               anbin=ta2(ibpp)
!           elseif(ipden == 1 .and. itden == 0)then   ! A+p
!               overlap function of A nucleus (1/fm^2)
!               anbin=ta1(ibpp)
!           else
!           endif
            pir    = part1(ibpp)
            tir    = part2(ibpp)
            evbin  = anbin*csnn
            pirr   = pir
            tirr   = tir
            vneump = pir
            vneumt = tir
            if(iii == 1)then
                write(19,*) "#! Random sampling of b from (bmin,bmax)"
                write(19,*) "iii, psno, b, N_part_p,N_part_t, N_bin=" // &
                            "(optical Glauber)", iii,psno,bp,vneump,vneumt,evbin
            end if
            averb = averb + bp
            psnon = psnon + anbin
            psnop = psnop + vneump
            psnot = psnot + vneumt
            ! For the error event.
            if( iii == iii_0 )then
                averb = averb - bp_0
                psnon = psnon - anbin_0
                psnop = psnop - vneump_0
                psnot = psnot - vneumt_0
            end if
            bp_0 = bp
            anbin_0 = anbin
            vneump_0 = vneump
            vneumt_0 = vneumt
        end if
        iii_0 = iii

!       Corrects output for i_mode=4. pathn: N_bin / N_part.
        pathn = evbin / ( pirr + tirr )
        if( .NOT.IS_NUCLEUS(KF_proj) .AND. .NOT.IS_NUCLEUS(KF_targ) ) &
            pathn = 1D0
!------------------------   Impact Parameter Sampling   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_initial_state_generation( i_error, i_pure_hadron )
!!      Generates the initial partonic/hadronic state.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa30/vneump,vneumt,mstptj
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        i_error = 0
        i_pure_hadron = 0
        i_stop_point  = INT( adj1(40) )


!-------------------------------------------------------------------------------
!------------------   hh, ll, lh(hl) & Angantyr Generating   -------------------
        if( ( .NOT.IS_NUCLEUS(KF_proj) .AND. .NOT.IS_NUCLEUS(KF_targ) ) &
            .OR. i_mode == 8 .OR. i_mode == 9 )then
            if( i_mode == 1 )then
                ! In parini.f90.
                call simulate_NN_in_framework_A
            else
                ! In parini.f90.
                call simulate_lh_and_angantyr_in_framework_BC
            end if
!------------------   hh, ll, lh(hl) & Angantyr Generating   -------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------   hA(Ah), lA(Al) & AB Generating   ----------------------
        else if( IS_NUCLEUS(KF_proj) .OR. IS_NUCLEUS(KF_targ) )then
            ! In parini.f90.
            call parini( time_ini, ijk )
            ! Infinite loop in parini.
            if( ijk == 1 )then
                i_error = 1
                return
            end if
        end if
!---------------------   hA(Ah), lA(Al) & AB Generating   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------   Diffractive & Hadron Remnants Event ---------------------
        if(N <= 0)then
            ! There are hadrons/leptons only. Tries the hadron rescattering.
            if( nbh > 0 .AND. i_stop_point == 4 )then
                i_pure_hadron = 1
                ! Moves leptons/hadrons from "sbh" to "PYJETS".
                N = 0
                do i=1,nbh,1
                    N = N + 1
                    do j=1,5,1
                        K(N,j) = kbh(i,j)
                        P(N,j) = pbh(i,j)
                        V(N,j) = vbh(i,j)
                    end do
                end do
            ! No particles produced. Regenerates an event.
            else if( nbh <= 0 )then
                i_error = 1
                iii = iii - 1
            end if
        end if
!-------------------   Diffractive & Hadron Remnants Event ---------------------
!-------------------------------------------------------------------------------


        if( i_mode == 4 .OR. mstptj == 1 ) i_pure_hadron = 1

!       full_events_history of OSC1999A
        if( i_pure_hadron == 0 )then
            ! Initial partonic state record.
            call oscar(1)
        else if( i_pure_hadron == 1 )then
            ! Initial hadronic state record.
            call oscar(3)
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_partonic_rescattering( i_error, i_pure_hadron )
!!      Stops the event evolution at the specified point.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       N.B.
!        For special Coal adj1(12) = 2 (i.e. adj1(12)=1 + kjp22=2 ), and
!           adj1(12)=1 + kjp22=2/4/7/9 :
!         The gluon splitting and/or the energetic quark deexcitation
!          before "parcas". It means that there might be no gluons
!          enter "parcas".
!        For special SFM adj1(12) = 3:
!         With the gluon splitting and the energetic quark deexcitation
!          before parcas and hadronized by string fragmentation model.


        i_error = 0
!       The pure hadronic state does not require the partonic rescattering.
        if( i_pure_hadron == 1 ) return
!       The hadronization model.
        i_had_model = INT( adj1(12) )
        mode_coal   = kjp22

!       Special Coal/SFM: gluon splitting and/or quark deexc. before parcas.
        ! In coales_40.f90.
        if( ( i_had_model == 1 .AND. ( mode_coal == 2 .OR. mode_coal == 7 ) ) &
            .OR. i_had_model == 3 ) call do_gluon_splitting
        ! In coales_40.f90.
        if( ( i_had_model == 1 .AND. ( mode_coal == 2 .OR. mode_coal == 4 &
                                  .OR. mode_coal == 7 .OR. mode_coal == 9 ) ) &
            .OR. i_had_model == 3 ) call do_quark_deexcitation
        ! In coales_40.f90.
            ! call compensate_conservation

!       Puts particles on-shell by hand to avoid potential errors.
        do i=1,N,1
            KS = K(i,1)
            KF = K(i,2)
            if( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 &
                .OR. IS_NUCLEUS(KF) ) cycle
            px0 = P(i,1)
            py0 = P(i,2)
            pz0 = P(i,3)
            E0  = P(i,4)
            dm  = P(i,5)
            P(i,4) = SQRT( px0**2 + py0**2 + pz0**2 + dm**2 )
            throe_p(4) = throe_p(4) + E0 - P(i,4)
        end do

!       Records the freeze-out information of particles.
        i_FZ_record = 0
        call update_freeze_out_information( i_FZ_record, 1, N )

        call parcas( time_par, iijk )

!       Infinite loop in parcas.
        if(iijk == 1)then
            i_error = 1
            return
        end if

!       Loads the freeze-out information of particles.
        ! i_FZ_load = -1
        ! call update_freeze_out_information( i_FZ_load, 1, N )

!       full_events_history of OSC1999A
        call oscar(2)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_hadronization( i_error, i_pure_hadron )
!!      Stops the event evolution at the specified point.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3


        i_error = 0
!       The pure hadronic state does not require the hadronization.
        if( i_pure_hadron == 1 ) return
!       The hadronization model.
        i_had_model = INT( adj1(12) )
        time_parcas_happen = time_par


!-------------------------------------------------------------------------------
!-----------------------   String Fragmentation Model   ------------------------
        if( i_had_model == 0 .OR. i_had_model == 3 )then
            time_parcas_happen = time_par
            call do_parton_string_fragmentation( i_error, time_parcas_happen )
!-----------------------   String Fragmentation Model   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Coalescence  Model  -----------------------------
        else if( i_had_model == 1 )then
            call do_parton_coalescence
        end if
!---------------------------   Coalescence  Model  -----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------   Rest Parton Re-hadronization  ------------------------
        if( i_error == 0 )  call rest_hadronization( i_error )
!----------------------   Rest Parton Re-hadronization  ------------------------
!-------------------------------------------------------------------------------


        if( i_error == 1 ) return
!       egam3: total gamma energy after hadronization
        call prt_sgam( ngam, egam3, 3 )


!-------------------------------------------------------------------------------
!-----------------------------   Spectator Moving  -----------------------------
        if(nbh > 0)then
!       Spectators in 'sbh' moved during the time of 'parcas'.
            if( time_parcas_happen > 1D-5 )then
                t_par = time_parcas_happen
                do l=1,nbh,1
                    KF  = kbh(l,2)
                    pT2 = pbh(l,1)**2 + pbh(l,2)**2
                    if( pT2 <= 1D-15 &
                        .AND. (kf == 2212 .OR. kf == 2112) )then
                        vbh(l,1) = vbh(l,1) + t_par * pbh(l,1) / pbh(l,4)
                        vbh(l,2) = vbh(l,2) + t_par * pbh(l,2) / pbh(l,4)
                        vbh(l,3) = vbh(l,3) + t_par * pbh(l,3) / pbh(l,4)
                        vbh(l,4) = t_par
                    end if
                end do
            end if
!       "sbh" -> "PYJTES".
            do l=1,nbh,1
                l1 = N + l
                do m=1,5,1
                    K(l1,m) = kbh(l,m)
                    P(l1,m) = pbh(l,m)
                    V(l1,m) = vbh(l,m)
                enddo
            enddo
            N = N + nbh
            nbh = 0
        end if
!-----------------------------   Spectator Moving  -----------------------------
!-------------------------------------------------------------------------------


!       full_events_history of OSC1999A
        call oscar(3)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_parton_string_fragmentation( i_error, time_parcas_happen )
!!      Performs the string fragmentation hadronization.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3


        i_error = 0

        ! In sfm_40.f90.
        call reconstruct_parton_string( time_parcas_happen )

!       full_events_history of OSC1999A (includes diquarks)
        ! call oscar(2)

        mode_SFM = kjp22
        if( mode_SFM == 2 .OR. mode_SFM == 4 )then
            ! In sfm_40.f90.
            call do_fragmentation_with_constant_tension( i_error )
        else if( mode_SFM == 1 .OR. mode_SFM == 3 )then
            ! In sfm_40.f90.
            call do_fragmentation_with_fluctuating_tension( i_error )
        else
            ! In sfm_40.f90.
            call do_fragmentation_pair_by_pair( i_error )
        end if

!       Throws away this event.
        if( i_error == 1 )then
            iii = iii - 1
            return
        end if

!       Removes photons from "PYJETS" to "sgam".
        KF_gamma_decay = 22
        call remo_gam( KF_gamma_decay )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_parton_coalescence
!       Performs the coalescence hadronization.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc


        mode_coal = kjp22

!       For Special Coal mode 2/4/7/9, the gluon splitting and quark
!        deexcitation have been done before "parcas".
        if( mode_coal == 1 .OR. mode_coal == 6 ) call do_gluon_splitting
        ! In coales_40.f90.
        if(      mode_coal == 1 .OR. mode_coal == 3 &
            .OR. mode_coal == 6 .OR. mode_coal == 8 ) &
            call do_quark_deexcitation

!       Puts particles on-shell by hand to avoid potential errors.
        do i=1,N,1
            KS = K(i,1)
            KF = K(i,2)
            if( .NOT.IS_EXIST(KS,i_mode) .OR. IS_NUCLEUS(KF) ) cycle
            px0 = P(i,1)
            py0 = P(i,2)
            pz0 = P(i,3)
            E0  = P(i,4)
            dm  = P(i,5)
            P(i,4) = SQRT( px0**2 + py0**2 + pz0**2 + dm**2 )
            throe_p(4) = throe_p(4) + E0 - P(i,4)
        end do

!       Performs the parton coalescence.
        call coales( mode_coal )

        ! call compensate_conservation


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_hadronic_rescattering( i_error )
!!      Performs the PACIAE built-in hadronic rescattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5, &
         disbe(100,100)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!       Gives status code KS=1 to avoid potential misidentification for PYTHIA8.
        do i=1,N,1
            K(i,1) = 1
        end do
        i_error = 0

!       Puts particles on-shell by hand to avoid potential errors.
        do i=1,N,1
            KF  = K(i,2)
            if( IS_NUCLEUS(KF) ) cycle
            px0 = P(i,1)
            py0 = P(i,2)
            pz0 = P(i,3)
            E0  = P(i,4)
            dm  = P(i,5)
            P(i,4) = SQRT( px0**2 + py0**2 + pz0**2 + dm**2 )
            throe_p(4) = throe_p(4) + E0 - P(i,4)
        end do

!       Records the freeze-out information of particles.
        i_FZ_record = 0
        call update_freeze_out_information( i_FZ_record, 1, N )

!       Moves (excludes) spectators with a large enough time 1D28.
        time_spectator = 1D28
        do i=1,N,1
            KF = K(i,2)
            pT2 = P(i,1)**2 + P(i,2)**2
            E = SQRT( P(i,1)**2 + P(i,2)**2 + P(i,3)**2 + P(i,5)**2 )
            if( ( KF == 2212 .AND. KF == 2112 .AND. pT2 <= 1D-15 ) &
                .OR. IS_NUCLEUS(KF) )then
                V(i,1) = V(i,1) + ( time_spectator + V(i,4) ) * P(i,1) / E
                V(i,2) = V(i,2) + ( time_spectator + V(i,4) ) * P(i,2) / E
                V(i,3) = V(i,3) + ( time_spectator + V(i,4) ) * P(i,3) / E
                V(i,4) = time_spectator + V(i,4)
            end if
        end do


!-------------------------------------------------------------------------------
!--------------------------   Hadronic Rescattering  ---------------------------
!       Hadronic cascade (rescattering, HRS).
        if(kjp21 == 1)then

!----------------------------   Hadron Specifying  -----------------------------
            ! In parini.f90.
            call filt
            nup = numbs(kfmax)
            nbh = N - nup
!           nup is the number of particles kept in 'PYJETS' (will joint HRS)
!           nbh is the number of particles storing in 'sbh' (wil not joint HRS)
!           Leptons will not rescatter with hadrons.
            do i = nup+1, N, 1
                i_nbh = i - nup
                do j=1,5,1
                    kbh(i_nbh,j) = K(i,j)
                    pbh(i_nbh,j) = P(i,j)
                    vbh(i_nbh,j) = V(i,j)
                end do
            end do
!           'PYJETS' to 'sa1_h'.
            nn = nup
            do i2=1,5
                do i1=1,nup,1
                    kn(i1,i2) = K(i1,i2)
                    pn(i1,i2) = P(i1,i2)
                    rn(i1,i2) = V(i1,i2)
                end do
            end do
!----------------------------   Hadron Specifying  -----------------------------

            time_had = 0D0
            ijkk = 0

            ! In hadcas.f90.
            call hadcas( time_had, ijkk )

!           Gives up the current event due to the infinite collision loop.
            if(ijkk == 1)then
                i_error = 1
                return
!           Empty initial collision list.
            else if( ijkk == 2 )then
                goto 9000
            end if

!           'sa1_h' to 'PYJETS'.
            N = nn
            do i2=1,5
                do i1=1,N,1
                    K(i1,i2) = kn(i1,i2)
                    P(i1,i2) = pn(i1,i2)
                    V(i1,i2) = rn(i1,i2)
                end do
            end do

!----------------------------   Gamma 66 Removing  -----------------------------
!           Moves photons "66" from 'PYJETS' to 'sgam'.
            n66 = 0
            do j=1,N,1
                kf = K(j,2)
                if(kf == 22)then
                    K(j,2) = 66
                    n66 = n66 + 1
                end if
            end do
            if(n66 > 0) call remo_gam(66)
!----------------------------   Gamma 66 Removing  -----------------------------

!-------------------------------   sbh Moving  ---------------------------------
!           'sbh' to 'PYJETS'
            do l=1,nbh,1
                l1 = N + l
                do m=1,5,1
                    K(l1,m) = kbh(l,m)
                    P(l1,m) = pbh(l,m)
                    V(l1,m) = vbh(l,m)
                end do
            end do
            N = N + nbh
!-------------------------------   sbh Moving  ---------------------------------

9000        continue
            nbh = 0
            nn  = 0

        endif   ! 1 241103


!-------------------------------------------------------------------------------
!--------------------------   Final Particles Moving  --------------------------
!       Moves hadrons to the "time_hadcas".
        do i=1,N,1
            E = SQRT( P(i,1)**2 + P(i,2)**2 + P(i,3)**2 + P(i,5)**2 )
            V(i,1) = V(i,1) + ( time_had - V(i,4) ) * P(i,1) / E
            V(i,2) = V(i,2) + ( time_had - V(i,4) ) * P(i,2) / E
            V(i,3) = V(i,3) + ( time_had - V(i,4) ) * P(i,3) / E
            V(i,4) = time_had
        end do
!--------------------------   Final Particles Moving  --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Kaon Mixing  --------------------------------
!       Changes K0, K0bar to K0L and K0S.
        do j=1,N,1
            kf = K(j,2)
            if(kf == 311 .or. kf == -311)then
                rrlu   = PYR(1)
                K(j,2) = 130
                if( rrlu > 0.5D0 ) K(j,2) = 310
            end if
        end do
!-------------------------------   Kaon Mixing  --------------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_unstable_hadron_decay
!!      Performs the decay of unstable hadrons after the hadronic rescattering.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/delt/ndel,nodel,kdel(kszj,5),pdel(kszj,5),vdel(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!-------------------------------------------------------------------------------
!----------------------------   A-framework Decay   ----------------------------
        if(i_mode == 1)then
!       decay the delta
!       moves delta from 'PYJETS' to 'delt' (In parini_40.f90)
            call remo_delt
!       'PYJETS' to 'saf'
            naf = n
            do i2=1,5,1
                do i1=1,n,1
                    kaf(i1,i2) = K(i1,i2)
                    paf(i1,i2) = P(i1,i2)
                    vaf(i1,i2) = V(i1,i2)
                enddo
            enddo
!       delta decay one by one
            do i1=1,ndel
                ! In parini.f90.
                call padecy_delt(i1)
!       move decayed particle from 'PYJETS' to 'saf'
                do i3=1,n
                    naf = naf + 1
                    do i4=1,5
                        kaf(naf,i4) = K(i3,i4)
                        paf(naf,i4) = P(i3,i4)
                        vaf(naf,i4) = V(i3,i4)
                    enddo
                enddo
            enddo
!       'saf' to 'PYJETS' (In parini_40.f90)
            call tran_saf

!       decay the pi0
            i11=0
            n0=n
907         continue
            do i1 = i11+1, n0
                kf = K(i1,2)
                if(kf == 111)then
                    call pydecy(i1)
!                   Use PYEDIT(1) not (2).
!                   call pyedit(2)
                    call pyedit(1)
                    i11 = i1
                    goto 907
                endif
            enddo
!       Removes photons from "PYJETS" to "sgam".
            call remo_gam(22)
        endif
!----------------------------   A-framework Decay   ----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------------   Hadron Decay   -------------------------------
        ! In Pythia8_fort_interface.f90.
        if( i_mode /= 1 ) call PADECY
!------------------------------   Hadron Decay   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Kaon Mixing  --------------------------------
!       Changes K0, K0bar to K0L and K0S.
        do j=1,N,1
            kf = K(j,2)
            if(kf == 311 .or. kf == -311)then
                rrlu   = PYR(1)
                K(j,2) = 130
                if( rrlu > 0.5D0 ) K(j,2) = 310
            end if
        end do
!-------------------------------   Kaon Mixing  --------------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine output_final_particle_information
!!      Final particle information and OSCAR output.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/work7/reac(20),crose(20)


!       "main.out" log.
        if( nout == 1 .or. iii == 1 .or. mod(iii,nout) == 0 &
            .or. iii == neve )then
            write(MSTU(11),*)
            write(MSTU(11),*)
            write(MSTU(11),*) "#!------------------------------------" // &
                              "-----------------------------------------"
            write(MSTU(11),*) "#! i-event =", iii
            ! In p_40.f.
            call PYLIST(1)
            ! In main.f90.
            call prt_final_info
        end if

!       OSCAR standard output
        call oscar(4)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine analyze_event
!!      Analyzes an event and outputs event-averaged statistics.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        LOGICAL IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYINT1/MINT(400),VINT(400)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa6_sum/sthroq,sthrob,sthroc,sthroe(4)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa7/ispmax,isdmax,iflmax,ispkf(100),n_bin_hist,asd(10), &
         afl(100,10,2),i_y_or_eta
        common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb, &
         pph,vnlep
        common/sa23/kpar,knn,kpp,knp,kep,NON_p,skpar,sknn,skpp,sknp,skep
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa30/vneump,vneumt,mstptj
        common/sa35/ncpart,ncpar(kszj)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/ctllist/npctl,npinel(600),npctl0,npctlm
        common/ctllist_p/nreac(20),nrel
        common/ctllist_h/nctl,noinel(600),nctl0,noel
        common/schuds/schun,schudn,schudsn,sfra
        common/work7/reac(20),crose(20)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For PACIAE statistical analysis. ( f: full phase space; else partial )
        common/anly1/an(100,10,100),bn(100,2),anf(100,10,100),bnf(100,2)
        common/anly1_sum1/ san(100,10,100),  sbn(100,2), &
                          sanf(100,10,100), sbnf(100,2), &
                           sao(100,10,100),  sbo(100,2), &
                          saof(100,10,100), sbof(100,2)
        common/anly1_sum2/ dnmin, dnminf, dncha, dnchaf
        common/anly4/ distr(100,10,5), dMult(5,2), distr_f(100,10,5), &
         dMult_f(5,2), san_distr(100,10,5), sbn_mult(5,2), &
         san_distr_f(100,10,5), sbn_mult_f(5,2), sao_distr(100,10,5), &
         sbo_mult(5,2), sao_distr_f(100,10,5), sbo_mult_f(5,2)
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
        common/string_kapa2/ snnc,sadiv, sgtime,sitime,sgpmax, skapa(6)
!       For the independent Optical Clauber calculation.
        common/Glauber1/ kjp23, kjp24
        common/Glauber2/ bpp(20),acoll(20),acollp(20),acollt(20),sbp(20)
        common/Glauber3/ avb, avneu, astbp, astbt, aanbin
        common/Glauber4/ averb, psnon, psnop, psnot, anbin
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the parini related statistics.
        common/anly_parini1/ woun, swoun, snpctl0, snpctlm, snpar
!       For the parcas and hadcas related statistics.
        common/anly_parcas1/ sreac(20), rinel
        common/anly_hadcas1/ sel, sinel, dinel(600), einel(600)
!       For the Angantyr related statistics.
        common/anly_angantyr1/ weight_event, sum_weight_event
        common/anly_angantyr2/ sum_bParam_ANG(3), sum_sigma_ANG(8,2), &
                               sum_Ncoll_ANG(8), sum_Npart_ANG(4,2)
!       For the evolution time.
        common/anly_time/ time_ini,  time_par,  time_had, &
                         stime_ini, stime_par, stime_had
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension snreac(20),  sreaco(20)
        dimension dineli(600), eineli(600)
        dimension skapao(6)
        dimension wthroe(4)
        real*8 nmin,nminf,ncha,nchaf
        dimension dmax_throe_p(4), dmax_throe_p0(4)
        data dmax_throe_p / 4*0D0 /

!       weight_event (sum_): weight (total weight) of the current (total) event.
!         Normally = 1.0 (Neve), but non-one (non-Neve) for the biased event,
!         mainly for Angantyr mode.


!-------------------------------------------------------------------------------
!-----------------------------   Event Analysis   ------------------------------
!       Records the weight of event generated.
        weight_event = 1D0
        if( i_mode == 1 ) VINT(100) = 1D0
        if( ABS( VINT(100) ) >= 1D-10 ) weight_event = VINT(100)
        sum_weight_event = sum_weight_event + weight_event

!       analyse the event on-line ( In analy.f90)
        call analy( nmin, nminf, ncha, nchaf, &
                    weight_event, sum_weight_event )
        call analy_parton( weight_event, sum_weight_event )

!---------------------------   Event Accumulation   ----------------------------
!       sum over events
        dnmin  = dnmin  + nmin  * weight_event
        dnminf = dnminf + nminf * weight_event
        dncha  = dncha  + ncha  * weight_event
        dnchaf = dnchaf + nchaf * weight_event

        sbn  = sbn  + bn  * weight_event
        sbnf = sbnf + bnf * weight_event
        do kk=1,100,1
            do jj=1,10,1
                real_weight = weight_event
                if( jj == 7 .OR. jj == 8 ) real_weight = 1D0
                do ii=1,100,1
                    san(ii,jj,kk)  = san(ii,jj,kk)  &
                                   + an(ii,jj,kk)  * real_weight
                    sanf(ii,jj,kk) = sanf(ii,jj,kk) &
                                   + anf(ii,jj,kk) * real_weight
                end do
            end do
        end do
        sbn_mult    = sbn_mult    + dMult   * weight_event
        sbn_mult_f  = sbn_mult_f  + dMult_f * weight_event
        do kk=1,5,1
            do jj=1,10,1
                real_weight = weight_event
                if( jj == 7 .OR. jj == 8 ) real_weight = 1D0
                do ii=1,100,1
                    san_distr(ii,jj,kk)   = san_distr(ii,jj,kk)   &
                                          + distr(ii,jj,kk)   * real_weight
                    san_distr_f(ii,jj,kk) = san_distr_f(ii,jj,kk) &
                                          + distr_f(ii,jj,kk) * real_weight
                end do
            end do
        end do

        stime_ini = stime_ini + time_ini * weight_event
        stime_par = stime_par + time_par * weight_event
        stime_had = stime_had + time_had * weight_event

        call prt_sgam(0,egam,4)
        segam  = segam  + egam  * weight_event
        segam1 = segam1 + egam1 * weight_event
        segam2 = segam2 + egam2 * weight_event
        segam3 = segam3 + egam3 * weight_event

        snnc   = snnc   + nnc   * weight_event
        sadiv  = sadiv  + adiv  * weight_event
        sgpmax = sgpmax + gpmax * weight_event
        skapa  = skapa  + akapa * weight_event
        sgtime = sgtime + gtime * weight_event
        sitime = sitime + itime * weight_event

        spathn = spathn + pathn * weight_event
        sevbin = sevbin + evbin * weight_event

!       npinel(592): # of nn collision calling 'PYTHIA' in parini_40.f90
!       npinel(593): # of nn collision not calling 'PYTHIA' in parini_40.f90
!       noel : statistics of elastic collisions in hadcas_40.f90
!       noinel(i): statistics of the i-th inelastic collisions in hadcas_40.f90
        sel = sel + noel * weight_event
        do i1=1,600
            sinel = sinel + noinel(i1) * weight_event
            dinel(i1) = dinel(i1) + noinel(i1) * weight_event
            einel(i1) = einel(i1) + npinel(i1) * weight_event
        enddo

        swoun   = swoun   + woun   * weight_event
        snpctl0 = snpctl0 + npctl0 * weight_event
        snpctlm = snpctlm + npctlm * weight_event
        snpar   = snpar   + ncpart * weight_event
!       swoun: # of wounded nucleons sumed up over enents
!       snpctl0: # of nn collision pairs (in parini.f90) sumed up over enents
!       snpar: # of collided nucleons (in parini.f90) sumed up over enents

!       nrel: statistics of blocked parton-parton scattering process in
!        parcas_40.f90
!       nreac(i): statistics of successful i-th parton-parton scattering
!        process in parcas_40.f90
!       nreac -> reac
!       rinel = 0
        do i1=1,20,1
!           rinel = rinel + nreac(i1)
            rinel = rinel + reac(i1) * weight_event
        enddo
        sreac = sreac + reac * weight_event

        skpar = skpar + kpar * weight_event
        sknn  = sknn  + knn  * weight_event
        skpp  = skpp  + kpp  * weight_event
        sknp  = sknp  + knp  * weight_event
        skep  = skep  + kep  * weight_event

        sthroq = sthroq + (ithroq + ithroq_p) * weight_event
        sthrob = sthrob + (ithrob + ithrob_p) * weight_event
        sthroc = sthroc + (ithroc + ich_p   ) * weight_event
        sthroe = sthroe + (throe  + throe_p ) * weight_event
        dmax_throe_p0 = throe  + throe_p

        do j=1,4,1
            do i=1,nbe,1
                KF = kbe(i,2)
                charge = PYCHGE( KF )
                if( charge > 0D0 ) sthroq = sthroq + 1D0 * weight_event
                if( charge < 0D0 ) sthrob = sthrob + 1D0 * weight_event
                sthroc = sthroc + charge         * weight_event
                sthroe(j) = sthroe(j) + pbe(i,j) * weight_event
                dmax_throe_p0(j) = dmax_throe_p0(j) + pbe(i,j)
            end do
        end do

!       Max lost 4-momentum.
        do i=1,4,1
            if( ABS( dmax_throe_p0(i) ) > ABS( dmax_throe_p(i) ) ) &
                dmax_throe_p(i) = dmax_throe_p0(i)
        end do
!---------------------------   Event Accumulation   ----------------------------

!       Counts total-event cross sections.
        call PASTAT(0,0)

!-----------------------------   Event Analysis   ------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------------   Event Averaged Output   --------------------------
!       internal and final printing and controlled return
        if( MOD(iii,nout) == 0 .OR. iii == neve )then

        open( 12, file = "rms.out", status = "unknown" )

!       flaa: the number of current events in type "float".
!       flaa = iii*1D0
!       Uses the sum of weights.
!       Normally sum_weight_event = Neve, but is different for the biased-event.
        ! flaa_wgt = sum_weight_event
        flaa = sum_weight_event

!-----------------------------   Event Averaging   -----------------------------
!       averaged over events
        averbo = averb / flaa
        psnono = psnon / flaa
        psnopo = psnop / flaa
        psnoto = psnot / flaa

        dnmino  = dnmin  / flaa
        dnminfo = dnminf / flaa
        dnchao  = dncha  / flaa
        dnchafo = dnchaf / flaa

        sbo  = sbn  / flaa
        sbof = sbnf / flaa
        sao  = san  / flaa
        saof = sanf / flaa
        sbo_mult    = sbn_mult    / flaa
        sbo_mult_f  = sbn_mult_f  / flaa
        sao_distr   = san_distr   / flaa
        sao_distr_f = san_distr_f / flaa

        segamo  = segam  / flaa
        segam1o = segam1 / flaa
        segam2o = segam2 / flaa
        segam3o = segam3 / flaa

        stime_inio = stime_ini / flaa
        stime_paro = stime_par / flaa
        stime_hado = stime_had / flaa

        snnco   = snnc   / flaa
        sgtimeo = sgtime / flaa
!       sgtimeo: average number of gluons in a string
        skapao  = skapa / flaa
        sadivo  = sadiv  / flaa
        sgpmaxo = sgpmax / flaa
        sitimeo = sitime / flaa

        spathni = spathn / flaa
        sevbini = sevbin / flaa

        reli   = nrel  / flaa
        rineli = rinel / flaa
        seli   = sel   / flaa
        sineli = sinel / flaa
        dineli = dinel / flaa
        eineli = einel / flaa
        swouni   = swoun   / flaa
        snpctl0i = snpctl0 / flaa
        snpctlmi = snpctlm / flaa
        snpari   = snpar   / flaa
!       swouni: # of wounded nucleons averaged over enents after parini.f90

!       snpctl0i: # of nn collision pairs averaged over enents
!       snpctlmi: maximum # of nn collision pairs averaged over enents
!       snpari: MC Glauber-like <N_part>
        skparo = skpar / flaa
        sknno  = sknn  / flaa
        skppo  = skpp  / flaa
        sknpo  = sknp  / flaa
        skepo  = skep  / flaa

!       Amounts of q, qbar, charge and 4-momentum thrown away.
        wthroq = sthroq / flaa
        wthrob = sthrob / flaa
        wthroc = sthroc / flaa
        wthroe = sthroe / flaa

        nrea = 0
        do i1=1,20,1
            nrea = nrea + nreac(i1)
        enddo
        ! Note iii here.
        snreac = nreac*1D0  / iii
        ! Note iii here.
        snrea  = DBLE(nrea) / iii

        srea = 0D0
        sreaco = sreac / flaa
        do i1=1,20,1
            srea = srea + sreaco(i1)
        enddo

!-------------------------------   User Output   -------------------------------
!       User output
!       Header. In analy.f90.
        call output_header

!       Corrects output.
        if( ( .NOT.IS_NUCLEUS(KF_proj) .AND. .NOT.IS_NUCLEUS(KF_targ) ) &
            .AND. i_mode /= 8 .AND. i_mode /= 9 )then
            snpctl0i    = 1.
            snpari      = 2.
            snpctlmi    = 1.
            eineli(592) = 1.
            swouni      = 2.
            skparo      = 2.
            if( nzp == 0 .AND. nzt == 0 )then
                sknno = 1.
            elseif( ABS(nzp) == 1 .AND. ABS(nzt) == 1 )then
                skppo = 1.
            elseif( (ABS(nzp) == 1 .AND. nzt == 0) .OR. &
                (nzp == 0 .AND. ABS(nzt) == 1) )then
                sknpo = 1.
            end if
        end if

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! Number of events ="
        write(12,*)  iii
        write(12,*) "#! NN(pp) total cross section in parini (mb) ="
        write(12,*)  para1_1
        write(12,*) "#! NN(pp) total cross section in hadcas (mb) ="
        write(12,*)  para1_2
        write(12,*) "#! MC Glauber-like <N_coll>, <N_part> in the " // &
                    "initial collision list of parini ="
        write(12,*) snpctl0i, snpari
        write(12,*) "#! largest ave. # of NN collision pairs in " // &
                    "initial collision list of parini ="
        write(12,*) snpctlmi

        write(12,*) "#! ave. # of NN collision pairs calling PYTHIA, "// &
                    "not calling PYTHIA in parini ="
        write(12,*) eineli(592), eineli(593)
        write(12,*) "#! ave. # of wounded nucleons ="
        write(12,*) swouni

        write(12,*) "#! N_bin/N_part (optical Glauber)"
        write(12,*) spathni
        write(12,*) "#! event averaged N_bin (optical Glauber)"
        write(12,*) sevbini

        write(12,*) "#! Npart with pT=0, Nnn, Npp in parini ="
        write(12,*) skparo, sknno, skppo
        write(12,*) "#! Nnp, Ntot, Nep in parini ="
        write(12,*) sknpo, sknno+skppo+sknpo, skepo

        if( INT(psno) == 0 .OR. INT(psno) == 4 )then
            write(12,*) "#! psno, ave. b, Npart_p, Npart_t and N_bin "// &
                        "(optical Glauber) ="
            write(12,*) INT(psno), bp, vneump, vneumt, evbin
        else if( INT(psno) == 1 .OR. INT(psno) == 5 )then
            write(12,*) "#! psno, ave. b, Npart_p, " // &
                        "Npart_t, T_pt (optical Glauber) ="
            write(12,*) INT(psno), avb, astbp, astbt, aanbin
        else if( INT(psno) == 3 )then
            write(12,*) "#! psno, ave. b, Npart_p, Npart_t and N_bin "
            write(12,*) INT(psno), sum_bParam_ANG(1)  / flaa, &
                                  sum_Npart_ANG(2,1) / flaa, &
                                  sum_Npart_ANG(2,2) / flaa, &
                        ( sum_Ncoll_ANG(3) + sum_Ncoll_ANG(4) &
                        + sum_Ncoll_ANG(5) + sum_Ncoll_ANG(6) ) / flaa
        else
            write(12,*) "#! psno, ave. b, Npart_p, Npart_t and N_bin "// &
                        "(optical Glauber) ="
            write(12,*) INT(psno), averbo, psnopo, psnoto, psnono*csnn
        end if
        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! Angantyr mode weighted output:"
        write(12,*) "#! b (fm), bPhi, bWeight ="
        write(12,*) ( sum_bParam_ANG(i) / flaa, i=1,3,1 )
        write(12,*) "#! glauberTot, glauberTotErr, " // &
                    "glauberND, glauberNDErr (mb) ="
        write(12,*) sum_sigma_ANG(1,1) / flaa, sum_sigma_ANG(1,2) / flaa, &
                    sum_sigma_ANG(2,1) / flaa, sum_sigma_ANG(2,2) / flaa
        write(12,*) "#! glauberINEL, glauberINELErr, " // &
                    "glauberEL, glauberELErr (mb) ="
        write(12,*) sum_sigma_ANG(3,1) / flaa, sum_sigma_ANG(3,2) / flaa, &
                    sum_sigma_ANG(4,1) / flaa, sum_sigma_ANG(4,2) / flaa
        write(12,*) "#! glauberDiffP, glauberDiffPErr, " // &
                    "glauberDiffT, glauberDiffTErr (mb) ="
        write(12,*) sum_sigma_ANG(5,1) / flaa, sum_sigma_ANG(5,2) / flaa, &
                    sum_sigma_ANG(6,1) / flaa, sum_sigma_ANG(6,2) / flaa
        write(12,*) "#! glauberDDiff, glauberDDiffErr (mb), " // &
                    "glauberBSlope, glauberBSlopeErr (GeV^-2) ="
        write(12,*) sum_sigma_ANG(7,1) / flaa, sum_sigma_ANG(7,2) / flaa, &
                    sum_sigma_ANG(8,1) / flaa, sum_sigma_ANG(8,2) / flaa

        write(12,*) "#! nCollTot, nCollND, nCollSDP, nCollSDT ="
        write(12,*) sum_Ncoll_ANG(1) / flaa, sum_Ncoll_ANG(2) / flaa, &
                    sum_Ncoll_ANG(4) / flaa, sum_Ncoll_ANG(5) / flaa
        write(12,*) "#! nCollDD, nCollCD, nCollEL ="
        write(12,*) sum_Ncoll_ANG(6) / flaa, sum_Ncoll_ANG(7) / flaa, &
                    sum_Ncoll_ANG(8) / flaa
        write(12,*) "#! nPartProj, nAbsProj, nDiffProj, nElProj ="
        write(12,*) ( sum_Npart_ANG(i,1) / flaa, i=1,4,1 )
        write(12,*) "#! nPartTarg, nAbsTarg, nDiffTarg, nElTarg ="
        write(12,*) ( sum_Npart_ANG(i,2) / flaa, i=1,4,1 )

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! event averaged energy of gamma after " // &
                    "partonic initiation, partonic cascade,"
        write(12,*) "#!  hadronization and end of event ="
        write(12,*) segam1o, segam2o, segam3o, segamo
        if(ipden >= 11.and.ipden <= 16) &
         write(12,*)'#! event average number of lepton studied ='
        if(ipden >= 11.and.ipden <= 16) write(12,*) vnlep/flaa

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! # of successful, blocked and all collision " // &
                    "in parton cascade ="
        write(12,*) rineli, reli, reli+rineli

        write(12,*) "#! average collision # in parton cascade ="
        write(12,*) srea

        write(12,*) "#! # of scaterring processes in parton cascade " // &
                    "(q: light u-d-s quark; Q: heavy c-b quark; " // &
                    "g: gluon)"
        write(12,*) "#! q1 + q2 -> q1 + q2" // &
                    "     q1 + q1 -> q1 + q1" // &
                    "      q1 + q2bar -> q1 + q2bar" // &
                    "    q1 + q1bar -> q2 + q2bar" // &
                    "    q1 + q1bar -> q1 + q1bar"
        write(12,*) (sreaco(i1),i1=1,5)
        write(12,*) "#! q1 + q1bar -> g + g" // &
                    "    g + g -> q1 + q1bar" // &
                    "     q1 + g -> q1 + g" // &
                    "            g + g -> g + g" // &
                    "              Q + q -> Q + q"
        write(12,*) (sreaco(i1),i1=6,10)
        write(12,*) "#! Q + g -> Q + g" // &
                    "         q + qbar -> Q + Qbar" // &
                    "    g + g -> Q + Qbar" // &
                    "           Q + Qbar -> q + qbar" // &
                    "        Q + Qbar -> g + g"
        write(12,*) (sreaco(i1),i1=11,15)

        write(12,*) "#! average frequency of the occurring of each " // &
                    "inela. in hadron cascade (at the end of the file)"
!       write(12,*) dineli
        write(12,*)"#! el. and inel. coll. # and sum in hadron cascade="
        write(12,*) seli, sineli, seli+sineli

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! default parj1, parj2, parj3, parj4, par21 ="
        write(12,*) parj1, parj2, parj3, parj4, parj21
        write(12,*) "#! Eff-parj1, parj2, parj3, parj4, parj21, keff ="
        write(12,*) skapao
        write(12,*)'#! averaged # of gluon in a string when kjp22=1,3'
        write(12,*) sgtimeo
        write(12,*)"#! event averaged value of the factor related to # "
        write(12,*)"#!  of gluons and hardest gluon in a string, event "
        write(12,*)"#!  averaged transverse momentum of hardest gluon,"
        write(12,*)"#!  event averaged # strings when kjp22=1,3 ="
        write(12,*) sadivo, sgpmaxo, sitimeo

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! times & sum ="
        write(12,*) stime_inio, stime_paro, stime_hado, &
                    stime_inio+stime_paro+stime_hado

        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*) "#! q, qbar, charge thrown away ="
        write(12,*) wthroq, wthrob, wthroc/3.
        write(12,*) "#! 3-momentum and energy thrown away ="
        write(12,*) wthroe
        write(12,*) "#! Max 3-momentum and energy thrown away ="
        write(12,*) dmax_throe_p

        write(12,*)
        write(12,*) "#!-------------------------------------" // &
                    "----------------------------------------"
        write(12,*)"#! multiplicity of negative, positive particles " // &
                   "and sums, partial & full ="
        write(12,*) dnmino, dnchao, dnmino+dnchao
        write(12,*) dnminfo, dnchafo, dnminfo+dnchafo

!       Outputs multiplicities, abscissa, 7-distributions of pi+K+p and ispmax
!        particles specified in usu.dat. In analy.f90.
        call output_hadron_distribution
!       Outputs multiplicities, abscissa, 7-distributions of g, u+d+s +anti-
!        and q/qbar. In analy.f90.
        call output_parton_distribution( sum_weight_event )

        ! Empty lines.
        write(12,"(2/)")
        write(12,*) "#! average frequency of the occurring of each " // &
                    "inela. in hadron cascade ="
        do i=1,60,1
            j=(i-1)*10
            write(12,*) (dineli(j+kt),kt=1,10,1)
        end do

        if(ipden >= 11.and.ipden <= 16)then
            do kk=1,10,1
                sbn1=sbn(1,1)
                sbn1=dmax1(sbn1,1.d-20)
                sbnf1=sbnf(1,1)
                sbnf1=dmax1(sbnf1,1.d-20)
                if(kk /= 1) sbo(kk,1)=sbn(kk,1)/sbn1
                if(kk /= 1) sbof(kk,1)=sbnf(kk,1)/sbnf1
                do i1=1,n_bin_hist,1
                    do i2=1,isdmax,1
                        san1=san(i1,i2,1)
                        san1=dmax1(san1,1.d-20)
                        sanf1=sanf(i1,i2,1)
                        sanf1=dmax1(sanf1,1.d-20)
                        if(kk /= 1) sao(i1,i2,kk)=san(i1,i2,kk)/san1
                        if(kk /= 1) saof(i1,i2,kk)=sanf(i1,i2,kk)/sanf1
                    enddo
                enddo
            enddo
        write(12,*) "#! relative multiplicity, p =",(sbo(ll,1),ll=1,10)
        write(12,*) "#! relative multiplicity, f =",(sbof(ll,1),ll=1,10)
            do m2=1,isdmax,1
                write(12,*) "#! ID of relative distribution m2 =",m2
                do m3=1,10
                    write(12,*) "#! distribution belong to m3 =",m3
                    write(12,*) (sao(m1,m2,m3),m1=1,n_bin_hist,1)
                    write(12,*) (saof(m1,m2,m3),m1=1,n_bin_hist,1)
                enddo
            enddo
        endif

        close(12)
!-------------------------------   User Output   -------------------------------

        endif
!--------------------------   Event Averaged Output   --------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine print_time_consumption( time_simulation_start )
!!      Prints the time consumption of the program running.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        call CPU_TIME( time_simulation_end )

!       Basic unit is second.
        time_consumption = time_simulation_end - time_simulation_start
        i_time_h = INT( time_consumption/60/60 )
        i_time_m = INT( (time_consumption - i_time_h*60D0*60D0) / 60D0 )
        i_time_s = INT( time_consumption  - i_time_h*60D0*60D0 - i_time_m*60D0 )
        write(22,*)
        write(22,*) "Time consumption =", INT( i_time_h ), "h", &
                    INT( i_time_m ), "m", INT( i_time_s ), "s"


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine stop_event_at_this_point( i_stop_point )
!!      Stops the event evolution at the specified point.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen


        iadj140 = i_stop_point


!-------------------------------------------------------------------------------
!------------------------   Simulation Stopping 1 & 3  -------------------------
        if( iadj140 == 1 .or. iadj140 == 3 )then
            call PAEDIT(1)
            eevp = 0D0
            do i1=1,N,1
!       Sets them "stable". For the ones from PYTHIA 8.
                K(i1,1) = 1
                eevp = eevp + P(i1,4)
            end do
            eevh = 0D0
            do i1=1,nbh,1
!       Sets them "stable". For the ones from PYTHIA 8.
                kbh(i1,1) = 1
                eevh = eevh + pbh(i1,4)
            end do
            if( iadj140 == 3 )then
                mode_coal = 1
                ! In coales.f90
                call coales( mode_coal )
            end if
            n44 = 0
            do j=1,nbh,1
                kf = kbh(j,2)
                if(kf == 22)then
                    ! '44': prompt direct photon
                    kbh(j,2) = 44
                    n44 = n44 + 1
                end if
            end do
!       Moves "44" from 'sbh' to 'sgam'.
            if(n44 > 0) call remo_gam_sbh(44)

!       'sbh' to 'PYJETS'.
            do l=1,nbh,1
                l1 = N + l
                do m=1,5,1
                    K(l1,m) = kbh(l,m)
                    P(l1,m) = pbh(l,m)
                    V(l1,m) = vbh(l,m)
                end do
            end do
            N = N + nbh
        end if
!------------------------   Simulation Stopping 1 & 3  -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------------   Simulation Stopping 2  ---------------------------
        if( iadj140 == 2 )then
            call PAEDIT(1)
            eevp = 0D0
            do i1=1,N,1
!       Sets them "stable". For the ones from PYTHIA 8.
                K(i1,1) = 1
                eevp = eevp + P(i1,4)
            end do
            eevh = 0D0
            do i1=1,nbh,1
!       Sets them "stable". For the ones from PYTHIA 8.
                kbh(i1,1) = 1
                eevh = eevh + pbh(i1,4)
            end do
!       'sbh' to 'PYJETS'.
            do l=1,nbh,1
                l1 = N + l
                do m=1,5,1
                    K(l1,m) = kbh(l,m)
                    P(l1,m) = pbh(l,m)
                    V(l1,m) = vbh(l,m)
                end do
            end do
            N = N + nbh
        end if
!--------------------------   Simulation Stopping 2  ---------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_freeze_out_information( i_FO, i_part1, i_part2 )
!!      Records and updates the freeze-out information of particles, i.e.
!!       the information fr the last collision, to /PYJETS_FO/.
!       i_FO: =  0, records the information from /PYJETS/ to /PYJETS_FO/
!             = -1, loads the information from /PYJETS_FO/ to /PYJETS/
!             =  1, updates the information of a particle "i_part1" in /PYJETS/
!                  to "i_part2" in /PYJETS_FO/
!             =  2, updates the information of a particle "i_part1" in
!                  /PYJETS_FO/ to "i_part2" in /PYJETS/
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        !       Arrays of freeze-out particle information.
        COMMON/PYJETS_FO/ N_FO, NPAD_FO, &
                          K_FO(KSZJ,8), P_FO(KSZJ,7), V_FO(KSZJ,5)


!       All entries in "PYJETS" -> "PYJETS_FO".
        if( i_FO == 0 )then
            do j=1,5,1
                do i=1,N,1
                    K_FO(i,j) = K(i,j)
                    P_FO(i,j) = P(i,j)
                    V_FO(i,j) = V(i,j)
                end do
            end do
            N_FO = N
!       All entries in "PYJETS_FO" -> "PYJETS".
        else if( i_FO == -1 )then
!           Only the position and time information.
            do j=1,5,1
                do i=1,N,1
                    ! K(i,j) = K_FO(i,j)
                    ! P(i,j) = P_FO(i,j)
                    V(i,j) = V_FO(i,j)
                end do
            end do
!       One entry in "PYJETS" -> "PYJETS_FO".
        else if( i_FO == 1 )then
            do i=1,5,1
                K_FO( i_part2, i ) = K( i_part1, i )
                P_FO( i_part2, i ) = P( i_part1, i )
                V_FO( i_part2, i ) = V( i_part1, i )
            end do
!       One entry in "PYJETS_FO" -> "PYJETS".
        else if( i_FO == 2 )then
            do i=1,5,1
                K( i_part2, i ) = K_FO( i_part1, i )
                P( i_part2, i ) = P_FO( i_part1, i )
                V( i_part2, i ) = V_FO( i_part1, i )
            end do
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remop
!!      Moves q, qbar, g, diquark, and anti-diquark from 'PYJETS' to 'sbe'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_PARTON
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        jb = 0
201     continue
        do i1 = jb+1, N, 1
            KS = K(i1,1)
            KF = K(i1,2)
            if( .NOT.IS_PARTON(KF) .OR. .NOT.IS_EXIST(KS,i_mode) )then
                jb = jb + 1
                cycle
            end if
            nbe = nbe + 1
            do i2=1,5,1
                kbe(nbe,i2) = K(i1,i2)
                pbe(nbe,i2) = P(i1,i2)
                vbe(nbe,i2) = V(i1,i2)
            end do
!       Moves particle list 'PYJETS' one step downward from i1+1 to N.
            do j = i1+1, N, 1
                do jj=1,5,1
                    K( j-1, jj ) = K(j,jj)
                    P( j-1, jj ) = P(j,jj)
                    V( j-1, jj ) = V(j,jj)
                end do
            end do
            N = N - 1
            goto 201
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rest_hadronization( i_error )
!!      Hadronizes the rest partons ("sbe") those failed in sfm/coal by
!!       calling two subroutines "rest_sfm" and/or "rest_coal".
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!       PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA.
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_c/ithroq,ithrob,ithroc,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)
        common/ctllist_p/nreac(20),nrel
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3


        i_error = 0
!       The hadronization model.
        i_had_model = INT( adj1(12) )
        i_rest_had_model = INT( adj1(18) )

!       Removes partons from 'PYJETS' to 'sbe'.
        call remop

!       Appends "trs" to "sbe". (partons of inel. coll. in parcas with sfm)
        if( ( i_had_model == 0 .OR. i_had_model == 3 ) &
            .AND. ( iparres == 1 .OR. iparres ==3 ) .AND. ntrs > 0 )then
            do ii1=1,ntrs,1
                ii3 = nbe + ii1
                do ii2=1,5,1
                    kbe(ii3,ii2) = ktrs(ii1,ii2)
                    pbe(ii3,ii2) = ptrs(ii1,ii2)
                    vbe(ii3,ii2) = vtrs(ii1,ii2)
                end do
            end do
            nbe  = nbe + ntrs
            ntrs = 0
        end if

!       Appends "sa37" to "sbe" (partons failed in coal).
        if( i_had_model == 1 .AND.  nth > 0 )then
            do ii1=1,nth,1
                ii3 = nbe + ii1
                do ii2=1,5,1
                    kbe(ii3,ii2) = kth(ii1,ii2)
                    pbe(ii3,ii2) = pth(ii1,ii2)
                    vbe(ii3,ii2) = vth(ii1,ii2)
                end do
            end do
            nbe = nbe + nth
            nth = 0
            ithroq = 0
            ithrob = 0
            ich    = 0
            throe  = 0D0
        end if

        if( nbe == 0 ) return
        ! No treatment will be performed for the rest partons.
        if( i_rest_had_model < 0 ) return

!       Dumps "PYJETS" into "aaff".
        naff = N
        do j=1,5,1
            do i=1,N,1
                kaff(i,j) = K(i,j)
                paff(i,j) = P(i,j)
                vaff(i,j) = V(i,j)
            end do
        end do
        N = 0

!       Breaks up potential diquarks directly to reduce the complexity for
!        "rest_sfm" and "rest_coal". (In main.f90)
        if( nbe > 0 ) call break_sbe

        kkii = 0
        ! In main.f90.
        if( i_rest_had_model == 0 .AND. nbe > 1 ) call rest_sfm
!       Something goes wrong. Throws away this event.
        if( kkii == 4 )then
            i_error = 1
            iii = iii - 1
            return
        end if

!       Performs the coalescence even though adj18=0. (In main.f90)
        if( nbe > 0 ) call rest_coal

!       Appends "aaff" to "PYJETS". (dumps "aaff" into "PYJETS" in fact)
        do ii1=1,naff,1
            ii3 = N + ii1
            do ii2=1,5,1
                K(ii3,ii2) = kaff(ii1,ii2)
                P(ii3,ii2) = paff(ii1,ii2)
                V(ii3,ii2) = vaff(ii1,ii2)
            enddo
        enddo
        N = N + naff
        naff = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rest_sfm
!!      Hadronizes/Fragments "sbe" with "PYEXEC".
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
!       Local variables.
        integer n_g, n_q, n_b   ! g, q, qbar
        dimension ps1(6), ps0(6), rc(3)
        dimension k_g(100,5), p_g(100,5), v_g(100,5)   ! n_g for gluon
        dimension k_q(100,5), p_q(100,5), v_q(100,5)   ! n_q for quark/di-
        dimension k_b(100,5), p_b(100,5), v_b(100,5)   ! n_b for anti-quark/di-
        data n_warning / 0 /
        logical sfm_elastic_parcas, coal_wo_diquark, coal_wo_gluon, normal_run


!       The hadronization model.
        i_had_model = INT( adj1(12) )
        if( nbe < 2 ) return
!       If the normal hadronization (sfm / coal) was selected and code
!        run normally. The rest (failed) partons "sbe" should nbe == 0 !
        sfm_elastic_parcas = (i_had_model == 0 .OR. i_had_model == 3) &
                             .AND. ( iparres == 0 .OR. iparres == 2 )
        coal_wo_diquark = i_had_model == 1 &
                          .AND. ( iparres == 0 .OR. iparres == 1 )
        coal_wo_gluon = i_had_model == 1 .AND. (kjp22 == 1 .OR. kjp22 == 2 &
                        .OR. kjp22 == 6 .OR. kjp22 == 7 )
        normal_run = sfm_elastic_parcas .OR. coal_wo_diquark .OR. coal_wo_gluon
!       We give a the tolerance of "nbe <= 50".
        if( nbe > 50 .AND. normal_run )then
            n_warning = n_warning + 1
            if( n_warning <= 10 )then
                write(22,*)
                write(22,*) "Warning! nbe >= 50 in rest_sfm! Gives up " // &
                "this event. Please check the code! Event No.:", iii
                write(22,*)
            end if
            kkii = 4
            return
        end if

!       Dumps ("sbe") g to "n_g", q to "n_q", and qbar to "n_b", respectively.
        n_g = 0
        k_g = 0
        p_g = 0.
        v_g = 0.
        do i=1,nbe,1
            if( kbe(i,2) == 21 )then   ! g
                n_g = n_g + 1
                do j=1,5,1
                    k_g(n_g,j) = kbe(i,j)
                    p_g(n_g,j) = pbe(i,j)
                    v_g(n_g,j) = vbe(i,j)
                end do
            end if
        end do

        n_q = 0
        k_q = 0
        p_q = 0D0
        v_q = 0D0
        do i=1,nbe,1
            if( kbe(i,2) > 0 .AND. kbe(i,2) /= 21 )then   ! q/diquark
                n_q = n_q + 1
                do j=1,5,1
                    k_q(n_q,j) = kbe(i,j)
                    p_q(n_q,j) = pbe(i,j)
                    v_q(n_q,j) = vbe(i,j)
                end do
            end if
        end do

        n_b = 0
        k_b = 0
        p_b = 0.
        v_b = 0.
        do i=1,nbe,1
            if( kbe(i,2) < 0 )then   ! qbar/anti-diquark
                n_b = n_b + 1
                do j=1,5,1
                    k_b(n_b,j) = kbe(i,j)
                    p_b(n_b,j) = pbe(i,j)
                    v_b(n_b,j) = vbe(i,j)
                end do
            end if
        end do

        nbe = 0

!       Fragments a single pure g-...-g string in "n_g".
        iteration = 0   ! No more than 50 times.
100     continue
        if( n_g > 1 )then
            do i=1,n_g,1
                do kk=1,5,1
                    P(i,kk) = p_g(i,kk)
                    V(i,kk) = v_g(i,kk)
                end do
                K(i,1) = 2   ! 'A'-'I'-
                K(i,2) = k_g(i,2)
                K(i,3) = 0
                K(i,4) = 0
                K(i,5) = 0
            end do
            K(n_g,1) = 1   ! 'V'
            N = n_g
            iikk = 0
            kkii = 0
!           Sums of px, py, pz, E, inv. m, and charge before the "PYEXEC".
            ps0  = 0.
            do i=1,6,1
                ps0(i) = PAPYP(0,i)
            end do
!           Sets fragmentation flag in PYEXEC.
            MSTP(111) = 1

            ! Consider px, py, pz, E. In p_40.f
            call PYEXEC

!           Sums of px, py, pz, E, inv. m, and charge after the "PYEXEC".
            ps1 = 0.
            do i=1,6,1
                ps1(i) = PAPYP(0,i)
            end do
!           Charge is not conserved or errors occur. Re-generate the event.
            if( ABS(ps0(6)-ps1(6)) > 1D-10 .OR. &
                iikk == 2 .OR. kkii == 2 )then
                N = 0
                iteration = iteration + 1
                if(iteration > 50) then   ! Throws away these gluons.
                    do kk=1,5,1
                        do i=1,n_g,1
                            kbe(i,kk) = k_g(i,kk)
                            pbe(i,kk) = p_g(i,kk)
                            vbe(i,kk) = v_g(i,kk)
                        end do
                    end do
                    nbe = n_g
                    n_g = 0
                    goto 200
                end if
                goto 100
            end if

            ! Succeeded.
            do i=1,4,1
                throe_p(i) = throe_p(i) + ps0(i) - ps1(i)
            end do
            call PYEDIT(1)
            ! Moves "22" from "PYJETS" to "sgam"
            call remo_gam( 22 )

!           Dumps produced hadrons (from PYEXEC, no gamma) into "sa1_h".
            do jj=1,5,1
                do ii=1,N,1
                    kn(ii,jj) = K(ii,jj)
                    pn(ii,jj) = P(ii,jj)
                    rn(ii,jj) = V(ii,jj)
                enddo
            enddo
            nn  = N
!           The first position in "sa1_h" is selected from random one of the
!            original gluons.
            i_g = INT( PYR(1)*n_g + 1 )
            do jj=1,5,1
                rn(1,jj) = v_g(i_g,jj)
            end do

!           Arranges produced particles on the surface of sphere (with
!            radius rrp and centred on parent position), first produced
!            particle is put on its parent position originally.
            rrp = 1.16
            rc(1)   = rn(1,1)
            rc(2)   = rn(1,2)
            rc(3)   = rn(1,3)
            !#TODO(Lei20230630): may need to be extended.
            rn(1,4) = 0D0
            time = rn(1,4)
!           In corresponding with the time set in "posi". (In sfm.f90)
            call posi(1,nn,rc,rrp,time)

!           Transfers four position messages from "sa1_h" to "PYJETS".
            do jj=1,4,1
                do ii=1,nn,1
                    v(ii,jj) = rn(ii,jj)
                enddo
            enddo

            n_g = 0
!           Appends "PYJETS" to "aaff". (In sfm.f90)
            call tran_pyjets
            N = 0
        end if

200     continue

        sn = (n_q + n_b)*1D0
        if( sn <= 0.9 ) goto 400   ! No q and qbar.

!       Fragments simple (di-)qbar-(di-)q string.
!#TODO(Lei20240220): need improvement for diquarks.
300     continue
        if( n_b > 0 .AND. n_q > 0 )then   ! 1
            do i=1,n_b,1   ! 2
                do kk=1,5,1
                    P(1,kk) = p_b(i,kk)
                    V(1,kk) = v_b(i,kk)
                end do
                K(1,1) = 2   ! 'A'
                K(1,2) = k_b(i,2)
                K(1,3) = 0
                K(1,4) = 0
                K(1,5) = 0
                do j=1,n_q,1   ! 3

                    do kk=1,5,1
                        P(2,kk) = p_q(j,kk)
                        V(2,kk) = v_q(j,kk)
                    end do
                    K(2,1) = 1   ! 'V'
                    K(2,2) = k_q(j,2)
                    K(2,3) = 0
                    K(2,4) = 0
                    K(2,5) = 0
                    N = 2
                    iikk = 0
                    kkii = 0
!                   Sums of px, py, pz, E, inv. m, and charge before "PYEXEC".
                    ps0=0.
                    do ii=1,6,1
                        ps0(ii) = PAPYP(0,ii)
                    end do
!                   Sets fragmentation flag in PYEXEC.
                    MSTP(111) = 1

                    ! Considers px, py, pz, E. In p_40.f
                    call PYEXEC

!                   Sums of px, py, pz, E, inv. m, and charge after "PYEXEC".
                    ps1=0.
                    do ii=1,6,1
                        ps1(ii) = PAPYP(0,ii)
                    end do
!                   Charge is not conserved or errors occur. Re-generate string.
                    if( ABS(ps0(6)-ps1(6)) > 1D-10 .OR. &
                        iikk == 2 .OR. kkii == 2 )then
                        N = 0
                        if( j == n_q )then
                            nbe = nbe + 1
                            do kk=1,5,1
                                ! Throws away i-qbar/anti-diq to "sbe".
                                kbe(nbe,kk) = k_b(i,kk)
                                pbe(nbe,kk) = p_b(i,kk)
                                vbe(nbe,kk) = v_b(i,kk)
                                ! Moves one step down.
                                do ii=i+1,n_b,1
                                    k_b(i,kk) = k_b(ii,kk)
                                    p_b(i,kk) = p_b(ii,kk)
                                    v_b(i,kk) = v_b(ii,kk)
                                end do
                            enddo
                            n_b = n_b - 1
                            if(n_b == 0) goto 400
                            goto 300
                        end if
                        cycle
                    end if

                    ! Succeeded.
                    do ii=1,4,1
                        throe_p(ii) = throe_p(ii) + ps0(ii) - ps1(ii)
                    end do
                    call PYEDIT(1)
                    ! Moves "22" from "PYJETS" to "sgam"
                    call remo_gam(22)

!                   Dumps produced hadrons (from PYEXEC, no gamma) into "sa1_h".
                    do jj=1,5
                        do ii=1,N,1
                            kn(ii,jj) = k(ii,jj)
                            pn(ii,jj) = p(ii,j)
                            rn(ii,jj) = v(ii,jj)
                        enddo
                    enddo
                    nn  = N
!                   The first position in "sa1_h" is selected from random one
!                    of partons.
                    i_p = 0
                    if( PYR(1) > 0.5 ) i_p = 1
                    do jj=1,5
                        rn(1,jj) = v_b(i,jj)
                        if(i_p == 1) rn(1,jj) = v_q(j,jj)
                    end do

!                   Arranges produced particles on the surface of sphere (with
!                    radius rrp and centred on parent position), first produced
!                    particle is put on its parent position originally.
                    rrp = 1.16
                    rc(1)   = rn(1,1)
                    rc(2)   = rn(1,2)
                    rc(3)   = rn(1,3)
                    !#TODO(Lei20230630): may need to be extended.
                    rn(1,4) = 0D0
                    time = rn(1,4)
!                   In corresponding with the time set in >posi". (In main.f90)
                    call posi(1,nn,rc,rrp,time)

!                   Transfers four position messages from "sa1_h" to "PYJETS".
                    do jj=1,4,1
                        do ii=1,nn,1
                            v(ii,jj) = rn(ii,jj)
                        enddo
                    enddo
!                   "PYJETS" to "aaff" (In main.f90)
                    call tran_pyjets
                    N = 0
                    ! Moves one step down. (j-q/diq)
                    do kk=1,5
                        do ii=j+1,n_q,1
                            k_q(j,kk) = k_q(ii,kk)
                            p_q(j,kk) = p_q(ii,kk)
                            v_q(j,kk) = v_q(ii,kk)
                        end do
                    end do
                    n_q = n_q - 1
                    ! Moves one step down. (i-qbar/anti-diq)
                    do kk=1,5
                        do ii=i+1,n_b,1
                            k_b(i,kk) = k_b(ii,kk)
                            p_b(i,kk) = p_b(ii,kk)
                            v_b(i,kk) = v_b(ii,kk)
                        end do
                    end do
                    n_b = n_b - 1
                    goto 300
                end do   ! 3
            enddo   ! 2
        endif   ! 1

400     continue

!       Moves rest qbar/anti-diqaruk to 'sbe'. (for the case only one gluon)
        if(n_g > 0)then
            do i=1,n_g,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_g(i,j)
                    pbe(nbe,j) = p_g(i,j)
                    vbe(nbe,j) = v_g(i,j)
                enddo
            enddo
        endif
        n_g = 0

!       Moves rest qbar/anti-diqaruk to 'sbe'.
        if(n_b > 0)then
            do i=1,n_b,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_b(i,j)
                    pbe(nbe,j) = p_b(i,j)
                    vbe(nbe,j) = v_b(i,j)
                enddo
            enddo
        endif
        n_b = 0

!       Moves rest q/diqaruk to 'sbe'.
        if(n_q > 0)then
            do i=1,n_q,1
                nbe = nbe + 1
                do j=1,5
                    kbe(nbe,j) = k_q(i,j)
                    pbe(nbe,j) = p_q(i,j)
                    vbe(nbe,j) = v_q(i,j)
                enddo
            enddo
        endif
        n_q = 0

!       Gives proper status code, etc.
        if(nbe > 0)then
            do i=1,nbe,1
                kbe(i,1) = 2
                kbe(i,3) = 0
                kbe(i,4) = 0
                kbe(i,5) = 0
            end do
        end if
        N = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rest_coal
!!      Hadronizes/Coalesces "sbe" with "coal".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)


!       The hadronization model.
        i_had_model = INT( adj1(12) )
!       The coalescence mode.
        mode_coal = 1
        if( i_had_model == 1 ) mode_coal = kjp22
        iteration = 0

!       Processes the rest partons until empty.
12345   continue
        ! No more than 50 times.
        iteration = iteration + 1
        ! No rest partons left.
        if( nbe == 0 ) goto 54231


!-------------------------------------------------------------------------------
!------------------------------   Data Dumping   -------------------------------
!       Dumps "sbe" to "PYJETS".
        N = nbe
        do i2=1,5,1
            do i1=1,nbe,1
                K(i1,i2) = kbe(i1,i2)
                P(i1,i2) = pbe(i1,i2)
                V(i1,i2) = vbe(i1,i2)
            end do
        end do
        nbe = 0
!------------------------------   Data Dumping   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Junctions Removing   ---------------------------
!       Removes potential junctions if sfm in main.
        if( ( i_had_model == 0 .OR. i_had_model == 3 ) .AND. &
            iteration == 1 )then
            jb = 0
100         continue
            do i1=jb+1,N,1
                kf   = K(i1,2)
                kfab = iabs(kf)
                if( kfab /= 88 )then
                    jb = jb + 1
                    cycle
                end if
                ! In sfm.f90.
                call updad_pyj(N,i1+1,1)
                N = N - 1
                goto 100
            end do
        end if
!----------------------------   Junctions Removing   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Parton Coalescence   ---------------------------
        call do_gluon_splitting

        call do_quark_deexcitation

        if( i_had_model == 3 ) nbe = 0

        call coales( mode_coal )

        ! call compensate_conservation
!----------------------------   Parton Coalescence   ---------------------------
!-------------------------------------------------------------------------------


!       Moves potential partons from "PYJETS" to "sbe".
        call remop

!       Appends "PYJETS" to "aaff". (In main.f90)
        call tran_pyjets
        N = 0

!       Appends "sa37" to "sbe".
        do ii1=1,nth,1
            ii3 = nbe + ii1
            do ii2=1,5,1
                kbe(ii3,ii2) = kth(ii1,ii2)
                pbe(ii3,ii2) = pth(ii1,ii2)
                vbe(ii3,ii2) = vth(ii1,ii2)
            end do
        end do
        nbe = nbe + nth
        nth = 0
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0D0

54231   if( iteration > 50 ) return

        if( nbe > 0 ) goto 12345

        RETURN


! 300623 Lei Do not use the following method.

!       Assume the lost 4-momentum (positive energy) excites one gluon
!        (off-shell), then try re-coalescence again.
!       The inv. mass should be greater than the mass of a pion (~140 MeV).
        p1  = throe_p(1)
        p2  = throe_p(2)
        p3  = throe_p(3)
        p4  = throe_p(4)
        sm2 = p4**2 - p1**2 - p2**2 - p3**2
        if( p4 > 0. .AND. nbe == 0 .AND. sm2 > 0.14D0**2 )then
            ! Status.
            kbe(1,1) = 2
            kbe(1,2) = 21
            kbe(1,3) = 0
            kbe(1,4) = 0
            kbe(1,5) = 0
            ! Momentum.
            pbe(1,1) = p1
            pbe(1,2) = p2
            pbe(1,3) = p3
            pbe(1,4) = p4
            pbe(1,5) = 0D0
            ! Position. Assume 10*10*10 fm^3 volume.
            !           (1~10 fm/c of QGP evolution)
            vbe(1,1) = PYR(1)*10D0
            vbe(1,2) = PYR(1)*10D0
            vbe(1,3) = PYR(1)*10D0
            !#TODO(Lei20230630): time of parcas or 0D0 ?
            vbe(1,4) = PYR(1)*10D0
            vbe(1,5) = 0D0
            nbe = 1
            throe_p = 0D0
            goto 12345
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_PYJETS
!!      Shares the lost 4-momentum in "throe_p" with partons in "PYJETS".
!       The critirion is inv. m^2 > 0 and E >= m, or original inv. m^2 < 0
!         but the new one is closer to 0 and E >= m, after sharing.
!       Sometimes there are junctions, hadrons (color neutral) and historical
!         entries, which should be excluded.   ! 300623 Lei 020224 Lei
!       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension p0(4),dp(4)


        if(N == 0) return
!       Potential number of junctions, hadrons and historical entries.
        n_junc = 0
        do i=1,N,1
!           Color neutral state ( PYK_12 =0, hadron).
            if( PYK(i,12) == 0 )then
                n_junc = n_junc + 1
!           Junction and historical entries in PYTHIA 6.
            else if( (i_mode == 2 .OR. i_mode == 3) .AND. &
                     (K(i,2) == 88 .OR. K(i,1) < 1 .OR. &
                                        K(i,1) > 10) )then
                n_junc = n_junc + 1
!           Junction (?) and historical entries in PYTHIA 8.
            else if( (i_mode == 5 .OR. i_mode == 6 .OR. &
                      i_mode == 8 .OR. i_mode == 9) .AND. &
                     (K(i,2) == 88 .OR. K(i,1) < 1) )then
                n_junc = n_junc + 1
            end if
        end do
        if(N == n_junc) return

        sn = (N - n_junc)*1D0
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
        do i1=1,N,1
            i_color = PYK(i1,12)
            ks = K(i1,1)
            kf = K(i1,2)
!           Exclude hadrons, junctions (PYTHIA 6) and historical entries.
            if( i_color == 0 .OR. ks < 1 .OR. kf == 88 ) cycle
!           Exclude extra historical entries in PYTHIA 6.
            if( (i_mode == 2 .OR. i_mode == 3) .AND. ks > 10 ) cycle
            p1 = P(i1,1)
            p2 = P(i1,2)
            p3 = P(i1,3)
            p4 = P(i1,4)
            p5 = P(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
!           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction,
!            hadrons or historical entries.
            if( (sm2 >= 0. .OR. (sm2 < 0. .AND. sm2 >= sm2_org)) .AND. &
                (p4+dp(4)) >= p5 )then
                do i2=1,4,1
                    P(i1,i2) = P(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
!       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing
!        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)) > 1D-15 .AND. &
                ABS(throe_p(i1)-p0(i1))  >  1D-15 ) goto 100
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_sbe
!!      Shares the lost 4-momentum in "throe_p" with particles in "sbe".
!       The critirion is inv. m^2 > 0 and E >= m, or original inv. m^2 < 0
!         but the new one is closer to 0 and E >= m, after sharing.
!       Sometimes there are junctions, which should be excluded.   ! 300623 Lei
!       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sbe/nbe,npad,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension p0(4),dp(4)

        if(nbe == 0) return
!       Potential number of junctions and spectators.
        n_junc = 0
        do i=1,nbe,1
            if( kbe(i,2) == 88 .OR. (pbe(i,1)**2+pbe(i,2)**2) <= 1D-15 ) &
                n_junc = n_junc + 1
        end do
        if(nbe == n_junc) return

        sn = (nbe - n_junc)*1D0
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
        do i1=1,nbe,1
            kf = kbe(i1,2)
            p1 = pbe(i1,1)
            p2 = pbe(i1,2)
            p3 = pbe(i1,3)
            p4 = pbe(i1,4)
            p5 = pbe(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
!           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2 >= 0. .OR. (sm2 < 0. .AND. sm2 >= sm2_org)) .AND. &
                (p4+dp(4)) >= p5 .AND. kf /= 88 .AND. pT2 > 1D-15 )then
    !             (p4+dp(4)) > 0. .AND. kf /= 88 .AND. pT2 > 1D-15 )then
                do i2=1,4,1
                    pbe(i1,i2) = pbe(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
!       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing
!        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)) > 1D-15 .AND. &
                ABS(throe_p(i1)-p0(i1))  >  1D-15 ) goto 100
        end do

        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine share_p_PYJETS_sa1h
!!      Shares the lost 4-momentum in "throe_p" with particles in "PYJETS".
!       The critirion is inv. m^2 > 0 and E >= m, or original inv. m^2 < 0
!         but the new one is closer to 0 and E >= m, after sharing.
!       Sometimes there are junctions, which should be excluded.   ! 300623 Lei
!       "throe_p" is 4-momentum accumulated before this calling.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension p0(4),dp(4)

        if( (N+nn) == 0 ) return
!       Potential number of junctions and spectators.
        n_junc = 0
    !     do i=1,N,1
    !         if( K(i,2) == 88 .OR. (P(i,1)**2+P(i,2)**2) <= 1D-15 )
    !  &          n_junc = n_junc + 1
    !     end do
    !     do i=1,nn,1
    !         if( kn(i,2) == 88 .OR. (pn(i,1)**2+pn(i,2)**2) <= 1D-15 )
    !  &          n_junc = n_junc + 1
    !     end do
    !     if( (N+nn) == n_junc ) return

        sn = (N + nn - n_junc)*1D0
100     continue
        dp = throe_p/sn
        p0 = throe_p   ! Old
!       Share with PYJETS.
        do i1=1,N,1
            kf = K(i1,2)
            p1 = P(i1,1)
            p2 = P(i1,2)
            p3 = P(i1,3)
            p4 = P(i1,4)
            p5 = P(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
!           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2 >= 0. .OR. (sm2 < 0. .AND. sm2 >= sm2_org)) .AND. &
                (p4+dp(4)) >= p5 .AND. kf /= 88 .AND. pT2 > 1D-15 )then
    !             (p4+dp(4)) > 0. .AND. kf /= 88 .AND. pT2 > 1D-15 )then
                do i2=1,4,1
                    P(i1,i2) = P(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
!       Share with sa1_h.
        do i1=1,nn,1
            kf = kn(i1,2)
            p1 = pn(i1,1)
            p2 = pn(i1,2)
            p3 = pn(i1,3)
            p4 = pn(i1,4)
            p5 = pn(i1,5)
            pT2= p1**2 + p2**2
            sm2_org = p4**2 - p1**2 - p2**2 - p3**2
            sm2=(p4+dp(4))**2-(p1+dp(1))**2-(p2+dp(2))**2-(p3+dp(3))**2
!           Inv. m^2 >= 0 or < 0 but closer to 0, E >= m, and not junction or spectator.
            if( (sm2 >= 0. .OR. (sm2 < 0. .AND. sm2 >= sm2_org)) .AND. &
                (p4+dp(4)) >= p5 .AND. kf /= 88 .AND. pT2 > 1D-15 )then
    !             (p4+dp(4)) > 0. .AND. kf /= 88 .AND. pT2 > 1D-15 )then
                do i2=1,4,1
                    pn(i1,i2) = pn(i1,i2) + dp(i2)
                end do
                throe_p = throe_p - dp   ! New
            end if
        end do
!       If remaining 4-momentum > sigma (accuracy, 1D-15 here) and sharing
!        suceeded, try another sharing (iteration).
        do i1=1,4,1   ! Compares the old with new one.
            if( ABS(throe_p(i1)) > 1D-15 .AND. &
                ABS(throe_p(i1)-p0(i1))  >  1D-15 ) goto 100
        end do

        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_final_info
!!      Prints the sum of momentum and energy for NN, NA(AN) and AA.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYINT1/MINT(400),VINT(400)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       For the charge and 4-momentum check.
        common/cp_check/ p_init(4), p_final(4)
        dimension p_h(4),p_l(4),p_rest(4),p_lost(4), p_tot(4)


        p_h    = 0.
        p_l    = 0.
        p_rest = 0.
        p_lost = 0.

!       Initial c & p.
        ic_init = 0
        if( .NOT.IS_NUCLEUS(KF_proj) )then
            ic_init = ic_init + PYCHGE( KF_proj )
        else
            ic_init = ic_init + nzp*PYCHGE(2212)
        end if
        if( .NOT.IS_NUCLEUS(KF_targ) )then
            ic_init = ic_init + PYCHGE( KF_targ )
        else
            ic_init = ic_init + nzt*PYCHGE(2212)
        end if
!       For Angantyr mode.
        if( i_mode == 8 .OR. i_mode == 9 )then
            p_init( 1 ) = VINT(397)
            p_init( 2 ) = VINT(398)
            p_init( 3 ) = VINT(399)
            p_init( 4 ) = VINT(400)
        end if

!       Final hadron ( + part of lepton e/mu/tau ), "PYJETS"
        ic_h = 0
        do i=1,N,1
            KS   = K(i,1)
            KF   = K(i,2)
!       Nuclui remnants from Angantyr.
            if( IS_EXIST(KS,i_mode) .AND. IS_NUCLEUS(KF) )then
                KF_ZZZ = MOD( KF, 10000000 ) / 10000
                ic_h = ic_h + KF_ZZZ*3
            else
                ic_h = ic_h + PYCHGE(kf)
            end if
            do j=1,4,1
                p_h(j) = p_h(j) + P(i,j)
            end do
        end do

!       Final lepton ( gamma ), "sgam"
        ic_l = 0
        do i=1,ngam,1
            kf   = kgam(i,2)
            ic_l = ic_l + PYCHGE(kf)
            do j=1,4,1
                p_l(j) = p_l(j) + pgam(i,j)
            end do
        end do

!       Rest parton, "sbe"
        ic_rest = 0
        do i=1,nbe,1
            kf   = kbe(i,2)
            ic_rest = ic_rest + PYCHGE(kf)
            do j=1,4,1
                p_rest(j) = p_rest(j) + pbe(i,j)
            end do
        end do

!       Lost charge and 4-momentum.
        ic_lost = ich + ich_p
        p_lost  = throe_p + throe

!       Total.
        ic_tot = ic_h + ic_l + ic_rest + ic_lost
        p_tot  = p_h + p_l + p_rest + p_lost

!       Prints information of 3-momentum and energy.
        write(22,*) "------------------------------------------------"// &
                    "------------------------------------------------"// &
                    "------------------------------------------------"
        write(22,*) "Initial              c & p =", &
                                REAL(ic_init) / 3., p_init
        write(22,*) "Final   ( h + l )    c & p =", &
                            REAL(ic_h + ic_l) / 3., (p_h + p_l)
        write(22,*) "Rest parton          c & p =", &
                                REAL(ic_rest) / 3., p_rest
        write(22,*) "Final + rest         c & p =", &
                REAL(ic_h + ic_l + ic_rest) / 3., (p_h + p_l + p_rest)
        write(22,*) "Lost                 c & p =", &
                                REAL(ic_lost) /3., p_lost
        write(22,*) "Final + rest + lost  c & p =", &
            REAL(ic_h + ic_l + ic_rest + ic_lost) / 3., &
                    (p_h + p_l + p_rest + p_lost)
        write(22,*) "------------------------------------------------"// &
                    "------------------------------------------------"// &
                    "------------------------------------------------"

!       Failure and gamma information.
        if( ABS(ppsa(5)) > 1D-10 ) write(22,*) "ppsa=", ppsa
        ! if( INT(adj1(40)) > 2 )then
!           Prints the partons those failed in coales.
            if(nth > 0)then
        write(22,*) "Summary and list of partons thrown away in coales"
                ! in coales.f90
                call prt_sa37(nth,cc)
            endif
!           Prints the rest partons those cannot be processed.
            if(nbe > 0)then
                write(22,*) "Summary and list of rest partons"
                ! in parini.f90
                call prt_sbe(nbe,cc)
            endif
!           Prints gamma.
            if(ngam > 0)then
                write(22,*)'summary and list of gammas'
                call prt_sgam(ngam,egam,8)
            endif
        ! endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tran_pyjets
!!      Appends 'PYJETS' to 'aaff'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)


        do l=1,N,1
            l1 = naff + l
            do m=1,5,1
                kaff(l1,m) = K(l,m)
                paff(l1,m) = P(l,m)
                vaff(l1,m) = V(l,m)
            end do
        end do
        naff = naff + N


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!++++++++++++++++++++++++++   OSCAR Standard Output   ++++++++++++++++++++++++++
        subroutine oscar(i_stage)
!!      Records final particle information or event history.
!!      Complies with the OSC1997A/1999A/2013A standard format.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYINT1/MINT(400),VINT(400)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
         iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/oscar0/ n0, npad0, k0(1000,5), p0(1000,5), v0(1000,5)
        common/oscar1/ n1, npad1, k1(KSZJ,5), p1(KSZJ,5), v1(KSZJ,5)
        common/oscar2/ n2, npad2, k2(KSZJ,5), p2(KSZJ,5), v2(KSZJ,5)
        common/oscar3/ n3, npad3, k3(KSZJ,5), p3(KSZJ,5), v3(KSZJ,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        character(8) PACIAE, VERSION, frame
        data PACIAE  / "PACIAE  "/
        data VERSION / "4.003   "/
!       win : energy
!       i_stage =-1 : prints file header.
!               = 0 : initial nucleons configuration.
!               = 1 : initial parton state (IPS), after parini calling.
!               = 2 : final parton state (FPS), after parcas calling.
!               = 3 : intermediate hadron state (IHS), after sfm/coal calling.
!               = 4 : final program output,after PACIAE running.(Final state FS)
!                     Normally final hadron state (FHS), after hadcas calling.
!                     May be final parton state if adj1(1)=0 & kjp21=0, etc.
!       nosc = 0 : no OSCAR output.
!            = 1 : OSCAR1997A (final_id_p_x, just PACIAE final output)
!            = 2 : OSCAR1999A (full_event_history)
!            = 3 : OSCAR2013A (full_event_history, dummy now)


!***********************************************************************
!***************************** File Header *****************************
        if(i_stage == -1)then
            ntest = 1
            if(ifram == 0)then
                frame = "lab"
                ebeam = win
            else if(ifram == 1)then
                frame = "nncm"
            else
                frame = "unknown"
                ! stop
            end if
            ebeam = win

!cccccccccccccccccccccccc   No OSCAR Output   cccccccccccccccccccccccccc
            if(nosc == 0)then

!cccccccccccccccccccccccccccc   OSC1997A   ccccccccccccccccccccccccccccc
            else if(nosc == 1)then
                write(34,"(A8)")  "OSC1997A"
                write(34,"(A12)") "final_id_p_x"
                write(34,100) PACIAE, VERSION, nap, nzp, nat, nzt, &
                              frame, ebeam, ntest
100             format(2(A8,2X),"(", I3, ",", I6, ")+(", I3, ",", &
                       I6, ")", 2X, A4, 2X, E10.4, 2X, I8)

!cccccccccccccccccccccccccccc   OSC1999A   ccccccccccccccccccccccccccccc
            else if(nosc == 2)then
                write(34,"(A10)") "# OSC1999A"
                write(34,"(A20)") "# full_event_history"
                write(34,"(A14)") "# PACIAE " //TRIM((ADJUSTL(VERSION)))
                write(34,200) nap, nzp, nat, nzt, frame, ebeam, ntest
200             format("# (", I3, ",", I6, ")+(", I3, ",", &
                       I6, ")", 2X, A4, 2X, E10.4, 2X, I8)

!cccccccccccccccccccccccccccc   OSC2013A   ccccccccccccccccccccccccccccc
            !#TODO(Lei20230214): need to extend.
            else if(nosc == 3)then
                write(34,*)"#!OSCAR2013 full_event_history ID "// &
                           "px py pz E m x y z t"
                write(34,*)"# Units: none "// &
                           "GeV GeV GeV GeV GeV fm fm fm fm "
                write(34,"(A14)") "# PACIAE " //TRIM((ADJUSTL(VERSION)))
            end if

            n0 = 0
            n1 = 0
            n2 = 0
            n3 = 0
            return
        endif
!***************************** File Header *****************************
!***********************************************************************


!***********************************************************************
!***************************** Data Dumping ****************************
        ! Initial nucleon configuration.
        if(i_stage == 0)then
            n0 = N
            do j=1,5,1
                do i=1,N,1
                    k0(i,j) = K(i,j)
                    p0(i,j) = P(i,j)
                    v0(i,j) = V(i,j)
                end do
            end do
        ! Initial parton state (IPS).
        else if(i_stage == 1)then
            n1 = 0
            do i=1,N,1
                KS = K(i,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                do j=1,5,1
                    k1(n1+1,j) = K(i,j)
                    p1(n1+1,j) = P(i,j)
                    v1(n1+1,j) = V(i,j)
                end do
                n1 = n1 + 1
            end do
            ! Sometime there might exist hadrons or leptons, which was generated
            !  from the spectators, the diffractive & the special sub-processes,
            !  e.g. NRQCD onia, or hadron beam remnants, or B-framework.
            !  We store them here, too.
            if(nbh > 0)then
                do i=1,nbh,1
                    KS = kbh(i,1)
                    if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                    do j=1,5,1
                        k1(n1+1,j) = kbh(i,j)
                        p1(n1+1,j) = pbh(i,j)
                        v1(n1+1,j) = vbh(i,j)
                    end do
                    n1 = n1 + 1
                end do
            end if
        ! Final parton state (FPS).
        else if(i_stage == 2)then
            n2 = 0
            do i=1,N,1
                KS = K(i,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                do j=1,5,1
                    k2(n2+1,j) = K(i,j)
                    p2(n2+1,j) = P(i,j)
                    v2(n2+1,j) = V(i,j)
                end do
                n2 = n2 + 1
            end do
            ! Sometime there might exist hadrons or leptons, which was generated
            !  from the spectators, the diffractive & the special sub-processes,
            !  e.g. NRQCD onia, or hadron beam remnants, or B-framework.
            !  We store them here, too.
            if(nbh > 0)then
                do i=1,nbh,1
                    ! Only for the particles that still exist.
                    KS = kbh(i,1)
                    if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                    do j=1,5,1
                        k2(n2+1,j) = kbh(i,j)
                        p2(n2+1,j) = pbh(i,j)
                        v2(n2+1,j) = vbh(i,j)
                    end do
                    n2 = n2 + 1
                end do
            end if
        ! Hadronization.
        else if(i_stage == 3)then
            n3 = N
            do j=1,5,1
                do i=1,N,1
                    k3(i,j) = K(i,j)
                    p3(i,j) = P(i,j)
                    v3(i,j) = V(i,j)
                end do
            end do
        end if
        if(i_stage < 4) return
!***************************** Data Dumping ****************************
!***********************************************************************


!***********************************************************************
!***************************** Event Block *****************************
!cccccccccccccccccccccccccccc   OSC1997A   ccccccccccccccccccccccccccccc
        phi = VINT(392)
        weight = 1D0
        if( ABS(VINT(100)) >= 1D-10 ) weight = VINT(100)
!       Prints final particles.
        if(nosc == 1)then
            write(34,101) iii, neve, N, bp, phi, weight
            do i=1,N,1
                KS = K(i,1)
                KF = K(i,2)
                if( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) cycle
                ! Puts particles on-shell.
                dm = P(i,5)
                if( .NOT.IS_NUCLEUS(KF) ) &
                E = SQRT( P(i,1)**2 + P(i,2)**2 + P(i,3)**2 + dm**2 )
                write(34,102) i, K(i,2), (P(i,j),j=1,3), E, dm, (V(i,j),j=1,4)
            end do
101         format(I10, 2(2X, I10), 3(2X, F8.3))
102         format(I10, 2X, I10, 2X, 9(E12.6, 2X))
!cccccccccccccccccccccccccccc   OSC1997A   ccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccc   OSC1999A   ccccccccccccccccccccccccccccc
!       Prints stage 0, 1, 2, 3 and 4 at once.
        else if(nosc == 2)then
!           Initial nucleons configuration.
            mark_stage = 0
            write(34,201) iii, neve, n0, bp, phi, mark_stage, weight
            do i=1,n0,1
              write(34,202) i, k0(i,2), (p0(i,j),j=1,5), (v0(i,j),j=1,4)
            end do
!           Initial parton state (IPS).
            mark_stage = 1
            write(34,201) iii, neve, n1, bp, phi, mark_stage, weight
            do i=1,n1,1
                KS = k1(i,1)
                KF = k1(i,2)
                if( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) cycle
                ! Puts particles on-shell.
                dm = p1(i,5)
                if( .NOT.IS_NUCLEUS(KF) ) &
                E = SQRT( p1(i,1)**2 + p1(i,2)**2 + p1(i,3)**2 + dm**2 )
              write(34,202) i, k1(i,2), (p1(i,j),j=1,3), E, dm, (v1(i,j),j=1,4)
            end do
!           Final parton state (FPS).
            mark_stage = 2
            write(34,201) iii, neve, n2, bp, phi, mark_stage, weight
            do i=1,n2,1
                KS = k2(i,1)
                KF = k2(i,2)
                if( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) cycle
                ! Puts particles on-shell.
                dm = p2(i,5)
                if( .NOT.IS_NUCLEUS(KF) ) &
                E = SQRT( p2(i,1)**2 + p2(i,2)**2 + p2(i,3)**2 + dm**2 )
              write(34,202) i, k2(i,2), (p2(i,j),j=1,3), E, dm, (v2(i,j),j=1,4)
            end do
!           Intermediate hadron state (IHS).
            mark_stage = 3
            write(34,201) iii, neve, n3, bp, phi, mark_stage, weight
            do i=1,n3,1
                KS = k3(i,1)
                KF = k3(i,2)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                ! Puts particles on-shell.
                dm = p3(i,5)
                if( .NOT.IS_NUCLEUS(KF) ) &
                E = SQRT( p3(i,1)**2 + p3(i,2)**2 + p3(i,3)**2 + dm**2 )
              write(34,202) i, k3(i,2), (p3(i,j),j=1,3), E, dm, (v3(i,j),j=1,4)
            end do
!           Final state / Final hadron state (FS/FHS).
            write(34,201) iii, neve, N, bp, phi, i_stage, weight
            do i=1,N,1
                KS = K(i,1)
                KF = K(i,2)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                ! Puts particles on-shell.
                dm = P(i,5)
                if( .NOT.IS_NUCLEUS(KF) ) &
                E = SQRT( P(i,1)**2 + P(i,2)**2 + P(i,3)**2 + dm**2 )
                write(34,202) i, K(i,2), (P(i,j),j=1,3), E, dm, (V(i,j),j=1,4)
            end do
!       Prints stage 0, 1, 2, 3 and 4 at once.
201         format(I10, 2(2X, I10), 2(2X, F8.3), 2X, I10, 2X, F8.3)
202         format(I10, 2X, I10, 2X, 9(E12.6, 2X))
!cccccccccccccccccccccccccccc   OSC1999A   ccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccc   OSC2013A   ccccccccccccccccccccccccccccc
        !#TODO(Lei20230214): need to extend.
        else if(nosc == 3)then
            ! Dummy.
        end if
!cccccccccccccccccccccccccccc   OSC2013A   ccccccccccccccccccccccccccccc


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam_sbh(ii)
!!      move particles with flavor code ii from  'sbh' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sbh/n,nonj,k(kszj,5),p(kszj,5),v(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        jb=0
201     do i1=jb+1,n
            kf=k(i1,2)
            if(kf /= ii)then
                jb=jb+1
                goto 202
            endif
            ngam=ngam+1
            do i2=1,5
                kgam(ngam,i2)=k(i1,i2)
                pgam(ngam,i2)=p(i1,i2)
                vgam(ngam,i2)=v(i1,i2)
            enddo
            if(i1 == n)then
                n=n-1
                goto 203
            endif
!       move particle list 'sbh' one step downward from i1+1 to n
            do j=i1+1,n
                j1=j-1
                do jj=1,5
                    k(j1,jj)=k(j,jj)
                    p(j1,jj)=p(j,jj)
                    v(j1,jj)=v(j,jj)
                enddo
            enddo
            n=n-1
            goto 201
202     enddo
203     continue
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam_par(ii)
!!      move particles with flavor code ii ('55') from 'PYJETS' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        jb=0
201     do i1=jb+1,n
            kf=k(i1,2)
            if(kf /= ii)then
                jb=jb+1
                goto 202
            endif
            ngam=ngam+1
            do i2=1,5
                kgam(ngam,i2)=k(i1,i2)
                pgam(ngam,i2)=p(i1,i2)
                vgam(ngam,i2)=v(i1,i2)
            enddo
            if(i1 == n)then
                n=n-1
                goto 203
            endif
!       move particle list 'PYJETS' one step downward from i1+1 to n
            do j=i1+1,n
                j1=j-1
                do jj=1,5
                    k(j1,jj)=k(j,jj)
                    p(j1,jj)=p(j,jj)
                    v(j1,jj)=v(j,jj)
                enddo
            enddo
            n=n-1
            goto 201
202     enddo
203     continue
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam(ii)
!!      move particles with flavor code ii from  'PYJETS' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        jb=0
201     do i1=jb+1,n
            kf=k(i1,2)
            if(kf /= ii)then
                jb=jb+1
                goto 202
            endif
            ngam=ngam+1
            do i2=1,5
                kgam(ngam,i2)=k(i1,i2)
                pgam(ngam,i2)=p(i1,i2)
                vgam(ngam,i2)=v(i1,i2)
            enddo
            if(i1 == n)then
                n=n-1
                goto 203
            endif
!       move particle list 'PYJETS' one step downward from i1+1 to n
            do j=i1+1,n
                j1=j-1
                do jj=1,5
                    k(j1,jj)=k(j,jj)
                    p(j1,jj)=p(j,jj)
                    v(j1,jj)=v(j,jj)
                enddo
            enddo
            n=n-1
            goto 201
202     enddo
203     continue
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_gam_hadro(ii)
!!      move particles with flavor code ii ('77') from  'PYJETS' to 'sgam'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        jb=0
201     do i1=jb+1,n
            kf=k(i1,2)
            if(kf /= ii)then
                jb=jb+1
                goto 202
            endif
            ngam=ngam+1
            do i2=1,5
                kgam(ngam,i2)=k(i1,i2)
                pgam(ngam,i2)=p(i1,i2)
                vgam(ngam,i2)=v(i1,i2)
            enddo
            if(i1 == n)then
                n=n-1
                goto 203
            endif
!       move particle list 'PYJETS' one step downward from i1+1 to n
            do j=i1+1,n
                j1=j-1
                do jj=1,5
                    k(j1,jj)=k(j,jj)
                    p(j1,jj)=p(j,jj)
                    v(j1,jj)=v(j,jj)
                enddo
            enddo
            n=n-1
            goto 201
202     enddo
203     continue
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sgam(nn,egam,iprt)
!!      Prints particle list "sgam" and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        dimension peo(4)

        ! in parini.f90
        call psum(pgam,1,ngam,peo)
        egam=peo(4)
        ch1=0.
        do i1=1,nn
            kf=kgam(i1,2)
            if( kf == 22 .or. kf == 44 .or. kf == 55 .or. kf == 66 &
             .or. kf == 77 ) goto 100
            ch1=ch1+pychge(kf)
100     enddo
        if( (nout == 1 .or. iii == 1 .or. mod(iii,nout) == 0 .or. iii &
          == neve) .and. iprt == 8 )then
            write(22,*)'c & p sum=                          ',ch1/3.,peo   !
            do i=1,nn
                write(22,*)i,kgam(i,1),kgam(i,2),(pgam(i,j),j=1,4)
            enddo
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_aaff(nn)
!!     print particle list and sum of momentum and energy
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        dimension peo(4)
!       do i=1,nn
!           write(22,*)i,kaff(i,2),(paff(i,j),j=1,4)
!       enddo
        call psum(paff,1,naff,peo)
        ich1=0.
        do i1=1,naff
            kf=kaff(i1,2)
            ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3D0
        write(22,*) 'aaff nn=',nn
        write(22,*) 'c & p sum=',cc,peo   !
        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_sbe   ! 030920
!!      Breaks up diquark (anti-diquark), gives four momenta
!!       and four positions to the broken partons
!       Note the names in "sbe".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_DIQUARK
        PARAMETER (KSZJ=80000)
        COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        N0 = N

        do i = 1, N0, 1
            KS = K(i,1)
            KF = K(i,2)
            if( .NOT.IS_DIQUARK(KF) .OR. .NOT.IS_EXIST(KS,i_mode) ) cycle
!       Converts KF code of diquark to correspoding ones of quarks.
            KF1 =  MOD( KF/1000, 10 )
            KF2 =  MOD( KF/100,  10 )
!       Gives properties to the broken quarks.
            K( i, 2 )   = KF1
            K( N+1, 1 ) = 1
            K( N+1, 2 ) = KF2
            K( N+1, 3 ) = i
            K( N+1, 4 ) = 0
            K( N+1, 5 ) = 0
!       Gives four momentum to the broken quarks.
            call bream_sbe( i, KF1, KF2 )
!       Gives four coordinate to the broken quarks.
            call coord_sbe(i)
            N = N + 1
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream_sbe( ii, kf1, kf2 )
!!      Give four momentum to the broken quarks.
!!      ii: line number of diquark in "sbe".
!!      kf1, kf2: flavor codes of broken quarks
!       Note the names in "sbe".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        dimension pi(4), ps(4), pp(20,5)
        integer decsuc


        am1 = amass(kf1)
        am2 = amass(kf2)
!       If zero mass diquark from PYTHIA.
        if( P(ii,5) <= 1D-15 )then
            am1 = 0D0
            am2 = 0D0
        end if
!       pp : four momenta & mass of broken quarks, local variable
        pp(1,5) = am1
        pp(2,5) = am2
!       ps : four momentum of diquark, local variable
        do i1=1,4,1
            ps(i1) = P(ii,i1)
        end do

!       Activates it for "decay method".
        goto 400

!       Broken quarks share out diquark three momentum randomly,
!        denoted as 'random three momentum method'
401     do i1=1,3,1
            pi(i1)   = PYR(1) * P(ii,i1)
            pp(2,i1) = ps(i1) - pi(i1)
            pp(1,i1) = pi(i1)
        end do
        pp11 = pp(1,1)
        pp12 = pp(1,2)
        pp13 = pp(1,3)
        pp14 = am1**2 + pp11**2 + pp12**2 + pp13**2
        if( pp14 <= 0D0 ) pp14 = 1D-20
        pp(1,4) = SQRT(pp14)
        pp21 = pp(2,1)
        pp22 = pp(2,2)
        pp23 = pp(2,3)
        pp24 = am2**2 + pp21**2 + pp22**2 + pp23**2
        if( pp24 <= 0D0 ) pp24 = 1D-20
        pp(2,4) = SQRT(pp24)

!       Activates it for 'random three momentum method'.
        goto 300

!       For 'decay method'.
400     continue
!       Decay method.
        decsuc = 1
        call decmom( ps, pp, am1, am2, decsuc )
        if( decsuc == 0 ) goto 401

300     continue
!       Adjusts four momentum conservation by iteration, no more than
!        4000 iterations
!       call conser(2,pp,ps)
        do i1=1,4,1
            p(ii,i1) = pp(1,i1)
        end do
        P(ii,5) = am1
        do i1=1,4,1
            P( N+1, i1 ) = pp(2,i1)
        enddo
        P( N+1, 5 ) = am2

!       Recalculates energies to ensure on-shell then collects the lost
!        4-momentum.
        P( ii,  4 ) = SQRT(P(ii,5)**2  + P(ii,1)**2  +P(ii,2)**2  + P(ii,3)**2)
        P( N+1, 4 ) = SQRT(P(N+1,5)**2 + P(N+1,1)**2 +P(N+1,2)**2 + P(N+1,3)**2)
        do i2=1,4,1
            throe_p(i2) = throe_p(i2) + ( ps(i2) - P(ii,i2) - P( N+1, i2 ) )
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord_sbe(ii)
!       Gives four positions to broken quarks.
!       First broken quark takes the four position of diquark.
!       Second broken quark is arranged around first ones within
!        0.5 fm randomly in each of three position coordinates and has
!        same fourth position coordinate as diquark.
!       Note the names in "sbe".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/sbe/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension rr(3)


        do i1=1,3,1
            rr(i1) = PYR(1) * 0.5D0
            V( N+1, i1 ) = V(ii,i1) + rr(i1)
            if( PYR(1) > 0.5D0 ) V( N+1, i1 ) = V(ii,i1) - rr(i1)
        end do
        V( N+1, 4 ) = V(ii,4)


        return
        end



!sa221123ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ratio_prob(x)
!!      calculates the ratio probability (probability ratio) of
!!       u(ubar):d(dbar):s(sbar):c(cbar):b(bbar):t(tbar) by the law of
!!       {m_q}^{-x}
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa38/ prob_ratio_q(6), am(6), amqq(6)
        dimension qmas(6)
        data qmas/0.33,0.33,0.50,1.50,4.80,175.0/


!       u,d,s,c,b,t kinematical mass
        dratio=qmas(1)**(-x)
        uratio=qmas(2)**(-x)
        sratio=qmas(3)**(-x)
        cratio=qmas(4)**(-x)
        bratio=qmas(5)**(-x)
        tratio=qmas(6)**(-x)
!       normalized by 'dratio'
        prob_ratio_q(1)=1.
        prob_ratio_q(2)=1.
        prob_ratio_q(3)=sratio/dratio
        prob_ratio_q(4)=cratio/dratio
        prob_ratio_q(5)=bratio/dratio
        prob_ratio_q(6)=tratio/dratio


        return
        end
!sa221123ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shanul(x,q2,rag,raq)
!!      calculate nuclear ratio R^A_i=f_{i/A}(x,Q2)/f_i(x,Q2) according
!!       to Xin-Nian Wang's paper (PLB 527 (2002) 85), multiply it to
!!       the parton distribution function in PYTHIA, resulted parton
!!       distribution function is including nuclear shadowing effect
!!      it was proved in Eur. Phys. J. C9(1999)61 that nuclear ratio does
!!       not depend strongly on the choice for the parton distribution
!!       function in nucleon f_i(x,Q2)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio


        q2_dummy = q2
        sq=0.1
        sg=0.26
        x35=x**0.35
        xsqr=sqrt(x)
        x2=x*x
        x3=x2*x
        xx=x3-1.2*x2+0.21*x
        a=nat
!       what is the definition of "a" for asymmetry reaction system ?
        a13=a**0.3333
        aa=(a13-1)**0.6
        coa=log(a)
        coa16=coa**0.16666
        bbq=1.-3.5*xsqr
        bbg=1.-1.5*x35
        eq=exp(-x2/0.01)
        eg=exp(-x2/0.004)
!       raq=a*(1.+1.19*coa16*xx-sq*aa*bbq*eq)
!       rag=a*(1.+1.19*coa16*xx-sg*aa*bbg*eg)
        raq=1.+1.19*coa16*xx-sq*aa*bbq*eq
        rag=1.+1.19*coa16*xx-sg*aa*bbg*eg


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!*******************************************************************************
!> This code calculates the nuclear overlap functions which are needed
!! to scale pp to pA and AB.
!! Units are fm
!! TAB is overlap function in 1/fm^2; divide the value by 10 to get it in mb^-1
!! scal1 is sigma(A1,A2)/sigma(p,p)
!! scal2 is n(A1,A2)/sigma(p,p)
!! The calculation is done in steps of 0.1 fm.
!! The cutoff R and b are 10 and 20 fm, respectively.
!! D.Miskowiec 1997, updated in 2001.
!*******************************************************************************
        subroutine xOVERLAP(A1,A2,rnp,rnt,sigma_NN,kjp23,denflag,density0,nshot)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        integer denflag,A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)


        rmax=10.0
!       calculate thickness functions for A1 and A2
        do i=1,200
        !longhy:change "i" type, b range is [0.05,19.95 ]
        bb(i)=dble(i)/10.0-0.05
        if(kjp23 == 2)then
        TA1(i)=TA(A1,bb(i),denflag,density0,nshot)
        TA2(i)=TA(A2,bb(i),denflag,density0,nshot)
        endif
        enddo
!       Calculate overlap function for A1+A2 collision using densities
        if(kjp23 == 2)then
        do i=1,200
        bbb=bb(i)
        TA1A2(i)=TAB(A1,A2,bbb,denflag,density0,nshot)
        enddo
        endif
!       Calculate Npart and Nbin using thickness functions
        do i=1,200
        bbb=bb(i)
        if(kjp23 == 2)then
        call PART(A1,A2,bbb,sigma_NN,nshot, &
         part1(i),part2(i))
        binn(i)=(sigma_NN/10.)*TA1A2(i)
        endif
        if(kjp23 == 1)then
        call irpt1(bbb,1,a1,a2,rou,rnp,rnt,1,0.d0,pir,tir,1.d0, &
         1.d0,1.d0,1.d0)
        part1(i)=pir
        part2(i)=tir
        endif
        enddo
!       Print data
!-----------
        write(19,*) "#! Results from the independent " // &
                    "Optical Glauber model"
        if(kjp23 == 2) write(19,901) 'i','b','TA','TB','TAB','Apart', &
                                     'Bpart','Nbin'
        if(kjp23 == 1) write(19,902) 'i','b','Apart','Bpart'
        do i=1,200
        if(kjp23 == 2) write(19,905) i,bb(i),TA1(i),TA2(i),TA1A2(i), &
                                     part1(i),part2(i),binn(i)
        if(kjp23 == 1) write(19,906) i,bb(i),part1(i),part2(i)
        enddo
 901    format(2a6,6a10)
 902    format(2a6,2a10)
 905    format(i5,2x,f6.2,6(f10.3))
 906    format(i5,2x,f6.2,2(f10.3))


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        function TA(A,b,denflag,density0,nshot)
!!      one dimension trapezoid integral
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        integer denflag,A
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)


        TA=0.0
!       do j=1,nshot
!         z=pyr(1)*rmax
!         TA=TA+xDENSITY(A,b,0.0d0,z,denflag,density0)
!       enddo
!       2: left and right symmetry, sa
!       ta=TA/nshot*2*rmax
        nshot2=2*nshot
        delta=rmax/nshot
        do j=0,nshot2
        if(j <= nshot)then
        z=-dfloat(nshot-j)*delta
        else
        z=dfloat(j-nshot)*delta
        endif
        density1=xDENSITY(A,b,0.0d0,z,denflag,density0)
        if(j == 0 .or. j == nshot2)density1=0.5*density1
        TA=TA+density1
        enddo
        TA=TA*delta


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> direct calculation of TAB (not via TA*TB)
!!      four dimensions trapezoid integral
        function TAB(A1,A2,b,denflag,density0,nshot)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        integer denflag,a1,a2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)


        TAB=0.0
!       do j=1,nshot
!         x=(pyr(1)-0.5)*2*rmax
!         y=(pyr(1)-0.5)*2*rmax
!         za=(pyr(1)-0.5)*2*rmax
!         zb=(pyr(1)-0.5)*2*rmax
        nshot2=2*nshot
        delta=rmax/nshot
        do 100 i=0,nshot2
        if(i <= nshot)then
        x=-dfloat(nshot-i)*delta
        else
        x=dfloat(i-nshot)*delta
        endif
        taby=0.0
        do 200 j=0,nshot2
        if(j <= nshot)then
        y=-dfloat(nshot-j)*delta
        else
        y=dfloat(j-nshot)*delta
        endif
        tabz=0.0
        do 300 k=0,nshot2
        if(k <= nshot)then
        z=-dfloat(nshot-k)*delta
        else
        z=dfloat(k-nshot)*delta
        endif
        tabz1=0.0
        do 400 k1=0,nshot2
        if(k1 <= nshot)then
        z1=-dfloat(nshot-k1)*delta
        else
        z1=dfloat(k1-nshot)*delta
        endif
        o=xDENSITY(A1,x+b,y,z,denflag,density0)* &
         xDENSITY(A2,x,y,z1,denflag,density0)
!       TAB=TAB+o
!       enddo
!       TAB=TAB/nshot*(2*rmax)**4
        if(k1 == 0 .or. k1 == nshot2)o=0.5*o
        tabz1=tabz1+o
400     enddo
        if(k == 0 .or. k == nshot2)tabz1=0.5*tabz1
        tabz=tabz+tabz1
300     enddo
        if(j == 0 .or. j == nshot2)tabz=0.5*tabz
        taby=taby+tabz
200     enddo
        if(i == 0 .or. i == nshot2)taby=0.5*taby
        tab=tab+taby
100     enddo
        tab=tab*delta**4


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> denflag=1 - sharp sphere
!! denflag=2 - Woods-Saxon
        function xDENSITY(A,x,y,z,denflag,density0)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        integer denflag,A


        r=sqrt(x*x+y*y+z*z)
        if (denflag == 1) then
        RA=(A/4.0D0*3.0/3.1415/density0)**0.3333333
        xDENSITY=0.0
        if (r <= RA) xDENSITY=density0
        elseif (denflag == 2) then
        ! 020511 070613 recovered
        RA=1.12*DBLE(A)**0.333333-0.86/DBLE(A)**0.333333
!       RA=1.19*A**0.333333-1.61/A**0.333333
!       density0=3./4.*A/3.1416/RA**3/(1+3.1416**2*0.54**2/RA**2)
        DR=0.54
        xDENSITY=density0/(1+exp((r-RA)/dr))
        else
        write(*,*)'wrong xDENSITY profile flag'
        stop
        endif


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!> npart calculation via TA and TB
!! part1 and part2 are the respective numbers of participants from A1 and A2
!! sigma_NN is in millibarn and TA1 and TA2 are in fm^-2
!! sigma_NN*TA/10 is mean number of NN collisions.
!! Vector b points from A to B.
!!      two dimensions trapezoid integral
        subroutine PART(A1,A2,b,sigma_NN,nshot,art1,art2)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        integer A1,A2
        common/sa31/rmax,bb(200),TA1(200),TA2(200),TA1A2(200), &
         part1(200),part2(200),binn(200)


        art1=0.0
        art2=0.0
!       do j=1,nshot
!       x=(pyr(1)-0.5)*2*rmax
!       y=(pyr(1)-0.5)*2*rmax
        nshot2=2*nshot
        delta=rmax/nshot
        do 100 i=0,nshot2
        if(i <= nshot)then
        x=-dfloat(nshot-i)*delta
        else
        x=dfloat(i-nshot)*delta
        endif
        art1y=0.0
        art2y=0.0
        do 200 j=0,nshot2
        if(j <= nshot)then
        y=-dfloat(nshot-j)*delta
        else
        y=dfloat(j-nshot)*delta
        endif
        b1=sqrt(x**2+y**2)
        b2=sqrt((x-b)**2+y**2)
        ib1=int(b1*10+1.0)
        ib2=int(b2*10+1.0)
        ib1=min(ib1,200)
        ib2=min(ib2,200)
        if(j == 0 .or. j == nshot2)then
        art1y=art1y+0.5*TA1(ib1)*(1-exp(-sigma_NN*TA2(ib2)/10))
        art2y=art2y+0.5*TA2(ib2)*(1-exp(-sigma_NN*TA1(ib1)/10))
        else
        art1y=art1y+TA1(ib1)*(1-exp(-sigma_NN*TA2(ib2)/10))
        art2y=art2y+TA2(ib2)*(1-exp(-sigma_NN*TA1(ib1)/10))
        endif
200     enddo
        if(i == 0 .or. i == nshot2)then
        art1y=0.5*art1y
        art2y=0.5*art2y
        endif
        art1=art1+art1y
        art2=art2+art2y
100     enddo
        A1_dummy=A1
        A2_dummy=A2
!       art1=art1+TA1*(1-exp(-sigma_NN*TA2/10))
!       art2=art2+TA2*(1-exp(-sigma_NN*TA1/10))
!       art1=art1+TA1(ib1)*(1-(1-sigma_NN*TA2(ib2)/10/A2)**A2)
!       art2=art2+TA2(ib2)*(1-(1-sigma_NN*TA1(ib1)/10/A1)**A1)
!       enddo
!       art1=art1/nshot*(2*rmax)**2
!       art2=art2/nshot*(2*rmax)**2
        art1=art1*delta**2
        art2=art2*delta**2


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IRPT1(BB1,IB1,IAP1,IAT1,ROU1,RP1,RT1,ISETA1,CSETA1, &
       PIR1,TIR1,AP1,BP1,AT1,BT1)
!!    'IRPT' IS A PROGRAM FOR COMPUTING THREE DIMENSIONS INTEGRAL BY GAUSSIAN
!!     TYPE QUADRATURE RULE.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      REAL*8 BB,ROU,RP,RT,CSETA,PIR,TIR,AP,BP,AT,BT
      BB=BB1
      IB=IB1
      IAP=IAP1
      IAT=IAT1
      ROU=ROU1
      RP=RP1
      RT=RT1
      ISETA=ISETA1
      CSETA=CSETA1
      PIR=PIR1
      TIR=TIR1
      AP=AP1
      BP=BP1
      AT=AT1
      BT=BT1
      CALL IRPT(BB,IB,IAP,IAT,ROU,RP,RT,ISETA,CSETA,PIR,TIR, &
      AP,BP,AT,BT)
      PIR1=PIR
      TIR1=TIR
      RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE IRPT(BB,IB,IAP,IAT,ROU,RP,RT,ISETA,CSETA,PIR,TIR,AP,BP,AT,BT)
!!    'IRPT' IS A PROGRAM FOR COMPUTING THREE DIMENSIONS INTEGRAL BY GAUSSIAN
!!     TYPE QUADRATURE RULE.
!!    HERE IT'S USED TO CALCULATE THE INTERSECT REGION OF PROJECTILE AND TARGET.
!!    IT IS WRITTEN BY SA BEN-HAO ON NOV. 6,1988.
!!    IAP AND IAT = MASS NUMBER OF PROJECTILE AND TARGET .
!!    RP AND RT = RADII OF PROJECTILE AND TARGET (IN SPHERE CASE).
!!    AP,BP AND AT,BT = THE LENGTH OF SINGLE, DOUBLE AXES OF PROJECTILE
!!     AND TARGET (IN SPHEROID CASE).
!!    AP=BP IF PROJECTILE IS A SPHERE,THE SAME FOR TARGET.
!!    IB = NUMBER OF IMPACT PARAMETERS NEEDED TO CALCULATE IN SINGLE RUN.
!!    BB = CURRENT IMPACT PARAMETER
!!    ISETA = NUMBER OF SETA (THE COSIN OF ANGLE OF SINGLE AXIS OF TARGET
!!     RELATIVE TO THE INCIDENT AXIS) NEEDED TO CALCULATE IN SINGLE RUN.
!!     FOR THE CASE OF SPHERE TARGET ISETA SHOULE BE EQUAL TO 1.
!!    CSETA = 0. FOR THE CASE OF PROJECTILE AND TARGET ARE SPHERES.
!!    CSETA = 1. FOR THE CASE OF TARGET IS SPHEROID BUT PROJECTILE IS SPHERE.
!!    CSETA = 2. FOR THE CASE OF TARGET AND PROJECTILE ARE SPHEROIDS.
!!    IRP OR PIR = THE NUMBER OF NUCLEONS IN INTERSECT REGION FROM PROJECTILE.
!!    IRT OR TIR = THE NUMBER OF NUCLEONS IN INTERSECT REGION FROM TARGET.
!!    IRPT = IRP*IRT.
!!    APIR OR AIRP = THE AVERAGE IRP OVER SETA.
!!    ATIR OR AIRT = THE AVERAGE IRT OVER SETA.
!!    RLARGE = THE LARGER ONE BETWEEN AP AND AT.
!!    D1,D2 AND D3 = LOWER LIMITS IN X,Y AND Z DIRECTIONS.
!!    U1,U2 AND U3 = UPPER LIMITS IN X,Y AND Z DIRECTIONS.
!!    N1,N2,N3 = THE NUMBER OF NODES IN X,Y AND Z DIRECTIONS.
!!    X11,X22 AND X33(200) = ARRAYS OF NODES IN X,Y AND Z DIRECTIONS.
!!    X1,X2 AND X3(200) = ARRAYS OF TRANSFORMED NODES IN X,Y AND Z DIRECTIONS.
!!    W1,W2 AND W3(200) = ARRAYS OF WEIGHTS IN X,Y AND Z DIRECTIONS.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/COM1/SINS,COSS
      COMMON/AA18/DNI(10),DPI(10),EDI(10),BMIN,BMAX
      DIMENSION ENDPTS(2),BSCR(200)
      DIMENSION X1(200),X2(200),X3(200),SUM(2),W1(200),W2(200),W3(200)
      DIMENSION X11(200),X22(200),X33(200)
      BP_dummy = BP
      N1=100
      N2=N1
      N3=N1
      JS=0
      IIB=0
      IF(ISETA == 1)GOTO 210
      RLARGE=AP
      IF(AT > AP)RLARGE=AT
      GOTO 211
210   RLARGE=RP
      IF(RT > RP)RLARGE=RT
211   D1=-RLARGE
      D2=D1
      D3=D1
      U1=RLARGE
      U2=U1
      U3=U1
!      COSS=SETA
!      SINS=SQRT(1-SETA*SETA)
!     THE SYMMETRY AROUND 3.1416/2 IS TAKEN INTO ACCOUNT.
      SSS=1./FLOAT(ISETA)
      BBB=(BMAX-BMIN)/FLOAT(IB)
      CALL GAUSSQ(1,N1,0D0,0D0,0,ENDPTS,BSCR,X11,W1)
      CALL GAUSSQ(1,N2,0D0,0D0,0,ENDPTS,BSCR,X22,W2)
      CALL GAUSSQ(1,N3,0D0,0D0,0,ENDPTS,BSCR,X33,W3)
!      WRITE(19,*)'B=',BB
!      WRITE(19,*)'AP,BP,AT,BT=',AP,BP,AT,BT
!      WRITE(19,*)'BMIN,BMAX,IB,ISETA,CSETA=',BMIN,BMAX,IB,ISETA,CSETA
!      WRITE(19,*)'N1,N2,N3=',N1,N2,N3
!      WRITE(19,*)'D1,D2,D3=',D1,D2,D3
!      WRITE(19,*)'U1,U2,U3=',U1,U2,U3
!      WRITE(19,110)(X11(I),I=1,N1)
!      WRITE(19,111)(W1(I),I=1,N1)
! 110  FORMAT(1X,'X1=',10E12.4/)
! 111  FORMAT(1X,'W1=',10E12.4/)
!     CALCULATE THE TRANSFORMED NODES.
!        F((U+D)/2+(U-D)/2*X).
!     X IS NODE.
      DO I=1,N1
      X1(I)=0.5*(U1+D1+(U1-D1)*X11(I))
      END DO
      DO I=1,N2
      X2(I)=0.5*(U2+D2+(U2-D2)*X22(I))
      END DO
      DO I=1,N3
      X3(I)=0.5*(U3+D3+(U3-D3)*X33(I))
      END DO
!     CALCULATE THE SUMS (INTEGRALS).
!       N
!      SUM  W  * F(X ) .
!      I=1   I      I
      IF(IB == 1)GOTO 201
      IIB=1
202   BB=BBB*IIB+BMIN
201   PIR=0.
      TIR=0.
      IF(ISETA == 1)GOTO 203
      JS=1
208   SETA=SSS*JS
      COSS=SETA
      SINS=SQRT(1-SETA*SETA)
203   ROUP=3./4./3.1416D0*IAP/(RP*RP*RP)
      ROUT=3./4./3.1416D0*IAT/(RT*RT*RT)
      DO 100 II=1,2
      ROU=ROUP
      IF(II == 2)ROU=ROUT
      SUM(II)=0.
      DO I1=1,N1
      XX1=X1(I1)
      SUM2=0.
      DO I2=1,N2
      XX2=X2(I2)
      SUM3=0.
      DO I3=1,N3
      XX3=X3(I3)
      SUM3=SUM3+W3(I3)* &
      FUNC(XX1,XX2,XX3,II,RP,RT,BB,AT,BT,CSETA)
      END DO
      SUM2=SUM2+W2(I2)*SUM3
      END DO
      SUM(II)=SUM(II)+W1(I1)*SUM2
      END DO
      SUM(II)=SUM(II)*(U1-D1)/2.*(U2-D2)/2.*(U3-D3)/2.*ROU
100   CONTINUE
      SUMM1=SUM(1)
      SUMM2=SUM(2)
      PIR=PIR+SUMM1
      TIR=TIR+SUMM2
      IF(ISETA == 1)GOTO 204
      IF(JS >= ISETA)GOTO 209
      JS=JS+1
      GOTO 208
209   APIR=PIR/ISETA
      ATIR=TIR/ISETA
      APTIR=APIR*ATIR
!      WRITE(19,205)
!205   FORMAT(/12X,'B',13X,'AIRP',13X,'AIRT',15X,'AIRPT')
!      WRITE(19,206)BB,APIR,ATIR,APTIR
!206   FORMAT(/3X,E12.4,4X,E12.4,4X,E12.4,10X,E12.4)
204   IF(IB == 1)GOTO 207
      IF(IIB >= IB)GOTO 207
      IIB=IIB+1
      GOTO 202
207   RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!>    SUBROUTINE FOR THE INTEGRANT.
      FUNCTION FUNC(X1,X2,X3,II,RP,RT,BB,AT,BT,CSETA)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/COM1/SINS,COSS
      X12=X1*X1
      X22=X2*X2
      X32=X3*X3
      BY=BB-X2
      BY2=BY*BY
      IF(INT(CSETA) /= 0)GOTO 400
      IF(II == 2)GOTO 200
!     CALCULATE THE VOLUME OF INTERSECTION FOR PROJECTILE
      ARG1=SQRT(X12+BY2+X32)
      F1=0.
      IF(RP >= ARG1)F1=1.
      ARG2=SQRT(X12+X22)
      F2=0.
      IF(RT >= ARG2)F2=1.
      GOTO 300
!     CALCULATE THE VOLUME OF INTERSECTION FOR TARGET
200   ARG1=SQRT(X12+X22+X32)
      F1=0.
      IF(RT >= ARG1)F1=1.
      ARG2=SQRT(X12+BY2)
      F2=0.
      IF(RP >= ARG2)F2=1.
300   FUNC=F1*F2
      GOTO 500
400   AT2=AT*AT
      BT2=BT*BT
      IF(II == 2)GOTO 600
      ARG1=SQRT(X12+BY2+X32)
      F1=0.
      IF(RP >= ARG1)F1=1.
      ARG2=SQRT(X12+X22)
      RR=AT*SINS
      RR=AMAX1(BT,RR)
      F2=0.
      IF(RR >= ARG2)F2=1.
      GOTO 700
600   YS=X2*COSS-X3*SINS
      ZS=X2*SINS+X3*COSS
      ARG1=(X12+YS*YS)/BT2+ZS*ZS/AT2
      F1=0.
      IF(1 >= ARG1)F1=1.
      ARG2=SQRT(X12+BY2)
      F2=0.
      IF(RP >= ARG2)F2=1.
700   FUNC=F1*F2
500   RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!>
!!       TITLE:  GAUSSQ
!!
!!                               1/20/75
!!
!!          THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS
!!       W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
!!       NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
!!
!!                INTEGRAL (FROM A TO B)  F(X) W(X) DX
!!
!!                             N
!!       BY                   SUM W  F(T )
!!                            J=1  J    J
!!
!!       (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)
!!       HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
!!       FUNCTIONS (LISTED BELOW), AND F(X) IS THE
!!       FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
!!       USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
!!       FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
!!
!!          ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
!!       ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES
!!       OF THE PROPER N-TH DEGREE POLYNOMIAL.
!!
!!    INPUT PARAMETERS (ALL REAL NUMBERS ARE IN DOUBLE PRECISION)
!!
!!       KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
!!                QUADRATURE RULE:
!!
!!       KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
!!       KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND
!!                  W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
!!       KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND
!!                  W(X) = SQRT(1 - X*X) ON (-1, 1)
!!       KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
!!                  (-INFINITY, +INFINITY)
!!       KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
!!                  BETA ON (-1, 1), ALPHA, BETA  >  -1.
!!                  NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
!!       KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
!!                  X**ALPHA ON (0, +INFINITY), ALPHA  >  -1
!!
!!       N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
!!       ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
!!       BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
!!                (OTHERWISE USE 0.D0)
!!       KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
!!                POINT (OR BOTH OF THE INTERVAL IS REQUIRED TO BE A
!!                NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
!!                QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
!!                ENDPOINTS (1 OR 2).
!!       ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
!!                ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
!!       B        REAL SCRATCH ARRAY OF LENGTH N
!!
!!    OUTPUT PARAMETERS (BOTH DOUBLE PRECISION ARRAYS OF LENGTH N)
!!
!!       T        WILL CONTAIN THE DESIRED NODES.
!!       W        WILL CONTAIN THE DESIRED WEIGHTS W(J).
!!
!!    SUBROUTINES REQUIRED
!!
!!       SLVE, CLASS, AND IMTQL2 ARE PROVIDED.  UNDERFLOW MAY SOMETIMES
!!       OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
!!       TURNED OFF.
!!
!!    ACCURACY
!!
!!       THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
!!       UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
!!       TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
!!       COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
!!       DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
!!       HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
!!       VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
!!       COMPLETELY HARMLESS.
!!
!!    METHOD
!!
!!          THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
!!       FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
!!       USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
!!       EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
!!       SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
!!       THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
!!       YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
!!       ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
!!       FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
!!       GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
!!
!!    REFERENCES
!!
!!       1.  GOLUB, G. H., AND WELSCH, J. H., "CALCULATION OF GAUSSIAN
!!           QUADRATURE RULES," MATHEMATICS OF COMPUTATION 23 (APRIL,
!!           1969), PP. 221-230.
!!       2.  GOLUB, G. H., "SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,"
!!           SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
!!       3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
!
      REAL*8 B(N), T(N), W(N), ENDPTS(2), MUZERO, T1, &
       GAM, SLVE, ALPHA, BETA
!
      CALL xCLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
!
!           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
!           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
!           B THE OFF-DIAGONAL ELEMENTS.
!           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
!           SUBMATRIX.
!
      IF (KPTS == 0)  GO TO 100
      IF (KPTS == 2)  GO TO  50
!
!           IF KPTS=1, ONLY T(N) MUST BE CHANGED
!
      T(N) = SLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
      GO TO 100
!
!           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
!
   50 GAM = SLVE(ENDPTS(1), N, T, B)
      T1 = ((ENDPTS(1) - ENDPTS(2))/(SLVE(ENDPTS(2), N, T, B) - GAM))
      B(N-1) = SQRT(T1)
      T(N) = ENDPTS(1) + GAM*T1
!
!           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
!           AND THUS THE VALUE OF B(N) IS ARBITRARY.
!           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
!           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
!           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
!
  100 W(1) = 1.0E0
      DO I = 2, N
         W(I) = 0.0E0
      END DO
!
      CALL IMTQL2 (N, T, B, W, IERR)
      DO I = 1, N
         W(I) = MUZERO * W(I) * W(I)
      END DO
!
      RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SLVE(SHIFT, N, A, B)
!!
!!      THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
!!      N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
!!
!!            (JN - SHIFT*IDENTITY) * DELTA  = EN,
!!
!!      WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
!!      THE N-TH POSITION.
!!
!!      THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
!!      ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
!!      MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
!!      2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
!!
!!
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      REAL*8 A(N), B(N), ALPHA
!
      ALPHA = A(1) - SHIFT
      NM1 = N - 1
      DO I = 2, NM1
         ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
      END DO
      SLVE = 1.0E0/ALPHA
      RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE xCLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
!!
!!          THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
!!       RECURRENCE RELATION
!!
!!            B P (X) = (X - A ) P   (X) - B   P   (X)
!!             J J            J   J-1       J-1 J-2
!!
!!       FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
!!       AND THE ZERO-TH MOMENT
!!
!!            MUZERO = INTEGRAL W(X) DX
!!
!!       OF THE GIVEN POLYNOMIAL'S WEIGHT FUNCTION W(X).  SINCE THE
!!       POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
!!       GUARANTEED TO BE SYMMETRIC.
!!
!!          THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
!!       JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
!!       JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
!!       REQUIRE THE xGAMMA FUNCTION.
!!
!!    ..................................................................
!!
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      REAL*8 A(N), B(N), MUZERO, ALPHA, BETA
      REAL*8 ABI, A2B2, xGAMMA, PI, AB
      DATA PI / 3.141592653589793D0 /
!
      NM1 = N - 1
      GO TO (10, 20, 30, 40, 50, 60), KIND
!
!              KIND = 1:  LEGENDRE POLYNOMIALS P(X)
!              ON (-1, +1), W(X) = 1.
!
   10 MUZERO = 2.0E0
      DO I = 1, NM1
         A(I) = 0.0E0
         ABI = I
         B(I) = ABI/SQRT(4*ABI*ABI - 1.0E0)
      END DO
      A(N) = 0.0E0
      RETURN
!
!              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
!              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
!
   20 MUZERO = PI
      DO I = 1, NM1
         A(I) = 0.0E0
         B(I) = 0.5E0
      END DO
      B(1) = SQRT(0.5E0)
      A(N) = 0.0E0
      RETURN
!
!              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
!              ON (-1, +1), W(X) = SQRT(1 - X*X)
!
   30 MUZERO = PI/2.0E0
      DO I = 1, NM1
         A(I) = 0.0E0
         B(I) = 0.5E0
      END DO
      A(N) = 0.0E0
      RETURN
!
!              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
!              +INFINITY), W(X) = EXP(-X**2)
!
   40 MUZERO = SQRT(PI)
      DO I = 1, NM1
         A(I) = 0.0E0
         B(I) = SQRT(I/2.0D0)
      END DO
      A(N) = 0.0E0
      RETURN
!
!              KIND = 5:  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
!              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
!              BETA GREATER THAN -1
!
   50 AB = ALPHA + BETA
      ABI = 2.0E0 + AB
      MUZERO = 2.0E0 ** (AB + 1.0E0) * xGAMMA(ALPHA + 1.0E0) * xGAMMA( &
       BETA + 1.0E0) / xGAMMA(ABI)
      A(1) = (BETA - ALPHA)/ABI
      B(1) = SQRT(4.0E0*(1.0E0 + ALPHA)*(1.0E0 + BETA)/((ABI + 1.0E0)* &
        ABI*ABI))
      A2B2 = BETA*BETA - ALPHA*ALPHA
      DO I = 2, NM1
         ABI = 2.0D0*I + AB
         A(I) = A2B2/((ABI - 2.0E0)*ABI)
         B(I) = SQRT (4.0D0*I*(I + ALPHA)*(I + BETA)*(I + AB)/ &
         ((ABI*ABI - 1)*ABI*ABI))
      END DO
      ABI = 2.0D0*N + AB
      A(N) = A2B2/((ABI - 2.0E0)*ABI)
      RETURN
!
!              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
!              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
!              THAN -1.
!
   60 MUZERO = xGAMMA(ALPHA + 1.0E0)
      DO I = 1, NM1
         A(I) = 2.0D0*I - 1.0E0 + ALPHA
         B(I) = SQRT(I*(I + ALPHA))
      END DO
      A(N) = 2.0D0*N - 1 + ALPHA
      RETURN
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!>    ------------------------------------------------------------------
!!
!!    TITLE:  IMTQL2
!!
!!
!!    THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
!!    NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
!!    AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
!!    HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
!!    THIS IS A MODIFIED VERSION OF THE 'EISPACK' ROUTINE IMTQL2.
!!
!!    THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
!!    EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
!!    METHOD.
!!
!!    ON INPUT:
!!
!!       N IS THE ORDER OF THE MATRIX;
!!
!!       D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
!!
!!       E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!!         IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
!!
!!       Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
!!
!!     ON OUTPUT:
!!
!!       D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!!         ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!!         UNORDERED FOR INDICES 1, 2, ..., IERR-1;
!!
!!       E HAS BEEN DESTROYED;
!!
!!       Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
!!         OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
!!         MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!!         EIGENVALUES;
!!
!!       IERR IS SET TO
!!         ZERO       FOR NORMAL RETURN,
!!         J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!!                    DETERMINED AFTER 30 ITERATIONS.
!!
      SUBROUTINE IMTQL2(N, D, E, Z, IERR)
!
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      REAL*8 D(N), E(N), Z(N), B, C, F, G, P, R, S, MACHEP
!
      MACHEP=2.2E-17
      IERR=0
      IF (N == 1)GOTO 1001
      E(N) = 0.0E0
      DO 240 L = 1, N
         J = 0
  105    DO 110 M = L, N
      IF(M == N)GOTO 120
      IF(ABS(E(M)) <= MACHEP*(ABS(D(M))+ABS(D(M+1))))GOTO 120
110   CONTINUE
!
120   P=D(L)
      IF(M == L)GOTO 240
      IF(J == 30)GOTO 1000
      J=J+1
!
      G=(D(L+1)-P)/(2.0E0*E(L))
      R=SQRT(G*G+1.0E0)
      G=D(M)-P+E(L)/(G+SIGN(R,G))
      S=1.0E0
      C=1.0E0
      P=0.0E0
      MML=M-L
!
      DO II=1,MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (ABS(F)  <  ABS(G)) GO TO 150
            C = G / F
            R = SQRT(C*C+1.0E0)
      E(I+1)=F*R
      S=1.0E0/R
      C=C*S
      GOTO 160
 150  S=F/G
      R=SQRT(S*S+1.0E0)
      E(I+1)=G*R
      C=1.0E0/R
      S=S*C
 160  G=D(I+1)-P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
!
            F = Z(I+1)
      Z(I+1)=S*Z(I)+C*F
      Z(I)=C*Z(I)-S*F
      END DO
!
      D(L)=D(L)-P
      E(L)=G
      E(M)=0.0E0
         GO TO 105
  240 CONTINUE
!
      DO 300 II = 2, N
      I=II-1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J)  >=  P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K  ==  I) GO TO 300
         D(K) = D(I)
         D(I) = P
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
  300 CONTINUE
!
      GO TO 1001
!                EIGENVALUE           AFTER 30 ITERATIONS
 1000 IERR = L
 1001 RETURN
!
      END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION xGAMMA(X)   ! 081010
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      C=1.0E0
      Y=X
      IF((Y-2.0E0) < 0E0) GOTO 1
      IF((Y-3.0E0) <= 0E0) GOTO 3
      IF((Y-10.0E0) > 0E0) GOTO 6
 5    Y=Y-1.0E0
      C=C*Y
      IF((Y-3.0E0) <= 0E0)THEN
         GOTO 3
      ELSE
         GOTO 5
      END IF
 1    C=Y*C
      Y=Y+1.0E0
      IF((Y-2.0E0) < 0E0) GOTO 1
      C=1.0E0/C
    3 B=Y-2.0E0
      G=(((((((-.5113262726698E-6*B+.51063592072582E-5)*B &
        -.248410053848712E-4)*B+.815530498066373E-4)*B    &
        -.2064476319159326E-3)*B+.4677678114964956E-3)*B  &
        -.9083465574200521E-3)*B+.002099759035077063E0)*B
      G=(((((((G-.002851501243034649E0)*B+.0111538196719067E0)*B &
        -.2669510287555266E-3)*B+.07424900794340127E0)*B         &
        +.08157691940138868E0)*B+.4118403304219815E0)*B          &
        +.4227843350985181E0)*B+.9999999999999999E0
      VALUE=G*C
      GO TO 7
    6 D=1.0E0/Y
      C=D*D
      G=(((((((((-1.392432216905901E0*C+.1796443723688306E0)*C &
        -.02955065359477124E0)*C+.00641025641025641E0)*C       &
        -.0019175269175269175E0)*C+.8417508417508418E-3)*C     &
        -.5952380952380952E-3)*C+.79365079365079365E-3)*C      &
        -.002777777777777778E0)*C+.08333333333333333E0)*D      &
        +.9189385332046727E0+(Y-.5E0)*LOG(Y)-Y
      VALUE=EXP(G)
    7 xGAMMA=VALUE
      RETURN
      END



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc