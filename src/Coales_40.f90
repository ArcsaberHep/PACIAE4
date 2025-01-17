!! coales_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2024 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, June 2004 - January 2025.

!> This is the program to perform the coalescence hadronization via the
!!  phenomenological coalescence model of PACIAE.

!!                                             By Ben-Hao at CIAE on 04/06/2004
!!                                  Last updated by An-Ke at UiO  on 17/01/2025


        subroutine coales( i_coal )
!!      A phenomenological coalescence model.
!
!       It was writen by Sa Ben-Hao on 04/06/2004.
!       Its input messages are in 'PYJETS'.
!       Its storing array is 'PYJETS'.
!       Its output message is in 'sa1_h' (in 'PYJETS' too).
!       iii: the run number
!
!       i_coal=1 or 0: whether to perform the gluon splitting and quark
!                      deexcitation or not (diffrent for the new coales).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sbh/nbh,nonh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         napp,natt,nzpp,nztt,pio
!       Arrays of freeze-out particle information.
        COMMON/PYJETS_FO/ N_FO, NPAD_FO, &
                          K_FO(KSZJ,8), P_FO(KSZJ,7), V_FO(KSZJ,5)



!===============================================================================
!       New method
        call coales_new( i_coal )
        return
!       The following old method is no longer used.
!===============================================================================



        iadj12=int(adj1(12))
!       Do gluon splitting and energetic quark deexcitation only, without
!        the real coalescence (i_coal=0).
        if( i_coal == 0 ) goto 888


!-------------------------------------------------------------------------------
!-------------------------   Variables Initialization   ------------------------
        rrp = 1.16
        nn  = 0
        nth = 0
        iphas = INT(adj1(21))
!-------------------------   Variables Initialization   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Junctions Removing   ---------------------------
        if( INT(adj1(40)) == 3 .or. iadj12 == 1 )then
!           Remove junctions.
            jb = 0
2010        do i1=jb+1,N,1
                kf   = K(i1,2)
                if(kf /= 88)then
                    jb = jb + 1
                    goto 2020
                end if
                ! "updad_pyj" is in sfm_30.f90.
                call updad_pyj(N,i1+1,1)
                N = N - 1
                goto 2010
2020        end do
        end if
!----------------------------   Junctions Removing   ---------------------------
!-------------------------------------------------------------------------------


!       g-splitting and q-deexcitation have been done before "parcas".
!       So jumps out them when do the real coalescence (i_coal=1).
888     if( iadj12 == 2 .AND. i_coal == 1 ) goto 1000


!-------------------------------------------------------------------------------
!-----------------------------   Gluon Splitting   -----------------------------
        n00 = N   ! Original total entries in PYJETS
!       Move gluons from 'pyjest' to 'sa36'.
        call remo_glu
!       Break-up gluon (with E_g>2E_u in 'sa36') -> qqbar string
!        (filling in 'PYJETS').
        call break_glu
!       So far, the parton list ('PYJETS') is composed of q and qbar only.
        i_deex_gen = INT( adj1(16) )   ! 180923 Lei
        adj17 = adj1(17)
!       adj17=max(4.0,adj17)

!       Shares 4-momentum in "throe_p" among partons.
        if( iadj12 == 1 ) call share_p_PYJETS
!-----------------------------   Gluon Splitting   -----------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Quark Deexcitation   ----------------------------
!       energetic q (qbar) de-excitation
        i_call_deex  = 0
        ! the #-th newly produced daughter qqbar
        i_daught_gen = 1
        ! the number of successful deexcitation
        n_deex = 0
        jb = 0
        ! Current total entries in PYJETS
        n0 = N
        if( i_deex_gen == 0 ) goto 900
700     continue
        do i1=jb+1,n0,1
            kf0   = K(i1,2)
            ee    = P(i1,4)
            iflav = 1
            if( kf0 < 0 ) iflav = -1
!           iflav = 1 : if source parton is quark
!                 =-1 : if source parton is anti-quark
            if( ee > adj17 )then
                if( i_deex == 1 ) call  deexcitation_W(  i1, nstep )
                if( i_deex == 2 ) call  deexcitation_E(  i1, nstep )
                if( i_deex == 3 ) call  deexcitation_PZ( i1, nstep )
                if( nstep > 0 ) n_deex = n_deex + 1
                ! times of 'call deexcitation'
                i_call_deex = i_call_deex + 1
            endif
!           nstep : number of deexcitation steps per source q (qbar)
!           Updates n0 and does deexcitation for newly produced qqbar pairs.
        if( i1 == n0 .AND. N > n0 .AND. i_daught_gen < i_deex_gen)then
            jb = i1
            i_daught_gen = i_daught_gen + 1
            n0 = N
            goto 700
        end if
        enddo
900     continue
!       energetic q (qbar) de-excitation, finished.

!       Shares the 4-momentum in 'throe_p' among partons.
        if( iadj12 == 1 ) call share_p_PYJETS

        if( iadj12 == 3 )then
!       Records the location (line numbers) of new strings (qqbar).
            do i1=n00+1,N,2
                nstr1 = nstr1 + 1
                nstr1a( nstr1 ) = i1
                nstr1v( nstr1 ) = i1 + 1
            end do
            nstr0 = nstr1
!       Updates "sbe".
            do ii=n00+1,N,1
                nbe = nbe + 1
                do jj=1,5,1
                    kbe( nbe, jj ) = K(ii,jj)
                    pbe( nbe, jj ) = P(ii,jj)
                    vbe( nbe, jj ) = V(ii,jj)
                end do
            end do
        end if
!---------------------------   Quark Deexcitation   ----------------------------
!-------------------------------------------------------------------------------


!       Just do the g-splitting and quark deexcitation, without real coalescence
1000    if( i_coal == 0 ) return


!-------------------------------------------------------------------------------
!-----------------------------   Parton Sorting   ------------------------------
!       Sorts out PYJETS randomly. Lower case only. (In coales.f90)
        call PASORT( 1, N, "pyjets", "null", "random" )
!-----------------------------   Parton Sorting   ------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Parton coalescence   ---------------------------
!       Load the table of hadron (meson: pseudoscalar-spin 0 & vector-spin 1
!        only, baryon: octet-spin 1/2 & decuplet-spin 3/2 only).
        if( iii == 1 ) call tabhb
!       iii is the event number.

!       Normal coalescence process.
        call hadpro(iphas)
!       ithroq : the total number of quarks thrown away
!       ithrob : the total number of anti-quarks thrown away
!       throe : total 4-momentum of the partons thrown away
!       ich : total charge of the partons thrown away

!       Re-coalesce failed parton after last 'call hadpro'.
!       Try coalescence without phase-space constraint agian for those partons
!        failed in the normal coalescence process.
        if( iphas /= 0 .AND. nth >= 2 )then
            call hadpro(0)
        endif
!----------------------------   Parton coalescence   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------------   Data dumping   -------------------------------
!       'sa1_h' to 'PYJETS'.
        N = nn
        do j2=1,5,1
            do j1=1,nn,1
                K(j1,j2) = kn(j1,j2)
                P(j1,j2) = pn(j1,j2)
                V(j1,j2) = rn(j1,j2)
            enddo
        enddo
!------------------------------   Data dumping   -------------------------------
!-------------------------------------------------------------------------------


!       Decay of unstable hadrons.
        call decayh(rrp)

!       Updates the freeze-out information in "PYJETS_FO" to be the
!        same as the ones in "PYJETS".
        do j=1,5,1
            do i=1,N,1
                K_FO(i,j) = K(i,j)
                P_FO(i,j) = P(i,j)
                V_FO(i,j) = V(i,j)
            end do
        end do
        N_FO = N


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coales_new( mode_coal )
!!      New coalescence model.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_QUARK
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


        nn  = 0
        nth = 0

!       Moves non-q/qbar (gluons, diquarks. etc.) from "PYJETS" to "aaff".
        call PASORT( 1, N, "pyjets", "|kf|", "min_to_max" )
        i_non_q_begin = N + 1
        do i=1,N,1
            KF = K(i,2)
            if( IS_QUARK(KF) ) cycle
            i_non_q_begin = i
            exit
        end do
        nsa = 0
        do i = i_non_q_begin, N, 1
            nsa = nsa + 1
            do j=1,5,1
                ksa( nsa, j ) = K(i,j)
                psa( nsa, j ) = P(i,j)
                vsa( nsa, j ) = V(i,j)
            end do
        end do
        N = i_non_q_begin - 1

!       Normal coalescence process.
        i_constraint = INT( adj1(21) )
        call do_coalescence( i_constraint, mode_coal )

!       Tries the coalescence without phase-space constraint agian for
!        partons failed in the normal coalescence process.
        if( i_constraint /= 0 .AND. nth >= 2 ) &
            call do_coalescence( 0, mode_coal )

!       "sa1_h" to "PYJETS".
        N = nn
        do j=1,5,1
            do i=1,nn,1
                K(i,j) = kn(i,j)
                P(i,j) = pn(i,j)
                V(i,j) = rn(i,j)
            end do
        end do
        nn = 0

!       Decay of unstable hadrons. (In Pythia8_fort_interface.f90)
        call PADECY

!       "sa2" to "PYJETS".
        do i = 1, nsa, 1
            N = N + 1
            do j=1,5,1
                K( N, j ) = ksa(i,j)
                P( N, j ) = psa(i,j)
                V( N, j ) = vsa(i,j)
            end do
        end do
        nsa = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_gluon_splitting
!!      Splits the gluons to qqbar in "PYJETS" and apends them to the end.
!!      According to the Altarelli-Parisi splitting function.
!!      Adjusts "sbe", the diquark list and the string list for thje specical
!!       version of SFM, too.
!       i_pT_split: transverse direction (pT sampling method)
!                   = 0, collinear.
!                   = 1, Gaussian.
!                   = 2, Exponential.
!       i_Z_split: longitudinal direction (W / pz splitting method).
!                  = 0,   random 3-momentum (ignore i_pT_split).
!                  = 1,   Altarelli-Parisi splitting function with light cone
!                          momentum scheme in the rotated frame.
!                  = 2,   Lund/FF/PS functions etc. with light cone momentum.
!                  = 3,   random Z fraction with light cone momentum.
!                  = 11, Altarelli-Parisi splitting function with
!                         energy scheme in the rotated frame.
!                  = 22, Lund/FF/PS functions etc. with energy.
!                  = 33, random Z fraction with energy.
!                  = 111,  Altarelli-Parisi splitting function with logitudinal
!                          momentum scheme in the rotated frame.
!                  = 222,  Lund/FF/PS functions etc. with logitudinal momentum.
!                  = 333,  random Z fraction with logitudinal momentum.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        real(kind=8) :: mass, mT2, P00(5), V00(5), rr(3)


!       Hadronization model.
        i_had_model = INT( adj1(12) )
!       For the debug.
        if( i_had_model == 1 .AND. i_deex > 99 ) return
!       Splitting threshold energy.
        E_threshold = adj1(39)
        dm_uubar    = 2D0 * amass(2)
        if( E_threshold < dm_uubar ) E_threshold = dm_uubar
!       Move gluons from "PYJEST" to "sa36".
        call remo_glu
!       Sorts gluons by energy from min to max.
        call PASORT( 1, nglu, "sa36", "e", "min_to_max" )
!       Width of pT sampling.
        i_pT_split = i_pT_coal
        ! i_pT_split = 0
        pT_width = adj1(35)
!       The spliting method (hardcode: Altarelli-Parisi splitting function).
        i_Z_split = 111
        if( i_deex == 1 )then
            i_Z_split = 1
        else if( i_deex == 2 )then
            i_Z_split = 11
        end if
!       Function for selecting z.
        i_function_z = INT( adj1(29) )
!       Backup.
        PARJ41 = PARJ(41)
        PARJ42 = PARJ(42)
        PARJ51 = PARJ(51)
        PARJ52 = PARJ(52)
        PARJ53 = PARJ(53)
        PARJ54 = PARJ(54)
        PARJ55 = PARJ(55)
        PARJ(41) = adj1(8)
        PARJ(42) = adj1(9)
        PARJ(51) = a_FF
        PARJ(52) = a_FF
        PARJ(53) = a_FF
        PARJ(54) = -aPS_c
        PARJ(55) = -aPS_b


!-------------------------------------------------------------------------------
!-----------------------------   Gluon Splitting   -----------------------------
        N0 = N
        nglu0 = nglu
        N_SHARE = 0
!       g -> qqbar
        loop_gluon: do i = 1, nglu0, 1
            eg = pglu(i,4)
            ! Too low energy to split.
            if( eg <= E_threshold )then
                ! Shares its 4-mom to the gluon with the next lowest energy.
                do j=1,4,1
                    pglu( i+1, j ) = pglu( i+1, j ) + pglu( i, j )
                end do
                ! Throws it away directly. (Better y midrapidity)
                ! do j=1,4,1
                !     throe_p(j) = throe_p(j) + pglu( i, j )
                ! end do
                ! Collects gluon with too low energy.
                ! nglu = nglu + 1
                ! do j=1,4,1
                !     pglu( nglu, j ) = pglu( nglu, j ) + pglu( i, j )
                ! end do
                N_SHARE = N_SHARE + 1
                cycle loop_gluon
            end if
            ! Records original information.
            do j=1,5,1
                P00(j) = pglu(i,j)
                V00(j) = vglu(i,j)
            end do
!       Rotates the gluon to the frame in which px* = py* = 0,
!        pz* = |p|, and the light-cone momentum W* = E+ pz* = E+ |p|.
            pxg = pglu(i,1)
            pyg = pglu(i,2)
            pzg = pglu(i,3)
            pg  = SQRT( pxg**2 + pyg**2 + pzg**2 )
            Wg  = eg + pg
            pxq = 0D0
            pyq = 0D0
            pzq = 0D0
            pq  = 0D0
            eq  = 0D0
            Wq  = 0D0
            px_qbar = 0D0
            py_qbar = 0D0
            pz_qbar = 0D0
            p_qbar  = 0D0
            e_qbar  = 0D0
            W_qbar  = 0D0
            KF   = 0D0
            mass = 0D0

!------------------------   Splitted Flavor Sampling   -------------------------
            N_TRY_KF = 0
100         continue
            N_TRY_KF = N_TRY_KF + 1
            ! No more than 50 times.
            if( N_TRY_KF > 50 )then
                write(22,*) "Warning, N_TRY_KF > 50 in g-splitting, " &
                         // "i-event, i-entry, KF, px py pz e m =", &
                            iii, i, 21, ( P00(j), j=1, 5, 1 )
                ! Shares its 4-mom to the gluon with the next lowest energy.
                do j=1,4,1
                    pglu( i+1, j ) = pglu( i+1, j ) + pglu( i, j )
                end do
                N_SHARE = N_SHARE + 1
                cycle loop_gluon
            end if
            ! Finds flavors of the splitted qqbar.
            call break_f( eg, KF, mass )
            if( PYR(1) > 0.5D0 ) KF = -KF
            mass = amass( KF )
!------------------------   Splitted Flavor Sampling   -------------------------

!--------------------------   Splitted pT Sampling   ---------------------------
            N_TRY_PT = 0
200         continue
            N_TRY_PT = N_TRY_PT + 1
            ! No more than 50 times.
            if( N_TRY_PT > 50 ) goto 100
            ! Generates compensated px py for splitted qqbar.
            call PAPTDI( i_pT_split, pT_width, pxq, pyq, pTq )
            px_qbar = -pxq
            py_qbar = -pyq
            mT2 = mass**2 + pxq**2 + pyq**2
!--------------------------   Splitted pT Sampling   ---------------------------

!---------------------------   Splitted Z Sampling   ---------------------------
            N_TRY_Z = 0
            loop_try_pz: do while(.true.)
                N_TRY_Z = N_TRY_Z + 1
                ! No more than 50 times.
                if( N_TRY_Z > 50 ) goto 200
                ! Altarelli-Parisi splitting function.
                if( i_Z_split ==1 .OR. i_Z_split ==11 .OR. i_Z_split ==111 )then
                    do while(.true.)
                        Z = PYR(1)
                        Zg2qqbar = Z**2 + ( 1D0 - Z )**2
                        if( PYR(1) <= Zg2qqbar ) exit
                    end do
                    ! Approximation.
                    ! Z = PYR(1)**( 1D0/3D0 )
                ! Lund/FF/PS functions etc.
                else if( i_Z_split==2 .OR.i_Z_split==22 .OR.i_Z_split==222 )then
                    dmT2_qqbar = (2D0 * mass)**2
                    call PAZDIS( i_function_z, 0D0, KF, -KF, dmT2_qqbar, Z )
                else
                    ! Z = 0.5D0
                    Z = PYR(1)
                end if
                ! Light cone momentum scheme.
                if( i_Z_split == 1 .OR. i_Z_split == 2 .OR. i_Z_split == 3 )then
                    Wq = Wg * Z
                    pzq = 0.5D0*( Wq - mT2 / MAX( 1D-4, Wq ) )
                    eq  = 0.5D0*( Wq + mT2 / MAX( 1D-4, Wq ) )
                    W_qbar = Wg - Wq
                    pz_qbar = 0.5D0*( W_qbar - mT2 / MAX( 1D-4, W_qbar ) )
                    e_qbar  = 0.5D0*( W_qbar + mT2 / MAX( 1D-4, W_qbar ) )
                    ! Checks if pz acceptable.
                    if( pzq > 0D0 .AND. pz_qbar > 0D0 )then
                        exit loop_try_pz
                    end if
                ! Energy scheme.
                else if(i_Z_split==11 .OR. i_Z_split==22 .OR. i_Z_split==33)then
                    eq = eg * Z
                    pzq2 = eq**2 - mass**2 - pxq**2 - pyq**2
                    e_qbar = eg - eq
                    pz_qbar2 = e_qbar**2 - mass**2 - px_qbar**2 - py_qbar**2
                    ! Checks if pz acceptable.
                    if( pzq2 > 0D0 .AND. pz_qbar2 > 0D0 )then
                        pzq = SQRT( pzq2 )
                        pz_qbar = SQRT( pz_qbar2 )
                        exit loop_try_pz
                    end if
                ! Logitudinal momentum (pz) scheme.
                else if(i_Z_split==111.OR.i_Z_split==222.OR.i_Z_split==333)then
                    pzq = pg * Z
                    pz_qbar = pg - pzq
                    exit loop_try_pz
                ! Random 3-momentum (badly sharp midrapidity).
                else
                    ! Same fraction.
                    ! rndm = PYR(1)
                    ! pxq = pxg * rndm
                    ! pyq = pyg * rndm
                    ! pzq = pzg * rndm
                    pxq = pxg * PYR(1)
                    pyq = pyg * PYR(1)
                    pzq = pzg * PYR(1)
                    px_qbar = pxg - pxq
                    py_qbar = pyg - pyq
                    pz_qbar = pzg - pzq
                    exit loop_try_pz
                end if
            end do loop_try_pz
            ! Ensures on-shell.
            eq     = SQRT( mass**2 + pxq**2 + pyq**2 + pzq**2 )
            e_qbar = SQRT( mass**2 + px_qbar**2 + py_qbar**2 + pz_qbar**2 )
!---------------------------   Splitted Z Sampling   ---------------------------

!       Gives splitted qqbar properties.
            K( N+1, 1 ) = 2
            K( N+1, 2 ) = KF
            K( N+1, 3 ) = 0
            K( N+1, 4 ) = 0
            K( N+1, 5 ) = 0

            K( N+2, 1 ) = 1
            K( N+2, 2 ) = -KF
            K( N+2, 3 ) = 0
            K( N+2, 4 ) = 0
            K( N+2, 5 ) = 0

            P( N+1, 1 ) = pxq
            P( N+1, 2 ) = pyq
            P( N+1, 3 ) = pzq
            P( N+1, 4 ) = eq
            P( N+1, 5 ) = mass

            P( N+2, 1 ) = px_qbar
            P( N+2, 2 ) = py_qbar
            P( N+2, 3 ) = pz_qbar
            P( N+2, 4 ) = e_qbar
            P( N+2, 5 ) = mass

            ! Rotates the qqbar back to the original frame.
            if( i_Z_split /= 0 )then
                theta = PYANGL( P00(3), SQRT( P00(1)**2 + P00(2)**2 ) )
                phi   = PYANGL( P00(1), P00(2) )
                beta_x = 0D0 ! + P00(1) / P00(4)
                beta_y = 0D0 ! + P00(2) / P00(4)
                beta_z = 0D0 ! + P00(3) / P00(4)
                MSTU(33) = 1
                call PYROBO( N+1, N+2, theta, phi, beta_x, beta_y, beta_z )
            end if

!       Collects the lost 4-momentum.
            do j=1,4,1
                throe_p(j) = throe_p(j) + ( P00(j) - P(N+1,j) - P(N+2,j) )
            end do

!       Gives four-positions to the splitted qqbar.
!       Arranged around the source gluon within 0-0.5 fm randomly
!        in each one of the three positions.
            do ii = N+1, N+2, 1
                do j=1,3
                    rr(j)  = PYR(1)*0.5D0
                    V(ii,j) = V00(j) + rr(j)
                    if( PYR(1) > 0.5D0 ) V(ii,j) = V00(j) - rr(j)
                end do
                V(ii,4) = V00(4)
                V(ii,5) = 0D0

!TODO(Lei20241016): gives formation time.
                ! IF( .true. )THEN
                IF( .false. )THEN
                    i_tau = 1
                    ! Formation time method 1.
                    ! Proper formation time of t0/10 was assumed for partons.
                    if( i_tau == 1 )then
                        tau0 = t0 / 10D0
                        V(ii,4) = V00(4) + tau0 * P(ii,4) / mass
                        V(ii,5) = 0D0
                    ! Formation time method 2, tau = E/mT^2.
                    else if( i_tau == 2 )then
                        V(ii,4) = V00(4) &
                       + P(ii,4) / ( P(ii,1)**2 + P(ii,2)**2 + P(ii,5)**2 ) &
                       * PARU(3)
                    end if
                END IF

!       For the special version of SFM.
                if( i_had_model == 3 )then
                    ! Updates them to "sbe".
                    nbe = nbe + 1
                    do j=1,5,1
                        kbe( nbe, j ) = K(ii,j)
                        pbe( nbe, j ) = P(ii,j)
                        vbe( nbe, j ) = V(ii,j)
                    end do
                end if
            end do

!       For the special version of SFM.
            if( i_had_model == 3 )then
                ! Updates the sting list.
                nstr1 = nstr1 + 1
                nstr1a( nstr1 ) = nbe - 1
                nstr1v( nstr1 ) = nbe
            end if

            N = N + 2
        end do loop_gluon
        nglu = 0
!-----------------------------   Gluon Splitting   -----------------------------
!-------------------------------------------------------------------------------


!       Recovers.
        PARJ(41) = PARJ41
        PARJ(42) = PARJ42
        PARJ(51) = PARJ51
        PARJ(52) = PARJ52
        PARJ(53) = PARJ53
        PARJ(54) = PARJ54
        PARJ(55) = PARJ55


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_quark_deexcitation
!!      Performs the deexcitation of exited quarks/anti-quarks.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_QUARK
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)


!       Hadronization model.
        i_had_model = INT( adj1(12) )
!       For the debug.
        if( i_had_model == 1 .AND. i_deex > 99 ) return
!       Allowed generations of q/qbar to be deexcited.
        i_deex_generation = INT( adj1(16) )
!       No deexcitation.
        if( i_deex_generation == 0 .OR. n_deex_step == 0 ) return

!       Deexcitation scheme.
        i_deex_scheme = i_deex
!       Generation counter of the newly produced daughter qqbars.
        i_daught_generation = 1
!       Original total entries in "PYJETS".
        N00 = N
        i_call_deex = 0
!       The number of successful deexcitations.
        n_deex_succ = 0
        i_begin = 0
!       Current total entries in PYJETS.
        N0 = N
!       Threshold energy.
        E_threshold = adj1(17)
        ! Kinemtical mass of uubar ~ 2*0.33 = 0.66 GeV.
        dm_uubar = 2D0 * amass(2)
        if( E_threshold <= dm_uubar ) E_threshold = dm_uubar

!       Performs the deexcitation.
700     continue
        do i = i_begin+1, N0, 1
            KF0 = K(i,2)
            ! Only quarks temporarily.
            if( .NOT.IS_QUARK(KF0) ) cycle
!           Light-cone momentum scheme (W).
            ee = SQRT( P(i,4)**2 + P(i,1)**2 + P(i,2)**2 + P(i,3)**2 )
!           Energy scheme (E).
            if( i_deex_scheme == 2 )then
                ee = P(i,4)
!           Momentum scheme (PZ).
            else if( i_deex_scheme == 3 )then
                ee = SQRT( P(i,1)**2 + P(i,2)**2 + P(i,3)**2 )
            end if
            if( ee > E_threshold )then
                select case( i_deex_scheme )
                case(1)
                    call deexcitation_W(  i, nstep )
                case(2)
                    call deexcitation_E(  i, nstep )
                case(3)
                    call deexcitation_PZ( i, nstep )
                case default
                    call deexcitation_W(  i, nstep )
                end select
!               nstep : number of deexcitation steps actually occurred.
                if( nstep > 0 ) n_deex_succ = n_deex_succ + 1
!               Number of "call deexcitation".
                i_call_deex = i_call_deex + 1
            end if
!           For the special version of SFM.
            if( N > N0 .AND. i_had_model == 3 )then
!               Updates them to "sbe" and the string list.
                do ii = N0+1, N, 2
                    nbe = nbe + 1
                    do j=1,5,1
                        kbe( nbe, j ) = K(ii,j)
                        pbe( nbe, j ) = P(ii,j)
                        vbe( nbe, j ) = V(ii,j)
                    end do
                    nbe = nbe + 1
                    do j=1,5,1
                        kbe( nbe, j ) = K( ii+1, j )
                        pbe( nbe, j ) = P( ii+1, j )
                        vbe( nbe, j ) = V( ii+1, j )
                    end do
                    nstr1 = nstr1 + 1
                    nstr1a( nstr1 ) = nbe - 1
                    nstr1v( nstr1 ) = nbe
                end do
            end if
!       Updates N0 and does the deexcitation for newly produced qqbar pairs.
            if( i == N0 .AND. N > N0 &
                .AND. i_daught_generation < i_deex_generation )then
                i_begin = i
                i_daught_generation = i_daught_generation + 1
                N0 = N
                goto 700
            end if
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine compensate_conservation
!!      Compensates the lost 4-momentum in "throe_p" for "PYJETS".
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        real(kind=8) :: ps0(6), ps1(6)


!       Counts (original) sums of 4-momentum.
        ps0 = 0D0
        do j=1,4,1
            do i=1,N,1
                ps0(j) = ps0(j) + P(i,j)
            end do
            ps0(j) = ps0(j) + throe_p(j)
        end do
        throe_p = 0D0

        n_column = KSZJ
        n_row = 5
        call conse_c( P, n_column, n_row, ps0, 1, N, 1 )

!       Counts sums of 4-momentum after adjustments and collects the lost ones.
        ps1 = 0D0
        do j=1,4,1
            do i=1,N,1
                ps1(j) = ps1(j) + P(i,j)
            end do
            throe_p(j) = ps0(j) - ps1(j)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_coalescence( i_constraint, mode_coal )
!!      Performs the parton coalescence hadronization.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
!       mode_coal: = 1-5, traversal coalescence.
!                  other, random coalescence.


!       Old scheme.
        if( mode_coal >= 1 .AND. mode_coal <= 5 )then
            call do_traversal_coalescence( i_constraint )

!       New scheme.
        else
            call do_random_coalescence( i_constraint )
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_traversal_coalescence( iphas )
!!      Performs the parton coalescence hadronization via the list traversal.
!       iphas: = 1, complete phase space constraint
!              = 2, position constraint only
!              = 3, momentum constraint only
!       ithroq : total number of quarks thrown away
!       ithrob : total number of anti-quarks thrown away
!       throe : total four momentum of the partons thrown away
!       ich : total charge of the partons thrown away
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/coal1/bmrat,i_mm
        logical is_fail_try


!-------------------------------------------------------------------------------
!-------------------------   Variables Initializing   --------------------------
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0D0
!       Appends last failed quarks to the list, i.e. PYJETS + sa37.
        do i=1,nth,1
            N = N + 1
            do j=1,5,1
                K(N,j) = kth(i,j)
                P(N,j) = pth(i,j)
                V(N,j) = vth(i,j)
            end do
        end do
        nth = 0
!-------------------------   Variables Initializing   --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------------   Parton Sorting   ------------------------------
!       Sorts out PYJETS randomly. Lower case only. (In coales.f90)
        call PASORT( 1, N, "pyjets", "null", "random" )
!-----------------------------   Parton Sorting   ------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Coalescence   -------------------------------
!       Numer of baryon generated (baryon plus).
        n_baryon_p = 0
!       Number of anti-baryon generated (baryon minus).
        n_baryon_m = 0
        n_baryon   = 0
!       Number of meson generated.
        n_meson    = 0

        i_iteration = 0
        i_fail_iteration = 0
        max_fail_iteration = 50
        if( iphas /= 0 ) max_fail_iteration = 10
        is_fail_try = .false.

100     continue
        i_iteration = i_iteration + 1

!       Re-coalesces the failed q & qbar. No more than 50 / 10 times.
        if( is_fail_try ) i_fail_iteration = i_fail_iteration + 1

        do i1 = 1, N-1, 1
            KF1 = K(i1,2)
            do i2 = i1+1, N, 1
                KF2 = K(i2,2)
                isucc = 0

!-----------------------------   Meson Producing   -----------------------------
!       Tries to produce a meason.
                if( KF1 * KF2 < 0 )then
                    rand = PYR(1)
                    ! bmrat: controls the ratio of baryon to meson
                    ! if( rand <= bmrat ) cycle
                    sume0 = P(i1,4) + P(i2,4)
                    sump1 = P(i1,1) + P(i2,1)
                    sump2 = P(i1,2) + P(i2,2)
                    sump3 = P(i1,3) + P(i2,3)
                    cm = sume0**2 - sump1**2 - sump2**2 - sump3**2
                    if( cm <= 0D0 )then
                        write(*,*) "Warning, cm <= 0 in coal=", cm
                        cycle
                    end if
                    cm = SQRT(cm)

                    ! Finds the meson state.
                    call findm( KF1, KF2, cm, KF_out, dmass_out, isucc )

                    ! Failed (meson).
                    if( isucc == 0 ) cycle

                    ! Phase space adjudgment.
                    if( iphas /= 0 )then
                        i_M_or_B = 2
                        call phas( i1, i2, 0, isucc, i_M_or_B, iphas )
                        ! Failed (meson).
                        if( isucc == 0 ) cycle
                    end if

                    ! Gives properties to the primary meson.
                    nn = nn + 1
                    kn(nn,1) = 1
                    kn(nn,2) = KF_out
                    kn(nn,3) = 0
                    kn(nn,4) = 0
                    kn(nn,5) = 0
                    pn(nn,5) = PYMASS( KF_out )
                    pn(nn,1) = sump1
                    pn(nn,2) = sump2
                    pn(nn,3) = sump3
                    e2 = sump1**2 + sump2**2 + sump3**2 + PYMASS( KF_out )**2
                    e = SQRT(e2)
                    pn(nn,4) = e
                    ! Collects the energy.
                    throe_p(4) = throe_p(4) + sume0 - e

                    ! Produced hadron is set in between contituent partons.
                    pyrx = PYR(1)
                    pyry = PYR(1)
                    rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1)
                    rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2)
                    rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3)

                    ! Succeeded.
                    n_meson  = n_meson  + 1

                    ! Moves parton list "PYJETS" one step downward since i2+1.
                    ! In sfm_40.f90.
                    call updad_pyj( N, i2+1, 1 )
                    N = N - 1
                    ! Moves parton list "PYJETS" one step downward since i1+1.
                    call updad_pyj( N, i1+1, 1 )
                    N = N - 1
!-----------------------------   Meson Producing   -----------------------------

!----------------------------   Baryon Producing   -----------------------------
!       Tries to produce a baryon.
                else if( KF1 > 0 .AND. KF2 > 0 )then
                    rand = PYR(1)
                    ! bmrat: controls the ratio of baryon to meson
                    reltive_prob_BM = bmrat
                    ! Only failed q left. Generates baryons definately.
                    if( is_fail_try )then
                        reltive_prob_BM = 999D0
                        do i_left = i2+1, N, 1
                            if( K( i_left, 2 ) < 0 )then
                                reltive_prob_BM = bmrat
                                exit
                            end if
                        end do
                    end if
                    if( rand > reltive_prob_BM ) cycle
                    do i3 = i2+1, N, 1
                        KF3 = K(i3,2)
                        if( KF3 < 0 ) cycle
                        sume0 = P(i1,4) + P(i2,4) + P(i3,4)
                        sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                        sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                        sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                        cm = sume0**2 - sump1**2 - sump2**2 - sump3**2
                        if( cm <= 0D0 )then
                            write(*,*) "Warning, cm <= 0 in coal=", cm
                            cycle
                        end if
                        cm = SQRT(cm)

                        ! Finds the baryon state.
                        call findb( KF1, KF2, KF3, cm, KF_out,dmass_out, isucc )

                        ! Failed (baryon).
                        if( isucc == 0 ) cycle

                        ! Phase space adjudgment.
                        if( iphas /= 0 )then
                            i_M_or_B = 3
                            call phas( i1,i2,i3, isucc,i_M_or_B,iphas )
                            ! Failed (baryon).
                            if( isucc == 0 ) cycle
                        end if

                        ! Gives properties to the baryon.
                        nnol = nn
                        nn   = nn + 1
                        kn(nn,1) = 1
                        kn(nn,2) = KF_out
                        kn(nn,3) = 0
                        kn(nn,4) = 0
                        kn(nn,5) = 0
                        pn(nn,5) = PYMASS( KF_out )
                        pn(nn,1) = sump1
                        pn(nn,2) = sump2
                        pn(nn,3) = sump3
                        e2 = sump1**2 + sump2**2 + sump3**2 + PYMASS( KF_out )**2
                        e = SQRT(e2)
                        pn(nn,4) = e
                        ! Collects the energy.
                        throe_p(4) = throe_p(4) + sume0 - e

                        ! Produced hadron is arranged among constituent partons.
                        pyrx = PYR(1)
                        pyry = PYR(1)
                        pyrz = PYR(1)
                        rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
                        rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
                        rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

                        ! Succeeded.
                        n_baryon_p = n_baryon_p + 1
                        n_baryon   = n_baryon   + 1

                        ! Moves parton list one step downward from i3+1 to n1.
                        call updad_pyj( N, i3+1, 1 )
                        N = N - 1
                        ! Moves parton list one step downward from i2+1 to n1.
                        call updad_pyj( N, i2+1, 1 )
                        N = N - 1
                        ! Moves parton list one step downward from i1+1 to n1.
                        call updad_pyj( N, i1+1, 1 )
                        N = N - 1
                        exit
                    end do
!----------------------------   Baryon Producing   -----------------------------

!--------------------------   Anti-Baryon Producing   --------------------------
!       Tries to produce an anti-baryon.
                else if( KF1 < 0 .AND. KF2 < 0)then
                    rand = PYR(1)
                    reltive_prob_BM = bmrat
                    ! Only failed q left. Generates baryons definately.
                    if( is_fail_try )then
                        reltive_prob_BM = 999D0
                        do i_left = i2+1, N, 1
                            if( K( i_left, 2 ) > 0 )then
                                reltive_prob_BM = bmrat
                                exit
                            end if
                        end do
                    end if
                    if( rand > reltive_prob_BM ) cycle
                    do i3 = i2+1, N, 1
                        KF3 = K(i3,2)
                        if( KF3 > 0 ) cycle
                        sume0 = P(i1,4) + P(i2,4) + P(i3,4)
                        sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                        sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                        sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                        cm = sume0**2 - sump1**2 - sump2**2 - sump3**2
                        if( cm <= 0D0 )then
                            write(*,*) "Warning, cm <= 0 in coal=", cm
                            cycle
                        end if
                        cm = SQRT(cm)

                        ! Finds the anti-baryon state.
                        call findb( KF1, KF2, KF3, cm, KF_out,dmass_out, isucc )

                        ! Failed (anti-baryon).
                        if( isucc == 0 ) cycle

                        ! Phase space adjudgment.
                        if( iphas /= 0 )then
                            i_M_or_B = 3
                            call phas( i1,i2,i3, isucc,i_M_or_B,iphas )
                            ! Failed (anti-baryon).
                            if( isucc == 0 ) cycle
                        end if

                        ! Gives properties to the anti-baryon.
                        nnol = nn
                        nn   = nn + 1
                        kn(nn,1) = 1
                        kn(nn,2) = KF_out
                        kn(nn,3) = 0
                        kn(nn,4) = 0
                        kn(nn,5) = 0
                        pn(nn,5) = PYMASS( KF_out )
                        pn(nn,1) = sump1
                        pn(nn,2) = sump2
                        pn(nn,3) = sump3
                        e2 = sump1**2 + sump2**2 + sump3**2 + PYMASS( KF_out )**2
                        e = SQRT(e2)
                        pn(nn,4) = e
                        ! Collects the energy.
                        throe_p(4) = throe_p(4) + sume0 - e

                        ! Produced hadron is arranged among contituent partons.
                        pyrx = PYR(1)
                        pyry = PYR(1)
                        pyrz = PYR(1)
                        rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
                        rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
                        rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

                        ! Succeeded.
                        n_baryon_m = n_baryon_m + 1
                        n_baryon   = n_baryon   + 1

                        ! Moves parton list one step downward from i3+1 to n1.
                        call updad_pyj( N, i3+1, 1 )
                        N = N - 1
                        ! Moves parton list one step downward from i2+1 to n1.
                        call updad_pyj( N, i2+1, 1 )
                        N = N - 1
                        ! Moves parton list one step downward from i1+1 to n1.
                        call updad_pyj( N, i1+1, 1 )
                        N = N - 1
                        exit
                    end do
                end if
!--------------------------   Anti-Baryon Producing   --------------------------

!       Generates the next hadron.
                if( isucc == 1 ) goto 100
            end do
        end do

!       Re-coalesces the failed q & qbar. No more than 50 / 10 times.
        if( N >= 2 ) is_fail_try = .true.
        if( N >= 2 .AND. i_fail_iteration < max_fail_iteration ) goto 100
!-------------------------------   Coalescence   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Failed Particle Collecting   ------------------------
!       Collects failed partons, i.e. PYJETS -> sa37.
!       "PYJETS" -> "sa37"
        do i=1,N,1
            do j=1,5,1
                kth(i,j) = K(i,j)
                pth(i,j) = P(i,j)
                vth(i,j) = V(i,j)
            end do
            if( K(i,2) > 0 )then
                ithroq = ithroq + 1
                ich    = ich    + PYCHGE( K(i,2) )
                throe  = throe  + P(i,4)
            elseif( K(i,2) < 0 )then
                ithrob = ithrob + 1
                ich    = ich    + PYCHGE( K(i,2) )
                throe  = throe  + P(i,4)
            end if
        end do
        nth = N
        N = 0
!-----------------------   Failed Particle Collecting   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_random_coalescence( iphas )
!!      Performs the parton coalescence hadronization via the random selection.
!       iphas: = 1, complete phase space constraint
!              = 2, position constraint only
!              = 3, momentum constraint only
!       ithroq : total number of quarks thrown away
!       ithrob : total number of anti-quarks thrown away
!       throe : total four momentum of the partons thrown away
!       ich : total charge of the partons thrown away
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/coal1/bmrat,i_mm


!-------------------------------------------------------------------------------
!-------------------------   Variables Initializing   --------------------------
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0D0
!       Appends last failed quarks to the list, i.e. PYJETS + sa37.
        do i=1,nth,1
            N = N + 1
            do j=1,5,1
                K(N,j) = kth(i,j)
                P(N,j) = pth(i,j)
                V(N,j) = vth(i,j)
            end do
        end do
        nth = 0
!-------------------------   Variables Initializing   --------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------------   Parton Sorting   ------------------------------
!       Sort the partons in order of qbar and q.
        call PASORT( 1, N, "pyjets", "kf", "min_to_max" )
        i_q_begin = 0
        do i=1,N,1
            if( K(i,2) > 0 )then
                i_q_begin = i
                exit
            end if
        end do
!       Then sorts them randomly, respectively.
        call PASORT( 1, i_q_begin-1, "pyjets", "null", "random" )
        call PASORT( i_q_begin, N, "pyjets", "null", "random" )
!       Numbers of qbar and q.
        n_qbar = i_q_begin - 1
        n_q = N - i_q_begin + 1
!-----------------------------   Parton Sorting   ------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Coalescence   -------------------------------
!       rel_prob_BM: controls the relative probability of baryon over meson.
!       rel_prob_s:  extra factor for strangeness suppression.
        ! rel_prob_BM = adj1(31) / ( 1. + adj1(31) )
        ! rel_prob_s  = adj1(31)*adj1(33) / ( 1D0 + adj1(31)*adj1(33) )
        rel_prob_BM = bmrat
        ! Number of generation tries.
        n_try = 0
        ! Number of successful generation tries.
        n_try_succ  = 0
        ! Number of failures before one successful generatiopn..
        n_fail = 0
        ! Max number of failures.
        n_fail_max = n_fail
        ! Max number of failures allowed.
        n_fail_allowed = 10*N
        if( N > 5000 ) n_fail_allowed = 2*N
        ! Counters of hadrons.
        n_baryon = 0
        n_baryon_anti = 0
        n_meson = 0
        N0 = N

!       Generates a meson or baryon/anti-baryon according to rel_prob_BM.
        do while( .true. )
            n_try = n_try + 1

!----------------------------   Failures Checking   ----------------------------
            if( N == N0 .AND. n_try > 1 )then
                n_fail = n_fail + 1
                if( n_fail > n_fail_max ) n_fail_max = n_fail
            else if( N < N0 )then
                n_try_succ = n_try_succ + 1
                ! Reset.
                n_fail = 0
            end if
            ! Throws away the partons remained.
            if( n_fail > n_fail_allowed )then
                write(22,*) "Warning, n_fail > n_fail_allowed in random coal!" &
                        //  " iii, iphas n_fail_allowed=", &
                              iii, n_fail_allowed, iphas
                exit
            end if
!----------------------------   Failures Checking   ----------------------------

            isucc = 0
            N0    = N
            rand  = PYR(1)

!----------------------------   Baryon Generating   ----------------------------
!           Generates a baryon/anti-baryon.
            if( (rand <= rel_prob_BM) .AND. ( n_q >= 3 .OR. n_qbar >= 3 ) )then
!               Determines baryon (i_flavor = 1) or anti-baryon (i_flavor = -1).
                i_flavor = 1
                if( n_q >= 3 .AND. n_qbar >= 3 )then
                    i_flavor = 1
                    if( PYR(1) < 0.5D0 ) i_flavor = -1
                else if( n_q >= 3 )then
                    i_flavor = 1
                else if( n_qbar >= 3 )then
                    i_flavor = -1
                end if
!               Baryon.
                if( i_flavor == 1 )then
                    i_begin = n_qbar + 1
                    i_end   = N
!               Anti-baryon.
                else if( i_flavor == -1 )then
                    i_begin = 1
                    i_end   = n_qbar
                end if
!               Selects 3 q/qbar randomly.
!               "iRandom" is a function to select an integer in [i_begin,i_end].
                i1 = iRandom( i_begin, i_end )
                do while(.true.)
                    i2 = iRandom( i_begin, i_end )
                    if( i1 /= i2 ) exit
                end do
                do while(.true.)
                    i3 = iRandom( i_begin, i_end )
                    if( i1 /= i3 .AND. i2 /= i3 ) exit
                end do
                KF1 = K(i1,2)
                KF2 = K(i2,2)
                KF3 = K(i3,2)
                rel_prob = rel_prob_BM
                ! Strangeness suppression.
                ! if( ABS(KF1) == 3 .OR. ABS(KF2) == 3 .OR. ABS(KF3) == 3) &
                !     rel_prob = rel_prob_s
                if( rand > rel_prob ) cycle
                sume0 = P(i1,4) + P(i2,4) + P(i3,4)
                sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                cm = sume0**2 - sump1**2 - sump2**2 - sump3**2
                cm = SQRT(cm)

                ! Finds the baryon state.
                call findb_new( KF1, KF2, KF3, cm, KF_out, isucc )

                ! Failed (baryon/anti-baryon).
                if( isucc == 0 ) cycle

                ! Phase space adjudgment.
                if( iphas /= 0 )then
                    i_M_or_B = 3
                    call phas( i1, i2, i3, isucc, i_M_or_B, iphas )
                    ! Failed (baryon/anti-baryon).
                    if( isucc == 0 ) cycle
                end if

                ! Gives properties to the baryon/anti-baryon.
                nnol = nn
                nn   = nn + 1
                kn(nn,1) = 1
                kn(nn,2) = KF_out
                kn(nn,3) = 0
                kn(nn,4) = 0
                kn(nn,5) = 0
                pn(nn,5) = PYMASS( KF_out )
                pn(nn,1) = sump1
                pn(nn,2) = sump2
                pn(nn,3) = sump3
                e2 = sump1**2 + sump2**2 + sump3**2 + PYMASS( KF_out )**2
                e = SQRT(e2)
                pn(nn,4) = e
                ! Collects the energy.
                throe_p(4) = throe_p(4) + sume0 - e

                ! Produced hadron is arranged among contituent partons.
                pyrx = PYR(1)
                pyry = PYR(1)
                pyrz = PYR(1)
                rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
                rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
                rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

                ! Succeeded.
                if( KF_out > 0 ) n_baryon      = n_baryon      + 1
                if( KF_out < 0 ) n_baryon_anti = n_baryon_anti + 1

                ! Sort i1, i2 and i3 from min to max.
                if( i1 > i2 )then
                    ix = i2
                    i2 = i1
                    i1 = ix
                end if
                if( i2 > i3 )then
                    ix = i3
                    i3 = i2
                    i2 = ix
                end if
                if( i1 > i2 )then
                    ix = i2
                    i2 = i1
                    i1 = ix
                end if

                ! Moves parton list one step downward from i3+1 to n1.
                ! In sfm_40.f90.
                call updad_pyj( N, i3+1, 1 )
                N = N - 1
                ! Moves parton list one step downward from i2+1 to n1.
                call updad_pyj( N, i2+1, 1 )
                N = N - 1
                ! Moves parton list one step downward from i1+1 to n1.
                call updad_pyj( N, i1+1, 1 )
                N = N - 1
                if( KF1*KF2*KF3 > 0 ) n_q    = n_q    - 3
                if( KF1*KF2*KF3 < 0 ) n_qbar = n_qbar - 3
!----------------------------   Baryon Generating   ----------------------------


!----------------------------   Meson Generating   -----------------------------
!           Generates a meson.
            else if( n_q > 0 .AND. n_qbar > 0 )then
!               Selects a q and a qbar randomly.
                i1 = iRandom( n_qbar+1, N )
                i2 = iRandom( 1, n_qbar )
                KF1 = K(i1,2)
                KF2 = K(i2,2)
                sume0 = P(i1,4) + P(i2,4)
                sump1 = P(i1,1) + P(i2,1)
                sump2 = P(i1,2) + P(i2,2)
                sump3 = P(i1,3) + P(i2,3)
                cm = sume0**2 - sump1**2 - sump2**2 - sump3**2
                cm = SQRT(cm)

                ! Finds the meson state.
                call findm_new( KF1, KF2, cm, KF_out, isucc )

                ! Failed (meson).
                if( isucc == 0 ) cycle

                ! Phase space adjudgment.
                if( iphas /= 0 )then
                    i_M_or_B = 2
                    call phas( i1, i2, 0, isucc, i_M_or_B, iphas )
                    ! Failed (meson).
                    if( isucc == 0 ) cycle
                end if

                ! Gives properties to the primary meson.
                nn = nn + 1
                kn(nn,1) = 1
                kn(nn,2) = KF_out
                kn(nn,3) = 0
                kn(nn,4) = 0
                kn(nn,5) = 0
                pn(nn,5) = PYMASS( KF_out )
                pn(nn,1) = sump1
                pn(nn,2) = sump2
                pn(nn,3) = sump3
                e2 = sump1**2 + sump2**2 + sump3**2 + PYMASS( KF_out )**2
                e = SQRT(e2)
                pn(nn,4) = e
                ! Collects the energy.
                throe_p(4) = throe_p(4) + sume0 - e

                ! Produced hadron is set in between contituent partons.
                pyrx = PYR(1)
                pyry = PYR(1)
                rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1)
                rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2)
                rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3)

                ! Succeeded.
                n_meson  = n_meson  + 1

                ! Moves parton list "PYJETS" one step downward since i1+1.
                call updad_pyj( N, i1+1, 1 )
                N = N - 1
                ! Moves parton list "PYJETS" one step downward since i2+1.
                call updad_pyj( N, i2+1, 1 )
                N = N - 1
                n_q    = n_q    - 1
                n_qbar = n_qbar - 1
            end if
!----------------------------   Meson Generating   -----------------------------

            if( N == 0 ) exit

        end do
!-------------------------------   Coalescence   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Failed Particle Collecting   ------------------------
!       Collects failed partons, i.e. PYJETS -> sa37.
!       "PYJETS" -> "sa37"
        if( N > 0 )then
            do i=1,N,1
                do j=1,5,1
                    kth(i,j) = K(i,j)
                    pth(i,j) = P(i,j)
                    vth(i,j) = V(i,j)
                end do
                if( K(i,2) > 0 )then
                    ithroq = ithroq + 1
                    ich    = ich    + PYCHGE( K(i,2) )
                    throe  = throe  + P(i,4)
                elseif( K(i,2) < 0 )then
                    ithrob = ithrob + 1
                    ich    = ich    + PYCHGE( K(i,2) )
                    throe  = throe  + P(i,4)
                end if
            end do
            nth = N
            N = 0
        else
            nth = 0
            N = 0
        end if
!-----------------------   Failed Particle Collecting   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tabhb
!!      Table of primary meson (pseudoscalar (spin 0) and vector (spin 1)
!!       only) and baryon (octet (spin 1/2) and decuplet (spin 3/2) only)
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
!       imc (ibc): dimension of meson (baryon) table considered
!       kqh(i,j): flavor code of j-th constituent q (or qbar) in i-th quark
!        configuration of meson
!       kfh(i,1), kfh(i,2): flavor code of pseudoscalar meson and vector meson
!        with i-th quark configuration of meson
!       proh(i,1), proh(i,2): proper probability of meson
!       amash(i,1), amash(i,2): mass of meson
!       kqb(i,j): flavor code of j-th constituent q (or qbar) in i-th quark
!        configuration of baryon
!       kfb(i,1), kfb(i,2): flavor code of octet baryon and decuplet baryon
!        with i-th quark configuration of baryon
!       prob(i,1), prob(i,2): probability of baryon
!       amasb(i,1), amasb(i,2): mass of baryon

!       New table.
!***********************************************************************
!       imc: 26 -> 30; ibc: 18 -> 45
!       Meson  : 20 + 20 + 6 + 4 = 50  (meson  + anti- + onium)
!       Baryon : 75 + 75         = 150 (baryon + anti-)
!       Total  : 50 + 150        = 200
!-----------------------------------------------------------------------
!       Onium (no anti-): 6 + 4
!       Light mixing onium: pi0, rho0, eta, omega, eta', phi
!       Heavy onium : eta_c, J/psi, eta_b, Upsilon
!***********************************************************************
!       Meson table.
!-----------------------------------------------------------------------
!       KF code of quark and anti-quark.
!       Note that kqh1 > 0 and kqh2 < 0.
        data (kqh(i,1),i=1,80) &
             /  1,  1,  1,  1,  2,  1,  3,  1,  4,  1, &   ! 10
                5,  2,  2,  2,  2,  3,  2,  4,  2,  5, &   ! 20
                3,  3,  3,  4,  3,  5,  4,  4,  5,  5, &   ! 30
                50*0 /                                     ! 80
        data (kqh(i,2),i=1,80) &
             / -1, -1, -1, -2, -1, -3, -1, -4, -1, -5, &   ! 10
               -1, -2, -2, -2, -3, -2, -4, -2, -5, -2, &   ! 20
               -3, -3, -4, -3, -5, -3, -4, -5, -4, -5, &   ! 30
                50*0 /                                     ! 80
!-----------------------------------------------------------------------
!       KF code of hadron (meson).
!       Pseudoscalar meson, s = 0.
        data (kfh(i,1),i=1,80) &
         /  111,  221,  331, -211,  211,  311, -311, -411,  411,  511, &   ! 10
           -511,  111,  221,  331,  321, -321, -421,  421,  521, -521, &   ! 20
            221,  331, -431,  431,  531, -531,  441,  541, -541,  551, &   ! 30
            50*0 /                                                         ! 80
!       Vector meson, s = 1.
        data (kfh(i,2),i=1,80) &
         /  113,  223,    0, -213,  213,  313, -313, -413,  413,  513, &   ! 10
           -513,  113,  223,    0,  323, -323, -423,  423,  523, -523, &   ! 20
            333,    0, -433,  433,  533, -533,  443,  543, -543,  553, &   ! 30
            50*0 /                                                         ! 80
!-----------------------------------------------------------------------
!       Probability of hadron (meson).
        data (proh(i,1),i=1,80) &
         /  0.5,0.167,0.333,   1.,   1.,   1.,   1.,   1.,   1.,   1., &   ! 10
             1.,  0.5,0.167,0.333,   1.,   1.,   1.,   1.,   1.,   1., &   ! 20
          0.667,0.333,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1., &   ! 30
            50*0. /                                                        ! 80
        data (proh(i,2),i=1,80) &
         /  0.5,  0.5,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   1., &   ! 10
             1.,  0.5,  0.5,   0.,   1.,   1.,   1.,   1.,   1.,   1., &   ! 20
             1.,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1., &   ! 30
            50*0. /                                                        ! 80
!***********************************************************************
!       Baryon table.
!-----------------------------------------------------------------------
!       KF code of 3-quarks.
! Note that kqb1 <= kqb2 <= kqb3.
        data (kqb(i,1),i=1,80) &
             /  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &   ! 10
                1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &   ! 20
                1,  2,  2,  2,  2,  2,  2,  2,  2,  2, &   ! 30
                2,  2,  2,  2,  3,  3,  3,  3,  3,  3, &   ! 40
                3,  4,  4,  4,  5,                     &   ! 45
                0,  0,  0,  0,  0,                     &   ! 50
                30*0 /                                     ! 80
        data (kqb(i,2),i=1,80) &
             /  1,  1,  1,  1,  1,  2,  2,  2,  2,  2, &   ! 10
                2,  2,  3,  3,  3,  3,  3,  4,  4,  4, &   ! 20
                5,  2,  2,  2,  2,  3,  3,  3,  3,  3, &   ! 30
                4,  4,  4,  5,  3,  3,  3,  4,  4,  4, &   ! 40
                5,  4,  4,  5,  5,                     &   ! 45
                0,  0,  0,  0,  0,                     &   ! 50
                30*0 /                                     ! 80
        data (kqb(i,3),i=1,80) &
             /  1,  2,  3,  4,  5,  2,  3,  3,  4,  4, &   ! 10
                5,  5,  3,  4,  4,  5,  5,  4,  5,  5, &   ! 20
                5,  2,  3,  4,  5,  3,  4,  4,  5,  5, &   ! 30
                4,  5,  5,  5,  3,  4,  5,  4,  5,  5, &   ! 40
                5,  4,  5,  5,  5,                     &   ! 45
                0,  0,  0,  0,  0,                     &   ! 50
                30*0 /                                     ! 80
!-----------------------------------------------------------------------
!       KF code of baryon.
!       Octet baryon, s = 1/2.
        data (kfb(i,1),i=1,80) &
         /    0, 2112, 3112, 4112, 5112, 2212, 3122, 3212, 4122, 4212, &   ! 10
           5122, 5212, 3312, 4132, 4312, 5132, 5312, 4412, 5142, 5412, &   ! 20
           5512,    0, 3222, 4222, 5222, 3322, 4232, 4322, 5232, 5322, &   ! 30
           4422, 5242, 5422, 5522,    0, 4332, 5332, 4432, 5342, 5432, &   ! 40
           5532,    0, 5442, 5542,    0,                               &   ! 45
               0,   0,    0,    0,    0,                               &   ! 50
           30*0 /                                                          ! 80
!       Decuplet baryon, s = 3/2.
        data (kfb(i,2),i=1,80) &
         / 1114, 2114, 3114, 4114, 5114, 2214,    0, 3214,    0, 4214, &   ! 10
              0, 5214, 3314,    0, 4314,    0, 5314, 4414,    0, 5414, &   ! 20
           5514, 2224, 3224, 4224, 5224, 3324,    0, 4324,    0, 5324, &   ! 30
           4424,    0, 5424, 5524, 3334, 4334, 5334, 4434,    0, 5434, &   ! 40
           5534, 4444, 5444, 5544, 5554,                               &   ! 45
               0,   0,    0,    0,    0,                               &   ! 50
           30*0 /                                                          ! 80
!-----------------------------------------------------------------------
!       Probability of baryon.
        data (prob(i,1),i=1,80) &
         /   0.,   1.,   1.,   1.,   1.,   1.,  0.5,  0.5,  0.5,  0.5, &   ! 10
            0.5,  0.5,   1.,  0.5,  0.5,  0.5,  0.5,   1.,  0.5,  0.5, &   ! 20
             1.,   0.,   1.,   1.,   1.,   1.,  0.5,  0.5,  0.5,  0.5, &   ! 30
             1.,  0.5,  0.5,   1.,   0.,   1.,   1.,   1.,  0.5,  0.5, &   ! 40
             1.,   0.,   1.,   1.,   0.,                               &   ! 45
             0.,   0.,   0.,   0.,   0.,                               &   ! 50
            30*0. /                                                        ! 80
        data (prob(i,2),i=1,80) &
          /  1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   1., &   ! 10
             0.,   1.,   1.,   0.,   1.,   0.,   1.,   1.,   0.,   1., &   ! 20
             1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   1., &   ! 30
             1.,   0.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   1., &   ! 40
             1.,   1.,   1.,   1.,   1.,                               &   ! 45
             0.,   0.,   0.,   0.,   0.,                               &   ! 50
            30*0. /                                                        ! 80
!***********************************************************************


        do i1=1,imc,1
            KF1 = kfh(i1,1)
            KF2 = kfh(i1,2)
            amash(i1,1) = PYMASS(KF1)
            amash(i1,2) = PYMASS(KF2)
        end do
        do i1=1,ibc,1
            KF1 = kfb(i1,1)
            KF2 = kfb(i1,2)
            amasb(i1,1) = PYMASS(KF1)
            amasb(i1,2) = PYMASS(KF2)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine hadpro(iphas)
!!      Parton coalescence (hadronization)
!       iphas: = 1, complete phase space constraint
!              = 2, position constraint only
!              = 3, momentum constraint only
!       ithroq : total number of quarks thrown away
!       ithrob : total number of anti-quarks thrown away
!       throe : total four momentum of the partons thrown away
!       ich : total charge of the partons thrown away
!       Double precision and integer declarations.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/coal1/bmrat,i_mm


!-------------------------------------------------------------------------------
!-------------------------   Variables initialization   ------------------------
!       Moved from 'coales' to here.
        ithroq = 0
        ithrob = 0
        ich    = 0
        throe  = 0D0
!       Appends the last failed quarks to list, i.e. PYJETS + sa37.
        if( nth > 0 )then
            do i=1,nth,1
                N = N + 1
                do j=1,5,1
                    K(N,j) = kth(i,j)
                    P(N,j) = pth(i,j)
                    V(N,j) = vth(i,j)
                end do
            end do
            nth = 0
        end if
!-------------------------   Variables initialization   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Coalescence   -------------------------------
        ibarp = 0   ! Numer of baryon generated. (baryon plus)
        ibarm = 0   ! Number of anti-baryon generated. (baryon minus)
        imes  = 0   ! Number of meson generated.
        nme   = 0
        nba   = 0

        i_fail_iteration = 0
        i_normal_coal = 0

390     continue

!       Re-coalesces the failed q & qbar. No more than 50 times.
        if( i_normal_coal == 1 ) i_fail_iteration = i_fail_iteration + 1

        do 400 i1=1,N-2,1
        kf1=K(i1,2)
        do 500 i2=i1+1,N-1
            kf2=K(i2,2)

!-----------------------------   Meson Producing   -----------------------------
!       Tries to produce a meason.
        if( (kf1 > 0.and.kf2 < 0) .or. (kf1 < 0.and.kf2 > 0) )then

            sume  = P(i1,4) + P(i2,4)
            sump1 = P(i1,1) + P(i2,1)
            sump2 = P(i1,2) + P(i2,2)
            sump3 = P(i1,3) + P(i2,3)
            cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
            if( cm > 1D20 ) cm = 1D20
            if( cm <= 0D0  ) goto 500
            cm = SQRT(cm)
            KF_in_1 = kf1
            KF_in_2 = kf2
!       Exchanges KFs of qbar and q to ensure the first one is q.
            if( kf1 < 0 )then
                KF_in_1 = kf2
                KF_in_2 = kf1
            end if
            call findm( KF_in_1, KF_in_2, cm, kfii, amasi, isucc )
            if(isucc == 0) goto 500

!       Phase space adjudgment.
            if( iphas /= 0 )then
            call phas(i1,i2,0,isucc,2,iphas)
            if( isucc == 0 ) goto 500
            endif

!       Proceed for success
            imes = imes + 1
            nme  = nme  + 1

!       Give proper variables to the primary meson.
            nnol = nn
            nn   = nn+1
            kn(nn,1) = 1
            kn(nn,2) = kfii
            kn(nn,3) = 0
            kn(nn,4) = 0
            kn(nn,5) = 0
            pn(nn,5) = amasi
            pn(nn,1) = sump1
            pn(nn,2) = sump2
            pn(nn,3) = sump3
            pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
            pnnmm = amasi*amasi + pnnm
            if( pnnmm > 1D20 ) pnnmm = 1D20
            if( pnnmm <= 0D0  ) pnnmm = 1D-20
            pnnn = SQRT(pnnmm)
            pn(nn,4)   = pnnn
            throe_p(4) = throe_p(4) + sume-pnnn

!       Produced hadron is set in between contituent partons randomly.
            pyrx = PYR(1)
            pyry = PYR(1)
            rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1)
            rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2)
            rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3)

!       Move parton list ('PYJETS') one step downward since i2+1.
            ! this subroutine is in 'sfm_40.f90'
            call updad_pyj(N,i2+1,1)
            N = N - 1

!       Move parton list ('PYJETS') one step downward since i1+1.
            call updad_pyj(N,i1+1,1)
            N = N - 1
!       Share the surplus 4-momentum in throe_p.
            ! call share_p_PYJETS
            call share_p_PYJETS_sa1h
            ! to construct three cycle again
            goto 390
        endif
!-----------------------------   Meson Producing   -----------------------------

!----------------------------   Baryon Producing   -----------------------------
!       Tries to produce a baryon.
            if(kf1 > 0.and.kf2 > 0)then
                rand=pyr(1)
                ! bmrat: ratio of baryon to meson
                if(rand > bmrat) goto 500
                do 600 i3=i2+1,N
                kf3=K(i3,2)
                if(kf3 < 0) goto 600
                sume  = P(i1,4) + P(i2,4) + P(i3,4)
                sump1 = P(i1,1) + P(i2,1) + P(i3,1)
                sump2 = P(i1,2) + P(i2,2) + P(i3,2)
                sump3 = P(i1,3) + P(i2,3) + P(i3,3)
                cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
                if( cm > 1D20 ) cm = 1D20
                if( cm <= 0D0  ) goto 600
                cm = SQRT(cm)

!       Find out the baryon from hadron table.
                call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc)
                if(isucc == 0) goto 600

!       Phase space adjudgment.
                if( iphas /= 0 )then
                call phas(i1,i2,i3,isucc,3,iphas)
                if( isucc == 0 ) goto 600
                endif

!       Proceed for success.
                ibarp = ibarp + 1
                nba   = nba   + 1

!       Give proper variables to the baryon.
                nnol  = nn
                nn    = nn+1
                kn(nn,1) = 1
                kn(nn,2) = kfii
                kn(nn,3) = 0
                kn(nn,4) = 0
                kn(nn,5) = 0
                pn(nn,5) = amasi
                pn(nn,1) = sump1
                pn(nn,2) = sump2
                pn(nn,3) = sump3
                pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
                pnnmm = amasi*amasi + pnnm
                if( pnnmm > 1D20 ) pnnmm = 1D20
                if( pnnmm <= 0D0  ) pnnmm = 1D-20
                pnnn = SQRT(pnnmm)
                pn(nn,4)   = pnnn
                throe_p(4) = throe_p(4) + sume-pnnn

!       Produced hadron is arranged among constituent partons randomly.
                pyrx = PYR(1)
                pyry = PYR(1)
                pyrz = PYR(1)
                rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
                rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
                rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

!       Move parton list one step downward from i3+1 to n1.
                call updad_pyj(N,i3+1,1)
                N = N - 1

!       Move parton list one step downward from i2+1 to n1.
                call updad_pyj(N,i2+1,1)
                N  = N  - 1

!       Move parton list one step downward from i1+1 to n1.
                call updad_pyj(N,i1+1,1)
                N = N - 1
!       Share the surplus 4-momentum in throe_p.
                ! call share_p_PYJETS
                call share_p_PYJETS_sa1h
                ! to construct three cycle again
                goto 390
600             continue
            endif
!----------------------------   Baryon Producing   -----------------------------

!--------------------------   Anti-Baryon Producing   --------------------------
!       Tries to produce an anti-baryon.
        if(kf1 < 0.and.kf2 < 0)then
            rand=pyr(1)
            if(rand > bmrat) goto 500
            do 700 i3=i2+1,N,1
            kf3=K(i3,2)
            if(kf3 > 0) goto 700
            sume  = P(i1,4) + P(i2,4) + P(i3,4)
            sump1 = P(i1,1) + P(i2,1) + P(i3,1)
            sump2 = P(i1,2) + P(i2,2) + P(i3,2)
            sump3 = P(i1,3) + P(i2,3) + P(i3,3)
            cm = sume*sume - sump1*sump1 - sump2*sump2 - sump3*sump3
            if( cm > 1D20 ) cm = 1D20
            if( cm <= 0D0  ) goto 700
            cm = SQRT(cm)

!       Find out the anti-baryon from hadron table.
            call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc)
            if(isucc == 0) goto 700

!       Phase space adjudgment.
            if( iphas /= 0 )then
            call phas(i1,i2,i3,isucc,3,iphas)
            if( isucc == 0 ) goto 700
            endif

            ibarm = ibarm + 1
            nba = nba + 1

!       Give proper variables to the anti-baryon.
            nnol = nn
            nn   = nn + 1
            kn(nn,1) = 1
            kn(nn,2) = kfii
            kn(nn,3) = 0
            kn(nn,4) = 0
            kn(nn,5) = 0
            pn(nn,5) = amasi
            pn(nn,1) = sump1
            pn(nn,2) = sump2
            pn(nn,3) = sump3
            pnnm  = sump1*sump1 + sump2*sump2 + sump3*sump3
            pnnmm = amasi*amasi + pnnm
            if( pnnmm > 1D20 ) pnnmm = 1D20
            if( pnnmm <= 0D0  ) pnnmm = 1D-20
            pnnn = SQRT(pnnmm)
            pn(nn,4)   = pnnn
            throe_p(4) = throe_p(4) + sume-pnnn

!       Produced hadron is arranged among contituent partons randomly.
            pyrx = PYR(1)
            pyry = PYR(1)
            pyrz = PYR(1)
            rn(nn,1) = pyrx*V(i1,1) + pyry*V(i2,1) + pyrz*V(i3,1)
            rn(nn,2) = pyrx*V(i1,2) + pyry*V(i2,2) + pyrz*V(i3,2)
            rn(nn,3) = pyrx*V(i1,3) + pyry*V(i2,3) + pyrz*V(i3,3)

!       Move parton list one step downward from i3+1 to n1.
            call updad_pyj(N,i3+1,1)
            N = N - 1
!       Move parton list one step downward from i2+1 to n1.
            call updad_pyj(N,i2+1,1)
            N = N - 1
!       Move parton list one step downward from i1+1 to n1.
            call updad_pyj(N,i1+1,1)
            N = N - 1
!       Share the surplus 4-momentum in throe_p.
            ! call share_p_PYJETS
            call share_p_PYJETS_sa1h
            ! to construct three cycle again
            goto 390
700         continue
        endif
!--------------------------   Anti-Baryon Producing   --------------------------

500     continue
400     continue
!       Re-coalesces the failed q & qbar. No more than 50 times.
        if( N >= 2 ) i_normal_coal = 1
        if( N >= 2 .AND. i_fail_iteration < 50 )then
        ! if( i_normal_coal == 1 .AND. i_fail_iteration < 10 )then
            goto 390
        end if
!-------------------------------   Coalescence   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Failed Particle Collecting   ------------------------
!       Collects failed quarks that may fail, i.e. PYJETS -> sa37.
!       'PYJETS' -> 'sa37'
        if( N > 0 )then
            do i=1,N,1
                do j=1,5,1
                    kth(i,j) = K(i,j)
                    pth(i,j) = P(i,j)
                    vth(i,j) = V(i,j)
                end do
                if( K(i,2) > 0 )then
                    ithroq = ithroq + 1
                    ich    = ich    + PYCHGE(K(i,2))
                    throe  = throe  + P(i,4)
                elseif( K(i,2) < 0 )then
                    ithrob = ithrob + 1
                    ich    = ich    + PYCHGE(K(i,2))
                    throe  = throe  + P(i,4)
                end if
            end do
            nth = N
            N = 0
        else
            nth = 0
            N = 0
        end if
!-----------------------   Failed Particle Collecting   ------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findm( kf1, kf2, cm_in, kfii, amasi, isucc )
!!      Finds the primary meson from mesonic table according to kf1 & kf2
!       cm : invariant mass of kf1 & kf2
!       kfii : flavor code of the primary meson
!       amasi : mass of the primary meson
!       isucc = 1 : success
!             = 0 : fail
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc


        kfii = 0
        amasi = 0D0
        isucc = 0
        cm = cm_in
        if( kf1*kf2 >=0 ) return
        if1 = kf1
        if2 = kf2
        ! Makes sure the first one is q.
        if( if1 < 0 )then
            if1 = kf2
            if2 = kf1
        end if

        n_succ_flavor_try = 0
        do i4=1,imc,1
            if( n_succ_flavor_try > 0 )Then
                ! uubar, ddbar light mixing onia. 3 combinations.
                if( (if1 == 1 .AND. if2 == -1) &
                    .OR. (if1 == 2 .AND. if2 == -2) )then
                        if( n_succ_flavor_try > 3 ) exit
                ! ssbar light mixing onia. 2 combinations.
                else if( (if1 == 3 .AND. if2 == -3) )then
                        if( n_succ_flavor_try > 2 ) exit
                ! Others. Only 1 combination. Failed.
                else
                    exit
                end if
            end if
            kfi = kqh(i4,1)
            kfj = kqh(i4,2)
            ! Flavor check.
            if( kfi /= if1 .OR. kfj /= if2 ) cycle
            n_succ_flavor_try = n_succ_flavor_try + 1
            amas1 = amash(i4,1)
            amas2 = amash(i4,2)
            ! The mass approaching method.
                ! Mass difference.
                ! amas1 = ABS( cm - amas1 )
                ! amas2 = ABS( cm - amas2 )
            ! The mass reciprocal method.
                ! Probability for pseudoscalar.
                xpseud = amas2 / ( amas1 + amas2 )
                ! Probability for vector.
                xvector = 1D0 - xpseud
            rand = PYR(1)
            i_pseud_or_vec = 1
            ! Vector wins.
            if( ABS( proh(i4,1) ) > 1D-5 .AND. ABS( proh(i4,2) ) > 1D-5 &
                .AND. rand > xpseud ) then
                ! .AND. amas2 <= amas1 )
                 i_pseud_or_vec = 2
            end if

            ! Succeeded.
            ran1  = PYR(1)
            proi  = proh( i4, i_pseud_or_vec )
            if( ran1 > proi ) cycle
            kfii  = kfh( i4, i_pseud_or_vec )
            amasi = amash( i4, i_pseud_or_vec )
            isucc = 1
            exit
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findb( kf1, kf2, kf3, cm_in, kfii, amasi, isucc )
!!      Finds the primary baryon (antibaryon) from baryonic table
!!       according to kf1, kf2 and kf3.
!       cm: invariant mass of kf1, kf2 & kf3
!       kfii : flavor code of the primary baryon
!       amasi : mass of the primary baryon
!       isucc = 1 : succeed
!       isucc = 0 : fail
!       iflav = 1 : if composing parton is quarks
!             =-1 : if composing parton is antiquarks
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc


        kfii = 0
        amasi = 0D0
        isucc = 0
        cm = cm_in
        if( kf1*kf2*kf3 == 0 &
            .OR. SIGN(1,kf1) /= SIGN(1,kf2) .OR. SIGN(1,kf1) /= SIGN(1,kf3) &
            .OR. SIGN(1,kf2) /= SIGN(1,kf3) ) return
        iflav = 1
        if( kf1 < 0 ) iflav = -1
        if1 = ABS( kf1 )
        if2 = ABS( kf2 )
        if3 = ABS( kf3 )
        ! Makes sure if1 <= if2 <= if3.
        if( if1 > if2 )then
            ifx = if1
            if1 = if2
            if2 = ifx
        end if
        if( if1 > if3 )then
            ifx = if1
            if1 = if3
            if3 = ifx
        end if
        if( if2 > if3 )then
            ifx = if2
            if2 = if3
            if3 = ifx
        end if

        i_baryon_type = 0
        do i4=1,ibc,1
            if( i_baryon_type == 3 ) exit
            kfi = kqb(i4,1)
            kfj = kqb(i4,2)
            kfk = kqb(i4,3)
            ! Flavor check.
            if( kfi /= if1 .OR. kfj /= if2 .OR. kfk /= if3 ) cycle
            i_spin = 1
            ! xxx type baryon. Only one state.
            if( if1 == if2 .AND. if1 == if3 )then
                i_baryon_type = 1
                i_spin = 2
            ! First xyz type baryon isospin-singlet. Only one state.
            else if( if1 /= if2 .AND. if1 /= if3 .AND. if2 /= if3 &
                .AND. ABS( prob(i4,2) ) < 1D-5 )then
                    i_baryon_type = 2
                    i_spin = 1
            ! Another xyz type and xxy(xyy) type baryon. Two states.
            else if( ABS( prob(i4,1) ) > 1D-5 &
               .AND. ABS( prob(i4,2) ) > 1D-5 )then
                i_baryon_type = 3
                amas1 = amasb(i4,1)
                amas2 = amasb(i4,2)
                ! The mass approaching method.
                ! Mass difference.
                    ! amas1 = ABS( cm - amas1 )
                    ! amas2 = ABS( cm - amas2 )
                ! The mass reciprocal method.
                    ! Probability for 1/2 octet.
                    xoct  = amas2 / ( amas1 + amas2 )
                    ! Probability for 3/2 decuplet.
                    xdec = 1D0 - xoct
                rand = PYR(1)
                ! 3/2 decuplet wins.
                if( rand > xoct ) i_spin = 2
                ! if( amas2 <= amas1 ) i_spin = 2
            end if

            ! Succeeded.
            ran1 = PYR(1)
            proi = prob( i4, i_spin )
            if( ran1 > proi ) cycle
            kfii = kfb( i4, i_spin )
            ! Anti-baryon.
            if( iflav == -1 ) kfii = -kfii
            amasi = amasb( i4, i_spin )
            isucc = 1
            exit
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findm_new( KF1, KF2, inv_mass_q, KF_h, i_succ )
!!       Finds the primary meson according to the rule of KF code directly.
!=======================================================================
!       Original method, i_find = 10 / 11 :
!       The state is determined by
!         10: the proper probability and the mass approaching method.
!         11: the proper probability and the mass reciprocal method.
!       The rule of KF code cf. [ JHEP 05 (2006) 026, P73 ].
!       The proper probability is defined as the normalization factor
!        of the hadronic wave function. Cf. "subroutine tabh" and
!        [ H. Georgi, Lie Algebras In Particle Physics : from Isospin
!          To Unified Theories (CRC Press, Taylor & Francis, Boca
!          Raton, 2000) ].
!=======================================================================
!       An alternative method, i_find = 20 / 21 :
!         20: the relative probability and the mass approaching method.
!         21: the relative probability and the mass reciprocal method.
!=======================================================================
!       Assume a normal meson with quark KF_i and anti-quark -KF_j,
!        KF_i /= KF_j (KF_i, KF_j > 0), total spin s.
!
!       KF of a normal meson not orbitally or radially excited:
!           KF_h = { MAX(KF_i, KF_j)*100 + MIN(KF_i, KF_j)*100 2*s + 1 }
!                  *SIGN( KF_i - KF_j ) *(-1)^MAX(KF_i, KF_j)
!
!       For onia:
!        (1) ddbar -- pi0 111, rho0 113; mixed with eta, omega and eta'.
!        (2) uubar -- eta 221, omega 223; mixed with pi0, rho0 and eta'.
!        (3) ssbar -- eta' 331, phi 333; mixed with eta.
!        (4) ccbar -- eta_c 441, J/psi 443;
!        (5) bbbar -- eta_b 551, Upsilon 553.
!
!       The meson are divied into 2 cases:
!        (1) x-ybar type, with 2 different flavors, e.g. pi+ -- udbar;
!        (2) x-xbar type, onia, with 2 same flavors, e.g. J/psi -- ccbar;
!=======================================================================
!       Only primary lowest-lying spin 0 pseudoscalar and spin 1 vector
!        states are considered. The orbitally or radially excitated
!        states are not included at present.
!        Spin 0, KF is **1.
!        Spin 1, KF is **3.
!       inv_mass_q: common invariant mass of kf1 and kf2.
!       KF_h : flavor code of the found primary meson.
!       i_succ = 1 : succeed,
!              = 0 : fail.
!=======================================================================
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension   KF_mix(2,3,5), n_mix_pair(5)
        dimension prob_mix(2,3,5)
        real(kind=8) inv_mass_q

!=======================================================================
!       For light mixing states and heavy onia.
        data n_mix_pair / 3, 3, 2, 1, 1 /
!       KF of onia.
!       ddbar
        data ( KF_mix(i,1,1), i=1,2,1 ) / 111, 113 /   ! pi0, rho0
        data ( KF_mix(i,2,1), i=1,2,1 ) / 221, 223 /   ! eta, omega
        data ( KF_mix(i,3,1), i=1,2,1 ) / 331, 331 /   ! eta'
!       uubar
        data ( KF_mix(i,1,2), i=1,2,1 ) / 111, 113 /   ! pi0, rho0
        data ( KF_mix(i,2,2), i=1,2,1 ) / 221, 223 /   ! eta, omega
        data ( KF_mix(i,3,2), i=1,2,1 ) / 331, 331 /   ! eta'
!       ssbar
        data ( KF_mix(i,1,3), i=1,2,1 ) / 221, 333 /   ! eta, phi
        data ( KF_mix(i,2,3), i=1,2,1 ) / 331, 331 /   ! eta'
!       ccbar
        data ( KF_mix(i,1,4), i=1,2,1 ) / 441, 443 /   ! eta_c, J/psi
!       bbbar
        data ( KF_mix(i,1,5), i=1,2,1 ) / 551, 553 /   ! eta_b, Upsilon

!       Proper probability of onia.
!       ddbar
        data ( prob_mix(i,1,1), i=1,2,1 ) / 0.5D0, 0.5D0 /
        data ( prob_mix(i,2,1), i=1,2,1 ) / 0.167D0, 0.5D0 /
        data ( prob_mix(i,3,1), i=1,2,1 ) / 0.333D0, 0.333D0 /
!       uubar
        data ( prob_mix(i,1,2), i=1,2,1 ) / 0.5D0, 0.5D0 /
        data ( prob_mix(i,2,2), i=1,2,1 ) / 0.167D0, 0.5D0 /
        data ( prob_mix(i,3,2), i=1,2,1 ) / 0.333D0, 0.333D0 /
!       ssbar
        data ( prob_mix(i,1,3), i=1,2,1 ) / 0.667D0, 1D0 /
        data ( prob_mix(i,2,3), i=1,2,1 ) / 0.333D0, 0.333D0 /
!       ccbar
        data ( prob_mix(i,1,4), i=1,2,1 ) / 1D0, 1D0 /
!       bbbar
        data ( prob_mix(i,1,5), i=1,2,1 ) / 1D0, 1D0 /
!=======================================================================


        i_succ = 0
        KF_h = 0

!       Incoming flavors check.
        if( KF1*KF2 >= 0 .OR. &
           ( ABS(KF1) > 5 .OR. ABS(KF2) > 5 ) )then
            write(*,*) "Wrong KF configuration in coales/findm_new:", &
                        KF1, KF2
            return
        end if

!       Assumes a meason comprises of KF_i and -KF_j ( KF_i, KF_j > 0 ).
        KF_i = ABS( KF1 )
        KF_j = ABS( KF2 )
        if( KF1 < 0 )then
            KF_i = ABS( KF2 )
            KF_j = ABS( KF1 )
        end if
        KF_max = MAX( KF_i, KF_j )
        KF_min = MIN( KF_i, KF_j )

        ! i_find = 10
        i_find = 11
        ! i_find = 20
        ! i_find = 22


!-------------------------------------------------------------------------------
!------------------------------   Normal Meson   -------------------------------
!       Case 1: x-ybar type meson. KF_i /= KF_j.
!       Proper probability: 1 & 1. Generates definitely.
        if( KF_i /= KF_j )then
            ! Spin 0: **1
            KF_h  = ( KF_max*100 + KF_min*10 + 1 )*SIGN( 1, KF_i - KF_j ) &
                  * (-1)**( KF_max )
            ! Spin 1: **3
            KF_h2 = ( KF_max*100 + KF_min*10 + 3 )*SIGN( 1, KF_i - KF_j ) &
                  * (-1)**( KF_max )
            ! The mass approaching method.
            if( i_find == 10 .OR. i_find == 20 )then
                ! Mass difference.
                d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                if( d_mass_2 < d_mass_1 ) KF_h = KF_h2
            ! The mass reciprocal method.
            else
                ! Mass ratio.
                d_mass_1 = PYMASS( KF_h  )
                d_mass_2 = PYMASS( KF_h2 )
                ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                if( PYR(1) < ratio_mass ) KF_h = KF_h2
            end if
            i_succ = 1
!------------------------------   Normal Meson   -------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------------   Onia   -----------------------------------
!       Case 2: x-xbar type meson, onium.
!       Mixing states: uubar, ddbar and ssbar. Heavy onia: ccbar, bbbar.
        else
!=======================================================================
!           ddbar: pi0, rho0.
!           uubar: eta, omega.
!           Mixing state      : pi0, rho0,   eta, omega,  eta'.
!           Proper probability: 0.5,  0.5, 0.167,   0.5, 0.333.
!           First pair : pi0 and rho.
!           Second pair: eta and omega.
!           Third one  : eta'.
!=======================================================================
!           ssbar: eta, phi; mixed with eta'.
!           Mixing state      :   eta, phi,  eta'.
!           Proper probability: 0.667,   1, 0.333.
!           First pair : eta and phi.
!           Second one : eta'.
!=======================================================================
!           ccbar: eta_c, J/psi.
!           Proper probability: 1 & 1. Generate definitely.
!=======================================================================
!           bbbar: eta_b, Upsilon.
!           Proper probability: 1 & 1. Generate definitely.
!=======================================================================


!------------------------   Proper Probability Method   ------------------------
!         Original method: proper probability.
          IF( i_find < 20 )THEN
            n_loop = n_mix_pair( KF_i )
            do j=1,n_loop,1
                i = 1
                KF_h  = KF_mix( 1, j, KF_i )
                KF_h2 = KF_mix( 2, j, KF_i )
                ! The mass approaching method.
                if( i_find == 10 )then
                    ! Mass difference.
                    d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                    d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                    if( d_mass_2 < d_mass_1 ) i = 2
                ! The mass reciprocal method.
                else
                    ! Mass ratio.
                    d_mass_1 = PYMASS( KF_h  )
                    d_mass_2 = PYMASS( KF_h2 )
                    ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                    if( PYR(1) < ratio_mass ) i = 2
                end if
                prob = prob_mix( i, j, KF_i )
                if( PYR(1) < prob )then
                    KF_h = KF_mix( i, j, KF_i )
                    i_succ = 1
                    return
                end if
                KF_h = 0
            end do
!------------------------   Proper Probability Method   ------------------------

!-----------------------   Relative Probability Method   -----------------------
!         Alternative method: relative probability.
          ELSE
            ! prob_tot = 2D0
            rand_numb = PYR(1)
            select case( KF_i )
            ! pi0, rho0, eta, emega and eta' mixing.
            case(1:2)
                ! prob_interval_0 = 0D0
                prob_pi_rho = (0.5D0 + 0.5D0) / 2D0
                prob_eta    = 0.167D0 / 2D0 + prob_pi_rho
                prob_omega  = 0.5D0 / 2D0 + prob_eta
                ! prob_eta_prim = 0.333D0 / 2D0 + prob_omega
                ! pio and rho0
                if( rand_numb <= prob_pi_rho )then
                    KF_h  = 111
                    KF_h2 = 113
                ! eta and omega
                else if( rand_numb <= prob_omega )then
                    KF_h  = 221
                    KF_h2 = 223
                ! eta'
                ! else if( rand_numb <= prob_eta_prim )then
                else
                    KF_h  = 331
                    KF_h2 = 331
                end if
            ! eta, phi and eta' mixing.
            case(3)
                prob_eta = 0.667D0 / 2D0
                prob_phi = 1D0 / 2D0 + prob_eta
                ! prob_eta_prim = 0.333 / 2. + prob_phi
                ! eta and phi
                if( rand_numb <= prob_phi )then
                    KF_h  = 221
                    KF_h2 = 333
                ! eta'
                ! else if( rand_numb <= prob_eta_prim )then
                else
                    KF_h  = 331
                    KF_h2 = 331
                end if
            ! eta_c and J/psi
            ! eta_b and Upsilon
            case(4:5)
                KF_h  = KF_i*100 + KF_i*10 + 1
                KF_h2 = KF_i*100 + KF_i*10 + 3
            end select
            ! The mass approaching method.
            if( i_find == 20 )then
                ! Mass difference.
                d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                if( d_mass_2 < d_mass_1 ) KF_h = KF_h2
            ! The mass reciprocal method.
            else
                ! Mass ratio.
                d_mass_1 = PYMASS( KF_h  )
                d_mass_2 = PYMASS( KF_h2 )
                ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                if( PYR(1) < ratio_mass ) KF_h = KF_h2
            end if
            i_succ = 1
          END IF
!-----------------------   Relative Probability Method   -----------------------

        end if
!----------------------------------   Onia   -----------------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findb_new( KF1, KF2, KF3, inv_mass_q, KF_h, i_succ )
!!      Finds the primary baryon (antibaryon) according to the rule
!        of KF code directly.
!=======================================================================
!       Original method, i_find = 10 / 11 :
!       The state is determined by
!         10: the proper probability and the mass approaching method.
!         11: the proper probability and the mass reciprocal method.
!       The rule of KF code cf. [ JHEP 05 (2006) 026, P73 ].
!       The proper probability is defined as the normalization factor 
!        of the hadronic wave function. Cf. "subroutine tabh" and 
!        [ H. Georgi, Lie Algebras In Particle Physics : from Isospin 
!          To Unified Theories (CRC Press, Taylor & Francis, Boca 
!          Raton, 2000) ].
!=======================================================================
!       An alternative method, i_find = 20 / 21 :
!         20: the relative probability and the mass approaching method.
!         21: the relative probability and the mass reciprocal method.
!=======================================================================
!       Assume KF_i >= KF_j >= KF_k > 0, total spin s.
!
!       KF of a baryon:
!           KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 2*s + 1
!
!       In particular, for spin 1/2 Lambda-like states with different
!        flavors of 3 quarks, reverse KF_j <-> KF_k, e.g. Lambda0 "3122"
!        while spin 1/2 Sigma-like states obey, e.g. Sigma0 "3212".
!
!       The baryons are divied into 3 cases:
!        (1) xyz type, with 3 different flavors, e.g. Lambda0 -- uds;
!        (2) xxy (xyy) type, with 2 same flavors, e.g. proton -- uud;
!        (2) xxx type, with 3 same flavors, e.g. Omega- -- sss.
!=======================================================================
!       Only primary lowest-lying spin 1/2 octets and spin 3/2 decuplets
!        are considered.
!        Spin 1/2, KF is ***2.
!        Spin 3/2, KF is ***4.
!       inv_mass_q: common invariant mass of incoming KF1, KF2 and KF3.
!       KF_h : flavor code of the found primary baryon.
!       i_succ = 1 : succeed,
!              = 0 : fail.
!=======================================================================
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real(kind=8) inv_mass_q


        i_succ = 0
        KF_h = 0

!       Incoming flavors check.
!       TODO(Lei20231125): extend to diquark + quark in the future?
        if( KF1*KF2*KF3 == 0 &
            .OR. ( ABS(KF1) > 5 .OR. ABS(KF2) > 5 .OR. ABS(KF3) > 5 ) &
            .OR. ( SIGN(1,KF1) /= SIGN(1,KF2) .OR. SIGN(1,KF1) /= SIGN(1,KF3) &
            .OR.   SIGN(1,KF2) /= SIGN(1,KF3) ) )then
            write(*,*) "Wrong KF configuration in coales/findb_new:", &
                        KF1, KF2, KF3
            return
        end if

!       Assume a baryon comprise of KF_i, KF_j and KF_k.
!       TODO(Lei20231125): extend to diquark + quark in the future?
        KF_i = ABS( KF1 )
        KF_j = ABS( KF2 )
        KF_k = ABS( KF3 )

        ! i_find = 10
        i_find = 11
        ! i_find = 20
        ! i_find = 22


!-------------------------------------------------------------------------------
!------------------------   Proper Probability Method   ------------------------
!   Original method: proper probability.
    IF( i_find < 20 )THEN

!       Case 1 & 2: xyz, xxy (xyy) type baryon.
        if( KF_i /= KF_j .OR. KF_i /= KF_k .OR. KF_j /= KF_k )then
            ! Sort KF_i, KF_j and KF_k from max to min (KF_i >= KF_j >= KF_k).
            if( KF_i < KF_j )then
                KF   = KF_j
                KF_j = KF_i
                KF_i = KF
            end if
            if( KF_j < KF_k )then
                KF   = KF_k
                KF_k = KF_j
                KF_j = KF
            end if
            if( KF_i < KF_j )then
                KF   = KF_j
                KF_j = KF_i
                KF_i = KF
            end if
!       Case 1: xyz type baryon.
            if( KF_i /= KF_j .AND. KF_i /= KF_k .AND. KF_j /= KF_k )then
                ! Spin 1/2 (Lambda, Xi_Q, Omega_Q): KF_i KF_k KF_j 2.
                ! Proper probability = 0.5.
                if( PYR(1) < 0.5D0 )then
                    KF_h = KF_i*1000 + KF_k*100 + KF_j*10 + 2
                    i_succ = 1
                    ! Anti-baryon.
                    if( KF1*KF2*KF3 < 0 ) KF_h = -KF_h
                    return
                end if
                ! Spin 3/2 (Sigma*, Xi*_Q, Omega*_Q): KF_i KF_j KF_k 4.
                ! Proper probability = 1.
                KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 4
                ! Spin 1/2 (Sigma, Xi'_Q, Omega'_Q): KF_i KF_j KF_k 2.
                ! Proper probability = 0.5.
                KF_h2 = KF_i*1000 + KF_j*100 + KF_k*10 + 2
                ! The mass approaching method.
                if( i_find == 10 )then
                    ! Mass difference.
                    d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                    d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                    if( d_mass_2 < d_mass_1 )then
                        KF_h = 0
                        ! Proper probability = 0.5.
                        if( PYR(1) > 0.5D0 ) return
                        KF_h = KF_h2
                    end if
                ! The mass reciprocal method.
                else
                    ! Mass ratio.
                    d_mass_1 = PYMASS( KF_h  )
                    d_mass_2 = PYMASS( KF_h2 )
                    ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                    if( PYR(1) < ratio_mass )then
                        KF_h = 0
                        ! Proper probability = 0.5.
                        if( PYR(1) > 0.5D0 ) return
                        KF_h = KF_h2
                    end if
                end if
!       Case 2: xxy (xyy) type baryon.
!       Proper probability = 1 & 1. Generates definitely.
            else
                ! Spin 1/2: KF_i KF_j KF_k 2.
                KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 2
                ! Spin 3/2: KF_i KF_j KF_k 4.
                KF_h2 = KF_i*1000 + KF_j*100 + KF_k*10 + 4
                ! The mass approaching method.
                if( i_find == 10 )then
                    ! Mass difference.
                    d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                    d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                    if( d_mass_2 < d_mass_1 ) KF_h = KF_h2
                ! The mass reciprocal method.
                else
                    ! Mass ratio.
                    d_mass_1 = PYMASS( KF_h  )
                    d_mass_2 = PYMASS( KF_h2 )
                    ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                    if( PYR(1) < ratio_mass ) KF_h = KF_h2
                end if
            end if
!       Case 3: xxx type baryon. Spin 3/2 only: xxx4.
!       Proper probability = 1. Generates definitely.
        else
            KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 4
        end if
!------------------------   Proper Probability Method   ------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Relative Probability Method   -----------------------
!   Alternative method: relative probability.
    ELSE

!       Case 1 & 2: xyz, xxy (xyy) type baryon.
        if( KF_i /= KF_j .OR. KF_i /= KF_k .OR. KF_j /= KF_k )then
            ! Sort KF_i, KF_j and KF_k from max to min (KF_i >= KF_j >= KF_k).
            if( KF_i < KF_j )then
                KF   = KF_j
                KF_j = KF_i
                KF_i = KF
            end if
            if( KF_j < KF_k )then
                KF   = KF_k
                KF_k = KF_j
                KF_j = KF
            end if
            if( KF_i < KF_j )then
                KF   = KF_j
                KF_j = KF_i
                KF_i = KF
            end if
!       Case 1: xyz type baryon.
            if( KF_i /= KF_j .AND. KF_i /= KF_k .AND. KF_j /= KF_k )then
                ! Spin 1/2 (Lambda, Xi_Q, Omega_Q): KF_i KF_k KF_j 2.
                ! Proper probability = 0.5. Only one state.
                if( PYR(1) <= 0.5D0 )then
                    KF_h = KF_i*1000 + KF_k*100 + KF_j*10 + 2
                ! Corresponding two states.
                else
                    ! Spin 1/2 (Sigma, Xi'_Q, Omega'_Q): KF_i KF_j KF_k 2.
                    KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 2
                    ! Spin 3/2 (Sigma*, Xi*_Q, Omega*_Q): KF_i KF_j KF_k 4.
                    KF_h2 = KF_i*1000 + KF_j*100 + KF_k*10 + 4
                    ! The mass approaching method.
                    if( i_find == 20 )then
                        ! Mass difference.
                        d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                        d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                        if( d_mass_2 < d_mass_1 ) KF_h = KF_h2
                    ! The mass reciprocal method.
                    else
                        ! Mass ratio.
                        d_mass_1 = PYMASS( KF_h  )
                        d_mass_2 = PYMASS( KF_h2 )
                        ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                        if( PYR(1) < ratio_mass ) KF_h = KF_h2
                    end if
                end if
!       Case 2: xxy (xyy) type baryon.
!       Proper probability = 1 & 1. Generates definitely.
            else
            ! Spin 1/2: KF_i KF_j KF_k 2.
                KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 2
            ! Spin 3/2: KF_i KF_j KF_k 4.
                KF_h2 = KF_i*1000 + KF_j*100 + KF_k*10 + 4
                ! The mass approaching method.
                if( i_find == 20 )then
                    ! Mass difference.
                    d_mass_1  = ABS( inv_mass_q - PYMASS( KF_h  ) )
                    d_mass_2  = ABS( inv_mass_q - PYMASS( KF_h2 ) )
                    if( d_mass_2 < d_mass_1 ) KF_h = KF_h2
                ! The mass reciprocal method.
                else
                    ! Mass ratio.
                    d_mass_1 = PYMASS( KF_h  )
                    d_mass_2 = PYMASS( KF_h2 )
                    ratio_mass = d_mass_1 / ( d_mass_1 + d_mass_2 )
                    if( PYR(1) < ratio_mass ) KF_h = KF_h2
                end if
            end if
!       Case 3: xxx type baryon. Spin 3/2 only: xxx4.
!       Proper probability = 1. Generate definitely.
        else
            KF_h = KF_i*1000 + KF_j*100 + KF_k*10 + 4
        end if

    END IF
!-----------------------   Relative Probability Method   -----------------------
!-------------------------------------------------------------------------------


        i_succ = 1
!       Anti-baryon.
        if( KF1*KF2*KF3 < 0 ) KF_h = -KF_h


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse_c( pp, n_column, n_raw, ps, npl, np, nstep )
!!      Adjusts four momentum conservation by iteration, no more than
!!       5000 iterations.
!       pp : four momentum of particles
!       ps : the four momentum should be conserved to
!       npl : order number of the first particle
!       np : order number of the last particle
!       nstep : interval of the step
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        dimension pp( n_column, n_raw ), ff(n_column)
        dimension ps(4), pxyz(3), arp(3)


        ps4  = ps(4)
        pxyz = 0D0

        jj = 0
100     es = 0D0
        do i = npl, np, nstep
            es = es + pp(i,4)
        end do
        fr = es / ps4

        do i = npl, np, nstep
            amas = pp(i,5)
            amas2 = amas * amas
            ppm = pp(i,4)
            ppf = ppm / fr
            den2 = ppm*ppm - amas2
            if( den2 <= 0D0 )then
                den2 = 1D-15
            end if
            ff(i) = SQRT( ABS( ppf*ppf - amas2 ) / den2 )
            do j=1,3,1
                ppp = ff(i)*pp(i,j)
                pp(i,j) = ppp
                pxyz(j) = pxyz(j) + ppp
            end do
        end do
        do i=1,3,1
            arp(i) = ABS( 1D0 - pxyz(i) / ps(i) )
        end do
        if( ABS( 1D0 - fr ) <= dep .AND. arp(1) <= dep .AND. arp(2) <= dep &
            .AND. arp(3) <= dep ) return
            do i=1,3,1
                pxyz(i) = pxyz(i) - ps(i)
                pxyz(i) = pxyz(i) / ( DBLE( np - npl ) / DBLE( nstep ) + 1D0 )
            end do
            do i = npl, np, nstep
            do j=1,3,1
                pp(i,j) = pp(i,j) - pxyz(j)
            end do
            pp5 = pp(i,5)
            pp52 = pp5*pp5
            pp(i,4) = SQRT( pp52 + pp(i,1)**2 + pp(i,2)**2 + pp(i,3)**2 )
        end do
        jj = jj+1
        if( jj < 5000 ) goto 100


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine deexcitation_W( i_parton_in, N_STEP )
!!      New deexcitation subroutine.
!!      Energetic q (qbar) de-exciatation.
!!      qqbar pair generation according to the light-cone momentum.
!
!       Sets the direction of the momentum of the jet (mother quark)
!        as +z (rotation, hence pz* = |p|, px* = py* = 0, star* means
!        the quantities after rotation) and excites qqbar.
!        After the deexciatation, the jet and qqbars will be rotated
!        back to the direction of the original jet. (For only one jet,
!        we don't need to boost it to its CMS.)
!
!       Assumes local pT compensation, i.e. px(-px) and py(-py) for
!        the exited q(qbar).
!
!       i_parton_in: the order number of the source quark (or anti-quark).
!       N_STEP: number of deexcitation steps per source q (qbar) occurred.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        real(kind=8) rr(3)
        real(kind=8) P00(5), V00(5)
        real(kind=8) M00, M0, M1, M2
        real(kind=8) lost_4mom(4)


!       Threshold energy.
        W_threshold = adj1(17)
        ! Kinemtical mass of uubar ~ 2*0.33 = 0.66 GeV.
        dm_uubar = 2D0 * amass( 2 )
        if( W_threshold < dm_uubar ) W_threshold = dm_uubar
!       Upper limit of deexcitation steps (times) per incoming q (qbar).
        ! n_deex_step
!       Width of pT sampling.
        i_pT_deex = i_pT_coal
        pT_width  = adj1(35)
!       Deexcitation function for selecting z.
        i_function_z = INT( adj1(29) )
!       Backup.
        PARJ41 = PARJ(41)
        PARJ42 = PARJ(42)
        PARJ51 = PARJ(51)
        PARJ52 = PARJ(52)
        PARJ53 = PARJ(53)
        PARJ54 = PARJ(54)
        PARJ55 = PARJ(55)
        PARJ(41) = adj1(8)
        PARJ(42) = adj1(9)
        PARJ(51) = a_FF
        PARJ(52) = a_FF
        PARJ(53) = a_FF
        PARJ(54) = -aPS_c
        PARJ(55) = -aPS_b


!-------------------------------------------------------------------------------
!---------------------   Original Information Recording   ----------------------
!       Original number of entries of "PYJETS".
        N00 = N
!       Original line number of the source q (qbar).
        i_parton00 = i_parton_in
        KF00 = K( i_parton00, 2 )
        do i=1,5,1
            P00( i ) = P( i_parton00, i )
            V00( i ) = V( i_parton00, i )
        end do
        px00  = P00(1)
        py00  = P00(2)
        pz00  = P00(3)
        E00   = P00(4)
        M00   = P00(5)
        vx00  = V00(1)
        vy00  = V00(2)
        vz00  = V00(3)
        t00   = V00(4)
        tau00 = V00(5)
!---------------------   Original Information Recording   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------   Source Parton Rotating   --------------------------
        N_TRY_PT0 = 0
100     continue
        N_TRY_PT0 = N_TRY_PT0 + 1
!       No more than 50 times.
        if( N_TRY_PT0 > 50 )then
            write(22,*) "Warning, N_TRY_PT0 > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
        i_parton0 = i_parton00
        KF0 = KF00
        px0 = 0D0
        py0 = 0D0
        pT0 = 0D0
!       If necessary, resamples px and py of the mother q (qbar).
        if( i_pT_endpoint == 1 ) &
            call PAPTDI( i_pT_deex, pT_width, px0, py0, pT0 )
!       Rotates the mother jet to +z, i.e. px* = py* = 0, pz* = |p|.
!        Hence W* = E + pz* = E + |p|. For only one particle, do not
!        boost it to the CM frame.
        pz0 = SQRT( px00**2 + py00**2 + pz00**2 )
        E0  = E00
        M0  = M00
        W0  = E00 + SQRT( px00**2 + py00**2 + pz00**2 )
        dmT2_0 = M0**2 + px0**2 + py0**2
!       Too low energy, exits.
        if( W0 <= W_threshold ) goto 999
!-------------------------   Source Parton Rotating   --------------------------
!-------------------------------------------------------------------------------


!       Excites qqbar from vacuum.


        N_STEP = 0
200     continue
!       Counts number of deexcitation steps (times).
        N_STEP = N_STEP + 1
!       Reached the maximum number of deexcitation steps. Stops deexcitation.
        if( N_STEP > n_deex_step ) goto 999
!       Too low energy to excite a qqbar.
        if( E0 <= dm_uubar ) goto 999


!-------------------------------------------------------------------------------
!-------------------------   Excited Flavor Sampling   -------------------------
        N_TRY_KF = 0
300     continue
        N_TRY_KF = N_TRY_KF + 1
!       No more than 50 times.
!       Special treatment for the source parton.
        if( N_TRY_KF >= 50 .AND. N_STEP == 1 .AND. i_pT_endpoint == 1 )then
            goto 100
!       Else failed, stops deexcitation.
        else if( N_TRY_KF > 50 )then
            write(22,*) "Warning, N_TRY_KF > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
!       Generates flavors of the excited qqbar.
        E = E0
        call break_f( E, KF, M1 )
        if( KF == 0 ) goto 999
!       Finds the flavors and masses of the excited qqbar.
        if( PYR(1) > 0.5D0 ) KF = -KF
        KF1 =  KF
        KF2 = -KF
        M1  = amass( KF )
        M2  = M1
!       Regards qqbar as an object.
        dm_qqbar = M1 + M2
!-------------------------   Excited Flavor Sampling   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pT Sampling   ---------------------------
        N_TRY_PT = 0
400     continue
        N_TRY_PT = N_TRY_PT + 1
!       No more than 50 times.
        if( N_TRY_PT > 50 ) goto 300
!       Generates compensated px py for excited qqbar.
        call PAPTDI( i_pT_deex, pT_width, px1, py1, pT1 )
!       Local pT compensation.
        px2 = -px1
        py2 = -py1
        dmT2_1 = M1**2 + px1**2 + py1**2
        dmT2_2 = M2**2 + px2**2 + py2**2
!       Regards qqbar as an object.
        px_qqbar   = px1 + px2
        py_qqbar   = py1 + py2
        dmT2_qqbar = dm_qqbar**2 + px_qqbar**2 + py_qqbar**2
!---------------------------   Excited pT Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pz Sampling   ---------------------------
        N_TRY_Z = 0
        loop_Z: do while(.true.)
            N_TRY_Z = N_TRY_Z + 1
!       No more than 50 times.
            if( N_TRY_Z > 50 ) goto 400
            ! In Coales_40.f90.
            call PAZDIS( i_function_z, 0D0, KF1, KF2, dmT2_qqbar, Z )
            W_qqbar  = W0 * Z
            pz_qqbar = 0.5D0*( W_qqbar - dmT2_qqbar / MAX( 1D-4, W_qqbar ) )
            W_0  = W0 - W_qqbar
            pz_0 =  0.5D0*( W_0 - dmT2_0 / MAX( 1D-4, W_0 ) )
            E_0  =  0.5D0*( W_0 + dmT2_0 / MAX( 1D-4, W_0 ) )
!       Checks if pz acceptable.
            if( pz_qqbar > 0D0 .AND. pz_0 > 0D0 )then
                N_TRY_Z1 = 0
                do while(.true.)
                    N_TRY_Z1 = N_TRY_Z1 + 1
!       No more than 50 times.
                    if( N_TRY_Z1 > 50 ) cycle loop_Z
                    Z1  = PYR(1)
                    W1  = W_qqbar * Z1
                    W2  = W_qqbar - W1
                    pz1 = 0.5D0*( W1 - dmT2_1 / MAX( 1D-4, W1 ) )
                    E1  = 0.5D0*( W1 + dmT2_1 / MAX( 1D-4, W1 ) )
                    pz2 = 0.5D0*( W2 - dmT2_2 / MAX( 1D-4, W2 ) )
                    E2  = 0.5D0*( W2 + dmT2_2 / MAX( 1D-4, W2 ) )
                    if( pz1 > 0D0 .AND. pz2 > 0D0 ) exit loop_Z
                end do
            end if
        end do loop_Z
        W0  = W_0
        pz0 = pz_0
        E0  = E_0
!---------------------------   Excited pz Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Properties Giving   ----------------------------
!       Gives four-momentum and mass to the mother q/qbar and the excited qqbar.
!       Mother.
        P( i_parton0, 1 ) = px0
        P( i_parton0, 2 ) = py0
        P( i_parton0, 3 ) = pz0
        P( i_parton0, 4 ) = E0
        P( i_parton0, 5 ) = M0
!       1st-excited.
        P( N+1, 1 ) = px1
        P( N+1, 2 ) = py1
        P( N+1, 3 ) = pz1
        P( N+1, 4 ) = E1
        P( N+1, 5 ) = M1
!       2nd-excited.
        P( N+2, 1 ) = px2
        P( N+2, 2 ) = py2
        P( N+2, 3 ) = pz2
        P( N+2, 4 ) = E2
        P( N+2, 5 ) = M2

!       Gives other status.
!       Mother.
        K( i_parton0, 4 ) = N + 1
        K( i_parton0, 5 ) = N + 2
!       1st-excited.
        K( N+1, 1 ) = 2
        K( N+1, 2 ) = KF1
        K( N+1, 3 ) = i_parton00
        K( N+1, 4 ) = 0
        K( N+1, 5 ) = 0
!       2nd-excited.
        K( N+2, 1 ) = 1
        K( N+2, 2 ) = KF2
        K( N+2, 3 ) = i_parton00
        K( N+2, 4 ) = 0
        K( N+2, 5 ) = 0
        N = N + 2
!----------------------------   Properties Giving   ----------------------------
!-------------------------------------------------------------------------------


!       Goes back for new qqbar if enough energy.
        if( W0 > W_threshold ) goto 200

999     continue
        if( N == N00 ) goto 1000
!       First daughter of the original parton.
        K( i_parton00, 4 ) = N00 + 1
!       The last daughter of the original parton.
        K( i_parton00, 5 ) = N


!-------------------------------------------------------------------------------
!---------------------   Rotating Back & Position Giving  ----------------------
!       Appends the momenta of the mother quark to "N+1" for the rotation.
        do i=1,5,1
            K( N+1, i ) = K( i_parton00, i )
            P( N+1, i ) = P( i_parton00, i )
            V( N+1, i ) = V( i_parton00, i )
        end do
        K( N+1, 1 ) = 1

!       Rotation.
        theta = PYANGL( P00(3), SQRT( P00(1)**2 + P00(2)**2 ) )
        phi   = PYANGL( P00(1), P00(2) )
!       The positions will be rotated too. Do not use them.
        ! MSTU(33) = 1
        call PYROBO( N00+1, N+1, theta, phi, 0D0, 0D0, 0D0 )

!       Gives three-position and time to excited qqbar.
!       Generated q and qbar are arranged around the source parton
!        within 0-0.5 fm randomly in each one of the three positions
!        and has the same time as source parton + formation time.
        do i = N00+1, N, 1
            do j=1,3
                rr(j) = PYR(1)*0.5D0
                V(i,j) = V00(j) + rr(j)
                if( PYR(1) > 0.5D0 ) V(i,j) = V00(j) - rr(j)
            end do
            V(i,4) = V00(4)
            V(i,5) = 0D0
!TODO(Lei20241016): gives formation time.
            ! IF( .true. )THEN
            IF( .false. )THEN
                i_tau = 1
!               Formation time method 1.
!               Proper formation time of t0/10 was assumed for partons.
                if( i_tau == 1 )then
                    tau0 = t0 / 10D0
                    V(i,4) = V00(4) + tau0 * P(i,4) / P(i,5)
                    V(i,5) = 0D0
!               Formation time method 2, tau = E/mT^2.
                else if( i_tau == 2 )then
                    V(i,4) = V00(4) &
                           + P(i,4) / (P(i,1)**2 + P(i,2)**2 + P(i,5)**2) &
                           * 0.197327D0
                end if
            END IF
        end do

!       Returns the momenta of the mother q/qbar from "N+1" back.
        do i=1,5,1
            P( i_parton00, i ) = P( N+1, i )
        end do
        do i=1,5,1
            V( i_parton00, i ) = V00( i )
        end do
!---------------------   Rotating Back & Position Giving  ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!       Collects the difference (lost) of four-momentum before and
!        after deexcitation, then tries to shares them.
        lost_4mom = 0D0
        do j=1,4,1
            do i= N00+1, N, 1
                lost_4mom(j) = lost_4mom(j) + P(i,j)
            end do
            lost_4mom(j) = lost_4mom(j) + P( i_parton00, j )
            lost_4mom(j) = P00(j) - lost_4mom(j)
        end do
!       Collects the rest four-momentum to "throe_p".
        throe_p = throe_p + lost_4mom
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!-------------------------------------------------------------------------------


1000    continue
!       Recovers.
        PARJ(41) = PARJ41
        PARJ(42) = PARJ42
        PARJ(51) = PARJ51
        PARJ(52) = PARJ52
        PARJ(53) = PARJ53
        PARJ(54) = PARJ54
        PARJ(55) = PARJ55


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine deexcitation_E( i_parton_in, N_STEP )
!!      New deexcitation subroutine.
!!      Energetic q (qbar) de-exciatation.
!!      qqbar pair generation according to the energy.
!
!       Sets the direction of the momentum of the jet (mother quark)
!        as +z (rotation, hence pz* = |p|, px* = py* = 0, star* means
!        the quantities after rotation) and excites qqbar.
!        After the deexciatation, the qqbars will be rotated
!        back to the direction of the original jet. (For only one jet,
!        we don't need to boost it to its CMS.)
!
!       Assumes local pT compensation, i.e. px(-px) and py(-py) for
!        the exited q(qbar).
!
!       i_parton_in: the order number of the source quark (or anti-quark).
!       N_STEP: number of deexcitation steps per source q (qbar) occurred.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        real(kind=8) rr(3)
        real(kind=8) P00(5), V00(5)
        real(kind=8) M00, M0, M1, M2
        real(kind=8) lost_4mom(4)


!       Threshold energy.
        E_threshold = adj1(17)
        ! Kinemtical mass of uubar ~ 2*0.33 = 0.66 GeV.
        dm_uubar = 2D0 * amass( 2 )
        if( E_threshold < dm_uubar ) E_threshold = dm_uubar
!       Upper limit of deexcitation steps (times) per incoming q (qbar).
        ! n_deex_step
!       Width of pT sampling.
        i_pT_deex = i_pT_coal
        pT_width  = adj1(35)
!       Deexcitation function for selecting z.
        i_function_z = INT( adj1(29) )
!       Backup.
        PARJ41 = PARJ(41)
        PARJ42 = PARJ(42)
        PARJ51 = PARJ(51)
        PARJ52 = PARJ(52)
        PARJ53 = PARJ(53)
        PARJ54 = PARJ(54)
        PARJ55 = PARJ(55)
        PARJ(41) = adj1(8)
        PARJ(42) = adj1(9)
        PARJ(51) = a_FF
        PARJ(52) = a_FF
        PARJ(53) = a_FF
        PARJ(54) = -aPS_c
        PARJ(55) = -aPS_b


!-------------------------------------------------------------------------------
!---------------------   Original Information Recording   ----------------------
!       Original number of entries of "PYJETS".
        N00 = N
!       Original line number of the source q (qbar).
        i_parton00 = i_parton_in
        KF00 = K( i_parton00, 2 )
        do i=1,5,1
            P00( i ) = P( i_parton00, i )
            V00( i ) = V( i_parton00, i )
        end do
        px00  = P00(1)
        py00  = P00(2)
        pz00  = P00(3)
        E00   = P00(4)
        M00   = P00(5)
        vx00  = V00(1)
        vy00  = V00(2)
        vz00  = V00(3)
        t00   = V00(4)
        tau00 = V00(5)
!---------------------   Original Information Recording   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------   Source Parton Rotating   --------------------------
        N_TRY_PT0 = 0
100     continue
        N_TRY_PT0 = N_TRY_PT0 + 1
!       No more than 50 times.
        if( N_TRY_PT0 > 50 )then
            write(22,*) "Warning, N_TRY_PT0 > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
        i_parton0 = i_parton00
        KF0 = KF00
        px0 = 0D0
        py0 = 0D0
        pT0 = 0D0
!       If necessary, resamples px and py of the mother q (qbar).
        if( i_pT_endpoint == 1 ) &
            call PAPTDI( i_pT_deex, pT_width, px0, py0, pT0 )
!       Rotates the mother jet to +z, i.e. px* = py* = 0, pz* = |p|.
!        Hence W* = E + pz* = E + |p|. For only one particle, do not
!        boost it to the CM frame.
        pz0 = SQRT( px00**2 + py00**2 + pz00**2 )
        E0  = E00
        M0  = M00
        ! W0  = E00 + SQRT( px00**2 + py00**2 + pz00**2 )
        dmT2_0 = M0**2 + px0**2 + py0**2
!       Too low energy, exits.
        if( E0 <= E_threshold ) goto 999
!-------------------------   Source Parton Rotating   --------------------------
!-------------------------------------------------------------------------------


!       Excites qqbar from vacuum.


        N_STEP = 0
200     continue
!       Counts number of deexcitation steps (times).
        N_STEP = N_STEP + 1
!       Reached the maximum number of deexcitation steps. Stops deexcitation.
        if( N_STEP > n_deex_step ) goto 999
!       Too low energy to excite a qqbar.
        if( E0 <= dm_uubar ) goto 999


!-------------------------------------------------------------------------------
!-------------------------   Excited Flavor Sampling   -------------------------
        N_TRY_KF = 0
300     continue
        N_TRY_KF = N_TRY_KF + 1
!       No more than 50 times.
!       Special treatment for the source parton.
        if( N_TRY_KF >= 50 .AND. N_STEP == 1 .AND. i_pT_endpoint == 1 )then
            goto 100
!       Else failed, stops deexcitation.
        else if( N_TRY_KF > 50 )then
            write(22,*) "Warning, N_TRY_KF > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
!       Generates flavors of the excited qqbar.
        E = E0
        call break_f( E, KF, M1 )
        if( KF == 0 ) goto 999
!       Finds the flavors and masses of the excited qqbar.
        if( PYR(1) > 0.5D0 ) KF = -KF
        KF1 =  KF
        KF2 = -KF
        M1  = amass( KF )
        M2  = M1
!       Regards qqbar as an object (alternative).
        dm_qqbar = M1 + M2
!-------------------------   Excited Flavor Sampling   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pT Sampling   ---------------------------
        N_TRY_PT = 0
400     continue
        N_TRY_PT = N_TRY_PT + 1
!       No more than 50 times.
        if( N_TRY_PT > 50 ) goto 300
!       Generates compensated px py for excited qqbar.
        call PAPTDI( i_pT_deex, pT_width, px1, py1, pT1 )
!       Local pT compensation.
        px2 = -px1
        py2 = -py1
        dmT2_1 = M1**2 + px1**2 + py1**2
        dmT2_2 = M2**2 + px2**2 + py2**2
!       Regards qqbar as an object.
        px_qqbar   = px1 + px2
        py_qqbar   = py1 + py2
        dmT2_qqbar = dm_qqbar**2 + px_qqbar**2 + py_qqbar**2
!---------------------------   Excited pT Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pz Sampling   ---------------------------
        N_TRY_Z = 0
        loop_Z: do while(.true.)
            N_TRY_Z = N_TRY_Z + 1
!       No more than 50 times.
            if( N_TRY_Z > 50 ) goto 400
            ! In Coales_40.f90.
            call PAZDIS( i_function_z, 0D0, KF1, KF2, dmT2_qqbar, Z )
            E_qqbar = E0 * Z
            E_0     = E0 - E_qqbar
!       Checks if E acceptable.
            if( E_qqbar**2 > dmT2_qqbar .AND. E_0**2 > dmT2_0 )then
                N_TRY_Z1  = 0
                do while(.true.)
                    N_TRY_Z1 = N_TRY_Z1 + 1
!       No more than 50 times.
                    if( N_TRY_Z1 > 50 ) cycle loop_Z
                    Z1 = PYR(1)
                    E1 = E_qqbar * Z1
                    E2 = E_qqbar - E1
                    if( E1**2 > dmT2_1 .AND. E2**2 > dmT2_2 )then
                        pz1 = SQRT( E1**2 - dmT2_1 )
                        pz2 = SQRT( E2**2 - dmT2_2 )
                        exit loop_Z
                    end if
                end do
            end if
        end do loop_Z
        pz0 = SQRT( E_0**2 - dmT2_0 )
        E0  = E_0
!---------------------------   Excited pz Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Properties Giving   ----------------------------
!       Gives four-momentum and mass to the mother q/qbar and the excited qqbar.
!       Mother.
        P( i_parton0, 1 ) = px0
        P( i_parton0, 2 ) = py0
        P( i_parton0, 3 ) = pz0
        P( i_parton0, 4 ) = E0
        P( i_parton0, 5 ) = M0
!       1st-excited.
        P( N+1, 1 ) = px1
        P( N+1, 2 ) = py1
        P( N+1, 3 ) = pz1
        P( N+1, 4 ) = E1
        P( N+1, 5 ) = M1
!       2nd-excited.
        P( N+2, 1 ) = px2
        P( N+2, 2 ) = py2
        P( N+2, 3 ) = pz2
        P( N+2, 4 ) = E2
        P( N+2, 5 ) = M2

!       Gives other status.
!       Mother.
        K( i_parton0, 4 ) = N + 1
        K( i_parton0, 5 ) = N + 2
!       1st-excited.
        K( N+1, 1 ) = 2
        K( N+1, 2 ) = KF1
        K( N+1, 3 ) = i_parton00
        K( N+1, 4 ) = 0
        K( N+1, 5 ) = 0
!       2nd-excited.
        K( N+2, 1 ) = 1
        K( N+2, 2 ) = KF2
        K( N+2, 3 ) = i_parton00
        K( N+2, 4 ) = 0
        K( N+2, 5 ) = 0
        N = N + 2
!----------------------------   Properties Giving   ----------------------------
!-------------------------------------------------------------------------------


!       Goes back for new qqbar if enough energy.
        if( E0 > E_threshold ) goto 200

999     continue
        if( N == N00 ) goto 1000
!       First daughter of the original parton.
        K( i_parton00, 4 ) = N00 + 1
!       The last daughter of the original parton.
        K( i_parton00, 5 ) = N


!-------------------------------------------------------------------------------
!---------------------   Rotating Back & Position Giving  ----------------------
!       Appends the momenta of the mother quark to "N+1" for the rotation.
        do i=1,5,1
            K( N+1, i ) = K( i_parton00, i )
            P( N+1, i ) = P( i_parton00, i )
            V( N+1, i ) = V( i_parton00, i )
        end do
        K( N+1, i ) = 1

!       Rotation.
        theta = PYANGL( P00(3), SQRT( P00(1)**2 + P00(2)**2 ) )
        phi   = PYANGL( P00(1), P00(2) )
!       The positions will be rotated too. Do not use them.
        ! MSTU(33) = 1
        call PYROBO( N00+1, N+1, theta, phi, 0D0, 0D0, 0D0 )

!       Gives three-position and time to excited qqbar.
!       Generated q and qbar are arranged around the source parton
!        within 0-0.5 fm randomly in each one of the three positions
!        and has the same time as source parton + formation time.
        do i = N00+1, N, 1
            do j=1,3
                rr(j) = PYR(1)*0.5D0
                V(i,j) = V00(j) + rr(j)
                if( PYR(1) > 0.5D0 ) V(i,j) = V00(j) - rr(j)
            end do
            V(i,4) = V00(4)
            V(i,5) = 0D0
!TODO(Lei20241016): gives formation time.
            ! IF( .true. )THEN
            IF( .false. )THEN
                i_tau = 1
!               Formation time method 1.
!               Proper formation time of t0/10 was assumed for partons.
                if( i_tau == 1 )then
                    tau0 = t0 / 10D0
                    V(i,4) = V00(4) + tau0 * P(i,4) / P(i,5)
                    V(i,5) = 0D0
!               Formation time method 2, tau = E/mT^2.
                else if( i_tau == 2 )then
                    V(i,4) = V00(4) &
                           + P(i,4) / (P(i,1)**2 + P(i,2)**2 + P(i,5)**2) &
                           * 0.197327D0
                end if
            END IF
        end do

!       Returns the momenta of the mother quark from "N+1" back.
        do i=1,5,1
            P( i_parton00, i ) = P( N+1, i )
        end do
        do i=1,5,1
            V( i_parton00, i ) = V00( i )
        end do
!---------------------   Rotating Back & Position Giving  ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!       Collects the difference (lost) of four-momentum before and
!        after deexcitation, then tries to shares them.
        lost_4mom = 0D0
        do j=1,4,1
            do i= N00+1, N, 1
                lost_4mom(j) = lost_4mom(j) + P(i,j)
            end do
            lost_4mom(j) = lost_4mom(j) + P( i_parton00, j )
            lost_4mom(j) = P00(j) - lost_4mom(j)
        end do
!       Collects the rest four-momentum to "throe_p".
        throe_p = throe_p + lost_4mom
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!-------------------------------------------------------------------------------


1000    continue
!       Recovers.
        PARJ(41) = PARJ41
        PARJ(42) = PARJ42
        PARJ(51) = PARJ51
        PARJ(52) = PARJ52
        PARJ(53) = PARJ53
        PARJ(54) = PARJ54
        PARJ(55) = PARJ55


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine deexcitation_PZ( i_parton_in, N_STEP )
!!      New deexcitation subroutine.
!!      Energetic q (qbar) de-exciatation.
!!      qqbar pair generation according to the pz (energy).
!
!       Sets the direction of the momentum of the jet (mother quark)
!        as +z (rotation, hence pz* = |p|, px* = py* = 0, star* means
!        the quantities after rotation) and excites qqbar.
!        After the deexciatation, the qqbars will be rotated
!        back to the direction of the original jet. (For only one jet,
!        we don't need to boost it to its CMS.)
!
!       Assumes local pT compensation, i.e. px(-px) and py(-py) for
!        the exited q(qbar).
!
!       i_parton_in: the order number of the source quark (or anti-quark).
!       N_STEP: number of deexcitation steps per source q (qbar) occurred.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        real(kind=8) rr(3)
        real(kind=8) P00(5), V00(5)
        real(kind=8) M00, M0, M1, M2
        real(kind=8) lost_4mom(4)


!       Threshold energy.
        PZ_threshold = adj1(17)
        ! Kinemtical mass of uubar ~ 2*0.33 = 0.66 GeV.
        dm_uubar = 2D0 * amass( 2 )
        if( PZ_threshold < dm_uubar ) PZ_threshold = dm_uubar
!       Upper limit of deexcitation steps (times) per incoming q (qbar).
        ! n_deex_step
!       Width of pT sampling.
        i_pT_deex = i_pT_coal
        pT_width  = adj1(35)
!       Deexcitation function for selecting z.
        i_function_z = INT( adj1(29) )
!       Backup.
        PARJ41 = PARJ(41)
        PARJ42 = PARJ(42)
        PARJ51 = PARJ(51)
        PARJ52 = PARJ(52)
        PARJ53 = PARJ(53)
        PARJ54 = PARJ(54)
        PARJ55 = PARJ(55)
        PARJ(41) = adj1(8)
        PARJ(42) = adj1(9)
        PARJ(51) = a_FF
        PARJ(52) = a_FF
        PARJ(53) = a_FF
        PARJ(54) = -aPS_c
        PARJ(55) = -aPS_b


!-------------------------------------------------------------------------------
!---------------------   Original Information Recording   ----------------------
!       Original number of entries of "PYJETS".
        N00 = N
!       Original line number of the source q (qbar).
        i_parton00 = i_parton_in
        KF00 = K( i_parton00, 2 )
        do i=1,5,1
            P00( i ) = P( i_parton00, i )
            V00( i ) = V( i_parton00, i )
        end do
        px00  = P00(1)
        py00  = P00(2)
        pz00  = P00(3)
        E00   = P00(4)
        M00   = P00(5)
        vx00  = V00(1)
        vy00  = V00(2)
        vz00  = V00(3)
        t00   = V00(4)
        tau00 = V00(5)
!---------------------   Original Information Recording   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------   Source Parton Rotating   --------------------------
        N_TRY_PT0 = 0
100     continue
        N_TRY_PT0 = N_TRY_PT0 + 1
!       No more than 50 times.
        if( N_TRY_PT0 > 50 )then
            write(22,*) "Warning, N_TRY_PT0 > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
        i_parton0 = i_parton00
        KF0 = KF00
        px0 = 0D0
        py0 = 0D0
        pT0 = 0D0
!       If necessary, resamples px and py of the mother q (qbar).
        if( i_pT_endpoint == 1 ) &
            call PAPTDI( i_pT_deex, pT_width, px0, py0, pT0 )
!       Rotates the mother jet to +z, i.e. px* = py* = 0, pz* = |p|.
!        Hence W* = E + pz* = E + |p|. For only one particle, do not
!        boost it to the CM frame.
        pz0 = SQRT( px00**2 + py00**2 + pz00**2 )
        E0  = E00
        M0  = M00
        ! W0  = E00 + SQRT( px00**2 + py00**2 + pz00**2 )
        dmT2_0 = M0**2 + px0**2 + py0**2
!       Too low energy, exits.
        if( pz0 <= PZ_threshold ) goto 999
!-------------------------   Source Parton Rotating   --------------------------
!-------------------------------------------------------------------------------


!       Excites qqbar from vacuum.


        N_STEP = 0
200     continue
!       Counts number of deexcitation steps (times).
        N_STEP = N_STEP + 1
!       Reached the maximum number of deexcitation steps. Stops deexcitation.
        if( N_STEP > n_deex_step ) goto 999
!       Too low energy to excite a qqbar.
        if( E0 < dm_uubar ) goto 999


!-------------------------------------------------------------------------------
!-------------------------   Excited Flavor Sampling   -------------------------
        N_TRY_KF = 0
300     continue
        N_TRY_KF = N_TRY_KF + 1
!       No more than 50 times.
!       Special treatment for the source parton.
        if( N_TRY_KF >= 50 .AND. N_STEP == 1 .AND. i_pT_endpoint == 1 )then
            goto 100
!       Else failed, stops deexcitation.
        else if( N_TRY_KF > 50 )then
            write(22,*) "Warning, N_TRY_KF > 50 in deexcitation, " &
                     // "i-event, i-entry, KF, px py pz e m =", &
                        iii, i_parton00, KF00, ( P00(i), i=1, 5, 1 )
!       Failed, exits.
            goto 999
        end if
!       Generates flavors of the excited qqbar.
        E = E0
        call break_f( E, KF, M1 )
        if( KF == 0 ) goto 999
!       Finds the flavors and masses of the excited qqbar.
        if( PYR(1) > 0.5D0 ) KF = -KF
        KF1 =  KF
        KF2 = -KF
        M1  = amass( KF )
        M2  = M1
!       Regards qqbar as an object (alternative).
        dm_qqbar = M1 + M2
!-------------------------   Excited Flavor Sampling   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pT Sampling   ---------------------------
        N_TRY_PT = 0
400     continue
        N_TRY_PT = N_TRY_PT + 1
!       No more than 50 times.
        if( N_TRY_PT > 50 ) goto 300
!       Generates compensated px py for excited qqbar.
        call PAPTDI( i_pT_deex, pT_width, px1, py1, pT1 )
!       Local pT compensation.
        px2 = -px1
        py2 = -py1
        dmT2_1 = M1**2 + px1**2 + py1**2
        dmT2_2 = M2**2 + px2**2 + py2**2
!       Regards qqbar as an object.
        px_qqbar   = px1 + px2
        py_qqbar   = py1 + py2
        dmT2_qqbar = dm_qqbar**2 + px_qqbar**2 + py_qqbar**2
!---------------------------   Excited pT Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!---------------------------   Excited pz Sampling   ---------------------------
        N_TRY_Z = 0
        loop_Z: do while(.true.)
            N_TRY_Z = N_TRY_Z + 1
!       No more than 50 times.
            if( N_TRY_Z > 50 ) goto 400
            ! In Coales_40.f90.
            call PAZDIS( i_function_z, 0D0, KF1, KF2, dmT2_qqbar, Z )
            pz_qqbar = pz0 * Z
            pz0      = pz0 - pz_qqbar
            Z1  = PYR(1)
            pz1 = pz_qqbar * Z1
            pz2 = pz_qqbar - pz1
            exit
        end do loop_Z
!       Recalculates energies of the mother parton and excited qqbar.
        E0 = SQRT( M0**2 + px0**2 + py0**2 + pz0**2 )
        E1 = SQRT( M1**2 + px1**2 + py1**2 + pz1**2 )
        E2 = SQRT( M2**2 + px2**2 + py2**2 + pz2**2 )
!---------------------------   Excited pz Sampling   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Properties Giving   ----------------------------
!       Gives four-momentum and mass to the mother q/qbar and the excited qqbar.
!       Mother.
        P( i_parton0, 1 ) = px0
        P( i_parton0, 2 ) = py0
        P( i_parton0, 3 ) = pz0
        P( i_parton0, 4 ) = E0
        P( i_parton0, 5 ) = M0
!       1st-excited.
        P( N+1, 1 ) = px1
        P( N+1, 2 ) = py1
        P( N+1, 3 ) = pz1
        P( N+1, 4 ) = E1
        P( N+1, 5 ) = M1
!       2nd-excited.
        P( N+2, 1 ) = px2
        P( N+2, 2 ) = py2
        P( N+2, 3 ) = pz2
        P( N+2, 4 ) = E2
        P( N+2, 5 ) = M2

!       Gives other status.
!       Mother.
        K( i_parton0, 4 ) = N + 1
        K( i_parton0, 5 ) = N + 2
!       1st-excited.
        K( N+1, 1 ) = 2
        K( N+1, 2 ) = KF1
        K( N+1, 3 ) = i_parton00
        K( N+1, 4 ) = 0
        K( N+1, 5 ) = 0
!       2nd-excited.
        K( N+2, 1 ) = 1
        K( N+2, 2 ) = KF2
        K( N+2, 3 ) = i_parton00
        K( N+2, 4 ) = 0
        K( N+2, 5 ) = 0
        N = N + 2
!----------------------------   Properties Giving   ----------------------------
!-------------------------------------------------------------------------------


!       Goes back for new qqbar if enough energy.
        if( pz0 > PZ_threshold ) goto 200

999     continue
        if( N == N00 ) goto 1000
!       First daughter of the original parton.
        K( i_parton00, 4 ) = N00 + 1
!       The last daughter of the original parton.
        K( i_parton00, 5 ) = N


!-------------------------------------------------------------------------------
!---------------------   Rotating Back & Position Giving  ----------------------
!       Appends the momenta of the mother quark to "N+1" for the rotation.
        do i=1,5,1
            K( N+1, i ) = K( i_parton00, i )
            P( N+1, i ) = P( i_parton00, i )
            V( N+1, i ) = V( i_parton00, i )
        end do
        K( N+1, i ) = 1

!       Rotation.
        theta = PYANGL( P00(3), SQRT( P00(1)**2 + P00(2)**2 ) )
        phi   = PYANGL( P00(1), P00(2) )
!       The positions will be rotated too. Do not use them.
        ! MSTU(33) = 1
        call PYROBO( N00+1, N+1, theta, phi, 0D0, 0D0, 0D0 )

!       Gives three-position and time to excited qqbar.
!       Generated q and qbar are arranged around the source parton
!        within 0-0.5 fm randomly in each one of the three positions
!        and has the same time as source parton + formation time.
        do i = N00+1, N, 1
            do j=1,3
                rr(j) = PYR(1)*0.5D0
                V(i,j) = V00(j) + rr(j)
                if( PYR(1) > 0.5D0 ) V(i,j) = V00(j) - rr(j)
            end do
            V(i,4) = V00(4)
            V(i,5) = 0D0
!TODO(Lei20241016): gives formation time.
            ! IF( .true. )THEN
            IF( .false. )THEN
                i_tau = 1
!               Formation time method 1.
!               Proper formation time of t0/10 was assumed for partons.
                if( i_tau == 1 )then
                    tau0 = t0 / 10D0
                    V(i,4) = V00(4) + tau0 * P(i,4) / P(i,5)
                    V(i,5) = 0D0
!               Formation time method 2, tau = E/mT^2.
                else if( i_tau == 2 )then
                    V(i,4) = V00(4) &
                           + P(i,4) / (P(i,1)**2 + P(i,2)**2 + P(i,5)**2) &
                           * 0.197327D0
                end if
            END IF
        end do

!       Returns the momenta of the mother quark from "N+1" back.
        do i=1,5,1
            P( i_parton00, i ) = P( N+1, i )
        end do
        do i=1,5,1
            V( i_parton00, i ) = V00( i )
        end do
!---------------------   Rotating Back & Position Giving  ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!       Collects the difference (lost) of four-momentum before and
!        after deexcitation, then tries to shares them.
        lost_4mom = 0D0
        do j=1,4,1
            do i= N00+1, N, 1
                lost_4mom(j) = lost_4mom(j) + P(i,j)
            end do
            lost_4mom(j) = lost_4mom(j) + P( i_parton00, j )
            lost_4mom(j) = P00(j) - lost_4mom(j)
        end do
!       Collects the rest four-momentum to "throe_p".
        throe_p = throe_p + lost_4mom
!-----------------------   Lost 4-Momentum Collecting   ------------------------
!-------------------------------------------------------------------------------


1000    continue
!       Recovers.
        PARJ(41) = PARJ41
        PARJ(42) = PARJ42
        PARJ(51) = PARJ51
        PARJ(52) = PARJ52
        PARJ(53) = PARJ53
        PARJ(54) = PARJ54
        PARJ(55) = PARJ55


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_debug( debug_mode )
!!      A simple debug treatment.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sin_local/nsin,non,ksin(kszj,5),psin(kszj,5),vsin(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension KFq(4), ps(4), vv(4), pp(20,5)


!       "sbh" will store spectators and non-hadorns.
        nbh = 0
        nsin = 0
        do i=1,N,1
            KS = K(i,1)
            KF = K(i,2)
            if( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) cycle
            KF_ABS = ABS( KF )
            KFq = 0
            i_debug = 0
            ! If KF of a particle = x-y-z-...-i4-i3-i2-i1.
            i4 = MOD( KF_ABS/1000, 10 )
            i3 = MOD( KF_ABS/100,  10 )
            i2 = MOD( KF_ABS/10,   10 )
            i1 = MOD( KF_ABS,      10 )
            pT2 = P(i,1)**2 + P(i,2)**2
            ! Excludes spectator nucleons.
            if( ( (KF == 2212 .OR. KF == 2112) .AND. pT2 < 1D-15 ) &
                .OR. IS_NUCLEUS(KF) )then
                i_debug = 0
            ! Excludes potential colored particles.
            ! Function "IPAPYK" is in Pythia8_fort_interface.f90.
            else if( IPAPYK(i,12) /= 0 )then
                i_debug = 0
            ! Hadrons.
            else if( KF_ABS > 100 .AND. KF_ABS < 1000000 )then
                rand_num = PYR(1)
                ! Baryon -> qqq / qbar-qbar-qbar.
                if( KF_ABS > 1000 .AND. i4*i3*i2*i1 > 0 )then
                    i_debug = 2
                    ! Light AAB type I (e.g. uud: 2212 p, 2214 Delta+)
                    !  or light AAA type (e.g. ddd: 1114 Delta-).
                    if( i4 == i3 .AND. i4 < 4 )then
                        ! Spin-1/2 baryon (e.g. 2212 p).
                        if( i1 == 2 )then
                            ! Quark + Spin-0 xy diquark (e.g. 2 u + 2101 ud_0).
                            if( rand_num <= 1D0/2D0 )then
                                KFq(1) = i4
                                KFq(2) = i3*1000 + i2*100 + 1
                            ! Quark + Spin-1 xy diquark (e.g. 2 u + 2103 ud_1).
                            else if( rand_num <= 2D0/3D0 )then
                                KFq(1) = i4
                                KFq(2) = i3*1000 + i2*100 + 3
                            ! Quark + Spin-1 xx diquark (e.g. 1 d + 2203 uu_1).
                            else
                                KFq(1) = i2
                                KFq(2) = i4*1000 + i3*100 + 3
                            end if
                        ! Spin-3/2 baryon (e.g. 2214 Delta+).
                        ! Or for light AAA type (e.g. ddd: 1114 Delta-),
                        !  Only to quark + spin-1 xx diquark.
                        else if( i1 == 4 )then
                            ! Quark + Spin-1 xy diquark (e.g. 2 u + 2103 ud_1).
                            if( rand_num <= 2D0/3D0 )then
                                KFq(1) = i4
                                KFq(2) = i3*1000 + i2*100 + 3
                            ! Quark + Spin-1 xx diquark (e.g. 1 d + 2203 uu_1).
                            else
                                KFq(1) = i2
                                KFq(2) = i4*1000 + i3*100 + 3
                            endif
                        end if
                    ! Light AAB type II (ABB, e.g. duu: 2112 n).
                    else if( i4 == 1 .OR. i4 == 2 )then
                        if( i1 == 2 ) then
                            if( rand_num <= 1D0/2D0 )then
                                KFq(1) = i2
                                KFq(2) = i4*1000 + i3*100 + 1
                            else if( rand_num <= 2D0/3D0 ) then
                                KFq(1) = i2
                                KFq(2) = i4*1000 + i3*100 + 3
                            else
                                KFq(1) = i4
                                KFq(2) = i3*1000 + i2*100 + 3
                            end if
                        else if( i1 == 4 )then
                            if( rand_num <= 2D0/3D0 )then
                                KFq(1) = i2
                                KFq(2) = i4*1000 + i3*100 + 3
                            else
                                KFq(1) = i4
                                KFq(2) = i3*1000 + i2*100 + 3
                            end if
                        end if
                    ! Light ABC type, Lambda-like and Sigma-like
                    !  (e.g. uds: 3122 Lambda0, 3212 Sigma0).
                    ! Or heavy baryon, in which the heavier quark will gain
                    !  "larger" momentum after decay.
                    else if( i4 > 2 )then
                        KFq(1) = i4
                        ! Lambda-like: quark + spin-0 diquark.
                        if( i3 < i2 )then
                            KFq(2) = i2*1000 + i3*100 + 1
                        ! Sigma-like:  quark + spin-1 diquark.
                        else
                            KFq(2) = i3*1000 + i2*100 + 3
                        end if
                    end if
                    ! Antibaryon.
                    if( KF < 0 )then
                        KFq(1) = -KFq(1)
                        KFq(2) = -KFq(2)
                    end if
                    ! Diquark -> qq / qbar-qbar.
                    KFq(3) = MOD( KFq(2) / 1000, 10 )
                    KFq(4) = MOD( KFq(2) / 100,  10 )
                ! Meson -> qqbar..
                else if( i4 == 0 .AND. i3*i2*i1 > 0 )then
                    i_debug = 1
                    ! Normal ABbar mesons(e.g. udbar: 211 pi+; dubar: -211 pi-).
                    if( i3 /= i2 )then
                        ! E.g. udbar: 211 pi+.
                        if( SIGN(1,KF)*(-1)**i3 == 1 )then
                            KFq(1) =  i3
                            KFq(2) = -i2
                        ! E.g. dubar: -211 pi-.
                        else
                            KFq(1) =  i2
                            KFq(2) = -i3
                        end if
                    ! Onia AAbar mesons.
                    else
                        ! eta, eta-eta' whit the mixing angle theta_P.
                        if( KF == 221 .OR. KF == 331 )then
                            ! theta_P = -14.1D0
                            theta_P = -15D0
                            ! eta
                            if( KF == 221 )then
                                prob_uubar_ddbar = ( COS(theta_P) / SQRT(6D0) &
                                                - SIN(theta_P) / SQRT(3D0) )**2
                                prob_ssbar = ( -2D0 * COS(theta_P) / SQRT(6D0) &
                                                - SIN(theta_P) / SQRT(3D0) )**2
                            ! eta'
                            else
                                prob_uubar_ddbar = ( SIN(theta_P) / SQRT(6D0) &
                                                + COS(theta_P) / SQRT(3D0) )**2
                                prob_ssbar = ( -2D0 * SIN(theta_P) / SQRT(6D0) &
                                                + COS(theta_P) / SQRT(3D0) )**2
                            end if
                            prob_tot = 2D0 * prob_uubar_ddbar + prob_ssbar
                            relative_uubar = prob_uubar_ddbar / prob_tot
                            relative_ddbar = 2D0 * relative_uubar
                            ! ddbar
                            if( rand_num <= relative_uubar )then
                                KFq(1) =  1
                                KFq(2) = -1
                            ! uubar
                            else if( rand_num <= relative_uubar )then
                                KFq(1) =  2
                                KFq(2) = -2
                            ! ssbar
                            else
                                KFq(1) =  3
                                KFq(2) = -3
                            end if
                        ! pio0, rho0, omega.
                        else if( i3 == 1 .OR. i3 == 2 )Then
                            ! ddbar
                            if( rand_num <= 0.5D0 )then
                                KFq(1) =  1
                                KFq(2) = -1
                            ! uubar
                            else
                                KFq(1) =  2
                                KFq(2) = -2
                            end if
                        ! phi and heavy onia.
                        else
                            KFq(1) =  i3
                            KFq(2) = -i3
                        end if
                    end if
                end if
            ! Something unrecognized.
            else
                i_debug = 0
            end if
            ! Moves particles that do not debug to "sbh".
            if( i_debug == 0 )then
                nbh = nbh + 1
                do j=1,5,1
                    kbh(nbh,j) = K(i,j)
                    pbh(nbh,j) = P(i,j)
                    vbh(nbh,j) = V(i,j)
                end do
            ! Debug hadrons to q/qbar.
            else
                ! Formation time.
                tau = PARU(3) * P(i,4) / ( P(i,1)**2 + P(i,2)**2 + P(i,5)**2 )
                ! Same time as parent hadron.
                ! tau = V(i,4)
                ! Propagates the hadron (debuged partons) to the formation time.
                do j=1,3,1
                    vv(j) = V(i,j) + tau * P(i,j) / P(i,4)
                end do
                vv(4) = V(i,4) + tau
                ! Gives the current algebra masses to quarks.
                dm1 = PARF( 90 + ABS( KFq(1) ) )
                dm2 = PYMASS( KFq(2) )
                if( i_debug == 1 ) dm2 = PARF( 90 + ABS( KFq(2) ) )
                ! Hadron -> qqbar or q-diquark (or qbar-diqbar).
                ! Gets momemta (will be stored to "pp" inside "decmom").
                do j=1,4,1
                    ps(j) = P(i,j)
                end do
                call decmom( ps, pp, dm1, dm2, i_succ )
                if( i_succ == 0 )then
                    write(22,*) "In debug decmom-1&2 i_succ=0! " &
                             // "KF1, KF2, m1, m2=", KFq(1), KFq(2), dm1, dm2
                    stop
                end if
                ! Gives debug quarks/diquark status.
                ! 1st.
                nsin = nsin + 1
                ksin(nsin,1) = 2
                ksin(nsin,2) = KFq(1)
                ksin(nsin,3) = 0
                ksin(nsin,4) = 0
                ksin(nsin,5) = 0
                do j=1,4,1
                    psin(nsin,j) = pp(1,j)
                    vsin(nsin,j) = vv(j)
                end do
                psin(nsin,5) = dm1
                vsin(nsin,5) = 0D0
                ! 2nd.
                nsin = nsin + 1
                ksin(nsin,1) = 1
                ksin(nsin,2) = KFq(2)
                ksin(nsin,3) = 0
                ksin(nsin,4) = 0
                ksin(nsin,5) = 0
                do j=1,4,1
                    psin(nsin,j) = pp(2,j)
                    vsin(nsin,j) = vv(j)
                end do
                psin(nsin,5) = dm2
                vsin(nsin,5) = 0D0
                ! Upto q + diquark.
                if( debug_mode > 101 ) cycle
                ! Diquark -> qq.
                if( i_debug == 2 )then
                    ! Gives the current algebra masses to quarks.
                    dm3 = PARF( 90 + ABS( KFq(3) ) )
                    dm4 = PARF( 90 + ABS( KFq(4) ) )
                    ! Gets momemta (will be stored to "pp" inside "decmom").
                    do j=1,4,1
                        ps(j) = pp(2,j)
                    end do
                    call decmom( ps, pp, dm3, dm4, i_succ )
                    if( i_succ == 0 )then
                        write(22,*) "In debug decmom-3&4 i_succ=0! " &
                                 // "KF3, KF4, m3, m4=",KFq(3),KFq(4),dm3,dm4
                        stop
                    end if
                    ! Gives debug quarks/diquark status.
                    ! 2nd.
                    ! nsin = nsin + 1
                    ksin(nsin,1) = 2
                    ksin(nsin,2) = KFq(3)
                    ksin(nsin,3) = 0
                    ksin(nsin,4) = 0
                    ksin(nsin,5) = 0
                    do j=1,4,1
                        psin(nsin,j) = pp(1,j)
                        vsin(nsin,j) = vv(j)
                    end do
                    psin(nsin,5) = dm3
                    vsin(nsin,5) = 0D0
                    ! 3rd.
                    nsin = nsin + 1
                    ksin(nsin,1) = 1
                    ksin(nsin,2) = KFq(4)
                    ksin(nsin,3) = 0
                    ksin(nsin,4) = 0
                    ksin(nsin,5) = 0
                    do j=1,4,1
                        psin(nsin,j) = pp(2,j)
                        vsin(nsin,j) = vv(j)
                    end do
                    psin(nsin,5) = dm4
                    vsin(nsin,5) = 0D0
                end if
            end if
        end do
        N = 0

!       "sin_local" -> "PYJETS".
        do j=1,5,1
            do i=1,nsin,1
                K(i,j) = ksin(i,j)
                P(i,j) = psin(i,j)
                V(i,j) = vsin(i,j)
            end do
        end do
        N    = nsin
        nsin = 0

!       Appends "sbh" to "PYJETS".
        ! do i=1,nbh,1
        !     N = N + 1
        !     do j=1,5,1
        !         K(N,j) = kbh(i,j)
        !         P(N,j) = pbh(i,j)
        !         V(N,j) = vbh(i,j)
        !     end do
        ! end do
        ! nbh = 0


        return
        end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine funcz( zz )
!!      Samples daughter energy fraction, zz, from mother according to
!!       simple fragmentation function.
!       adj1(29) = 11 : Lund string fragmentation function
!                = 12 : Field-Feymman fragmentation function (FF)
!                = 13 : Perterson/SLAC fragmentation function (P/S)
!       Lund string fragmentation function:
!        f(z) = (1/z) * (1-z)^a * exp( - b * mT^2 / z )
!             ~ (1/z) * (1-z)^a * exp( - b * 2.36 / z )
!        assume: mT^2 ~ <pT>^2 + m_u^2 ~ 1.5^2 + 0.333^2 = 2.36
!       FF fragmentation function:
!        f(z) = 1 - a + 3 * a * (1-z)**2, 0 <= z <= 1,
!        and its largest value: fmax = f(0) = 1 - a + 3*a = 1 + 2*a.
!       P/S fragmentation function:
!        f(z) = 1 / ( z * ( 1 - 1/z - a/(1-z) )^2 ).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/local_fmax/ fmax_value(5)
        common/sa24/adj1(40),nnstop,non24,zstop
        data i_call / 0 /


        i_call = i_call + 1
!       Sets parameters.
        ! Lund a parameter
        a   = adj1(6)
        ! Lund b parameter
        b   = adj1(7)
        ! Field-Feynman parameter
        aFF = a_FF
        ! Peterson/SLAC parameter
        c   = aPS_c
        ! 11, 12, 13
        i_z_samp_func = INT( adj1(29) )
        if( i_z_samp_func < 1 .OR. i_z_samp_func > 3 ) i_z_samp_func = 12

!       Finds max values.
        if( i_call == 1 )then
            fmax_value = 0.
!           1 = Lund
            fmax = 0.
            z1   = 0.025
            do i1=1,20
                z   = z1 + (i1-1)*0.05D0
                fzz = (1./z) * (1.-z)**a * dexp(-2.36*b/z)
                if(fzz > fmax) fmax = fzz
            end do
            fmax_value(1) = fmax
!           2 = FF
            fmax = 1. + 2.*aFF
            if(aFF < 0) fmax = 1. - aFF
            fmax_value(2) = fmax
!           3 = P/S
            fmax = 0.
            z1   = 0.025
            do i1=1,20
                z   = z1 + (i1-1)*0.05D0
                fzz = 1./( z * (1.-1./z-c/(1.-z)) * (1.-1./z-c/(1.-z)) )
                if(fzz > fmax) fmax = fzz
            end do
            fmax_value(3) = fmax
!           4 = dummy, but set FF as default
            ! fmax =
            fmax = fmax_value(2)
            fmax_value(4) = fmax
!           5 = dummy, but set FF as default
            ! fmax =
            fmax = fmax_value(2)
            fmax_value(5) = fmax
        end if

!       Samples z.
        fmax = fmax_value( i_z_samp_func )

100     ran1 = PYR(1)
        ran2 = PYR(1)
        fm = ran1 * fmax

        select case( i_z_samp_func )
        ! Lund
        case(11)
            fran2 = (1./ran2) * (1.-ran2)**a * dexp(-2.36*b/ran2)
        ! FF
        case(12)
            fran2 = 1. - aFF + 3.*aFF*(1.-ran2)**2
            if(aFF < 0) fran2 = aFF + 3.*aFF*(1.-ran2)**2
        ! P/S
        case(13)
            fran2 = 1./( ran2 * (1.-1./ran2-c/(1.-ran2))**2 )
        ! Sets FF as default.
        case default
            fran2 = 1. - aFF + 3.*aFF*(1.-ran2)**2.
            if(aFF < 0) fran2 = aFF + 3.*aFF*(1.-ran2)**2
        end select

        if(fm > fran2) goto 100
        zz = ran2


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine phas( i1, i2, i3, isucc, i_M_or_B, iphas )
!!      The phase space judgement.
!       i_M_or_B = 2 for meson
!       i_M_or_B = 3 for baryon
!       iphas: = 1, complete phase space constraint
!              = 2, position constraint only
!              = 3, momentum constraint only
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension ri1(3),ri2(3),ri3(3),pi1(3),pi2(3),pi3(3)


        isucc = 1
!       Critical value.
        delc = adj1(22)
        if( iphas == 2 .OR. iphas == 3 ) delc = 0.5D0 * delc

!       For baryon.
        if( i_M_or_B == 3 )then
            ri1(1) = V(i1,1)
            ri1(2) = V(i1,2)
            ri1(3) = V(i1,3)
            ri2(1) = V(i2,1)
            ri2(2) = V(i2,2)
            ri2(3) = V(i2,3)
            ri3(1) = V(i3,1)
            ri3(2) = V(i3,2)
            ri3(3) = V(i3,3)
            ri121  = ri1(1) - ri2(1)
            ri122  = ri1(2) - ri2(2)
            ri123  = ri1(3) - ri2(3)
            ri131  = ri1(1) - ri3(1)
            ri132  = ri1(2) - ri3(2)
            ri133  = ri1(3) - ri3(3)
            ri231  = ri2(1) - ri3(1)
            ri232  = ri2(2) - ri3(2)
            ri233  = ri2(3) - ri3(3)
            delr12 = SQRT( ri121*ri121 + ri122*ri122 + ri123*ri123 )
            delr13 = SQRT( ri131*ri131 + ri132*ri132 + ri133*ri133 )
            delr23 = SQRT( ri231*ri231 + ri232*ri232 + ri233*ri233 )
            pi1(1) = P(i1,1)
            pi1(2) = P(i1,2)
            pi1(3) = P(i1,3)
            pi2(1) = P(i2,1)
            pi2(2) = P(i2,2)
            pi2(3) = P(i2,3)
            pi3(1) = P(i3,1)
            pi3(2) = P(i3,2)
            pi3(3) = P(i3,3)
            pi121  = pi1(1) - pi2(1)
            pi122  = pi1(2) - pi2(2)
            pi123  = pi1(3) - pi2(3)
            pi131  = pi1(1) - pi3(1)
            pi132  = pi1(2) - pi3(2)
            pi133  = pi1(3) - pi3(3)
            pi231  = pi2(1) - pi3(1)
            pi232  = pi2(2) - pi3(2)
            pi233  = pi2(3) - pi3(3)
            delp12 = SQRT( pi121*pi121 + pi122*pi122 + pi123*pi123 )
            delp13 = SQRT( pi131*pi131 + pi132*pi132 + pi133*pi133 )
            delp23 = SQRT( pi231*pi231 + pi232*pi232 + pi233*pi233 )
            ! Complete phase space constaint.
            if( iphas == 1 )then
                del12 = delr12 * delp12
                del13 = delr13 * delp13
                del23 = delr23 * delp23
            ! Position constraint.
            else if( iphas == 2 )then
                del12 = delr12
                del13 = delr13
                del23 = delr23
            ! Momentum constraint.
            elseif( iphas == 3 )then
                del12 = delp12
                del13 = delp13
                del23 = delp23
            ! No constraint.
            else
                del12 = 0D0
                del13 = 0D0
                del23 = 0D0
            end if
            if( del12 <= delc .AND. del13 <= delc .AND. del23 <= delc )then
                isucc = 1
            else
                isucc = 0
            end if

        else if( i_M_or_B == 2 )then
            ri1(1) = V(i1,1)
            ri1(2) = V(i1,2)
            ri1(3) = V(i1,3)
            ri2(1) = V(i2,1)
            ri2(2) = V(i2,2)
            ri2(3) = V(i2,3)
            ri121  = ri1(1) - ri2(1)
            ri122  = ri1(2) - ri2(2)
            ri123  = ri1(3) - ri2(3)
            delr   = SQRT( ri121*ri121 + ri122*ri122 + ri123*ri123 )
            pi1(1) = P(i1,1)
            pi1(2) = P(i1,2)
            pi1(3) = P(i1,3)
            pi2(1) = P(i2,1)
            pi2(2) = P(i2,2)
            pi2(3) = P(i2,3)
            pi121  = pi1(1) - pi2(1)
            pi122  = pi1(2) - pi2(2)
            pi123  = pi1(3) - pi2(3)
            delp   = SQRT( pi121*pi121 + pi122*pi122 + pi123*pi123 )
            ! Complete phase space constaint.
            if( iphas == 1 )then
                delrp = delr * delp
            ! Position constraint.
            else if( iphas == 2 )then
                delrp = delr
            ! Momentum constraint.
            else if( iphas == 3 )then
                delrp = delp
            ! No constraint.
            else
                delrp = 0D0
            end if
            if( delrp <= delc )then
                isucc = 1
            else
                isucc = 0
            end if
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_f( energy, kf, amq )
!!      Samples the flavor (mass) of the generated qqbar pair.
!       energy: energy of original q (qbar) or g.
!       kf (amq): flavor code (mass) of the generated quark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa38/ prob_ratio_q(6), am(6), amqq(6)


        kf  = 0
        amq = 0D0
!       Quark masses.
        do i=1,6,1
            am(i) = amass(i)
        end do
!       Mass of qqbar pair.
        amqq = 2*am
!       Mass of uubar pair.
        amuu = amqq(1)
!       Too low energy to excite qqbar pair.
        if( energy  <  amuu ) return

!       Samples the flavor of generated quark pair according to prob_ratio_q.
        prob_ratio_tot = 0D0
        do i=1,6,1
!       Sum of ratio upto i.
            if( energy >=  amqq(i) )then
                n_flavor = i
                prob_ratio_tot = prob_ratio_tot + prob_ratio_q(i)
            end if
        end do
        prob_interval_low = 0D0
        prob_interval_upp = 0D0
!       Calculate probability value at upper end point of each probability
!        interval and determine simulated flavor code.
        rand_num = PYR(1)
        do j=1,n_flavor,1
            ratio_interval = prob_ratio_q(j) / prob_ratio_tot
            prob_interval_upp = prob_interval_upp + ratio_interval
            if( rand_num > prob_interval_low .AND. &
                rand_num <= prob_interval_upp )then
                kf = j
                amq = am(j)
                return
            end if
            prob_interval_low = prob_interval_upp
        end do

!       For example:
!        in the case of ratio probability is u:d:s:c:b:t=1:1:0.5:*:*:*
!        and has 'energy >= amqq(3)', the possible candidates are u, d
!        and s, the corresponding probability intervals are (0,0.4),
!        (0.4,0.8], and (0.8,1).
!        Ratio: u:d:s=          1          :         1          :   0.5
!                     |____________________|____________________|__________|
!        Interval:    0                   0.4                  0.8         1
!        0.4=1/sum (sum=1+1+0.5); 0.8=(1+1)/sum.

!       Mass defined in PYTHIA, in GeV.
!       Kinematical mass:
!           amd = 0.33D0
!           amu = 0.33D0   ! amu=amd
!           ams = 0.5D0
!           amc = 1.5D0
!           amb = 4.8D0
!           amt = 175D0
!       Current algebra mass
!           amd = 0.0099D0
!           amu = 0.0056D0   ! amu < amd
!           ams = 0.199D0
!           amc = 1.23D0
!           amb = 4.17D0
!           amt = 165D0
!       Constituent mass (Cf. PYTHIA6 manual P465 and program)
!           amd = 0.325D0
!           amu = 0.325D0   ! amu=amd
!           ams = 0.5D0
!           amc = 1.6D0
!           amb = 5.0D0
!           amt = 0D0 (null)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tdgaus(v,pmax,np,pp)
!!..... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
!!..... set pmax < 0 if pmax should be infinity.
!!..... np : the total number of particles wanted to sample their pT
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa33/smadel,ecce,secce,parecc,iparres
        dimension pp(50,2)


        do 30 i=1,np
            p2 = 0
            if(v <= 1.e-8)return
            if(pmax < 1.E-9)return
            if(pmax < 0)then
                a = 1.
                goto 10
            endif
            aa=-pmax/v
            if(aa < -70)then
                a=1.
                goto 10
            endif
            a = 1. - exp(aa)
10          p2 = -v*log(max(1.e-20,1. - a*PYR(1)))
            if(p2 < 0.) goto 10
            ps=sqrt(p2)
            fi=2.*3.1415926*PYR(1)
!       randomly sample [px,py] on circle of ellipsoid with half major axis
!       of ps*(1+smadel) and half minor axis of ps*(1-smadel)
          pp(i,1)=ps*cos(fi)*(1+smadel)
          pp(i,2)=ps*sin(fi)*(1-smadel)
!       note: ps is not in the dimension list
30      continue


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_glu
!!      Removes gluons from "PYJETS" to "sa36".
!!      Adjusts "sbe", the diquark list and the string list too.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)


!       Hadronization model.
        i_had_model = INT( adj1(12) )

        nglu = 0
        jb = 0
201     continue
        do i1 = jb+1, N, 1
            KF = K(i1,2)
            ! Keeps.
            if( KF /= 21 )then
                jb = jb + 1
                cycle
            end if
            nglu = nglu + 1
            do i2=1,5,1
                kglu( nglu, i2 ) = K( i1, i2 )
                pglu( nglu, i2 ) = P( i1, i2 )
                vglu( nglu, i2 ) = V( i1, i2 )
            end do
!           Moves the particle list one step downward from i1+1 to N.
            do jj=1,5,1
                do j = i1+1, N, 1
                    K( j-1, jj ) = K(j,jj)
                    P( j-1, jj ) = P(j,jj)
                    V( j-1, jj ) = V(j,jj)
                end do
            end do
            N = N - 1
!           For the special version of SFM.
            if( i_had_model == 3 )then
                ! Adjusts line numbers of components of broken diquarks.
                do i_diq = 1, idi, 1
                    if( i1 < ifcom(i_diq) ) ifcom(i_diq) = ifcom(i_diq) - 1
                    if( i1 <   npt(i_diq) )   npt(i_diq) =   npt(i_diq) - 1
                end do
                do j = i1+1, N+1, 1
                    ndiq( j-1 ) = ndiq(j)
                end do
                ! Adjusts line numbers of the string location.
                do i_string = 1, nstr1, 1
                    if( i1 <= nstr1a(i_string) ) &
                        nstr1a(i_string) = nstr1a(i_string) - 1
                    if( i1 <= nstr1v(i_string) ) &
                        nstr1v(i_string) = nstr1v(i_string) - 1
                end do
            end if
            goto 201
        end do

!       Removes gluons from "sbe" for the special version of SFM.
        if( i_had_model == 3 )then
            jb = 0
301         continue
            do i1 = jb+1, nbe, 1
                KF = kbe(i1,2)
                ! keep
                if( KF /= 21 )then
                    jb = jb + 1
                    goto 302
                end if
                if( i1 == nbe )then
                    nbe = nbe - 1
                    goto 303
                end if
!               Moves the particle list one step downward from i1+1 to nbe.
                do jj=1,5,1
                    do j = i1+1, nbe, 1
                        kbe( j-1, jj ) = kbe(j,jj)
                        pbe( j-1, jj ) = pbe(j,jj)
                        vbe( j-1, jj ) = vbe(j,jj)
                    end do
                end do
                nbe = nbe - 1
                goto 301
302         end do
303         continue
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_glu
!!      Breaks up gluons (in 'sa36') and gives flavor (mass) and four momentum
!!       (position) to broken objects, which are assumed to be a string and
!!       filling in 'PYJETS'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)


!       Kinematical mass.
        amu = amass(2)
        amuu=2*amu
!       throw away g with energy<amuu
        jb=1
2010    continue
        do i1=jb,nglu,1
        eg=pglu(i1,4)
        if(eg < amuu)then
        do i2=1,4
        throe_p(i2)=throe_p(i2)+pglu(i1,i2)
        enddo
        if(i1 == nglu)then
            nglu = nglu - 1
            goto 2020
        endif
!       move particle list ('sa36') one step downward from i1+1 to nglu
        do jj=1,5
        do j=i1+1,nglu
        kglu(j-1,jj)=kglu(j,jj)
        pglu(j-1,jj)=pglu(j,jj)
        vglu(j-1,jj)=vglu(j,jj)
        enddo
        enddo
        nglu=nglu-1
        jb=i1
        goto 2010
        endif
        enddo
2020    continue

!       g (in 'sa36') -> qq_{bar} (as a string filling in 'PYJETS')
        do i1=1,nglu   ! do loop over gluons
        eg=pglu(i1,4)
        call break_f(eg,kf,amq)
        kf1=kf
        kf2=-kf
        am1=amq
        am2=amq

        K(N+1,1)=2   ! 'A'
        K(N+2,1)=1   ! 'V'
        K(N+1,2)=kf1
        K(N+2,2)=kf2
        K(N+1,3)=0
        K(N+2,3)=0
        K(N+1,4)=0
        K(N+2,4)=0
        K(N+1,5)=0
        K(N+2,5)=0
!       P(N+1,5)=am1
!       P(N+2,5)=am2

!       give four momentum to the breaked quarks
        call bream_glu(i1,kf1,kf2)
!       give four coordinate to the breaked quarks
        call coord_glu(i1)
        if(i1 == nglu)then
        N=N+2
        goto 200
        endif
        N=N+2
        enddo


200     return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream_glu(ii,kf1,kf2)
!!      Gives four momentum to the broken quarks.
!       ii: line number of initial gluon in 'sa36'
!       kf1,kf2: flavor codes of broken quarks
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension pi(4),ps(4),pp(20,5)
        integer decsuc


        am1=amass(kf1)
        am2=amass(kf2)
        pp(1,5)=am1   ! mass of first broken quark
        pp(2,5)=am2   ! mass of second broken quark
!       pp : four momenta & mass of broken quarks, local variable
        do i1=1,4
        ps(i1)=pglu(ii,i1)
        enddo
!       ps : four momentum of initial gluon, local variable
        ! for 'decay method'
        goto 400

        ! for 'random three momentum method'
401     do i1=1,3
        pi(i1)=PYR(1)*pglu(ii,i1)
        pp(1,i1)=pi(i1)
        pp(2,i1)=ps(i1)-pi(i1)
        enddo
        pp11=pp(1,1)
        pp12=pp(1,2)
        pp13=pp(1,3)
        pp14=am1*am1+pp11*pp11+pp12*pp12+pp13*pp13
        if(pp14 <= 0.)pp14=1.e-20
        pp(1,4)=dsqrt(pp14)
        pp21=pp(2,1)
        pp22=pp(2,2)
        pp23=pp(2,3)
        pp24=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
        if(pp24 <= 0.)pp24=1.e-20
        pp(2,4)=dsqrt(pp24)
        ! activate it for 'random three momentum method'
        goto 300

        ! for 'decay method'
400     continue
!       Decay method.
        decsuc=1
        call decmom(ps,pp,am1,am2,decsuc)
        ! return to random three momentum method
        if(INT(decsuc) == 0) goto 401
300     continue
!       adjust four momentum conservation by iteration,no more than
!        4000 iterations
!       call conser(2,ps)
!       do i1=1,2
!       do i2=1,4
!       ppi1i2=pp(i1,i2)
!       pp12=abs(ppi1i2)
!       if(pp12 > 1.d6)pp12=1.d6
!       pp(i1,i2)=sign(pp12,ppi1i2)
!       enddo
!       enddo
        do i1=1,4
        P(N+1,i1)=pp(1,i1)
        enddo
        P(N+1,5)=am1
        do i1=1,4
        P(N+2,i1)=pp(2,i1)
        enddo
        P(N+2,5)=am2

        ! Collects lost 4-momentum.
        do i2=1,4
        throe_p(i2) = throe_p(i2) + ( ps(i2) - P(N+1,i2) - P(N+2,i2) )
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord_glu(ii)   ! 160822
!!      Gives four positions to broken quarks.
!       first broken quark takes the four position of gluon
!       second broken quark is arranged around first ones within
!        0.5 fm randumly in each of three positions and has same
!        fourth position as gluon
!       ii: order number of gluon in 'sa36'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension rr(3)


        do i1=1,3
        V(N+1,i1)=vglu(ii,i1)
        rr(i1)=PYR(1)*0.5   ! 261002
        V(N+2,i1)=vglu(ii,i1)+rr(i1)
        if(PYR(1) > 0.5d0)V(N+2,i1)=vglu(ii,i1)-rr(i1)
        enddo
        V(N+1,4)=vglu(ii,4)
        V(N+2,4)=vglu(ii,4)
        return


        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa36(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension peo(4)


        do i=1,nn
            write(mstu(11),*)i,kglu(i,2),(pglu(i,j),j=1,4)
            write(9,*)i,kglu(i,2),(pglu(i,j),j=1,4)
        enddo
        call psum(pglu,1,nglu,peo)
        ich1=0
        do i1=1,nn
        kf=kglu(i1,2)
        ich1=ich1+PYCHGE(kf)
        enddo
        cc=ich1/3D0
        write(22,*)'sa36 nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   !
!       write(9,*)peo,ich1/3   !


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa37(nn,cc)
!!      Prints particle list and sum of momentum and energy.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa37/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)


        call psum(pbh,1,nbh,peo)
        ich1=0
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+PYCHGE(kf)
        enddo
        cc=ich1/3D0
        ! write(22,*)'sa37 nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        do i=1,nn
        write(22,*)i,kbh(i,1),kbh(i,2),(pbh(i,j),j=1,4)
        enddo
        write(22,*) "------------------------------------------------"// &
                    "------------------------------------------------"// &
                    "------------------------------------------------"


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine thptdi(pt,px,py)
!!      Generates transverse momentum according to thermodynamic distribution
!!       (exponential distribution)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa33/smadel,ecce,secce,parecc,iparres


!       Generate p_T and azimuthal angle, gives p_x and p_y.
        pt=-parj(21)*log(max(1d-10,PYR(1)))
        phi=paru(2)*PYR(1)
!       randomly sample [px,py] on circle of sphere with radius pt
        px=pt*cos(phi)*(1+smadel)
        py=pt*sin(phi)*(1-smadel)
        pt=SQRT(px*px+py*py)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine PAPTDI( i_pT, width, px, py, pT )
!!      Generates transverse momentum.
!       i_pT: = 0, px = py = pT = 0.
!             = 1, Gaussian pT (px and py) with width.
!             = 2, Exponential pT with width.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa33/smadel,ecce,secce,parecc,iparres


        pT = 0D0
        px = 0d0
        py = 0D0
        if( i_pT == 0 ) return

!       Generates p_T and azimuthal angle, gives Gaussian px and Gaussian py.
        if( i_pT == 1 ) pT = width * SQRT( -LOG( MAX( 1D-10, PYR(1) ) ) )

!       Generates p_T according to the Exponential distribution (thermodynamic).
        if( i_pT == 2 ) pT = - width * LOG( MAX( 1D-10, PYR(1) ) )

!       Non-uniform tail.
        if( PARJ(23) > PYR(1) ) pT = PARJ(24) * pT

!       Generates azimuthal angle.
        phi = PARU(2) * PYR(1)

!       Generates px and py.
        ! px = pT * COS(phi)
        ! py = pT * SIN(phi)
!       Randomly samples [px,py] on the ellipsoid with the half major axis
!        of pT*(1+smadel) and the half minor axis of pT*(1-smadel).
        px = pT * COS(phi) * (1 + smadel )
        py = pT * SIN(phi) * (1 - smadel )

!       pT might deviate little bit from original one.
        pT = SQRT( px**2 + py**2 )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine PAZDIS( i_Z, Z_width, KF1_in, KF2_in, mT2, Z )
!!      Generates the longitudinal splitting variable z.
!       i_z: = 0, random z.
!            = 1, Lund symmetric fragmentation function.
!            = 2, Field–Feynman(FF) + Peterson/SLAC(PS).
!            = 3, Lund + Peterson/SLAC.
!            = 4, Lund + Bowler.
!            = 5, as = 4, but interpolate for c and b.
!            = 6, Gaussian.
!            = 7, Gaussian (uds) + Peterson/SLAC (cb).
!            = 8, exponential.
!            = 9, exponential (uds) + Peterson/SLAC (cb).
!            = 10, Z = Z_width if 0 < Z_width < 1, or random Z.
!            = 11, simple Lund from "funcz"
!            = 12, simple Field–Feynman from "funcz"
!            = 13, simple Peterson/SLAC from "funcz"
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        real(kind=8) :: mT2


        Z = PYR(1)
        if( i_Z == 10 .AND. Z_width > 0D0 .AND. Z_width < 1D0 ) Z = Z_width
        if( i_Z == 0 .OR. i_Z == 10 )  return
!       Checks if heavy flavor.
        KF1 = KF1_in
        KF2 = KF2_in
        KFA = ABS( KF1_in )
        KFB = ABS( KF2_in )
!       For diquark.
        if( KFA >= 10 ) KFA = MOD( KFA / 1000, 10 )
        if( KFB >= 10 ) KFB = MOD( KFB / 1000, 10 )
        KFH = MAX( KFA, KFB )
!       Ensures KF1 has the larger heavy flavor.
        if( ( KFB == 4 .OR. KFB == 5 ) .AND. KFB > KFA )then
            KF1 = KF2_in
            KF2 = KF1_in
        end if

!       Samples z.
        select case( i_Z )
        case( 1:5 )
            MSTJ11 = MSTJ(11)
            MSTJ(11) = i_Z
            call PYZDIS( KF1, KF2, mT2, Z )
            MSTJ(11) = MSTJ11
        case( 6:9 )
            i_sample_z = 1
            if( i_Z == 8 .OR. i_Z == 9 ) i_sample_z = 2
            do while(.true.)
                call PAPTDI( i_sample_z, Z_width, ZX, ZY, Z )
                if( Z > 0D0 .OR. Z < 1D0 ) exit
            end do
            if( (KFH == 4 .OR. KFH == 5) .AND. &
                ( i_Z == 7 .OR. i_Z == 9 ) )then
                MSTJ11 = MSTJ(11)
                MSTJ(11) = 2
                call PYZDIS( KF1, KF2, mT2, Z )
                MSTJ(11) = MSTJ11
            end if
        case( 11:13 )
            call funcz( Z )
        end select


        RETURN
        END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine PASORT(i_begin,i_end, name_common,quantity,operation)
!!      Sorts the arrays of corresponding COMMON blocks by
!!       Indexing Quicksort algorithm (from max to min or min to max) or
!!       Knuth shuffle algorithm (random sorting).
!
!       i_begin: the beginning index of arrays in COMMON to be sorted.
!       i_end: the ending index of arrays in COMMON to be sorted.
!       Note: the following CHARACTERs are accepted with lower case only.
!       name_common: name of COMMON block to be sorted, CHARACTER,
!                    could be "pyjets", "sa2", "sbe", etc.
!       quantity: name of quantity to be brought into sorted order.
!                 CHARACTER, could be "kf",
!                                     "px", "py", "pz", "e", "m",
!                                     "rx", "ry", "rz", "t", "tau",
!                                     "pt", "y", "eta", "charge", etc.
!                 Or "1, 2, 3, ... , 25", corresponds to the FUNCTION
!                   PAPYP(I,J) (PYP) in PYTHIA 6. Cf. PYTHIA 6.4 manual.
!       operation: name of the operation, CHARACTER, could be
!                  "max_to_min", "min_to_max" or "random".
!
!       Usage examples:
!       (Note: the CHARACTERs are accepted with the lower case only.)
!        call PASORT( 5, N, "pyjets", "null",  "random" )
!        call PASORT( 1, N, "pyjets", "e",     "max_to_min" )
!        call PASORT( 1, N, "pyjets", "e",     "min_to_max" )
!        call PASORT( 1, N, "pyjets", "eta",   "min_to_max" )
!        call PASORT( 1, N, "pyjets", "|eta|", "min_to_max" )
!        call PASORT( 1, N, "sbh",    "|eta|", "min_to_max" )
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)
        common/delt/ndel,nodel,kdel(kszj,5),pdel(kszj,5),vdel(kszj,5)
        common/sin_local/nsin,non,ksin(kszj,5),psin(kszj,5),vsin(kszj,5)
!       Arrays of particle information of one NN pair for PYTHIA 8.
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
!       Arrays of particle information of a AA collision (total NN pairs).
        COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
               K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
!       Arrays of freeze-out particle information.
        COMMON/PYJETS_FO/ N_FO, NPAD_FO, &
                          K_FO(KSZJ,8), P_FO(KSZJ,7), V_FO(KSZJ,5)
!       Local temporary variables.
        character*(*) name_common, quantity, operation
        character*1 char_1
        allocatable k_tmp(:,:)
        allocatable p_tmp(:,:), v_tmp(:,:)
        allocatable quantity_tmp(:), i_order(:)


!       Firstly, saves /PYJETS/ to temporary arrays and dumps data from
!        "name_commom" to /PYJETS/.

!       Allocates local memory and saves /PYJETS/.
        allocate( k_tmp(N,8),  p_tmp(N,7),  v_tmp(N,5) )
        do j=1,5,1
            do i=1,N,1
                k_tmp(i,j) = K(i,j)
                p_tmp(i,j) = P(i,j)
                v_tmp(i,j) = V(i,j)
            end do
        end do

!       Identifies the operation.
        i_operation = 0
        if(      TRIM(ADJUSTL(operation))  ==  "min_to_max" )then
            i_operation =  1
            N_sort = i_end
        else if( TRIM(ADJUSTL(operation))  ==  "max_to_min" )then
            i_operation = -1
            N_sort = i_end
        else if( TRIM(ADJUSTL(operation))  ==  "random" )then
            i_operation =  0
            N_sort = i_end
            goto 666
        end if

!       Identifies "name_common" and dumps data.
        ! if(      TRIM(ADJUSTL(name_common))  ==  "pyjets"   )then
        if(      TRIM(ADJUSTL(name_common))  ==  "sa2"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = ksa(i,j)
                    P(i,j) = psa(i,j)
                    V(i,j) = vsa(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa36"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kglu(i,j)
                    P(i,j) = pglu(i,j)
                    V(i,j) = vglu(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa37"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kth(i,j)
                    P(i,j) = pth(i,j)
                    V(i,j) = vth(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sbe"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kbe(i,j)
                    P(i,j) = pbe(i,j)
                    V(i,j) = vbe(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "saf"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kaf(i,j)
                    P(i,j) = paf(i,j)
                    V(i,j) = vaf(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "aaff"  )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kaff(i,j)
                    P(i,j) = paff(i,j)
                    V(i,j) = vaff(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sbh"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kbh(i,j)
                    P(i,j) = pbh(i,j)
                    V(i,j) = vbh(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa1_h" )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kn(i,j)
                    P(i,j) = pn(i,j)
                    V(i,j) = rn(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sgam"  )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kgam(i,j)
                    P(i,j) = pgam(i,j)
                    V(i,j) = vgam(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "trs"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = ktrs(i,j)
                    P(i,j) = ptrs(i,j)
                    V(i,j) = vtrs(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "delt"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = kdel(i,j)
                    P(i,j) = pdel(i,j)
                    V(i,j) = vdel(i,j)
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sin_local"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = ksin(i,j)
                    P(i,j) = psin(i,j)
                    V(i,j) = vsin(i,j)
                end do
            end do
        ! Potential  overflow error.
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_py8"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = K_PY8(i,j)
                    P(i,j) = P_PY8(i,j)
                    V(i,j) = V_PY8(i,j)
                end do
            end do
            do j=6,7,1
                do i = i_begin, i_end, 1
                    k_tmp(i,j) = K_PY8(i,j)
                    p_tmp(i,j) = P_PY8(i,j)
                end do
            end do
            do i = i_begin, i_end, 1
                k_tmp(i,8) = K_PY8(i,8)
            end do
        ! Potential  overflow error.
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_aa"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = K_AA(i,j)
                    P(i,j) = P_AA(i,j)
                    V(i,j) = V_AA(i,j)
                end do
            end do
            do j=6,7,1
                do i = i_begin, i_end, 1
                    k_tmp(i,j) = K_AA(i,j)
                    p_tmp(i,j) = P_AA(i,j)
                end do
            end do
            do i = i_begin, i_end, 1
                k_tmp(i,8) = K_AA(i,8)
            end do
        ! Potential  overflow error.
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_fo"   )then
            do j=1,5,1
                do i = i_begin, i_end, 1
                    K(i,j) = K_FO(i,j)
                    P(i,j) = P_FO(i,j)
                    V(i,j) = V_FO(i,j)
                end do
            end do
            do j=6,7,1
                do i = i_begin, i_end, 1
                    k_tmp(i,j) = K_FO(i,j)
                    p_tmp(i,j) = P_FO(i,j)
                end do
            end do
            do i = i_begin, i_end, 1
                k_tmp(i,8) = K_FO(i,8)
            end do
        end if

!       Secondly, identifies the sorting-quantity.
        do i=1,10,1
            char_1 = quantity(1:1)
            if( char_1 /= " " ) exit
        end do

!       = 0, KF code.
        if(      TRIM(ADJUSTL(quantity))  ==  "kf" )then
            i_quantity = 0
!       = [1,25], PAPYP(i,j).
        else if( TRIM(ADJUSTL(quantity))  ==  "px" )then
            i_quantity = 1
        else if( TRIM(ADJUSTL(quantity))  ==  "py" )then
            i_quantity = 2
        else if( TRIM(ADJUSTL(quantity))  ==  "pz" )then
            i_quantity = 3
        else if( TRIM(ADJUSTL(quantity))  ==  "e"  )then
            i_quantity = 4
        else if( TRIM(ADJUSTL(quantity))  ==  "m"  )then
            i_quantity = 5
        else if( TRIM(ADJUSTL(quantity))  ==  "charge" )then
            i_quantity = 6
        else if( TRIM(ADJUSTL(quantity))  ==  "|p|" )then
            i_quantity = 8
        else if( TRIM(ADJUSTL(quantity))  ==  "pt" )then
            i_quantity = 10
        else if( TRIM(ADJUSTL(quantity))  ==  "y" )then
            i_quantity = 17
        else if( TRIM(ADJUSTL(quantity))  ==  "eta" )then
            i_quantity = 19
!       = [26,29], Velocity.
        else if( TRIM(ADJUSTL(quantity))  ==  "vx" )then
            i_quantity = 26
        else if( TRIM(ADJUSTL(quantity))  ==  "vy" )then
            i_quantity = 27
        else if( TRIM(ADJUSTL(quantity))  ==  "vz" )then
            i_quantity = 28
        else if( TRIM(ADJUSTL(quantity))  ==  "|v|" )then
            i_quantity = 29
!       = [-6, -1], space-time coordinates and lifetime.
        else if( TRIM(ADJUSTL(quantity))  ==  "rx" )then
            i_quantity = -1
        else if( TRIM(ADJUSTL(quantity))  ==  "ry" )then
            i_quantity = -2
        else if( TRIM(ADJUSTL(quantity))  ==  "rz" )then
            i_quantity = -3
        else if( TRIM(ADJUSTL(quantity))  ==  "t"  )then
            i_quantity = -4
        else if( TRIM(ADJUSTL(quantity))  ==  "tau" )then
            i_quantity = -5
        else if( TRIM(ADJUSTL(quantity))  ==  "|r|" )then
            i_quantity = -6
!       ABS( x + 100 ), Absolute values.
        else if( TRIM(ADJUSTL(quantity))  ==  "|kf|" )then
            i_quantity = 100
        else if( TRIM(ADJUSTL(quantity))  ==  "|px|" )then
            i_quantity = 101
        else if( TRIM(ADJUSTL(quantity))  ==  "|py|" )then
            i_quantity = 102
        else if( TRIM(ADJUSTL(quantity))  ==  "|pz|" )then
            i_quantity = 103
        else if( TRIM(ADJUSTL(quantity))  ==  "|y|" )then
            i_quantity = 117
        else if( TRIM(ADJUSTL(quantity))  ==  "|eta|" )then
            i_quantity = 119
        else if( TRIM(ADJUSTL(quantity))  ==  "|vx|" )then
            i_quantity = 126
        else if( TRIM(ADJUSTL(quantity))  ==  "|vy|" )then
            i_quantity = 127
        else if( TRIM(ADJUSTL(quantity))  ==  "|vz|" )then
            i_quantity = 128
        else if( TRIM(ADJUSTL(quantity))  ==  "|rx|" )then
            i_quantity = -101
        else if( TRIM(ADJUSTL(quantity))  ==  "|ry|" )then
            i_quantity = -102
        else if( TRIM(ADJUSTL(quantity))  ==  "|rz|" )then
            i_quantity = -103
!       = 999, null.
        else if( TRIM(ADJUSTL(quantity))  ==  "null" )then
            i_quantity = 999
!       "Number", uses ASCII.
        else if(  IACHAR(char_1) >= 48 .AND. IACHAR(char_1) <= 57 )then
            read( quantity, "(I2)" ) i_quantity
        end if

!       Gets sorting-quantity array.

!       Allocates memory for sorting-quantity.
        allocate( quantity_tmp( N_sort ) )

!       Sets small enough values for i < i_begin.
        if( i_begin > 1 )then
            small_value = -999999999D0
            do i = 1, i_begin-1, 1
                quantity_tmp(i) = small_value
                small_value = small_value + 1D0
            end do
        end if

!       = 0, KF code.
        if( i_quantity == 0 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = K(i,2) * 1D0
            end do
!       = [1,25], PAPYP(i,j).
        else if( i_quantity >= 1 .AND. i_quantity <= 25 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = PAPYP( i, i_quantity )
            end do
!       = [26,28], velocity components.
        else if( i_quantity >= 26 .AND. i_quantity <= 28 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = P( i, i_quantity - 25 ) / P(i,4)
            end do
!       = 29, absolute velocity |v|.
        else if( i_quantity == 29 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = SQRT( (P(i,1)/P(i,4))**2 + &
                     (P(i,2)/P(i,4))**2 + (P(i,3)/P(i,4))**2 )
            end do
!       = [-5, -1], space-time coordinate components and lifetime.
        else if( i_quantity >= -5 .AND. i_quantity <= -1 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = V( i, -i_quantity )
            end do
!       = -6, absolute "space coordinate", displacement, |r|.
        else if( i_quantity == -6 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = SQRT( V(i,1)**2 +  &
                                          V(i,2)**2 + V(i,3)**2 )
            end do
!       = 100, absolute KF.
        else if( i_quantity == 100 )then
            do i= i_begin, i_end, 1
                quantity_tmp( i ) = ABS( K(i,2) )
            end do
!       = [101,125], ABS( PAPYP(i,j) ).
        else if( i_quantity >= 101 .AND. i_quantity <= 125 )then
            do i=i_begin,i_end,1
                quantity_tmp( i ) = ABS( PAPYP( i, i_quantity-100 ) )
            end do
!       = [126,128], absolute velocity components.
        else if( i_quantity >= 126 .AND. i_quantity <= 128 )then
            do i=i_begin,i_end,1
                quantity_tmp( i ) = ABS( P(i, i_quantity-125 ) /P(i,4) )
            end do
!       = [-103,-101], absolute space-time coordinate components.
        else if( i_quantity >= -103 .AND. i_quantity <= -101 )then
            do i=i_begin,i_end,1
                quantity_tmp( i ) = ABS( V(i, -(i_quantity+100) ) )
            end do
        end if

666     continue

!       Allocates memory for the array of the order index.
        allocate( i_order(N_sort) )
        j_begin = 1
        do i = 1, i_end, 1
            i_order( i ) = i
        end do

!       Indexing Quicksort algorithm. "max_to_min" or "min_to_max".
        if( ABS(i_operation) == 1 )then
!       Gets the index in ascending order. (In coales.f90)
            call PASORT_indexx( N_sort, quantity_tmp, i_order )

!       Knuth shuffle algorithm. "random".
        else if( i_operation == 0 )then
            do i = i_end, i_begin+1, -1
                i_exchange = i_begin + FLOOR( (i + 1 - i_begin)*PYR(1) )
                i_tmp = i_order( i_exchange )
                i_order( i_exchange ) = i_order( i )
                i_order( i ) = i_tmp
            end do
        end if

!       Feeds the sorted rearrangement back.
        j_begin = i_begin
        j_end = i_end
        ! if( i_operation == -1 )then
        !     j_begin = i_end
        !     j_end = i_begin
        ! end if
        if(      TRIM(ADJUSTL(name_common))  ==  "pyjets"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K(i_insert,j) = k_tmp( i_new, j )
                    P(i_insert,j) = p_tmp( i_new, j )
                    V(i_insert,j) = v_tmp( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa2"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    ksa(i_insert,j) = K( i_new, j )
                    psa(i_insert,j) = P( i_new, j )
                    vsa(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa36"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kglu(i_insert,j) = K( i_new, j )
                    pglu(i_insert,j) = P( i_new, j )
                    vglu(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa37"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kth(i_insert,j) = K( i_new, j )
                    pth(i_insert,j) = P( i_new, j )
                    vth(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sbe"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kbe(i_insert,j) = K( i_new, j )
                    pbe(i_insert,j) = P( i_new, j )
                    vbe(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "saf"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kaf(i_insert,j) = K( i_new, j )
                    paf(i_insert,j) = P( i_new, j )
                    vaf(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "aaff"  )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kaff(i_insert,j) = K( i_new, j )
                    paff(i_insert,j) = P( i_new, j )
                    vaff(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sbh"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kbh(i_insert,j) = K( i_new, j )
                    pbh(i_insert,j) = P( i_new, j )
                    vbh(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sa1_h" )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kn(i_insert,j) = K( i_new, j )
                    pn(i_insert,j) = P( i_new, j )
                    rn(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sgam"  )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kgam(i_insert,j) = K( i_new, j )
                    pgam(i_insert,j) = P( i_new, j )
                    vgam(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "trs"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    ktrs(i_insert,j) = K( i_new, j )
                    ptrs(i_insert,j) = P( i_new, j )
                    vtrs(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "delt"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    kdel(i_insert,j) = K( i_new, j )
                    pdel(i_insert,j) = P( i_new, j )
                    vdel(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "sin_local"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    ksin(i_insert,j) = K( i_new, j )
                    psin(i_insert,j) = P( i_new, j )
                    vsin(i_insert,j) = V( i_new, j )
                end do
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_py8"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_PY8(i_insert,j) = K( i_new, j )
                    P_PY8(i_insert,j) = P( i_new, j )
                    V_PY8(i_insert,j) = V( i_new, j )
                end do
            end do
            do j=6,7,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_PY8(i_insert,j) = k_tmp( i_new, j )
                    P_PY8(i_insert,j) = p_tmp( i_new, j )
                end do
            end do
            do i=i_begin,i_end,1
                i_new = i_order( i )
                i_insert = i
                if( i_operation == -1 )  i_insert= i_end -i_insert+1
                K_PY8(i_insert,8) = k_tmp( i_new, 8 )
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_aa"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_AA(i_insert,j) = K( i_new, j )
                    P_AA(i_insert,j) = P( i_new, j )
                    V_AA(i_insert,j) = V( i_new, j )
                end do
            end do
            do j=6,7,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_AA(i_insert,j) = k_tmp( i_new, j )
                    P_AA(i_insert,j) = p_tmp( i_new, j )
                end do
            end do
            do i=i_begin,i_end,1
                i_new = i_order( i )
                i_insert = i
                if( i_operation == -1 )  i_insert= i_end -i_insert+1
                K_AA(i_insert,8) = k_tmp( i_new, 8 )
            end do
        else if( TRIM(ADJUSTL(name_common))  ==  "pyjets_fo"   )then
            do j=1,5,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_FO(i_insert,j) = K( i_new, j )
                    P_FO(i_insert,j) = P( i_new, j )
                    V_FO(i_insert,j) = V( i_new, j )
                end do
            end do
            do j=6,7,1
                do i=i_begin,i_end,1
                    i_new = i_order( i )
                    i_insert = i
                    if( i_operation == -1 )  i_insert= i_end -i_insert+1
                    K_FO(i_insert,j) = k_tmp( i_new, j )
                    P_FO(i_insert,j) = p_tmp( i_new, j )
                end do
            end do
            do i=i_begin,i_end,1
                i_new = i_order( i )
                i_insert = i
                if( i_operation == -1 )  i_insert= i_end -i_insert+1
                K_FO(i_insert,8) = k_tmp( i_new, 8 )
            end do
        end if

!       Recovers /PYJETS/ if the incoming COMMOM was not /PYJETS/.
        if( TRIM(ADJUSTL(name_common))  /=  "pyjets" )then
            do j=1,5,1
                do i=1,N,1
                    K(i,j) = k_tmp(i,j)
                    P(i,j) = p_tmp(i,j)
                    V(i,j) = v_tmp(i,j)
                end do
            end do
        end if

!       Releases memory.
        deallocate( k_tmp,  p_tmp,  v_tmp, i_order )
        if( ABS(i_operation) == 1 ) deallocate( quantity_tmp )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE PASORT_indexx(n,arr,indx)
!!      Generates the sorted index using the Insertion sort and
!!       the Quicksort algorithms.
!
!       Indexes an array arr(1:n), i.e., outputs the array indx(1:n)
!        such that arr(indx(j)) is in ascending order for j = 1,2,...,N.
!       The input quantities n and arr are not changed.
!       Parameters: M is the size of subarrays sorted by straight
!        insertion and NSTACK is the required auxiliary storage.
!
!       Cf. NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC
!        COMPUTING (ISBN 0-521-43064-X)
!
        INTEGER n,indx(n),M,NSTACK
        REAL*8 arr(n)
        PARAMETER (M=7,NSTACK=50)
        INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
        REAL*8 a
        do j=1,n
            indx(j)=j
        enddo
        jstack=0
        l=1
        ir=n
!       Insertion sort when subarray small enough.
 1      if(ir-l < M)then
            do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,l,-1
                    if(arr(indx(i)) <= a) goto 2
                    indx(i+1)=indx(i)
                enddo
                i=l-1
 2              indx(i+1)=indxt
            enddo
            if(jstack == 0)then
                return
            end if
            ir=istack(jstack)    ! Pop stack and begin a new round of
            l=istack(jstack-1)   !  partitioning.
            jstack=jstack-2
        else
            k=(l+ir)/2          ! Choose median of left, center and right
            itemp=indx(k)       ! elements as partitioning element a. Also
            indx(k)=indx(l+1)   ! rearrange so that a(l) <= a(l+1) <= a(ir).
            indx(l+1)=itemp
            if(arr(indx(l)) > arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
            endif
            if(arr(indx(l+1)) > arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
            endif
            if(arr(indx(l)) > arr(indx(l+1)))then
                itemp=indx(l)
                indx(l)=indx(l+1)
                indx(l+1)=itemp
            endif
            i=l+1   ! Initialize pointers for partitioning.
            j=ir
            indxt=indx(l+1)
            a=arr(indxt)   ! Partitioning element.
 3          continue       ! Beginning of innermost loop.
            i=i+1          ! Scan up to find element > a
            if(arr(indx(i)) < a) goto 3
 4          continue
            j=j-1          ! Scan down to find element < a.
            if(arr(indx(j)) > a) goto 4
            if(j < i) goto 5! Pointers crossed. Exit with partitioning complete.
            itemp=indx(i)   ! Exchange elements of both arrays.
            indx(i)=indx(j)
            indx(j)=itemp
            goto 3              ! End of innermost loop.
 5          indx(l+1)=indx(j)   ! Insert partitioning element in both arrays.
            indx(j)=indxt
            jstack=jstack+2
!       Push pointers to larger subarray on stack, process smaller
!        subarray immediately.
            if(jstack > NSTACK) write(*,*) 'NSTACK too small in indexx'
            if(ir-i+1 >= j-l)then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
            else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
            endif
        endif
        goto 1


        END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function iRandom( i_begin, i_end )
!       Generates an integer random number between [ i_begin, i_end ].
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        iRandom = i_begin + FLOOR( (i_end + 1 - i_begin)*PYR(1) )


        return
        end



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc