!! sfm_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2024 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, July 2002 - November 2024.

!> This is the program to perform the string fragmentation hadronization by
!!  interfacing with PYTHIA 6/8.

!!                                             By Ben-Hao at CIAE on 31/07/2002
!!                                  Last updated by An-Ke at UiO  on 23/11/2024


        subroutine sfm
!!      Performs the string fragmentation hadronization by calling "PAEXEC".
!
!       It was written by Ben-Hao Sa on 31/07/2002.
!       Its input messages are in 'PYJETS'.
!       Its internal working block is 'sa1_h'.
!       Its output messages are in 'PYJETS' ('sa1_h' the same).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_EXIST
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
        dimension rc(3)


!       Backup.
        MSTP111   = MSTP(111)
        MSTP(111) = 1
        nn = 0
        iikk = 0
        kkii = 0

!       Calls PAEXEC repalacing PYEXEC. (In Pythia8_fort_interface.f90)
        call PAEXEC

!       Recovers.
        MSTP(111) = MSTP111

        if(kkii == 2) return

!       Removes unnecessary entries.
        if( .NOT.IS_PYTHIA8(i_mode) )then
            call PAEDIT(1)
!       Moves non-historical from /PYJETS_PY8/ to /PYJETS/.
        else if( IS_PYTHIA8(i_mode) )then
            N = 0
            do i=1,N_PY8,1
                KS = K_PY8(i,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                N = N + 1
                do j=1,5,1
                    K(N,j) = K_PY8(i,j)
                    P(N,j) = P_PY8(i,j)
                    V(N,j) = V_PY8(i,j)
                end do
            end do
        end if

!       Angantyr has given 4-positions.
        if( i_mode == 8 .OR. i_mode == 9 ) return

!       "PYJETS" to "sa1_h".
        nn = N
        do j2=1,5,1
            do j1=1,N
                kn(j1,j2) = K(j1,j2)
                pn(j1,j2) = P(j1,j2)
                rn(j1,j2) = V(j1,j2)
            end do
        end do
        rrp = 1.16D0

!       Arranges produced particles on the surface of sphere (with radius rrp
!        and centred on parent position), produced particle is put on its
!        parent position originally.
        ipp = 1
        do while(.true.)
            r1  = rn(ipp,1)
            r2  = rn(ipp,2)
            r3  = rn(ipp,3)
            rc(1) = r1
            rc(2) = r2
            rc(3) = r3
!       in corresponding with the time set in "posi"
            rn(ipp,4) = 0D0
            time = rn(ipp,4)
!       Finds out the hadrons with same parent (rc).
!       Note: produced particles with same position are arranged continuously
!        in PYJETS.
            ip = 0
            do j2 = ipp+1, N, 1
                s1 = rn(j2,1)
                s2 = rn(j2,2)
                s3 = rn(j2,3)
                if( ABS(s1-r1) <= 1D-6 .AND. ABS(s2-r2) <= 1D-6 &
                    .AND. ABS(s3-r3) <= 1D-6 ) ip = ip + 1
            end do
            ipp1 = ipp + ip
            call posi(ipp,ipp1,rc,rrp,time)
            ipp = ipp1 + 1
            if( ipp >= N ) exit
        end do

!       Transfers four position messages from "sa1_h" to "PYJETS".
        do j=1,4,1
            do i=1,nn,1
                V(i,j) = rn(i,j)
            end do
        end do
        nn = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine reconstruct_parton_string( time_coord_recover )
!!      Reconstructs the partons as strings.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_DIQUARK
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       Local arrays.
        common/diq_local/ kdiq(kszj,5), dgmas(kszj)


!       Hadronization model.
        i_had_model = INT( adj1(12) )


!-------------------------------------------------------------------------------
!---------------------------   Diquark Recovering   ----------------------------
!       Recovers parton configuration in "sbe" if the diquarks have been broken.

!       Without diquarks breaking-up.
        if( iparres == 2 .OR. iparres == 3 ) goto 1800

!       Loops over "sbe".
        idii = 0
        do i=1,nbe,1
            KS = kbe(i,1)
            kf = kbe(i,2)
            if( IS_DIQUARK(KF) .AND. IS_EXIST(KS,i_mode) )then
                idii = idii + 1
                do j=1,5,1
                    kdiqq = kbe(i,j)
                    kdiq(idii,j) = kdiqq
                end do
                dgmas(idii) = pbe(i,5)
            end if
        end do

!       Loops over "PYJETS"
        idij = 0
        jb = 0
880     continue
        do i = jb+1, N, 1
            jb = jb + 1
            ndiqi = ndiq(i)
            if( ndiqi <= 0 ) cycle
            ! Diquark (anti-diquark).
            idij = idij + 1
            j = npt( ndiqi )
            do i1=1,5
                K(i,i1) = kdiq(idij,i1)
            end do
            do i1=1,3,1
                P(i,i1) = P(i,i1) + P(j,i1)
            end do
            dimass = dgmas(idij)
            pi1 = P(i,1)
            pi2 = P(i,2)
            pi3 = P(i,3)
            pi4 = dsqrt( pi1*pi1 + pi2*pi2 + pi3*pi3 + dimass*dimass )
            ! Collects the lost energy.
            throe_p(4) = throe_p(4) + P(i,4) + P(j,4) - pi4
            P(i,4) = pi4
            P(i,5) = dimass
            ! Sets 3-coordinate of the diquark as one of the
            !  corresponding quarks if it has been broken.
            if( PYR(1) > 0.5D0 )then
                do i1=1,3,1
                    V(i,i1) = V(j,i1)
                end do
            end if
            ! Moves particle list "PYJETS" ("ndiq") one step downward from
            !  j+1 to N.
            do j1 = j+1,N , 1
                ndiq(j1-1) = ndiq(j1)
                do jj=1,5,1
                    K( j1-1, jj ) = K(j1,jj)
                    P( j1-1, jj ) = P(j1,jj)
                    V( j1-1, jj ) = V(j1,jj)
                end do
            end do
            N = N - 1
            ! Subtracts "npt" by one from idij+1 to idi.
            do j1 = idij+1, idi, 1
                npt(j1) = npt(j1) - 1
            end do
            ! Adjusts line numbers of the string list.
            do i_string = 1, nstr1, 1
                if( j < nstr1a(i_string) ) &
                    nstr1a(i_string) = nstr1a(i_string) - 1
                if( j < nstr1v(i_string) ) &
                    nstr1v(i_string) = nstr1v(i_string) - 1
            end do
            goto 880
        end do

1800    continue
!       Assigns 4-coordinate of random one of the rest partons (or one end?)
!        to the first parton in the string.
!#TODO(Lei20240219): not work well for PYTHIA 8, need improvement.
        if( .NOT.IS_PYTHIA8(i_mode) )then
            do i1=1,nstr1,1
                i_string_A = nstr1a(i1)
                i_string_V = nstr1v(i1)
                do while(.true.)
                    i_A = INT( PYR(1)*( i_string_V -i_string_A +1) +i_string_A )
                    ! Excludes junctions.
                    if( K( i_A, 2 ) /= 88 ) exit
                end do
                ! if( PYR(1) > 0.5D0 ) i_A = i_string_V
                do j1=1,5,1
                    V( i_string_A, j1 ) = V( i_A, j1 )
                end do
            end do
        end if
!       Recovers the 4-coordinate from 'parini'. (not ?) This treatment is
!        equivalent to giving medium correction in momentum space.
        i_coord_recover = 1
        if( ( i_had_model == 0 .OR. i_had_model == 3 ) .AND. &
            i_coord_recover == 1 )then
            do i2=1,5,1
                do i1=1,N,1
                    V(i1,i2) = vbe(i1,i2)
                end do
            end do
            time_coord_recover = 0D0
        end if
        nbe = 0
!---------------------------   Diquark Recovering   ----------------------------
!-------------------------------------------------------------------------------


!       Puts particles on-shell by hand to avoid potential errors.
        do i=1,N,1
            KS = K(i,1)
            KF = K(i,2)
            ! Erases the bad color flow information for ion-induced collisions.
            if( kjp22 < 5 .AND. ( IS_NUCLEUS(KF_proj) &
                             .OR. IS_NUCLEUS(KF_targ) ) )then
                do j=3,5,1
                    K(i,j) = 0
                end do
                cycle
            end if
            if( .NOT.IS_EXIST(KS,i_mode) .OR. IS_NUCLEUS(KF) &
                .OR. KF == 88 ) cycle
            px0 = P(i,1)
            py0 = P(i,2)
            pz0 = P(i,3)
            E0  = P(i,4)
            dm  = P(i,5)
            P(i,4) = SQRT( px0**2 + py0**2 + pz0**2 + dm**2 )
            throe_p(4) = throe_p(4) + E0 - P(i,4)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine calculate_string_tension( i_begin, i_end, i_tension )
!!      Calculates the effective string tension and takes the popcorn mechanism
!!       correction into account.
!       i_tension: = 1, effective single string tension with the gluon
!                       wrinkle correction
!                  = 2, effective multiple string tension considering the
!                       string overlapping.
!                  = 3, gluon wrinkle correction + string overlapping
!       i_begin, e_end: the endpoints of a single string.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/ctllist/npctl,npinel(600),npctl0,npctlm
!       For the calculation of the single, multiple string tension.
        common/string_kapa1/ parj10, ww0, seco
!       For the independent Optical Clauber calculation.
        common/Glauber5/ pathn, evbin, pirr, tirr, spathn, sevbin
!       For the parini related statistics.
        common/anly_parini1/ woun, swoun, snpctl0, snpctlm, snpar
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       Local arrays.
        common/sin_local/nsin,non,ksin(kszj,5),psin(kszj,5),vsin(kszj,5)
!       pathn: N_bin / N_part
!       itime: number of strings in current event, when i_tension = 1 or 3
!       gtime: number of gluons in current event, when i_tension = 1 or 3
!       akapa(1): parj(1)
!       akapa(2): parj(2)
!       akapa(3): parj(3)
!       akapa(4): parj(4)
!       akapa(5): parj(21)
!       akapa(6): effective string tension


        i1 = i_begin
        i2 = i_end
        ckapa = 1D0


!-------------------------------------------------------------------------------
!-------------------------   Multiple String Tension   -------------------------
        if( i_tension == 2 .OR. i_tension == 3 )then
!       calculate the multiple effective string tension and PARJ(1) etc.

!----------------------   String Overlapping Correction   ----------------------
            ampi = MINT(31)
!#TODO(Lei20241003): consider the different n_MPI for the diffenrent NN pair or
!                    the avarage n_MPI?
!           note MINT(31) = 0 for a low_pT event
!           pathn = ( npinel(1) + npinel(592) ) / pirr
!           pathn = ( npinel(1) + npinel(592) ) / ( 0.5*woun )
!           numerator: NN collision # calculated
!           denominator: N_part of projectile nucleus calculated (pirr:
!            from Glauber model)
            if(ampi <= 0)then
                ckapa = 1D0
            elseif(ampi > 0)then
!               ampi = ampi*evbin
                if( (ipden == 0 .and. itden == 0) .or.  &
                    (ipden == 2 .and. itden == 2) )then
                    pathn = 1.
                else
                    pathn = evbin / ( pirr + tirr )
                end if
!               evbin: N_bin of collision system (A+B)
!               pirr: N_part of projectile nucleus (Glauber)
!               tirr: N_part of target nucleus (Glauber)
                ampi  = ampi*pathn
                ckapa = ( 1D0 + (ampi - 1D0) / (1D0 + 1D0/(cp0**2) ) )**cr0
!               ckapa = 1D0 + (ampi - 1D0)
!               denoc = 1D0 + 1D0/(cp0**2)
!               if(denoc < 1D-20) denoc = 1D-20
!               ckapa = ckapa/denoc
!               if(ckapa < 1D-20) ckapa = 1D-20
!               ckapa = (ckapa)**cr0
            end if
!           ckapa is multiple string tension
!           string tension of the pure qqbar string, kapa0, is assumed to be 1
!----------------------   String Overlapping Correction   ----------------------

            if(i_tension == 2)then
                parj_2   = parj2**( 1D0/ckapa )
                parj_21  = parj21*( ( ckapa/1D0 )**(0.5D0) )
                parj_1   = parj1**( 1D0/ckapa )
                parj_3   = parj3**( 1D0/ckapa )
                parj_4   = parj4**( 1D0/ckapa )
                PARJ(1)  = parj_1
                PARJ(2)  = parj_2
                PARJ(3)  = parj_3
                PARJ(4)  = parj_4
                PARJ(21) = parj_21
!               recalculate PARJ(1) with popcorn mechanism correction
                wxef     = parj_3
                wyef     = parj_4
                wrhoef   = parj_2
                wnumer   = 1D0 + 2D0*wxef*wrhoef + 9D0*wyef &
                            + 6D0 * wxef * wyef * wrhoef       &
                            + 3D0 * wxef*wxef * wyef * wrhoef*wrhoef
                wdenom   = 2D0 + wrhoef
                wweff    = wnumer / wdenom
                PARJ(1)  = seco * wweff * ( parj10/seco/ww0 )**( 1D0/ckapa )
                akapa(1) = PARJ(1)
                akapa(2) = PARJ(2)
                akapa(3) = PARJ(3)
                akapa(4) = PARJ(4)
                akapa(5) = PARJ(21)
                akapa(6) = ckapa
            end if
        end if
!-------------------------   Multiple String Tension   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-----------------   Single / Single+Multiple String Tension   -----------------
        if( i_tension == 1 .OR. i_tension == 3 )then

!------------------------   Gluon Wrinkle Correction  --------------------------
            ! parameter alpha
            vfr24  = 3.5D0
            ! sqrt( s_0 ) in GeV
            vfr25  = 0.8D0
            vfr252 = vfr25*vfr25
            toteng = 0D0
            toten  = 0D0
            totglukt = 0D0
            pmax = 0D0
            ggg  = 0D0
            do i3=i1,i2
                ! root_s, string total energy
                toten = toten + psin(i3,4)
                ! pT*pT
                pp2 = psin(i3,1)**2 + psin(i3,2)**2
                ! pT (kT)
                ppp = dsqrt(pp2)
                if( ksin(i3,2) == 21 )then
                    ! k_{Tmax}
                    if( ppp > pmax ) pmax = ppp
                    if( pp2 >= vfr252 )then
                        ! sum over gluons in a string
                        toteng = toteng + dlog( pp2/vfr252 )
                        ggg = ggg + 1D0
                    endif
                endif
            enddo
            pmax2 = pmax*pmax
            ! numerator
            if( pmax2 >= vfr252 ) totglukt = totglukt + dlog(pmax2/vfr252)
            ! denominator
            sss = dlog( toten*toten/vfr252 ) + toteng
            div = totglukt/sss
!           div: factor related to number of gluons and hardest gluon in
!            current string
!           pmax: transverse momentum of hardest gluon in current string
            adiv  = adiv  + div
            gpmax = gpmax + pmax
!           string tension of the pure qqbar string, kapa0, is assumed to be 1
!           calculate kapa and PARJ(1) etc. of current string
            effk2 = (1D0 - div)**( -vfr24 )
            itime = itime + 1
            gtime = gtime + ggg
!------------------------   Gluon Wrinkle Correction  --------------------------

!       single
            if(i_tension == 1)then
                parj_2   = parj2**( 1D0/effk2 )
                parj_21  = parj21*( ( effk2/1D0 )**(0.5D0) )
                parj_1   = parj1**( 1D0/effk2 )
                parj_3   = parj3**( 1D0/effk2 )
                parj_4   = parj4**( 1D0/effk2 )
                PARJ(1)  = parj_1
                PARJ(2)  = parj_2
                PARJ(3)  = parj_3
                PARJ(4)  = parj_4
                PARJ(21) = parj_21
!               recalculate PARJ(1) with popcorn mechanism correction
                wxef     = parj_3
                wyef     = parj_4
                wrhoef   = parj_2
                wnumer   = 1D0 + 2D0*wxef*wrhoef + 9D0*wyef &
                         + 6D0 * wxef * wyef * wrhoef       &
                         + 3D0 * wxef*wxef * wyef * wrhoef*wrhoef
                wdenom   = 2D0 + wrhoef
                wweff    = wnumer / wdenom
                PARJ(1)  = seco*wweff * ( parj10/seco/ww0 )**( 1D0/effk2 )
!               sum over strings
                akapa(1) = akapa(1) + PARJ(1)
                akapa(2) = akapa(2) + parj_2
                akapa(3) = akapa(3) + parj_3
                akapa(4) = akapa(4) + parj_4
                akapa(5) = akapa(5) + parj_21
                akapa(6) = akapa(6) + effk2
            end if

!       single + multiple
            if(i_tension == 3)then
!           ckapa is multiple string tension
                parj_2   = parj2**( 1D0/( ckapa*effk2 ) )
                parj_21  = parj21*( ( (effk2*ckapa) / 1D0 )**(0.5D0) )
                parj_1   = parj1**( 1D0 / ( ckapa*effk2 ) )
                parj_3   = parj3**( 1D0 / ( ckapa*effk2 ) )
                parj_4   = parj4**( 1D0 / ( ckapa*effk2 ) )
                PARJ(1)  = parj_1
                PARJ(2)  = parj_2
                PARJ(3)  = parj_3
                PARJ(4)  = parj_4
                PARJ(21) = parj_21
!               recalculate parj(1) with popcorn mechanism correction
                wxef     = parj_3
                wyef     = parj_4
                wrhoef   = parj_2
                wnumer   = 1D0 + 2D0*wxef*wrhoef + 9D0*wyef &
                         + 6D0 * wxef * wyef * wrhoef       &
                         + 3D0 * wxef*wxef * wyef * wrhoef*wrhoef
                wdenom   = 2D0 + wrhoef
                wweff    = wnumer / wdenom
                PARJ(1)  = seco * wweff &
                         * ( parj10/seco/ww0 )**( 1D0 / ( ckapa + effk2 ) )
!               sum over strings
                akapa(1) = akapa(1) + parj(1)
                akapa(2) = akapa(2) + parj_2
                akapa(3) = akapa(3) + parj_3
                akapa(4) = akapa(4) + parj_4
                akapa(5) = akapa(5) + parj_21
                akapa(6) = akapa(6) + ckapa*effk2
            end if
        end if
!-----------------   Single / Single+Multiple String Tension   -----------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_fragmentation_with_constant_tension( i_error )
!!      Fragments the string system with the constant string tension.
!       i_error: 0/1, errors occured and throws away the current event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa34/itorw,iikk,cp0,cr0,kkii
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension ps0(6), ps1(6)


        i_error = 0
        ps0 = 0D0
        ps1 = 0D0
!       Total 4-momentum before "sfm".
        do jj=1,4,1
            do ii=1,N,1
                KS = K(ii,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                ps0(jj) = ps0(jj) + P(ii,jj)
            end do
        end do

        kkii = 0

!       fragments strings all at once (In sfm.f90.)
        call calculate_string_tension( 1, N, kjp22 )
        call sfm

!        Throw away the current event if any errors exist.
        ! if( MSTU(23) > 0 .OR. MSTU(30) > 0 )then
        ! if( MSTU(24) > 0 )then
        !     iii = iii - 1
        !     goto 300
        ! end if
        if(kkii == 2)then
            ! throw away current event
            i_error = 1
            return
        end if

!       Total 4-momentum after "sfm".
        do jj=1,4,1
            do ii=1,N,1
                KS = K(ii,1)
                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                ps1(jj) = ps1(jj) + P(ii,jj)
            end do
            ! Collects lost 4-momentum.
            throe_p(jj) = throe_p(jj) + ps0(jj) - ps1(jj)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_fragmentation_with_fluctuating_tension( i_error )
!!      Fragments the string system string-by-string with the fluctuated
!!       string tension.
!       i_error: 0/1, errors occured and throws away the current event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
         parj21,parj4,adiv,gpmax,nnc
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       Local arrays.
        common/sin_local/nsin,non,ksin(kszj,5),psin(kszj,5),vsin(kszj,5)
        dimension ps0(6), ps1(6)
!       itime: number of strings in current event, when kjp22 = 1 or 3
!       astr: not used
!       gtime: number of gluons in current event, when kjp22 = 1 or 3
!       akapa(1): parj(1)
!       akapa(2): parj(2)
!       akapa(3): parj(3)
!       akapa(4): parj(4)
!       akapa(5): parj(21)
!       akapa(6): effective string tension


!-------------------------------------------------------------------------------
!----------------------   Tension Counters Initializing   ----------------------
        nstr  = 0
        itime = 0
        gtime = 0D0
        adiv  = 0D0
        gpmax = 0D0
        akapa = 0D0
!----------------------   Tension Counters Initializing   ----------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------------   String Fragmentation   ---------------------------
!       find string and line number of its first and last components,
!        calculate PARJ(1) etc. if kjp22=1 or 3, and then hadronize string

!       'PYJETS' to 'sin'
        nsin = N
        do i2=1,5,1
            do i1=1,N,1
                ksin(i1,i2) = K(i1,i2)
                psin(i1,i2) = P(i1,i2)
                vsin(i1,i2) = V(i1,i2)
            end do
        end do
        N = 0
        nbe  = 0
        naff = 0
        i_error = 0


!       loop over string (begin)
10001   do i1=1,nsin,1
!       find a string
            ! i1 is 'A'
            if( ksin(i1,1) == 2 )then
                do i2=i1+1,nsin
                    ! i2 is 'V'
                    if( ksin(i2,1) == 1 )then

!---------------------   Single String Locating   ----------------------
                        nstr = nstr + 1
                        ! line number of first component of nstr-th string
                        nstra(nstr) = i1
                        ! line number of last component of nstr-th string
                        nstrv(nstr) = i2
                        call calculate_string_tension( i1, i2, kjp22 )
!       "sin" to "PYJETS" (the part of current string in "sin" to "PYJETS").
                        n = i2 - i1 + 1
                        do ii1=i1,i2,1
                            ii3=ii1-i1+1
                            do ii2=1,5,1
                                K(ii3,ii2) = ksin(ii1,ii2)
                                P(ii3,ii2) = psin(ii1,ii2)
                                V(ii3,ii2) = vsin(ii1,ii2)
                            end do
                        end do
!---------------------   Single String Locating   ----------------------

!------------------------   Single String Fragmentation ------------------------
!       hadronization of current string
                        iikk = 0
                        kkii = 0
                        ps0  = 0D0
                        ps1  = 0D0
!       Total 4-momentum before "sfm".
                        do jj=1,4,1
                            do ii=1,N,1
                                KS = K(ii,1)
                                if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                                ps0(jj) = ps0(jj) + P(ii,jj)
                            end do
                        end do

                        call sfm

!                       Throw away the current event if any errors exist.
                        ! if( MSTU(23) > 0 .OR. MSTU(30) > 0 )then
                        ! if( MSTU(24) > 0 )then
                        !     iii = iii - 1
                        !     goto 300
                        ! end if
!------------------------   Single String Fragmentation ------------------------

!-------------------------   Failed String Treating   --------------------------
                        if( iikk == 2 .and. &
                          ( (ipden == 0 .and. itden == 0) .or. &
                            (ipden == 2 .and. itden == 2) .or. &
                             ipden >= 11 ) )then
                            ! Throws away this event.
                            i_error = 1
                            return
                        else if(iikk == 2 .OR. kkii == 2)then
!                       moves the part of current string in "sin" to "sbe"
                            do ii1=i1,i2,1
                                nbe = nbe + 1
                                do ii2=1,5,1
                                    kbe(nbe,ii2) = ksin(ii1,ii2)
                                    pbe(nbe,ii2) = psin(ii1,ii2)
                                    vbe(nbe,ii2) = vsin(ii1,ii2)
                                end do
                            end do
                            N = 0
                        end if
!-------------------------   Failed String Treating   --------------------------

!-----------------------   Lost 4-momentum Collecting   ------------------------
                        if(n > 0)then
!       Total 4-momentum after 'sfm'
                            do jj=1,4,1
                                do ii=1,N,1
                                    KS = K(ii,1)
                                    if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                                    ps1(jj) = ps1(jj) + P(ii,jj)
                                end do
                                ! Collects lost 4-momentum.
                                throe_p(jj) = throe_p(jj) + ps0(jj) - ps1(jj)
                            end do
                        end if
!-----------------------   Lost 4-momentum Collecting   ------------------------

!----------------------------   Geberated Hadrons   ----------------------------
!       'PYJETS' to 'aaff'
                        if(n > 0)then
                            do ii1=1,N,1
                                ii3 = naff + ii1
                                do ii2=1,5,1
                                    kaff(ii3,ii2) = K(ii1,ii2)
                                    paff(ii3,ii2) = P(ii1,ii2)
                                    vaff(ii3,ii2) = V(ii1,ii2)
                                end do
                            end do
                            naff = ii3
                        end if
!----------------------------   Geberated Hadrons   ----------------------------

!--------------------------   String List Updating   ---------------------------
                        if(i2 < nsin)then
!       revamps 'sin', i.e. moves parton list 'sin' ii (= i2 - i1 + 1) steps
!        downward from i2+1 to nsin
                            ii = i2 - i1 + 1
                            do jj=1,5,1
                                do j = i2+1, nsin, 1
                                    ksin(j-ii,jj) = ksin(j,jj)
                                    psin(j-ii,jj) = psin(j,jj)
                                    vsin(j-ii,jj) = vsin(j,jj)
                                end do
                            end do
                            nsin = nsin - ii
                            goto 10001
                        end if
!--------------------------   String List Updating   ---------------------------

                        ! without rest partons
                        if(i2 == nsin) goto 20001
                    end if
                end do
            end if
        end do

!--------------------------   Rest Parton Dumping   ----------------------------
!       rest partons which can not compose a string
        if(nsin >= 1)then
!       "sin" to "sbe"
            do ii1=1,nsin,1
                nbe = nbe + 1
                do ii2=1,5,1
                    kbe(nbe,ii2) = ksin(ii1,ii2)
                    pbe(nbe,ii2) = psin(ii1,ii2)
                    vbe(nbe,ii2) = vsin(ii1,ii2)
                end do
            end do
        end if
!--------------------------   Rest Parton Dumping   ----------------------------

!       loop over string endded
!       fragments string by string endded

20001   continue
!--------------------------   String Fragmentation   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!----------------------------   Tension Counting   -----------------------------
!       average over strings in current hh collision
        atime = dfloat( itime )
!       itime, # of strings in current hh collision
        if(atime > 0D0)then
            do i1=1,6,1
                akapa(i1) = akapa(i1) / atime
            end do
            gtime = gtime / atime
            adiv  = adiv  / atime
            gpmax = gpmax / atime
!       gtime: averaged # of gluons in a string in current hh collision
        end if
!----------------------------   Tension Counting   -----------------------------
!-------------------------------------------------------------------------------


!       "aaff" to "PYJETS"
        N = naff
        do ii2=1,5,1
            do ii1=1,N,1
                K(ii1,ii2) = kaff(ii1,ii2)
                P(ii1,ii2) = paff(ii1,ii2)
                V(ii1,ii2) = vaff(ii1,ii2)
            end do
        end do
        naff = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine do_fragmentation_pair_by_pair( i_error )
!!      Fragments strings one NN (ll, lh/hl) pair by pair for the complete CR.
!       NOTE: one MUST NOT remove any entries before while this fragmentation
!             mode on.
!       /PYJETS_AA/ is the full input (strings of all NN pairs).
!       /aaff/ is the full output (results of all NN pairs).
!       /PYJETS/ or "PYJETS_PY8" is single input and output.
!       i_error: 0/1, errors occured and throws away the current event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_NUCLEUS,IS_PARTON
        PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sa34/itorw,iikk,cp0,cr0,kkii
        common/sbe/nbe,non_be,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/aaff/naff,nonff,kaff(kszj,5),paff(kszj,5),vaff(kszj,5)
!       Arrays of particle information of one NN pair for PYTHIA 8.
        COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
               K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
!       Arrays of particle information of a AA collision (total NN pairs).
        COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
               K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension ps0(6), ps1(6)


!       "PYJETS" to "PYJETS_AA", uodates the coordinate and momentum information
!        after the partonic rescattering using the "pointers" V(i,5).
        ! Only to (4).
        do j=1,4,1
            do i=1,N,1
                i_org = INT( V(i,5) )
                P_AA( i_org, j ) = P(i,j)
                V_AA( i_org, j ) = V(i,j)
            end do
        end do
        N = 0
        nbe = 0
        naff = 0
        i_error = 0

!       Loops over NN pairs.
!       Inputs the strings of one NN pair and performs the fragmentation.
        n_NN_pair = N_COLL_NN_AA
        i_begin = 0
        i_end   = 0
        do i_NN = 1, n_NN_pair, 1
            iikk = 0
            kkii = 0
            ps0  = 0D0
            ps1  = 0D0
            ! Updaes i_begin.
            i_begin = i_end + 1
            ! Required for both PYTHIA 6 and 8.
            N = N_AA( i_NN )
            ! Required for PYTHIA 8.
            N_PY8 = N_AA( i_NN )
            ! Updaes i_end.
            i_end   = i_end + N
!           Feeds data of one NN pair from /PYJETS_AA/ into /PYJETS/ , etc.
            ! In Pythia8_fort_interface.f90.
            call PAFEED( i_NN, i_begin, i_end, ps0 )

!           Does string fragmentation.
            call sfm

!           Throws away the current event if any errors exist.
            ! if( MSTU(23) > 0 .OR. MSTU(30) > 0 )then
            ! if( MSTU(24) > 0 )then
            !     iii = iii - 1
            !     goto 300
            ! end if

!--------------------------   Failed Pair Treating   ---------------------------
!           If non-culeus-induced collisions, throws away the event.
            if( (iikk == 2 .OR. kkii == 2) .AND. .NOT.IS_NUCLEUS(KF_proj) &
                                           .AND. .NOT.IS_NUCLEUS(KF_targ) )then
                    i_error = 1
                    return
!           Other cases, appends the existed partons, non-partons of current
!            NN pair in "PYJETS_AA" to "sbe" and "sbh", respectively.
            else if( iikk == 2 .OR. kkii == 2 )then
                do ii1 = i_begin, i_end, 1
!                   Only for the particles that still exist.
                    KS = K_AA(ii1,1)
                    KF = K_AA(ii1,2)
                    if( .NOT.IS_EXIST(KS,i_mode) ) cycle
                    if( IS_PARTON(KF) )then
                        nbe = nbe + 1
                        do ii2=1,5,1
                            kbe(nbe,ii2) = K_AA(ii1,ii2)
                            pbe(nbe,ii2) = P_AA(ii1,ii2)
                            vbe(nbe,ii2) = V_AA(ii1,ii2)
                        end do
                    else
                        nbh = nbh + 1
                        do ii2=1,5,1
                            kbh(nbh,ii2) = K_AA(ii1,ii2)
                            pbh(nbh,ii2) = P_AA(ii1,ii2)
                            vbh(nbh,ii2) = V_AA(ii1,ii2)
                        end do
                    end if
                end do
                N = 0
                N_PY8 = 0
!--------------------------   Failed Pair Treating   ---------------------------

!------------------------   Successful Pair Treating   -------------------------
            else if( N > 0 )then
!               Calculates the 4-momentum after sfm.
                do i1 = 1, N, 1
                    KS = K(i1,1)
                    if( IS_EXIST(KS,i_mode) )then
                        do i2=1,5,1
                            ps1(i2) = ps1(i2) + P(i1,i2)
                            ! Sets it as 1.
                            ! if( kjp21 == 1 ) K(i1,1) = 1
                        end do
                    end if
                end do
                do i1=1,4,1
                    throe_p(i1) = throe_p(i1) + ps0(i1) - ps1(i1)
                end do
!               Appends "PYJETS" to "aaff".
                do i1=1,N,1
                    naff = naff + 1
                    do i2=1,5
                        kaff(naff,i2) = K(i1,i2)
                        paff(naff,i2) = P(i1,i2)
                        vaff(naff,i2) = V(i1,i2)
                    end do
                end do
            end if
!------------------------   Successful Pair Treating   -------------------------

        end do

!       "aaff" to "PYJETS"
        N = naff
        do i2=1,5,1
            do i1=1,N,1
                K(i1,i2) = kaff(i1,i2)
                P(i1,i2) = paff(i1,i2)
                V(i1,i2) = vaff(i1,i2)
            end do
        end do
        naff = 0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine decayh(rrp)
!!      Decay of unstable hadrons.
!       Its input messages are in 'PYJETS'
!       Its internal wooking block is 'sa1_h'
!       Its output message is in 'PYJETS' ('sa1_h' the same)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST, IS_NUCLEUS
        PARAMETER (KSZJ=80000)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/sa6_p/ithroq_p,ithrob_p,ich_p,non6_p,throe_p(4)
        common/sgam/ngam,nongam,kgam(kszj,5),pgam(kszj,5),vgam(kszj,5)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
        dimension rc(3), ps0(4), ps1(4)


!       Backup.
        MSTJ1   = MSTJ(1)
        MSTJ21  = MSTJ(21)
        MSTP111 = MSTP(111)
        MSTJ(1)   = 1
        MSTJ(21)  = 2
        MSTP(111) = 1

!       Puts particles on-shell by hand to avoid potential errors.
        do i=1,N,1
            KS0 = K(i,1)
            KF0 = K(i,2)
            if( .NOT.IS_EXIST(KS0,i_mode) .OR. IS_NUCLEUS(KF0) ) cycle
            px0 = P(i,1)
            py0 = P(i,2)
            pz0 = P(i,3)
            E0  = P(i,4)
            dm  = P(i,5)
            P(i,4) = SQRT( px0**2 + py0**2 + pz0**2 + dm**2 )
            throe_p(4) = throe_p(4) + E0 - P(i,4)
        end do
!       Sums of px, py, pz, and E before decay.
        ps0 = 0D0
        do i=1,4,1
            ps0(i) = PAPYP(0,i)
        end do

!       Performs decay. Call PYEXEC to consider the shower off the decay
!        production and the fragmentation of partons from the hadronic decay.
        call PYEXEC
        call PYEDIT(1)

!       Sums of px, py, pz, and E after decay.
        ps1 = 0D0
        do i=1,4,1
            ps1(i) = PAPYP(0,i)
        end do
!       Collects lost 4-momentum.
        do i=1,4,1
            throe_p(i) = throe_p(i) + ps0(i) - ps1(i)
        end do
!       Moves gammas from "PYJETS" to "sgam".
        call remo_gam(22)
!       "PYJETS" to "sa1_h".
        nn = N
        do j2=1,5,1
            do j1=1,N,1
                kn(j1,j2) = K(j1,j2)
                pn(j1,j2) = P(j1,j2)
                rn(j1,j2) = V(j1,j2)
            end do
        end do

!       Arranges decayed particles on the surface of sphere (with radius rrp
!        and centred on parent position), produced particle is put on its
!        parent position originally.
        i_begin = 1
        do while(.true.)
            r1  = rn( i_begin, 1 )
            r2  = rn( i_begin, 2 )
            r3  = rn( i_begin, 3 )
!       In corresponding with the time set in "posi".
            ! rn( i_begin, 4 ) = 0D0
            time  = rn( i_begin, 4 )
            rc(1) = r1
            rc(2) = r2
            rc(3) = r3
!       Finds out the hadrons with same parent (rc).
!       Note: produced particles with same position are arranged continuously
!        in PYJETS.
            n_same_parent = 0
            do j2 = i_begin+1, N, 1
                s1 = rn(j2,1)
                s2 = rn(j2,2)
                s3 = rn(j2,3)
                if( ABS(s1-r1) <= 1D-6 .AND. ABS(s2-r2) <= 1D-6 &
                    .AND. ABS(s3-r3) <= 1D-6 ) n_same_parent = n_same_parent + 1
            end do
            i_end = i_begin + n_same_parent
            call posi( i_begin, i_end, rc, rrp, time )
            i_begin = i_end + 1
            if( i_begin >= N ) exit
        end do

!       Transfers four position messages from "sa1_h" to "PYJETS".
        do j=1,4,1
            do i=1,nn,1
                V(i,j) = rn(i,j)
            end do
        end do
        nn = 0

!       Recovers.
        MSTJ(1)   = MSTJ1
        MSTJ(21)  = MSTJ21
        MSTP(111) = MSTP111


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updah(j2,jc,i)
!!      Moves the hadron list i steps downward from jc till j2.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)


        do jj=1,5,1
            do j=jc,j2,1
                kn( j-i, jj ) = kn(j,jj)
                pn( j-i, jj ) = pn(j,jj)
                rn( j-i, jj ) = rn(j,jj)
            end do
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine posi( nn1, nn2, rc, rrp, time )
!!      Arranges produced particles (from nn1+1 to nn2) on the surface of
!!       sphere with radius rrp and centred on parent position
!       rc : the coordinate of center of the parent
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        common/posi_local/ rr(3,kszj)
        dimension rc(3)


        iii = 0
        do ii = nn1+1, nn2, 1
200         call samp1( rrp, ii )
            if(ii == 1) cycle
            if(  ABS( kn(ii,2) ) < 1000 .OR. ( ABS( kn(ii,2) ) > 10000 &
                .AND. kn(ii,2) /= 100443 .AND. kn(ii,2) /= 100553      &
                .AND. kn(ii,2) /= 10441  .AND. kn(ii,2) /= 10551       &
                .AND. kn(ii,2) /= 20443  .AND. kn(ii,2) /= 20553 ) )    cycle
            do jj = 1, ii-1, 1
                if(  ABS( kn(jj,2) ) < 1000 .OR. ( ABS( kn(jj,2) ) > 10000 &
                    .AND. kn(jj,2) /= 100443 .AND. kn(jj,2) /= 100553      &
                    .AND. kn(jj,2) /= 10441  .AND. kn(jj,2) /= 10551       &
                    .AND. kn(jj,2) /= 20443  .AND. kn(jj,2) /= 20553 ) )  cycle
                dcc = dsqrt( ( rr(1,jj) - rr(1,ii) )**2 &
                    +        ( rr(2,jj) - rr(2,ii) )**2 &
                    +        ( rr(3,jj) - rr(3,ii) )**2 )
!       Non-overlapping demand among baryons. 0.8 is the hard core of a baryon.
                if( dcc < 0.8D0 )then
                    iii = iii + 1
                    if( iii < 10000 ) goto 200
                    exit
                end if
            end do
        end do

        do ii = nn1+1, nn2, 1
!       Gives zero to the time of hadronized and decayed particles
!           rn(ii,4) = 0D0
!       Gives time.
            rn(ii,4) = time
            do jj=1,3
                rn(ii,jj) = rr(jj,ii) + rc(jj)
            end do
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine samp1( xf, i )
!!      Gives the position on the surface of sphere with
!!       radius xf, to particle i.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/posi_local/ rr(3,kszj)


        cita = 2D0 * PYR(1) - 1D0
        fi   = 2D0 * 3.1416 * PYR(1)
        sita = SQRT( 1D0 - cita**2 )
        rr(1,i) = xf * sita * COS(fi)
        rr(2,i) = xf * sita * SIN(fi)
        rr(3,i) = xf * cita


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa1_h(nn1)
!!      Prints particle list and the sum of four momentum and charge.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYCHGE
        PARAMETER (KSZJ=80000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
        dimension peo(4)


        do i=1,nn1
            write(MSTU(11),*) i, kn(i,2), ( pn(i,j), j=1,4,1 )
        end do
        call psum(pn,1,nn1,peo)
        ich1 = 0D0
        do i1=1,nn1,1
            kf = kn(i1,2)
            ich1 = ich1 + PYCHGE(kf)
        end do
        write(MSTU(11),*) peo, ich1 / 3D0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updad_pyj(j2,j1,i)
!!      Moves the parton list (PYJETS) i steps downward from j1 to j2.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)


        do jj=1,5,1
            do j=j1,j2,1
                K( j-i, jj ) = K(j,jj)
                P( j-i, jj ) = P(j,jj)
                V( j-i, jj ) = V(j,jj)
            end do
        end do
!       N = N - i


        return
        end



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc