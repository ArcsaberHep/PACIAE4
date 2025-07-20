!! Parcas_40.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2025 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: Ben-Hao Sa, November 2002 - July 2025.

!> This is the program to deal with the parton cascade (partonic rescattering).

!!                                             By Ben-Hao at CIAE on 19/11/2002
!!                                  Last updated by An-Ke at UiO  on 20/07/2025


        subroutine parcas( time_par, iijk )
!!      Deals with parton cascade (partonic rescattering).
!
!       It was written by Ben-Hao Sa on 19/11/2002.
!       Its input messages are in 'PYJETS'.
!       Its working block is 'PYJETS'.
!       Its output messages are in 'PYJETS'.
!       jjj: jjj-th loop interaction in a event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_EXIST,IS_PARTON,IS_DIQUARK
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/work7/reac(20),crose(20)
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
         nap,nat,nzp,nzt,pio
!       For the gamma related statistics.
        common/anly_gamma1/ egam,  egam1,  egam2,  egam3, &
                           segam, segam1, segam2, segam3
!       For the simulation control.
        COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
               KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!       Arrays of freeze-out particle information.
        COMMON/PYJETS_FO/ N_FO, NPAD_FO, &
                          K_FO(KSZJ,8), P_FO(KSZJ,7), V_FO(KSZJ,5)
        common/e_loss1/ rpo(kszj,5),ppo(kszj,5)
!       icol : current total number of collision pairs in collision time list
!       lc   : line number (in [1,n]) of colliding pair, eg.
!       lc(1,100): line number of first particle of 100-th colliding pair
!       tc   : the collision time of colliding pair
!       ic,jc: line number of colliding particles
!       n0: the n before current collision
!       pi,pj: four momentum of colliding particles
!       taup(i) : formation time of particle i.
!       ishp(i)=1 if i-th particle inside the simulated volume
!              =0 if i-th particle outside the simulated volume
!       rpo,ppo: parton four coordinate before Newton motion, four momentum
!        before energy loss


        iijk = 0
        time_par = 0D0
        if( ABS( adj1(1) ) <= 1D-15 ) return
        if( N < 2 ) return
        time = time_par

!       Initialization.
        call reset_eve
        do i1=1,N,1
            V(i1,4)  = 0D0
            taup(i1) = 0D0
        end do


!-------------------------------------------------------------------------------
!-----------------------------   Parton Checking   -----------------------------
!       Checks and excludes non-parton, historical entries and the
!        partons outside the simulation range.
        dpmax = adj1(27)
        drmax = PARAM(10)*dmax1(rnt,rnp)
        do i1=1,N,1
            ishp(i1) = 0
            KS = K(i1,1)
            KF = K(i1,2)
!#TODO(Lei20241018): diquark collision in the future?
            if( .NOT.IS_EXIST(KS,i_mode) .OR. .NOT.IS_PARTON(KF) &
                .OR. IS_DIQUARK(KF) ) cycle
            ppp1 = P(i1,1)
            ppp2 = P(i1,2)
            ppp3 = P(i1,3)
            ppp4 = P(i1,4)
            rrr1 = V(i1,1)
            rrr2 = V(i1,2)
            rrr3 = V(i1,3)
            pppm = ppp1*ppp1 + ppp2*ppp2 + ppp3*ppp3
            if( pppm < 1D-28 ) pppm = 1D-28
            if( pppm > 1D28  ) cycle
            pppm = SQRT(pppm)
            rrrm = rrr1*rrr1 + rrr2*rrr2 + rrr3*rrr3
            if( rrrm < 1D-28 ) rrrm = 1D-28
            if( rrrm > 1D28 ) cycle
            rrrm = SQRT(rrrm)
            if( ( pppm <= dpmax .AND. ppp4 <= dpmax ) &
                .AND. rrrm <= drmax )then
                ishp(i1) = 1
            end if
        end do
!-----------------------------   Parton Checking   -----------------------------
!-------------------------------------------------------------------------------


!       Calculates the position for the center of mass of the
!        non-freeze-out system. The distance of a particle, when checking
!        is it freezing out or not, is measured with respect to this center
        call copl_p(time)

!       Step 1
!       Create the parton-parton (initial) collision time list.
        call ctlcre_par(iijk)
        ! Initial collis. list is empty.
        if(iijk == 2) return


!-------------------------------------------------------------------------------
!----------------------------   Parton Collision   -----------------------------
!       The loop over parton-parton collisions within an event.
        jjj = 0
        icolo = icol
!       Statistic of the number of loops in parton cascade within an event.
24      jjj = jjj + 1
        if( jjj > 100*icolo )then
            write(19,*) "Warning, infinite loop may have happened in" &
                    // "parcas iii, jjj, icolo=", iii, jjj, icolo
            iii = iii - 1
            iijk = 1
            return
        end if
        n0 = N

        if(jjj > 1) call copl_p(time)

!       Step 2
!       Finds the binary collision (icp) with the minimum colli. time.
        call find_par( icp, tcp )
!       The collision time list is empty if icp=0.
        if(icp == 0) goto 25
        ic = lc(1,icp)
        jc = lc(2,icp)
!       Parton rescattering is assumed no longer than 10000 fm/c.
        if( tcp > 10000D0 )then
            do i1 = icp+1, icol, 1
                lc( 1, i1-1 ) = lc(1,i1)
                lc( 2, i1-1 ) = lc(2,i1)
                tc( 1, i1-1 ) = tc(1,i1)
                tc( 2, i1-1 ) = tc(2,i1)
            end do
            icol = icol - 1
            jjj = jjj - 1
            goto 24
        end if
        time = tcp

!       Step 3
!       Perform classical Neuton motion for ic & jc partons
        do i=1,N,1
            do j=1,3,1
                rpo(i,j) = V(i,j)
            end do
            rpo(i,4) = V(i,4)
        end do
        call his_p( time, rnt, rnp, istop )
!       istop=1 means all particles have get out of considered volume.
        if(istop == 1) goto 25

!       Considers parton energy loss phenomenologically.
!       Do not use this option with the paronic rescattering openning
!        altogether. If so, double-counting!
        iadj136 = INT( adj1(36) )
        if( iadj136 == 1 )then
            bmax  = rnt + rnp
            bmax2 = bmax * bmax
            bp2   = bp * bp
!       Energy loss per unit distance.
            adj137  = adj1(37)
            dell = adj137 * ( nap / 197D0 )**0.75D0 * ( 1D0 - bp2 / bmax2 ) &
                 * ( win * win / 40000D0 )**0.25D0
            call eloss(dell)
        end if

!       Step 4
!       Performs parton-parton collision & updates particle list.
!       If lmn = 4, 6, & 7 updates 'PYJETS', 'sbe', diquark list, string
!        list, & lc(1-2,m) too.
        kkk  = 0
        iway = 0
        call collis( ic,jc,kf3,kf4, tcp, jjj,kkk, iway, icnew,jcnew, lmn, time )

!       Updates the information of the parton pair that actually collided.
        if( iway == 0 .AND. kkk == 0 )then
            i_FZ_update = 1
            call update_freeze_out_information( i_FZ_update, ic, ic )
            call update_freeze_out_information( i_FZ_update, jc, jc )
!       Updates the information of the radiated partons (shower) if any.
            if( N  >  n0 )then
                do i = n0+1, N, 1
                    N_FO = N_FO + 1
                    call update_freeze_out_information( i_FZ_update, i, N_FO )
                end do
            end if
        end if

!       Step 5
!       Updates the collision list.
!#TODO(Lei20241018): need to be extended for heavy quarks.
        if( lmn == 4 .OR. lmn == 6 .OR. lmn == 7 ) goto 26
!       In above case 'update collision list' has been executed in "collis'.
        call update_ctl( iway, tcp )

!       'new' method
!1      if( iadj136 == 0 ) call update_ctl( iway, tcp )
!1      if( iadj136 == 1 .AND. iway == 1 ) call update_ctl( iway, tcp )
!       If collision does not happen remove the collision pairs (one of
!        company is colliding particle) from collision time list only.
!1      if( iadj136 == 1 .AND. iway == 0)then
!       Creates the collision time list completely
!1      icol = 0
!1      do i=1,mclis,1
!1          lc(1,i) = 0
!1          lc(2,i) = 0
!1          tc(1,i) = 0D0
!1          tc(2,i) = 0D0
!1      end do
!1      call ctlcre_para( time )
!1      end if
26      continue
!       The loop over collisions within an event.
        goto 24
!----------------------------   Parton Collision   -----------------------------
!-------------------------------------------------------------------------------


25      continue
        time_par = time
!       time_par: is the time lasted in parton cascade hereafter
!       Cross sections statistics.
        do i=1,20,1
            reaci = reac(i)
            if( reaci > 0D0 ) crose(i) = crose(i) / reaci
        end do


!-------------------------------------------------------------------------------
!----------------------------   Gamma 55 Removing   ----------------------------
!       Removes gamma from "PYJETS" to "sgam".
!#TODO(Lei20240218): The id/KF "55" is not so safe for PYTHIA 8, which
!                    represents "DMmed(s=1)" (Dark Matter vector mediator).
!       The hadronization model.
        i_had_model = INT(adj1(12))
        if( i_had_model /= 0 )then
            n55 = 0
            do i1=1,N,1
                KF = K(i1,2)
                if(KF == 22)then
                    K(i1,2) = 55
                    n55 = n55 + 1
                end if
            end do
!           Moves "55" from "PYJETS" to "sgam".
            if(n55 > 0) call remo_gam(55)
!           egam2: energy of gammas after parton cascade.
            call prt_sgam( n55, egam2, 2 )
        end if
!----------------------------   Gamma 55 Removing   ----------------------------
!-------------------------------------------------------------------------------


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine reset_eve
!!      Initializes the collision time list.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MCLIS=280000)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/work7/reac(20),crose(20)
!       reac and crose: the arrays to account for the number and
!        the value of cross section for 2->2 partonic processes


        ic = 0
        jc = 0
        icol  = 0
        reac  = 0D0
        crose = 0D0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_par( iijk )
!!      Creates the initial collision list.
!       Every process (i.e. parton casecade) starts from time equal to 0.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        LOGICAL IS_PARTON
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
        icol = 1
!       Every process (i.e. parton casecade) starts from time equal to 0.
        ijk = 0
!       Considers d, u, s, c, b, their antiquarks, and gluon only.
        do i = 1, N-1, 1
            kfi = K(i,2)
            if( .NOT.IS_PARTON( kfi ) ) cycle
            do j = i+1, N, 1
                kfj = K(j,2)
                if( .NOT.IS_PARTON( kfj ) ) cycle
                iflag = 0
                call rsfilt_p( i, j, iflag )
                if( iflag == 0 ) cycle
                tc( 1, icol ) = 0D0
                tc( 2, icol ) = 0D0
                kji = 0
                call coij_p( i, j, icol, kji )
                if( kji == 1 )then
                    ijk = ijk + 1
                    cycle
                end if
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if( tc1 > tc2 ) tcicol = tc2
                if( tcicol > 0D0 )then
                    tci = tcicol
                    do j1 = 1, icol-1, 1
                        tcj = tc(1,j1)
                        tck = tc(2,j1)
                        if( tcj > tck ) tcj = tck
                        if( ABS( tcj - tci ) < dddt )then
                            icol = icol - 1
                            exit
                        end if
                    end do
                    icol = icol + 1
                end if
            end do
        end do
        icol = icol - 1
        if( icol == 0 ) iijk = 2


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine find_par(icp,tcp)
!!      Finds out the binary collision (icp) with minimum colli. time (tcp).
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MCLIS=280000)
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/collist/lc(2,mclis),tc(2,mclis),icol


        icp = 0
        tcp = 10000D0
        do i=1,icol,1
            ! Plays role after first colli.
            ia = ( ABS( lc(1,i) - ic ) + ABS( lc(2,i) - jc ) ) &
               * ( ABS( lc(1,i) - jc ) + ABS( lc(2,i) - ic ) )
            if( ia == 0 ) cycle
            tc1 = tc(1,i)
            tc2 = tc(2,i)
            tci = tc1
            if( tc1 > tc2 ) tci = tc2
            ! Alternative choice.
            ! tci = MAX( tc(1,i), tc(2,i) )
            if( tci <= 0D0 ) cycle
            if( tcp < tci )  cycle
            icp = i
            tcp = tci
        enddo


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctl(iway,time)
!!      Updates collision time list for both of w/o & w/ inelastic
!!       parton-parton scattering.
!       if iway=1 the collision does not happen
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
        j = 0
!       Loops over old colliding pairs.
        do i=1,icol,1
            i1  = lc(1,i)
            j1  = lc(2,i)
            iiii = ( i1 - ic ) * ( j1 - ic )
            jjjj = ( i1 - jc ) * ( j1 - jc )
            ! Avoids the particle collideing with itself.
            if(i1 == j1) cycle
!       Throws away the pairs composed of ic and/or jc
            if( iiii == 0 .OR. jjjj == 0 ) cycle
            tc1 = tc(1,i)
            tc2 = tc(2,i)
            tci = tc1
            if( tc1 > tc2 ) tci = tc2
!           if( tci <= time ) cycle
!           if( ddt == 0D0 ) dddt = ddt + 0.03D0
!           if( ddt /= 0D0 ) dddt = ddt * 300D0
!       Throws away the pairs with tc-time <= dddt (time accuracy).
            if( (tci - time) <= dddt ) cycle
!       Proceeds for survivor.
            j = j + 1
            tc(1,j) = tc(1,i)
            tc(2,j) = tc(2,i)
            lc(1,j) = lc(1,i)
            lc(2,j) = lc(2,i)
        end do
        icol = j
        if( iway == 1 ) return

!       Loops over ic,jc (new) and old partons (i.e. construct colli. pair
!        by partons, one of which is ic or jc and another one is in parton list)
        icol = j + 1
        do i=1,n0,1
            i1 = i
            if(i1 == ic) cycle
            if(i1 == jc) cycle
            kfi = ABS(K(i,2))
!           if( kfi > 3 .and. kfi /= 21 ) cycle
!       consider d,u,s, their antiquarks, and gluon only
            if( kfi > 6 .and. kfi /= 21 ) cycle
!       consider d,u,s,c,b, their antiquarks, and gluon only
            do k1=1,2,1
                if(k1 == 1) j1 = ic
                if(k1 == 2) j1 = jc
                call rsfilt_p(j1,i1,iflag)
                if(iflag == 0) cycle
                tc(1,icol) = 0D0
                tc(2,icol) = 0D0
                call tcolij_par(i1,j1,time,icol)
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if(tc1 > tc2) tcicol = tc2
                if( tcicol > 0D0 )then
                    tci = tcicol
                    do jj = 1, icol-1, 1
                        tcj = tc(1,jj)
                        tck = tc(2,jj)
                        if(tcj > tck) tcj = tck
                        if( ABS(tcj - tci) < dddt )then
                            icol = icol - 1
                            exit
                        end if
                    end do
                    icol = icol + 1
                end if
            end do
        end do
        icol = icol - 1
        n0 = N


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctlm(time)
!!      A part of updating collision time list (throw away old collission
!!       pairs only) for inela. parton-parton scattering 7.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MCLIS=280000)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
        j = 0
!       loop over old colliding pairs
        do i=1,icol,1
            i1 = lc(1,i)
            j1 = lc(2,i)
            iiii = ( i1 - ic ) * ( j1 - ic )
            jjjj = ( i1 - jc ) * ( j1 - jc )
            ! Avoid the particle collideing with itself.
            if( i1 == j1 ) cycle
!       throw away the pairs composed of ic and/or jc
            if( iiii == 0 .OR. jjjj == 0 ) cycle
            tc1=tc(1,i)
            tc2=tc(2,i)
            tci=tc1
            if(tc1 > tc2)  tci=tc2
!           if( tci <= time ) cycle
!           if( ddt == 0D0 ) dddt = ddt + 0.03D0
!           if( ddt /= 0D0 ) dddt = ddt * 300D0
!       throw away the pairs with tc-time <= dddt (time accuracy)
            if( (tci-time) <= dddt ) cycle
            j = j + 1
            tc(1,j) = tc(1,i)
            tc(2,j) = tc(2,i)
            lc(1,j) = lc(1,i)
            lc(2,j) = lc(2,i)
        end do
        icol = j


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_ctln(time)
!!      Updates the collision time list (a part) for inela. parton-parton
!!       scattering 7.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
!       loop over ic (jc) and old 'PYJETS' (i.e. construct colli. pair
!        composed of partons one of which is ic (jc) and another one
!        in old 'PYJETS')
        icol = icol + 1
        do i=1,n0,1
            i1 = i
            if(i1 == ic) cycle
            if(i1 == jc) cycle
            kfi = ABS(k(i,2))
!          if(kfi > 3 .and. kfi /= 21) cycle
!       consider d,u,s, their antiquarks, and gluon only
            if( kfi > 6 .and. kfi /= 21 ) cycle
!       consider d,u,s,c,b, their antiquarks, and gluon only
            do k1=1,2,1
                if(k1 == 1) j1 = ic
                if(k1 == 2) j1 = jc
                call rsfilt_p(j1,i1,iflag)
                if(iflag == 0) cycle
                tc(1,icol) = 0D0
                tc(2,icol) = 0D0
                call tcolij_par(i1,j1,time,icol)
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if( tc1 > tc2 ) tcicol = tc2
                if( tcicol > 0D0 )then
                    tci = tcicol
                    do jj = 1, icol-1, 1
                        tcj = tc(1,jj)
                        tck = tc(2,jj)
                        if(tcj > tck) tcj = tck
                        if( ABS(tcj - tci) < dddt )then
                            icol = icol - 1
                            exit
                        end if
                    end do
                    icol = icol + 1
                end if
            end do
        end do
        icol = icol - 1
        n0 = N


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coij_p( i, j, icp, kji )
!!      Calculates the collision time for construction
!!       of initial collision time list.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/syspar_p/rsig1,pio,tcut
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0
        dimension px(4),py(4),pij(4)
        dimension dr(3),db(3),vi(3),vj(3),pic(4),pjc(4),pxc(4),pyc(4)
        double precision b(3)


        kji = 0
        pio = 3.141592653589793D0
        pi(4) = P(i,4)
        pj(4) = P(j,4)
        if( pi(4) < 1D-10 ) pi(4) = 1D0-10
        if( pj(4) < 1D-10 ) pj(4) = 1D0-10
        pij(4) = pi(4) + pj(4)
        pic(4) = pi(4)
        pjc(4) = pj(4)
        do k1=1,3,1
            pi(k1)  = P(i,k1)
            pj(k1)  = P(j,k1)
            pij(k1) = pi(k1) + pj(k1)
            pic(k1) = pi(k1)
            pjc(k1) = pj(k1)
            b(k1)   = ( pi(k1) + pj(k1) ) / pij(4)
        end do
        rmi = p(i,5)
        rmj = p(j,5)
!       Invariant mass.
        eiej2 = dot(pij,pij)
!       Inserts the energy cut
        if( eiej2 < 0D0 ) then
            kji = 1
            return
        end if
        do n1=1,4,1
            px(n1)  = V(i,n1)
            py(n1)  = V(j,n1)
            pxc(n1) = px(n1)
            pyc(n1) = py(n1)
        enddo
        ilo  = 0
        kf1  = K(i,2)
        kf2  = K(j,2)
        ikf1 = ABS(kf1)
        ikf2 = ABS(kf2)
        call fsig( ilo, kf1, kf2, kf3, kf4, eiej2, sig, tsmp, lmn, jjj )
        if( sig <= 0D0 ) return
        if( ilo == -2 ) return
        rsig1 = SQRT( sig / pio )
        ! Momentum.
        call lorntz( 0, b, pic, pjc )
        ! Position.
        call lorntz( 0, b, pxc, pyc )
        rb = 0D0
        bb = 0D0
        rr = 0D0
        rtai = 0D0
        do k1=1,3,1
            vi(k1) = pic(k1) / pic(4)
            vj(k1) = pjc(k1) / pjc(4)
        end do
        do k1=1,3,1
            dr(k1) = pxc(k1) - pyc(k1) - ( vi(k1)*pxc(4) - vj(k1)*pyc(4) )
            db(k1) = vi(k1) - vj(k1)
            rb = rb + dr(k1)*db(k1)
            bb = db(k1)**2 + bb
            rr = rr + dr(k1)*dr(k1)
        end do
        if( bb <= 1D-10 ) return
        tcol = 0D0 - rb / bb
        do ik=1,3,1
            dr(ik) = dr(ik) + tcol * db(ik)
            rtai = rtai + dr(ik) * dr(ik)
        end do
        sg = rtai
        dmin = SQRT(sg)
        if( dmin > rsig1 ) return
        do ik=1,3,1
            pxc(ik) = pxc(ik) + vi(ik) * ( tcol - pxc(4) )
            pyc(ik) = pyc(ik) + vj(ik) * ( tcol - pyc(4) )
        end do
!       Moves along Newton trajectory in CMS,
        pxc(4) = tcol
        pyc(4) = tcol
!       Transform back to Lab ((causality violation).
        call lorntz(1,b,pxc,pyc)
!       Max time 10000 fm/c.
        if( pxc(4) <= 10000D0 .AND. pyc(4) <= 10000D0 )then
            lc(1,icp) = i
            lc(2,icp) = j
            tc(1,icp) = pxc(4)
            tc(2,icp) = pyc(4)
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine tcolij_par(i,j,time,icp)
!!      Calculates the collision time of i and j.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/syspar_p/rsig1,pio,tcut
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/scatt/pi(4),pj(4),ic,jc,n0
        dimension px(4),py(4),pij(4),pic(4),pjc(4),pxc(4),pyc(4)
        dimension dr(3),db(3),vi(3),vj(3),rfi(4),rfj(4)
        double precision b(3)


        pi(4) = P(i,4)
        pj(4) = P(j,4)
        if( pi(4) < 1D-10 ) pi(4) = 1D-10
        if( pj(4) < 1D-10 ) pj(4) = 1D-10
        pij(4) = pi(4) + pj(4)
        pic(4) = pi(4)
        pjc(4) = pj(4)
        do k1=1,3,1
            pi(k1)  = P(i,k1)
            pj(k1)  = P(j,k1)
            pij(k1) = pi(k1) + pj(k1)
            pic(k1) = pi(k1)
            pjc(k1) = pj(k1)
            b(k1) = ( pi(k1) + pj(k1) ) / ( pi(4) + pj(4) )
        end do
        rmi = P(i,5)
        rmj = P(j,5)
        eiej2 = dot(pij,pij)
!       squared invariant mass
!       insert! energy cut
        if( eiej2 < 0D0 ) return
!       ecut  = SQRT(eiej2) - rmi - rmj
!       ecut0 = 0.02D0    ! P.Yang
!       if(ecut <= ecut0) return
! Note! energy cut can be taken into account HERE to decide coll. pairs
! etc.  According to Y. Pang's opinions, May,1994, CCAST.(05/24/95)
        do n1=1,4,1
            px(n1)  = V(i,n1)
            py(n1)  = V(j,n1)
            pxc(n1) = px(n1)
            pyc(n1) = py(n1)
        end do
        ilo  = 0
        kf1  = K(i,2)
        kf2  = K(j,2)
        ikf1 = ABS(kf1)
        ikf2 = ABS(kf2)
        if( (ikf1 <= 6 .OR. ikf1==21 ) .AND. ( ikf2 <= 6 .OR. ikf2 == 21 ) )then
!       d, u, s, c, b quarks, their anti quarks, and gluon only
            call fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)
        else
            sig = 0D0
        end if
!       if(il == 2) return
        if(ilo == -2) return
        if(sig <= 0D0) return
        rsig1 = SQRT( sig / pio )

        do i1=1,3
            rfi(i1) = px(i1) + ( taup(i) - time ) * pi(i1) / pi(4)
            rfj(i1) = py(i1) + ( taup(j) - time ) * pj(i1) / pj(4)
        enddo
        rfi(4) = taup(i)
        rfj(4) = taup(j)
!       spatial coordinates of colliding particles at formation
!        moment in Lab. frame
        call lorntz(0,b,rfi,rfj)
        ctaui = rfi(4)
        ctauj = rfj(4)
        tcol = ctaui
        if(ctaui < ctauj) tcol = ctauj
!       for back to back collision the collision time is equal to
!        the lager one of their formation times
        call lorntz(0,b,pic,pjc)
        call lorntz(0,b,pxc,pyc)
        rb = 0D0
        bb = 0D0
        rr = 0D0
        rtai  = 0D0
        kflag = 0
        do ik=1,3,1
            vi(ik) = pic(ik) / pic(4)
            vj(ik) = pjc(ik) / pjc(4)
        end do
        do ik=1,3,1
            db(ik) = vi(ik)  - vj(ik)
            dr(ik) = pxc(ik) - pyc(ik) &
                   - ( vi(ik)*pxc(4) - vj(ik)*pyc(4) ) + db(ik)*tcol
            rtai   = rtai + dr(ik)*dr(ik)
        end do
        dott = 0D0
        do ik=1,3,1
            dott = dr(ik)*pic(ik) + dott
        end do

        if( dott >= 0D0 )then
            kflag = 1
            if(tcol <= pxc(4)) return
            if(tcol <= pyc(4)) return
!       for the back to back collisions (collision happens in future)
        else
            rtai = 0D0
            do ik=1,3,1
                dr(ik) = pxc(ik) - pyc(ik) - ( vi(ik)*pxc(4) - vj(ik)*pyc(4) )
                rb = rb + dr(ik) * db(ik)
                bb = bb + db(ik) * db(ik)
                rr = rr + dr(ik) * dr(ik)
            end do
            if(bb  <=  1D-10) return
            tcol = 0D0 - rb/bb
            if(tcol <= pxc(4)) return
            if(tcol <= pyc(4)) return
!       collision happens in future
            if( tcol-ctaui <= 0D0 ) return
            if( tcol-ctauj <= 0D0 ) return
!       collision must be after formation
            do ik=1,3,1
                dr(ik) = dr(ik) + db(ik) * tcol
                rtai = rtai + dr(ik) * dr(ik)
            enddo
!       for the face to face collisions
        end if

        sg = rtai
        dmin = SQRT(sg)
        if(dmin > rsig1) return
        do ik=1,3,1
            pxc(ik) = pxc(ik) + vi(ik) * ( tcol - pxc(4) )
            pyc(ik) = pyc(ik) + vj(ik) * ( tcol - pyc(4) )
        enddo
!       move along Newton trajectory in CMS
        pxc(4) = tcol
        pyc(4) = tcol
        call lorntz(1,b,pxc,pyc)
!       transform back to Lab.
        tcol1 = pxc(4)
        tcol2 = pyc(4)
        if(kflag == 0)then
            if( tcol1 - taup(i) < 0D0 ) return
            if( tcol2 - taup(j) < 0D0 ) return
        else
            if( tcol1-taup(i) < -1D-4 ) return
            if( tcol2-taup(j) < -1D-4 ) return
        end if
        if( tcol1 <= time ) return
        if( tcol2 <= time ) return
!       collision happens in the future
        if( tcol1 <= 10000D0 .AND. tcol2 <= 10000D0 )then
            tc(1,icp) = tcol1
            tc(2,icp) = tcol2
            lc(1,icp) = i
            lc(2,icp) = j
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function dot(v1,v2)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension v1(4),v2(4)


        dot = v1(4)*v2(4) - v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function ichge(kf)
!       Calculates the charge (in unit of 1/3) of parton.
!       kf: parton flaver
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        ichge = 0
        if(kf == 21 .or. kf == -21)then
            ichge=0   ! gluon charge
        else if(kf > 0)then
            if(kf == 1) ichge = -1   ! d
            if(kf == 2) ichge =  2   ! u
            if(kf == 3) ichge = -1   ! s
            if(kf == 4) ichge =  2   ! c
            if(kf == 5) ichge = -1   ! b
            if(kf == 6) ichge =  2   ! t
        else if(kf < 0)then
            if(kf == -1) ichge =  1   ! dbar
            if(kf == -2) ichge = -2   ! ubar
            if(kf == -3) ichge =  1   ! sbar
            if(kf == -4) ichge = -2   ! cbar
            if(kf == -5) ichge =  1   ! bbar
            if(kf == -6) ichge = -2   ! tbar
        else
            write(22,*) "Abort! Non-parton in ichge of parcas.f90."
            stop
        endif


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)
!!      Calculates the total cross section and decide the type of reaction
!!       (when ilo=0) or sample the t value as well (when ilo=1).
!       tsmp: t value sampled.
!       lmn: order number happened process
!       sig: the total cross section of parton kf1 collides with kf2
!       kf3 and kf4: kf code of the colliding pair after collision
!       eiej2: squared invariant mass of colliding pair
!       jjj: jjj-th loop within a event
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
!       leadding order pQCD differential cross section of 2->2 processes
        external fs11_0,fs11_1,fs11_2,fs12_0,fs12_1,fsqq,fsgg_0,fsgg_1,fsqg
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/syspar_p/rsig1,pio,tcut
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        dimension ssig(0:3), eee(0:INT(adj1(4))), dd(0:INT(adj1(4))), &
                             eee1(0:INT(adj1(4))), eee2(0:INT(adj1(4))), &
                             eee3(0:INT(adj1(4)))
        logical lgg,lqg,lqq,lqaq,lqqp



!===============================================================================
!       New method
        call fsig_new(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)
        return
!       The following old method is no longer used.
!===============================================================================



!       classify the incident channel
        kf0  = ABS(kf1)
        ikf2 = ABS(kf2)
!       lqq: two incoming partons are quarks (or antiquarks) with same flavor
        lqq  = (kf1 == kf2) .AND. (kf0 /= 21)
!       lgg: two gluons
        lgg = (kf1 == 21) .AND. (kf2 == 21)
!       lqg: one quark (or antiquark) and one gluon
        lqg = ( (kf0 == 21)  .AND. (kf2 /= 21)) .OR. &
              ( (ikf2 == 21) .AND. (kf1 /= 21))
!       lqaq: quark (antiquark) and its corresponding antiquark (quark)
        lqaq = (kf1 == -kf2) .AND. (kf0 /= 21)
!       lqqp: two quarks,or two antiquarks,or one quark and one antiquark
!        with different flavor.
        lqqp = (kf0 /= ikf2) .AND. (kf0 /= 21) .AND. (ikf2 /= 21)

!       calculate the total cross section for each incedent channel
        sig = 0D0
        lmn = 0
        do i=0,3,1
            ssig(i) = 0D0
        end do
        idw = INT( adj1(4) )
        iadj120 = INT( adj1(20) )
        ! mass of qqbar (kinematical mass)
        dmass = amass(1)*2D0
        umass = amass(2)*2D0
        smass = amass(3)*2D0
        cmass = amass(4)*2D0
        bmass = amass(5)*2D0
        tmass = amass(6)*2D0
        xs = eiej2
!       xs: squared invariant mass of colliding pair
        if(eiej2 <= 0.) eiej2 = 1D-16
!       invariant mass of colliding pair
        eiej = SQRT(eiej2)
        sm = 0D0
!       the sum of squared mass of partons in initial and final states, it
!        is set to zero (zero mass appoxi. in kinematics treatment)
        nf = 3
!       alpha_s
!       as = 12D0 * pio / ( 33D0 - 2D0*nf )
!       as = 12D0 * pio / ( 33D0 - 2D0*nf ) / xs
!       as = 0.3D0
        as = adj1(2)
!       tmax
!       uxt = 0D0
!       tmin
!       dxt = -xs
!       cf. my note book p.15
!       above uxt and dxt are used in limited and regularized cross section
!
!       tcut = 0.4D0   ! 5.86
!sa's   tcut: the cut for avoidding divergence
!song's tcut: the cut value of t, a process is considered to be 'hard'
!song's  when t > tcut, 'soft' when t < tcut.
        tcut = adj1(3)
        uxt  = -tcut
!       uxt: the upper limit (max.) of t
!       note,the integral variable t is negative
!       dxt = ( sm - xs ) + tcut
        dxt = - xs - tcut
!       dxt = - xs + tcut
!       dxt: the lower limit (min.) of t
!       above uxt and dxt are used in LO pQCD cross section
        if(uxt <= dxt) then
!       the process can not happen
            ilo=-2
            return
        end if

        eee = 0D0
        if( lqaq .OR. lgg ) then
            ! gg->
            if( kf0 == 21 )then
                ! ->gg process 9
                call integ(fsgg_0,uxt,dxt,idw,eee,sum)
                ! fsgg_0: differential cross section gg->gg
                sum1 = sum
                do i1=0,idw,1
                    eee1(i1) = eee(i1)
                end do
                ! ->qqbar process 7 (inelastic)
                iparreso = iparres
                if( ( iparres == 1 .OR. iparres == 3 ) .AND. eiej < umass )then
                    iparres = 0
                    ! treates as elastic (iparres=0), process 9 happens
                    ! process 7 doesn't happen
                    call integ(fsgg_1,uxt,dxt,idw,eee,sum)
                    iparres = iparreso
                else
                    call integ(fsgg_1,uxt,dxt,idw,eee,sum)
                end if
                sum2 = sum
                do i1=0,idw,1
                    eee2(i1) = eee(i1)
                end do
                sums = sum1 + sum2
                sum3 = sum1 / sums
                sele = PYR(1)
                ! gg->gg process 9
                if( sele <= sum3 )then
                    do i1=0,idw,1
                        eee(i1) = eee1(i1)
                    end do
                    kf3 = kf1
                    kf4 = kf2
                    lmn = 9
                ! gg->qqbar process 7 (inelastic)
                else
                    do i1=0,idw,1
                        eee(i1) = eee2(i1)
                    end do
                    ! ->qqbar branch
                    ! which is in coales_40.f90
                    call break_f(eiej,kf7,amq)
                    kf3 =  kf7
                    kf4 = -kf3
                    lmn = 7
                end if
                sum = sums
                sig = sums
            ! qqbar->
            else
                ! ->qqbar process 5
                call integ(fs11_2,uxt,dxt,idw,eee,sum)
                ! fs11_2: differential cross section qqbar->qqbar
                sum1 = sum
                do i1=0,idw,1
                    eee1(i1) = eee(i1)
                end do
                ! ->gg process 6 (inelastic)
                call integ(fsqq,uxt,dxt,idw,eee,sum)
                sum2 = sum
                do i1=0,idw,1
                    eee2(i1) = eee(i1)
                end do
                ! ->q'qbar' process 4 (inelastic)
                iparreso = iparres
                if( ( iparres == 1 .OR. iparres == 3 ) .AND. eiej < umass )then
                    iparres = 0
                    ! treates as elastic (iparres=0), process 5 or 6 happens
                    ! process 4 (inelastic) doesn't happen.
                    call integ(fs11_1,uxt,dxt,idw,eee,sum)
                    iparres = iparreso
                else
                    call integ(fs11_1,uxt,dxt,idw,eee,sum)
                end if
            sum3 = sum
            do i1=0,idw,1
                eee3(i1) = eee(i1)
            end do
            sums = sum1 + sum2 + sum3
            sum4 = sum1 / sums
            sum5 = ( sum1 + sum2 ) / sums
            sele = PYR(1)
            ! ->qqbar process 5
            if( sele <= sum4 )then
                do i1=0,idw,1
                    eee(i1) = eee1(i1)
                end do
                kf3 = kf1
                kf4 = kf2
                lmn = 5
            ! ->gg process 6 (inelastic)
            else if( sele > sum4 .and. sele <= sum5 )then
                do i1=0,idw,1
                    eee(i1) = eee2(i1)
                end do
                kf3 = 21
                kf4 = 21
                lmn = 6
            ! ->q'qbar', process 4 (inelastic)
            else
                do i1=0,idw,1
                    eee(i1) = eee3(i1)
                end do
                aa = PYR(1)
                ! ->qqbar branch
                do while(.true.)
                    ! which is in coales_40.f90
                    call break_f(eiej,kf7,amq)
                    if(kf7 /= kf0) exit
                end do
                kf3 = kf7
                kf4 = -kf3
                lmn = 4
                end if
                sum = sums
                sig = sums
            end if
        else if( lqq )then
            !  process 2
            call integ(fs11_0,uxt,dxt,idw,eee,sum)
            ! fs11_0 : differential x section of qq->qq
            sig = sum
            kf3 = kf1
            kf4 = kf2
            lmn = 2
        else if( lqqp )then
            ! process 1 and 3
            call integ(fs12_0,uxt,dxt,idw,eee,sum)
            ! process 1
            if( kf1*kf2 > 0 )then
                sig = sum
                kf3 = kf1
                kf4 = kf2
                lmn = 1
            end if
            ! process 3
            if( kf1*kf2 < 0 )then
                sig = sum
                kf3 = kf1
                kf4 = kf2
                lmn = 3
            end if
        else if( lqg )then
            ! process 8
            call integ(fsqg,uxt,dxt,idw,eee,sum)
            ! fsqg: differential x section of qg->qg
            sig = sum
            kf3 = kf1
            kf4 = kf2
            lmn = 8
        end if

        sig = sig * 0.04D0
!       0.04 is the transformation factor from (GeV)^-2 to (fm)^2.
        if( iadj120 == 0 .or. iadj120 == 1 )then
            do i=0,idw,1
                if( eee(i) < 0D0 )then
                    sig = 0D0
                    return
                end if
            end do
        end if

!       sample the t value
        if( ilo == 1 )then
            if( iadj120 == 2 .OR. iadj120 == 3 )then
                ! tsmp = -PYR(1) * xs
                seta = PYR(1) * pio
                tsmp = xs * ( COS(seta) - 1D0 ) / 2D0 + tcut
                return
            end if
            call eee_dd(idw,eee,dd)
            call samp_integ(idw,uxt,dxt,tt,dd)
            tsmp = tt + tcut
!       tsmp: the t value sampled.
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fsig_new(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)
!!      Calculates the total cross section, decides the type of reaction
!!       and/or samples the t value as well.
!
!       The masses of the heavy c and b are retained (massive), and
!        other masses of light quarks are set to zero (zero mass
!        approximation in kinematics treatment, massless).
!
!       ilo: flag for sampling t or not, or the special failure (-2).
!       kf1 and kf2: kf code of the incoming pair after collision
!       kf3 and kf4: kf code of the outgoing pair after collision
!       eiej2: squared invariant mass of colliding pair
!       sig: the total cross section of parton kf1 collideing with kf2
!       tsmp: t value sampled.
!       lmn: order number happened process.
!
!       q: light (anti-)quark; hQ/Q: heavy (anti-)quark; g: gluon.
!       We ignore the collision between two heavy Q.
!       Massless:
!             1: q1 + q2 -> q1 + q2
!             2: q1 + q1 -> q1 + q1
!             3: q1 + q2bar -> q1 + q2bar
!             4: q1 + q1bar -> q2 + q2bar   (inelastic)
!             5: q1 + q1bar -> q1 + q1bar
!             6: q1 + q1bar -> g + g        (inelastic)
!             7: g  + g -> q1 + q1bar       (inelastic)
!             8: q  + g -> q + g
!             9: g  + g -> g + g
!       Massive:
!            10: hQ + q -> hQ + q
!            11: hQ + g -> hQ + g
!       Not yet used:
!            12: q + qbar -> hQ + hQbar   (inelastic)
!            13: g + g    -> hQ + hQbar   (inelastic)
!       Not yet included:
!            14: hQ  + hQbar  -> q + qbar (inelastic)
!            15: hQ  + hQbar  -> g + g    (inelastic)
!       Not taken into accout:
!            16: hQ  + hQbar  -> hQ + hQbar     (? elastic）
!            17: hQ1 + hQ1bar -> hQ2 + hQ2bar   (? inelastic)
!
!       jjj: jjj-th collision in the collisions list.
!
!       References:
!        Massless:
!            B.L. Combridge et al., Phys. Lett. B 70 (1977) 234;
!            Bin Zhang, Comput.Phys.Commun. 109 (1998) 193-206;
!            Jussi Auvinen et al., Phys.Rev.C 82 (2010) 024906;
!            Ben-Hao Sa, Comput.Phys.Commun. 183 (2012) 333–346.
!       Massive:
!            B.L. Combridge Nucl. Phys. B 151 (1979) 429-456;
!            Hamza Berrehrah et al., Phys.Rev.C 89 (2014) 5, 054901.
!       t sampling:
!            Klaus Geiger et al., Nucl.Phys.B 369 (1992) 600-654;
!            Ben-Hao Sa et al., Comput.Phys.Commun. 183 (2012) 333–346.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER( pi=3.141592653589793D0 )
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa33/smadel,ecce,secce,parecc,iparres
!       Note the following logical type.
        logical l_q1q1bar, l_gg, l_hQhQ
        logical l_q1q2_to_q1q2,       l_q1q1_to_q1q1, &
                l_q1q2bar_to_q1q2bar, l_qg_to_qg, &
                l_hQq_to_hQq,         l_hQg_to_hQg
        real(kind=8) :: sigma_relative(0:6), K_factor
        real(kind=8) :: Lambda_QCD, Lambda_QCD2
        integer :: KF_outgoing( 2, 0:6 )


!-------------------------------------------------------------------------------
!---------------------------   Process Classifying   ---------------------------
!       Classifies the incident partons.
!       kf: 1 ~ 5 -- d, u, s, c, b; 21 -- g.
!       kh: 0, 1, 20 -- q, hQ, g
        ikf1 = ABS(kf1)
        kh1 = 0
        if( kf1 == 21 ) kh1 = 20
        if( ikf1 == 4 .OR. ikf1 == 5 ) kh1 = 1
        ikf2 = ABS(kf2)
        kh2 = 0
        if( kf2 == 21 ) kh2 = 20
        if( ikf2 == 4 .OR. ikf2 == 5 ) kh2 = 1
!       i_heavy: 0, light q + light q; 20, light q + g; 40, g + g;
!                1, heavy Q + light q; 21, heavy Q + g;
!                2, heavy Q + heavy Q.
!          i.e. > 0, with heavy Q;
!               = 0/20/40 ( MOD(*,1)=0 ), without heavy Q.
        i_heavy = kh1 + kh2

!       Classifies the process.
!       Process 1: q1 + q2 -> q1 + q2.
        l_q1q2_to_q1q2 = ( kf1 /= kf2 ) .AND. ( kf1*kf2 > 0 ) &
                         .AND. ( i_heavy == 0 )
!       Process 2: q1 + q1 -> q1 + q1.
        l_q1q1_to_q1q1 = ( kf1 == kf2 ) .AND. ( i_heavy == 0 )
!       Process 3: q1 + q2bar -> q1 + q2bar.
        l_q1q2bar_to_q1q2bar = ( ikf1 /= ikf2 ) .AND. ( kf1*kf2 < 0 ) &
                              .AND. ( i_heavy == 0 )
!       Process 4, 5, 6 and 12: q1 + q1bar.
        l_q1q1bar = ( kf1 == -kf2 ) .AND. ( i_heavy == 0 )
!       Process 7, 9 and 13: g + g.
        l_gg = i_heavy == 40
!       Process 8: q + g -> q + g.
        l_qg_to_qg = i_heavy == 20
!       Process 10: hQ + q -> hQ + q.
        l_hQq_to_hQq = i_heavy == 1
!       Process 11: hQ + g -> hQ + g.
        l_hQg_to_hQg = i_heavy == 21
!       We will ignore hQ + hQ processes.
        l_hQhQ = i_heavy == 2
!---------------------------   Process Classifying   ---------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------   Preparation   -------------------------------
!       Initialization of the total cross section calculated, t sampled
!        and order number of process to be returned.
        sig  = 0D0
        tsmp = 0D0
        lmn  = 0
!       jjj-th collision.
        i_coll = jjj
!       Relative cross sections among different outgoing partons in the
!        inelastic processed.
        sigma_relative = 0D0
!       Switch of inelastic processes.
        i_inelastic = iparres
!       Selection of inelastic processes when i_inelastic = 1/3.
!       i_inel_proc = 6, all of 4, 6, 7 (12, 13); = 7, only 7 (13).
        i_inelastic_process = i_inel_proc
!       Selection of the form of differential cross section in LO-pQCD
!        and the scattering angle distribution to be used.
!       i_LOpQCD:
!        =0, full LO-pQCD parton-parton cross section.
!        =1, keeping only leading divergent terms in the LO-pQCD.
!        =2, same as 0 but flat scattering angle distribution.
!        =3, same as 1 but flat scattering angle distribution.
        i_LOpQCD = INT( adj1(20) )
!       Phenomenological high-order correction factor K.
        K_factor = adj1(1)
!       Squared invariant mass of the colliding pair.
        s = eiej2
        if( s <= 0D0 ) return
!       Invariant mass of the colliding pair.
        sqrt_s = SQRT( s )
!       Effective coupling constant.
!       For massless processes, use the constant adj1(2).
!       For massive processes, it will be calculated automatically.
        Lambda_QCD  = adj1(25)
        Lambda_QCD2 = Lambda_QCD**2
        alpha_s     = adj1(2)
!       t_cut: the cutoff used to regulate 't' divergence, AKA mu^2.
        t_cut = adj1(3)
!       t_max: the upper limit (max.) of t.
        t_max = -t_cut
!       t_min: the lower limit (min.) of t.
        t_min = -s - t_cut
!       For full LO-pQCD, regulates 'u' divergence in some cases, too.
        if( i_LOpQCD == 0 .OR. i_LOpQCD == 2 ) t_min = -s + t_cut
!       Above t_max and t_min will be used in LO-pQCD cross section.
!       Note that the integral variable t is negative.

!       Wrong integration limits.
        if( t_max <= t_min )then
            ilo = -2
            return
        end if
!       Takes the "threshold energy" into account for massless processes
!        due to the introduction of the cutoff. For the cuoffs of the
!        extra "u" divergence, it becomes 2*t_cut (see below).
        if( MOD(i_heavy,10) == 0 .AND. s < t_cut )then
            ilo = -2
            return
        end if

!       Mass of heavy quark for massive processes. It will be given
!        automatically.
        dM = 0D0
!-------------------------------   Preparation   -------------------------------
!-------------------------------------------------------------------------------


!       Calculates the total cross section for each colliding pair.


!-------------------------------------------------------------------------------
!------------------------   Cross Section Calculating   ------------------------
!       Process 1: q1 + q2 -> q1 + q2.
        if( l_q1q2_to_q1q2 )then
!       't'-divergence in both full and simplified LO-pQCD.
            t_min = -s - t_cut
            sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                     K_factor, 1, i_LOpQCD )
            lmn = 1
            sig = sig_out
            kf3 = kf1
            kf4 = kf2

!-------------------------------------------------------------------------------
!       Process 2: q1 + q1 -> q1 + q1.
        else if( l_q1q1_to_q1q1 )then
!       't'-divergence in both full and simplified LO-pQCD.
!       Extra 'u'-divergence in full one. "threshold energy" becomes 2*t_cut.
            if( (i_LOpQCD == 0 .OR. i_LOpQCD == 2) &
                .AND. s < 2D0*t_cut )then
                ilo = -2
                return
            end if
            sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                     K_factor, 2, i_LOpQCD )
            lmn = 2
            sig = sig_out
            kf3 = kf1
            kf4 = kf2

!-------------------------------------------------------------------------------
!       Process 3: q1 + q2bar -> q1 + q2bar.
        else if( l_q1q2bar_to_q1q2bar )then
!       "t"-divergence in both full and simplified LO-pQCD.
            t_min = -s - t_cut
            sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                     K_factor, 3, i_LOpQCD )
            lmn = 3
            sig = sig_out
            kf3 = kf1
            kf4 = kf2

!-------------------------------------------------------------------------------
!       Process 4, 5, 6 and 12: q1 + q1bar ->
        else if( l_q1q1bar )then
!       Elastic process 5, q1 + q1bar -> q1 + q1bar.
!       't'-divergence in both full and simplified LO-pQCD.
            t_min = -s - t_cut
            sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                     K_factor, 5, i_LOpQCD )
            sigma_relative(0)   = sig_out
            KF_outgoing( 1, 0 ) = kf1
            KF_outgoing( 2, 0 ) = kf2
!       Inelastic process 4, 12 and 6, q1 + q1bar -> q2+q2bar / hQ+hQbar / g+g.
!       Takes the threshold energy (2*mass) into account.
            if( ( i_inelastic == 1 .OR. i_inelastic == 3 ) &
                .AND. i_inelastic_process == 6 )then
!       Inelastic process 4, q1 + q1bar -> u + ubar / d + dbar / s + sbar.
!       With the massless approximation, their cross sections are equal.
!       No divergence (?).
                t_min = - s
                t_max = 0D0
!       Inelastic Process 4, q1 + q1bar -> u + ubar.
                if( sqrt_s  >=  2D0*amass(2) )then
                    sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                    alpha_s, K_factor, 4, i_LOpQCD )
                    sigma_relative(2)   = sig_out
                    KF_outgoing( 1, 2 ) =  2
                    KF_outgoing( 2, 2 ) = -2
!       Inelastic Process 4, q1 + q1bar -> d + dbar.
                    if( sqrt_s  >=  2D0*amass(1) )then
                        sigma_relative(1)   = sig_out
                        KF_outgoing( 1, 1 ) =  1
                        KF_outgoing( 2, 1 ) = -1
                    end if
!       Inelastic Process 4, q1 + q1bar -> s + sbar.
                    if( sqrt_s  >=  2D0*amass(3) )then
                        sigma_relative(3)   = sig_out
                        KF_outgoing( 1, 3 ) =  3
                        KF_outgoing( 2, 3 ) = -3
                    end if
                end if
!===============================================================================
!       DO NOT USE THE MASSIVE PROCESSES OF q1 + q1bar -> hQ + hQbar NOW.
!       BECAUSE THEIR REVERSE PROCESSES HAVE NOT YET BEEN INTRODUCED.
                IF( .FALSE. )THEN
!       Inelastic process 12, q1 + q1bar -> c + cbar.
!       Modifies the integral limits and alpha_s due to the heavy mass.
                if( sqrt_s  >=  2D0*amass(4) )then
                    dM  = amass(4)
                    dM2 = dM**2
                    t_min = dM2 - s*( 1D0 + SQRT(1D0 - 4D0*dM2/s) ) /2D0
                    t_max = dM2 - s*( 1D0 - SQRT(1D0 - 4D0*dM2/s) ) /2D0
                    Q2 = s
                    nf = 4
                    alpha_s = 12D0 * pi / ( 33D0 - 2D0*nf ) &
                              / LOG( Q2 / Lambda_QCD2 )
                    sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                    alpha_s, K_factor, 12, i_LOpQCD )
                    sigma_relative(4)   = sig_out
                    KF_outgoing( 1, 4 ) =  4
                    KF_outgoing( 2, 4 ) = -4
!       Inelastic process 12, q1 + q1bar -> b + bbar.
!       Modifies the integral limits and alpha_s due to the heavy mass.
                    if( sqrt_s  >=  2D0*amass(5) )then
                        dM  = amass(5)
                        dM2 = dM**2
                        t_min = dM2 -s*( 1D0 +SQRT(1D0 -4D0*dM2/s) )/2D0
                        t_max = dM2 -s*( 1D0 -SQRT(1D0 -4D0*dM2/s) )/2D0
                        Q2 = s
                        nf = 5
                        alpha_s = 12D0 * pi / ( 33D0 - 2D0*nf ) &
                                  / LOG( Q2 / Lambda_QCD2 )
                        sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                        alpha_s, K_factor, 12, i_LOpQCD)
                        sigma_relative(5)   = sig_out
                        KF_outgoing( 1, 5 ) =  5
                        KF_outgoing( 2, 5 ) = -5
                    end if
                end if
                END IF
!===============================================================================
!       Inelastic process 6, q1 + q1bar -> g + g.
!       't'-divergence in both full and simplified LO-pQCD.
!       Extra 'u'-divergence in full one. "threshold energy" becomes 2*t_cut.
                if( (i_LOpQCD == 0 .OR. i_LOpQCD == 2) &
                    .AND. s < 2D0*t_cut )then
                    ilo = -2
                    return
                end if
                sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                alpha_s, K_factor, 6, i_LOpQCD )
                sigma_relative(6)   = sig_out
                KF_outgoing( 1, 6 ) = 21
                KF_outgoing( 2, 6 ) = 21
            end if
!       Determines the process to be happened according to the relative
!        cross section.
            sig = 0D0
            do i=0,6,1
                sig = sig + sigma_relative(i)
            end do
            rand_num = PYR(1)
            prob_low = 0D0
            prob_upp = 0D0
            do i=0,6,1
                prob_upp = prob_upp + sigma_relative(i) / sig
                if(rand_num > prob_low .AND. rand_num < prob_upp) exit
                prob_low = prob_upp
            end do
            lmn = 5
            if( i >= 1 .AND. i <= 3 )then
                lmn = 4
            else if( i == 4 .AND. i == 5 )then
                lmn = 12
            else if( i == 6 )then
                lmn = 6
            end if
            ! sig  = sig
            kf3  = KF_outgoing( 1, i )
            kf4  = KF_outgoing( 2, i )

!-------------------------------------------------------------------------------
!       Process 7 + 9 + 13, g + g ->
        else if( l_gg )then
!       't'-divergence in both full and simplified LO-pQCD.
!       Extra 'u'-divergence in full one. "threshold energy" becomes 2*t_cut.
            if( (i_LOpQCD == 0 .OR. i_LOpQCD == 2) &
                .AND. s < 2D0*t_cut )then
                ilo = -2
                return
            end if
!       Elastic process 9, g + g -> g + g.
            sig_out = sigma_2parton( t_min, t_max, s, dM, &
                            alpha_s, K_factor, 9, i_LOpQCD )
            sigma_relative(0)   = sig_out
            KF_outgoing( 1, 0 ) = kf1
            KF_outgoing( 2, 0 ) = kf2

!       Inelastic process 7 and 13, g + g -> q + qbar / hQ + hQbar.
!       Takes the threshold energy (2*mass) into account.
            if( i_inelastic == 1 .OR. i_inelastic == 3 )then

!       Inelastic Process 7, g + g -> u + ubar / d + dbar / s + sbar.
!       With the massless approximation, their cross sections are equal.
!       Inelastic Process 7, g + g -> u + ubar.
                if( sqrt_s  >=  2D0*amass(2) )then
                    sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                    alpha_s, K_factor, 7, i_LOpQCD )
                    sigma_relative(2)   = sig_out
                    KF_outgoing( 1, 2 ) =  2
                    KF_outgoing( 2, 2 ) = -2
!       Inelastic Process 7, g + g -> d + dbar.
                    if( sqrt_s  >=  2D0*amass(1) )then
                        sigma_relative(1)   = sig_out
                        KF_outgoing( 1, 1 ) =  1
                        KF_outgoing( 2, 1 ) = -1
                    end if
!       Inelastic Process 7, g + g -> s + sbar.
                    if( sqrt_s  >=  2D0*amass(3) )then
                        sigma_relative(3)   = sig_out
                        KF_outgoing( 1, 3 ) =  3
                        KF_outgoing( 2, 3 ) = -3
                    end if
                end if
!===============================================================================
!       DO NOT USE THE MASSIVE PROCESSES OF g + g -> hQ + hQbar NOW.
!       BECAUSE THEIR REVERSE PROCESSES HAVE NOT YET BEEN INTRODUCED.
                IF( .FALSE. )THEN
!       Inelastic process 13, g + g -> c + cbar.
!       Modifies the integral limits and alpha_s due to the heavy mass.
                if( sqrt_s  >=  2D0*amass(4) )then
                    dM  = amass(4)
                    dM2 = dM**2
                    t_min = dM2 - s*( 1D0 + SQRT(1D0 - 4D0*dM2/s) ) /2D0
                    t_max = dM2 - s*( 1D0 - SQRT(1D0 - 4D0*dM2/s) ) /2D0
                    Q2 = 1D0/2D0 * s
                    nf = 4
                    alpha_s = 12D0 * pi / ( 33D0 - 2D0*nf ) &
                              / LOG( Q2 / Lambda_QCD2 )
                    sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                    alpha_s, K_factor, 13, i_LOpQCD )
                    sigma_relative(4)   = sig_out
                    KF_outgoing( 1, 4 ) =  4
                    KF_outgoing( 2, 4 ) = -4
!       Inelastic process 13, g + g -> b + bbar.
!       Modifies the integral limits and alpha_s due to the heavy mass.
                    if( sqrt_s  >=  2D0*amass(5) )then
                        dM  = amass(5)
                        dM2 = dM**2
                        t_min = dM2 -s*( 1D0 +SQRT(1D0 -4D0*dM2/s) )/2D0
                        t_max = dM2 -s*( 1D0 -SQRT(1D0 -4D0*dM2/s) )/2D0
                        Q2 = 1D0/2D0 * s
                        nf = 5
                        alpha_s = 12D0 * pi / (33D0 - 2D0*nf) &
                                  / LOG( Q2 / Lambda_QCD2 )
                        sig_out = sigma_2parton( t_min, t_max, s, dM, &
                                        alpha_s, K_factor, 13, i_LOpQCD)
                        sigma_relative(5)   = sig_out
                        KF_outgoing( 1, 5 ) =  5
                        KF_outgoing( 2, 5 ) = -5
                    end if
                end if
                END IF
!===============================================================================
            end if
!       Determines the process to be happened according to the relative
!        cross section.
            sig = 0D0
            do i=0,5,1
                sig = sig + sigma_relative(i)
            end do
            rand_num = PYR(1)
            prob_low = 0D0
            prob_upp = 0D0
            do i=0,5,1
                prob_upp = prob_upp + sigma_relative(i) / sig
                if(rand_num > prob_low .AND. rand_num < prob_upp) exit
                prob_low = prob_upp
            end do
            lmn = 9
            if( i >= 1 .AND. i <= 3 ) lmn = 7
            if( i == 4 .AND. i == 5 ) lmn = 13
            ! sig  = sig
            kf3  = KF_outgoing( 1, i )
            kf4  = KF_outgoing( 2, i )

!-------------------------------------------------------------------------------
!       Process 8: q + g -> q + g.
        else if( l_qg_to_qg )then
!       't'-divergence in both full and simplified LO-pQCD.
!       Extra 'u'-divergence in full one. "threshold energy" becomes 2*t_cut.
            if( (i_LOpQCD == 0 .OR. i_LOpQCD == 2) &
                .AND. s < 2D0*t_cut )then
                ilo = -2
                return
            end if
            sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                     K_factor, 8, i_LOpQCD )
            lmn  = 8
            sig  = sig_out
            kf3  = kf1
            kf4  = kf2

!-------------------------------------------------------------------------------
!       Process 10: hQ + q -> hQ + q.
        else if( l_hQq_to_hQq )then
!       't'-divergence.
!       Modifies the integral limits and alpha_s due to the heavy mass.
            nf = MAX( ikf1, ikf2 )
            dM = amass( nf )
            dM2 = dM**2
!       How to determine the Q02 ?
            ! Q02  = 1D0
            Q02 = 0.5D0 * dM2
            ! Q02 = dM2
            Q04 = Q02**2
            Q2  = 1D0/2D0 * ( Q02 + (s - dM2)**2 / s )
            t_cut = Q02
            ! t_min = - ( s - dM2 )**2 / s - t_cut
            t_min = - ( s - dM2 )**2 / s
            t_max = -t_cut
            alpha_s = 12D0 * pi / ( 33D0 - 2D0*nf ) &
                      / LOG( Q2 / Lambda_QCD2 )
            E_threshold = dM2 + 1D0/2D0*Q02 &
                         + SQRT( dM2*Q02 + 1D0/4D0*Q04 )
            if( s > E_threshold )then
                sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                         K_factor, 10, i_LOpQCD )
                lmn  = 10
                sig  = sig_out
                kf3  = kf1
                kf4  = kf2
            end if

!-------------------------------------------------------------------------------
!       Process 11: hQ + g -> hQ + g.
        else if( l_hQg_to_hQg )then
!       't'- and 'u'-divergence.
!       Modifies the integral limits and alpha_s due to the heavy mass.
            nf = MIN( ikf1, ikf2 )
            dM = amass( nf )
            dM2 = dM**2
!       How to determine the Q02 ?
            ! Q02 = 1D0
            Q02 = 0.5D0 * dM2
            ! Q02 = dM2
            Q04 = Q02**2
            Q2  = 1D0/2D0 * ( s - dM2 )
            t_cut = Q02
            t_min = - MIN( (s - dM2 - Q02) , ( ( s - dM2 )**2 / s ) )
            t_max = -t_cut
            alpha_s = 12D0 * pi / ( 33D0 - 2D0*nf ) &
                      / LOG( Q2 / Lambda_QCD2 )
            if( Q02  <  1D0/2D0*dM2 )then
                E_threshold = dM2 + 1D0/2D0*Q02 &
                              + SQRT( dM2*Q02 + 1D0/4D0*Q04 )
            else
                E_threshold = dM2 + 2D0*Q02
            end if
            if( s >= E_threshold )then
                sig_out = sigma_2parton( t_min, t_max, s, dM, alpha_s, &
                                         K_factor, 11, i_LOpQCD )
                lmn  = 11
                sig  = sig_out
                kf3  = kf1
                kf4  = kf2
            end if

!-------------------------------------------------------------------------------
!       Ignores hQ + hQ processes.
        else if( l_hQhQ )then
            return
        end if
!------------------------   Cross Section Calculating   ------------------------
!-------------------------------------------------------------------------------


        if( sig <= 0D0 ) return


!-------------------------------------------------------------------------------
!-------------------------------   t Sampling   --------------------------------
!       tsmp: the t value sampled.
        i_sample_t = ilo
        sigma_tot  = sig
        if( i_sample_t == 1 )then
!       Optionly, samples the t value with flat scattering angle distribution.
            if( i_LOpQCD == 2 .or. i_LOpQCD == 3 )then
!               tsmp = -PYR(1) * s
                theta = PYR(1) * pi
!#TODO(Lei20240608): the following formula needs to be modified.
                tsmp  = s * ( COS(theta) - 1D0 ) / 2D0 + t_cut
            else
                t = sample_t( t_min, t_max, s, dM, lmn, i_LOpQCD, t_cut)

                IF(.false.)THEN
!       Samples t from the probability density function
!           f(t) = 1/sigma_tot * sigma(t) ,
!        in ( t_min, t_max ).
                do while(.true.)
                    rand_num_1 = PYR(1)
                    t = t_max + rand_num_1 * (t_min - t_max)
                    sigma_t = sigma_2parton( t_min, t, s, dM, alpha_s, &
                                             K_factor, lmn, i_LOpQCD )
                    rand_num_2 = PYR(2)
                    if( rand_num_2  <=  sigma_t/sigma_tot )then
                        if( ABS(t_min) > s )then
                            if( (t-t_max) > t_max ) cycle
                            t = t - t_max
                        end if
                        exit
                    end if
                end do
                END IF

                tsmp = t
            end if
        end if
!-------------------------------   t Sampling   --------------------------------
!-------------------------------------------------------------------------------


!       0.0389379 (~0.04) is the transformation factor from GeV^-2 to fm^2.
        sig = sig*0.0389379


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function sigma_2parton( t_min_in, t_max_in, s_in, mass_in, &
                 alpha_s_in, K_factor_in, i_process_in, i_LOpQCD_in )
!!      Total cross section of q(bar) / Q(bar) /g + q(bar) / Q(bar) /g .
!       The masses of the heavy c and b are retained ( M_in, M ), and
!        other masses of light quarks are approximately to 0.
!       i_process: order number of the process.
!       Massless:
!             1: q1 + q2 -> q1 + q2
!             2: q1 + q1 -> q1 + q1
!             3: q1 + q2bar -> q1 + q2bar
!             4: q1 + q1bar -> q2 + q2bar   (inelastic)
!             5: q1 + q1bar -> q1 + q1bar
!             6: q1 + q1bar -> g + g        (inelastic)
!             7: g  + g -> q1 + q1bar       (inelastic)
!             8: q  + g -> q + g
!             9: g  + g -> g + g
!       Massive:
!            10: hQ + q -> hQ + q
!            11: hQ + g -> hQ + g
!       Not yet used:
!            12: q + qbar -> hQ + hQbar   (inelastic)
!            13: g + g    -> hQ + hQbar   (inelastic)
!       Not yet included:
!            13: hQ  + hQbar  -> q + qbar       (inelastic)
!            15: hQ  + hQbar  -> g + g          (inelastic)
!            16: hQ  + hQbar  -> hQ + hQbar
!             ?: hQ1 + hQ1bar -> hQ2 + hQ2bar   (? inelastic)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER( pi=3.141592653589793D0 )
        real(kind=8) :: mass_in, M, K_factor_in, K_factor


!       The lower and upper limits of the integratoin.
        t_min = t_min_in
        t_max = t_max_in
        s = s_in
        M = mass_in
        alpha_s   = alpha_s_in
        K_factor  = K_factor_in
        i_process = i_process_in
        i_LOpQCD  = i_LOpQCD_in

        sigma_2parton = 0D0
        sigma = 0D0

!       LO-pQCD keeping leading divergent terms.
        if( i_LOpQCD == 1 .OR. i_LOpQCD == 3 )then
!           1: q1 + q2 -> q1 + q2
!           2: q1 + q1 -> q1 + q1
!           3: q1 + q2bar -> q1 + q2bar
!           5: q1 + q1bar -> q1 + q1bar
            if( i_process == 1 .OR. i_process == 2 .OR. &
                i_process == 3 .OR. i_process == 5       )then
                sigma = 8D0/9D0 * pi * alpha_s**2 &
                        * ( 1D0/t_min - 1D0/t_max )
!           4: q1 + q1bar -> q2 + q2bar
            else if( i_process == 4 )then
                sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                    * ( t_max-t_min + (t_max**3-t_min**3)/3D0/s**2 )
                ! sigma = 0D0
!           6: q1 + q1bar -> g + g
            else if( i_process == 6 )then
                sigma = 32D0/27D0 * pi * alpha_s**2 / s &
                        * LOG( ABS(t_min/t_max) )
!           7: g  + g -> q1 + q1bar
            else if( i_process == 7 )then
                sigma = 1D0/3D0 * pi * alpha_s**2 / s &
                        * LOG( ABS(t_min/t_max) )
!           8: q  + g -> q + g
            else if( i_process == 8 )then
                sigma = 2D0 * pi * alpha_s**2 &
                        * ( 1D0/t_min - 1D0/t_max )
!           9: g  + g -> g + g
            else if( i_process == 9 )then
                sigma = 9D0/2D0 * pi * alpha_s**2 &
                        * ( 1D0/t_min - 1D0/t_max )
            end if

!       Full LO-pQCD.
        else
!           1: q1 + q2    -> q1 + q2
!           3: q1 + q2bar -> q1 + q2bar
            if( i_process == 1 .OR. i_process == 3 )then
                sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                        *( &
                           primitive_func( t_max, s, M, 1 ) &
                         - primitive_func( t_min, s, M, 1 ) &
                         )
!           2: q1 + q1 -> q1 + q1
            else if( i_process == 2 )then
                sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                        *( &
                           primitive_func( t_max, s, M, 1 ) &
                         - primitive_func( t_min, s, M, 1 ) &
                         + primitive_func( t_max, s, M, 2 ) &
                         - primitive_func( t_min, s, M, 2 ) &
                         + 2D0/3D0 * ( &
                           s*LOG( ABS(t_max / ( t_max + s )) ) &
                         - s*LOG( ABS(t_min / ( t_min + s )) ) ) &
                         )
!           4: q1 + q1bar -> q2 + q2bar
            else if( i_process == 4 )then
                sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                        *( &
                           primitive_func( t_max, s, M, 4 ) &
                         - primitive_func( t_min, s, M, 4 ) &
                         )
                ! sigma = 0D0
!           5: q1 + q1bar -> q1 + q1bar
            else if( i_process == 5 )then
                sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                        *( &
                           primitive_func( t_max, s, M, 1 ) &
                         - primitive_func( t_min, s, M, 1 ) &
                         + primitive_func( t_max, s, M, 4 ) &
                         - primitive_func( t_min, s, M, 4 ) &
                         - 2D0/3D0 * ( &
                           ( s*LOG( ABS(t_max) ) &
                                    + 2D0*t_max + t_max**2/2D0/s ) &
                         - ( s*LOG( ABS(t_min) ) &
                                    + 2D0*t_min + t_min**2/2D0/s ) ) &
                         )
!           6: q1 + q1bar -> g + g
            else if( i_process == 6 )then
                sigma = pi * alpha_s**2 / s**2 &
                        *( &
                           32D0/27D0 * ( &
                           primitive_func( t_max, s, M, 6 ) &
                         - primitive_func( t_min, s, M, 6 ) ) &
                         - 8D0/3D0 * ( &
                           primitive_func( t_max, s, M, 4 ) &
                         - primitive_func( t_min, s, M, 4 ) ) &
                         )
!           7: g  + g -> q1 + q1bar
            else if( i_process == 7 )then
                sigma = pi * alpha_s**2 / s**2 &
                        *( &
                           1D0/6D0 * ( &
                           primitive_func( t_max, s, M, 6 ) &
                         - primitive_func( t_min, s, M, 6 ) ) &
                         - 3D0/8D0 * ( &
                           primitive_func( t_max, s, M, 4 ) &
                         - primitive_func( t_min, s, M, 4 ) ) &
                         )
!           8: q  + g -> q + g
            else if( i_process == 8 )then
                sigma = pi * alpha_s**2 / s**2 &
                        *( &
                         - 4D0/9D0 * ( &
                           primitive_func( t_max, s, M, 8 ) &
                         - primitive_func( t_min, s, M, 8 ) ) &
                         + primitive_func( t_max, s, M, 1 ) &
                         - primitive_func( t_min, s, M, 1 ) &
                         )
!           9: g  + g -> g + g
            else if( i_process == 9 )then
                sigma = 9D0/2D0 * pi * alpha_s**2 / s**2 &
                        *( &
                           primitive_func( t_max, s, M, 9 ) &
                         - primitive_func( t_min, s, M, 9 ) &
                         )
            end if

        end if

!       Full LO-pQCD for massive quarks.
!       Process 10, Q + q -> Q + q.
        if( i_process == 10 )then
            sigma = 4D0/9D0 * pi * alpha_s**2 / (s-M**2)**2 &
                    *( &
                       primitive_func( t_max, s, M, 10 ) &
                     - primitive_func( t_min, s, M, 10 ) &
                     )
!       Process 11, Q + g -> Q + g.
        else if( i_process == 11 )then
            sigma = 1D0/16D0 * pi * alpha_s**2 / (s-M**2)**2 &
                    *( &
                       primitive_func( t_max, s, M, 11 ) &
                     - primitive_func( t_min, s, M, 11 ) &
                     )
!       Process 12, q + qbar -> Q + Qbar.
        else if( i_process == 12 )then
            sigma = 4D0/9D0 * pi * alpha_s**2 / s**2 &
                    *( &
                       primitive_func( t_max, s, M, 12 ) &
                     - primitive_func( t_min, s, M, 12 ) &
                     )
!       Process 13, g + g -> Q + Qbar.
        else if( i_process == 13 )then
            sigma = 1D0/16D0 * pi * alpha_s**2 / s**2 &
                    *( &
                       primitive_func( t_max, s, M, 13 ) &
                     - primitive_func( t_min, s, M, 13 ) &
                     )
        end if

!       K factor enhancement.
        sigma_2parton = sigma * K_factor


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function primitive_func( t_in, s_in, mass_in, i_func_in )
!!      The primitive function F(t) of the differential x-section f(t).
!!      Note that there are no case 3, 5 and 7. They are redundant.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real(kind=8) :: mass_in, M


        t = t_in
        s = s_in
        M = mass_in
        i_func = i_func_in

        primitive_func = 0D0

        select case (i_func)
!-------------------------------------------------------------------------------
!       Massless.
!       Cf. B.L. Combridge et al., Phys. Lett. B 70 (1977) 234.

!       f(t) = ( s^2 + u^2 ) / t^2
        case(1)
            primitive_func = -2D0*s**2/t + 2D0*s*LOG( ABS(t) ) +t
!       f(t) = ( s^2 + t^2 ) / u^2
        case(2)
            primitive_func = -2D0*s**2/(s + t) + (s + t) &
                             - 2D0*s*LOG( ABS(s + t) )
!       f(t) = ( t^2 + u^2 ) / s^2
        case(4)
            primitive_func = 1D0/3D0/s**2 * ( t**3 + (s+t)**3 )
!       f(t) = ( u^2 + t^2 ) / ut
        case(6)
            primitive_func = s*LOG( ABS( (s+t)/t ) ) -2D0*(s+t)
!       f(t) = ( u^2 + s^2 ) / us
        case(8)
            primitive_func = -(s + t)**2/2D0/s -s*LOG( ABS(s + t) )
!       f(t) = 3 - u*t/s^2 - u*s/t^2 - s*t/u^2
        case(9)
            primitive_func = 3D0*t + t**2/2D0/s + t**3/3D0/s**2 - s**2/t &
                             + s*LOG( ABS( t / (s + t) ) ) - s**2/(s+t)

!-------------------------------------------------------------------------------
!       Massive.
!       Cf. B.L. Combridge Nucl. Phys. B 151 (1979) 429-456.

!       f(t) = [ ( M^2 - u )^2 + ( s - M^2 )^2 + 2 * M^2 * t ] / t^2
        case(10)
            primitive_func = t + 2D0*s*LOG( ABS(t) ) &
                             - 2D0*( s - M**2 )**2 / t
!       f(t) = 32*[ (s-M^2)*(M^2-u) ] / t^2
!              + 64/9*[ (s-M^2)*(M^2-u) + 2*M^2*(s+M^2) ] / (s-M^2)^2
!              + 64/9*[ (s-M^2)*(M^2-u) + 2*M^2*(M^2+u) ] / (M^2-u)^2
!              + 16/9*[ M^2*(4*M^2-t) ] / [ (s-M^2)*(M^2-u) ]
!              + 16*[ (s-M^2)*(M^2-u) + M^2*(s-u) ] / [ t*(s-M^2) ]
!              - 16*[ (s-M^2)*(M^2-u) - M^2*(s-u) ] / [ t*(M^2-u) ]
        case(11)
            T1 = 32D0 *(s - M**2)*( -(s - M**2)/t + LOG(ABS(t)) )
            T2 = 64D0/9D0 / (s - M**2)**2 &
                 *( (s -M**2)*(s +t -M**2)**2/2D0 + 2*M**2*(s +M**2)*t )
            T3 = 64D0/9D0 *( (s -3*M**3)*LOG( ABS(s + t - M**2) ) &
                 - 4D0*M**4 / (s + t - M**2) )
            T4 = 16D0/9D0 * M**2/(s - M**2) * ( - (s + t - M**2) &
                 + (3D0*M**2 + s)*LOG( ABS(s + t - M**2) ) )
            T5 = 16D0 * ( s/(s - M**2)*t &
                 + ( (s - M**2) + M**2*(2D0*s - 2D0*M**2)/(s - M**2) ) &
                 * LOG( ABS(t) ) )
            T6 = 16D0 * ( (s - 3D0*M**2)*LOG( ABS(t) ) &
                 + M**2*LOG( ABS(s+t-M**2) ) )
            primitive_func = T1 + T2 + T3 + T4 + T5 - T6
!       f(t) = [ (M^2-t)^2 + (M^2-u)^2 + 2*M^2*s ] / s^2
        case(12)
            primitive_func = 1D0/s**2 * ( -( M**2 - t )**3/3D0 &
                             + ( s + t - M**2 )**3/3D0 + 2D0*M**2*s*t )
!       f(t) = 12/s^2 * (M^2-t) * (M^2-u)
!              + 8/3 * [ (M^2-t)*(M^2-u) - 2*M^2*(M^2+t) ] / (M^2-t)^2
!              + 8/3 * [ (M^2-t)*(M^2-u) - 2*M^2*(M^2+u) ] / (M^2-u)^2
!              - 2/3*M^2 * (s-4*M^2) / [ (M^2-t)*(M^2-u) ]
!              - 6 * [ (M^2-t)*(M^2-u) + M^2*(u-t) ] / [ s*(M^2-t)]
!              - 6 * [ (M^2-t)*(M^2-u) + M^2*(t-u) ] / [ s*(M^2-u)]
        case(13)
            T1 = 12D0/s**2 * ( (M**2*s - M**4) &
                 + (2D0*M**2 - s)/2D0*t**2 -1D0/3D0*t**3 )
            T2 = 8D0/3D0 * ( M**2 - t + 4D0*M**4/(t - M**2) &
                 - (s + 2D0*M**2)*LOG( ABS(t - M**2) ) )
            T3 = 8D0/3D0 * ( M**2 - s - t + 4D0*M**4/(s + t - M**2) &
                 + (s + 2D0*M**2)*LOG( ABS(s + t - M**2) ) )
            T4 = 2D0/3D0 * M**2 * (s - 4D0*M**2) / s &
                 * LOG( ABS( (s + t - M**2)/(t - M**2) ) )
            T5 = 6D0/s * ( (s + M**2)*t + t**2/2D0 &
                 + M**2*s*LOG( ABS(t - M**2) ) )
            T6 = 6D0/s * ( 2D0*(s - M**2) + (2D0 + M**2)*t -1D0/2D0*t**2 &
                 + (M**2*s - 2D0*s + 2D0*M**2 -2D0*M**4) &
                 *LOG( ABS(s +t -M**2) ) )
            primitive_func = T1 + T2 + T3 - T4 - T5 - T6
        case default
            primitive_func = 0D0
        end select


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function sample_t( t_min_in, t_max_in, s_in, mass_in, &
                 i_process_in, i_LOpQCD_in, t_cut_in )
!!      Samples t from the regularized differential cross section
!!       distribution dsigma/dt.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real(kind=8) :: mass_in, M, mu2


!       The lower and upper limits.
        t_min = t_min_in
        t_max = t_max_in
        s = s_in
        M = mass_in
        i_process = i_process_in
        i_LOpQCD  = i_LOpQCD_in
        t_cut = t_cut_in
        mu2   = t_cut_in
        Q02   = t_cut_in

        sample_t = 0E0
        t = 0D0
!       Maximum number of sapmling.
        n_sample_max = 100000

!       LO-pQCD keeping leading divergent terms.
        if( i_LOpQCD == 1 .OR. i_LOpQCD == 3 )then
!           1: q1 + q2 -> q1 + q2
!           2: q1 + q1 -> q1 + q1
!           3: q1 + q2bar -> q1 + q2bar
!           5: q1 + q1bar -> q1 + q1bar
!           8: q  + g -> q + g
!           9: g  + g -> g + g
!           With the common form of f(t) = 1/(t - mu2)^2. Direct sampling.
            if( i_process == 1 .OR. i_process == 2 .OR. &
                i_process == 3 .OR. i_process == 5 .OR. &
                i_process == 8 .OR. i_process == 9      )then
                ! do while(.true.)
                do i_sample = 1, n_sample_max, 1
                    rand_num = PYR(1)
                    t = mu2*( 1D0 - (s + mu2)/(s*rand_num + mu2) )
!       Limits t in [ -s, -mu2 ] ?
                    ! if( t <= -mu2 .AND. t >= -s ) exit
                    exit
                end do
!           4: q1 + q1bar -> q2 + q2bar
!           With the form of f(t) = t^2 + s^2. Rejection sampling.
            else if( i_process == 4 )then
                ! do while(.true.)
                do i_sample = 1, n_sample_max, 1
                    rand_num_1 = PYR(1)
                    t = t_max + rand_num_1 * (t_min - t_max)
                    dsigma_dt     = t**2 + s**2
                    dsigma_dt_max = s**2 + s**2
                    rand_num_2 = PYR(2)
                    if( rand_num_2  <=  dsigma_dt/dsigma_dt_max ) exit
                end do
!           6: q1 + q1bar -> g + g
!           7: g  + g -> q1 + q1bar
!           With the common form of f(t) = -1/(t - mu2). Direct sampling.
            else if( i_process == 6 .OR. i_process == 7 )then
                ! do while(.true.)
                do i_sample = 1, n_sample_max, 1
                    rand_num = PYR(1)
                    t = mu2 - (mu2 + s)*( mu2/(mu2 + s) )**rand_num
!       Limits t in [ -s, -mu2 ] ?
                    ! if( t <= -mu2 .AND. t >= -s ) exit
                    exit
                end do
            end if

!       Full LO-pQCD for massless processes.
!       Uses the rejection sampling method.
        else if( i_process >= 1 .AND. i_process <= 9 )then
            ! do while(.true.)
            do i_sample = 1, n_sample_max, 1
                rand_num_1 = PYR(1)
                t = t_max + rand_num_1 * (t_min - t_max)
                dsigma_dt     = dsigma_dt_func( t, s, M, i_process )
                dsigma_dt_max = MAX( &
                              dsigma_dt_func( t_min, s, M, i_process ) , &
                              dsigma_dt_func( t_max, s, M, i_process ) )
                rand_num_2 = PYR(2)
                if( rand_num_2  <=  dsigma_dt/dsigma_dt_max ) exit
            end do
        end if

!       Full LO-pQCD for massive quarks.
!       Uses the rejection sampling method.
!       Process 10, Q + q -> Q + q.
!       Process 11, Q + g -> Q + g.
!       Process 12, q + qbar -> Q + Qbar.
!       Process 13, g + g -> Q + Qbar.
        if( i_process >= 10 )then
            ! do while(.true.)
            do i_sample = 1, n_sample_max, 1
                rand_num_1 = PYR(1)
                t = t_max + rand_num_1 * (t_min - t_max)
                dsigma_dt     = dsigma_dt_func( t, s, M, i_process )
                dsigma_dt_max = MAX( &
                              dsigma_dt_func( t_min, s, M, i_process ) , &
                              dsigma_dt_func( t_max, s, M, i_process ) )
                rand_num_2 = PYR(2)
                if( rand_num_2  <=  dsigma_dt/dsigma_dt_max ) exit
            end do
        end if

!       Hard to sample t. Gives up it.
        ! if( i_sample == n_sample_max ) t = 0D0

        sample_t = t


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function dsigma_dt_func( t_in, s_in, mass_in, i_process_in )
!!      The differential cross section f(t) = dsigma/dt.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        real(kind=8) :: mass_in, M


        t = t_in
        s = s_in
        M = mass_in
        u = 2*M**2 - s - t
        i_process = i_process_in

        dsigma_dt_func = 0D0

        select case (i_process)
!-------------------------------------------------------------------------------
!       Massless.
!       Cf. B.L. Combridge et al., Phys. Lett. B 70 (1977) 234.

!       f(t) = 4/9*( s^2 + u^2 ) / t^2
        case(1)
            dsigma_dt_func = 4D0/9D0 * ( s**2 + u**2 ) / t**2
!       f(t) = 4/9*[ (s^2 + u^2)/t^2 + (s^2 + t^2)/u^2 ] - 8/27 * s^2/(ut)
        case(2)
            dsigma_dt_func = 4D0/9D0*( (s**2 + u**2)/t**2 &
                             + (s**2 + t**2)/u**2 ) - 8D0/27D0 *s**2/u/t
!       f(t) = 4/9*( s^2 + u^2 ) / t^2
        case(3)
            dsigma_dt_func = 4D0/9D0 * ( s**2 + u**2 ) / t**2
!       f(t) = 4/9*( t^2 + u^2 ) / s^2
        case(4)
            dsigma_dt_func = 4D0/9D0 * ( t**2 + u**2 ) / s**2
!       f(t) = 4/9*[ (s^2 + u^2)/t^2 + (t^2 + u^2)/s^2 ] - 8/27 * u^2/(st)
        case(5)
            dsigma_dt_func = 4D0/9D0*( (s**2 + u**2)/t**2 &
                             + (t**2 + u**2)/s**2 ) - 8D0/27D0 *u**2/s/t
!       f(t) = 32/27*( u^2 + t^2 ) / (ut) - 8/3*( u^2 + t^2 ) / s^2
        case(6)
            dsigma_dt_func = 32D0/27D0*(u**2 + t**2)/u/t &
                             - 8D0/3D0*(u**2 + t**2)/s**2
!       f(t) = 1/6*( u^2 + t^2 ) / (ut) - 3/8*( u^2 + t^2 ) / s^2
        case(7)
            dsigma_dt_func = 1D0/6D0*(u**2 + t**2)/u/t &
                             - 3D0/8D0*(u**2 + t**2)/s**2
!       f(t) = -4/9*( u^2 + s^2 ) / (us) + ( u^2 + s^2 ) / t^2
        case(8)
            dsigma_dt_func = -4D0/9D0*(u**2 + s**2)/u/s &
                             + (u**2 + s**2)/t**2
!       f(t) = 3 - u*t/s^2 - u*s/t^2 - s*t/u^2
        case(9)
            dsigma_dt_func = 9D0/2D0*(3D0 -u*t/s**2 -u*s/t**2 -s*t/u**2)

!-------------------------------------------------------------------------------
!       Massive.
!       Cf. B.L. Combridge Nucl. Phys. B 151 (1979) 429-456.

!       f(t) = [ ( M^2 - u )^2 + ( s - M^2 )^2 + 2 * M^2 * t ] / t^2
        case(10)
            dsigma_dt_func = ( (M**2 - u)**2 + (s - M**2)**2 &
                             + 2D0*M**2*t ) / t**2
!       f(t) = 32*[ (s-M^2)*(M^2-u) ] / t^2
!              + 64/9*[ (s-M^2)*(M^2-u) + 2*M^2*(s+M^2) ] / (s-M^2)^2
!              + 64/9*[ (s-M^2)*(M^2-u) + 2*M^2*(M^2+u) ] / (M^2-u)^2
!              + 16/9*[ M^2*(4*M^2-t) ] / [ (s-M^2)*(M^2-u) ]
!              + 16*[ (s-M^2)*(M^2-u) + M^2*(s-u) ] / [ t*(s-M^2) ]
!              - 16*[ (s-M^2)*(M^2-u) - M^2*(s-u) ] / [ t*(M^2-u) ]
        case(11)
            T1 = 32D0 * (s - M**2) * (M**2 - u) / t**2
            T2 = 64D0/9D0 * ( (s - M**2)*(M**2 - u) &
                 + 2D0*M**2*(s + M**2) ) / (s - M**2)**2
            T3 = 64D0/9D0 * ( (s - M**2)*(M**2 - u) &
                 + 2D0*M**2*(M**2 + u) ) / (M**2 - u)**2
            T4 = 16D0/9D0 * M**2*(4D0*M**2 - t) &
                 / (s - M**2) /(M**2 - u)
            T5 = 16D0 * ( (s - M**2)*(M**2 - u) &
                 + M**2*(s - u) ) / t / (s - M**2)
            T6 = 16D0 * ( (s - M**2)*(M**2 - u) &
                 - M**2*(s - u) ) / t / (M**2 - u)
            dsigma_dt_func = T1 + T2 + T3 + T4 + T5 - T6
!       f(t) = [ (M^2-t)^2 + (M^2-u)^2 + 2*M^2*s ] / s^2
        case(12)
            dsigma_dt_func = ( (M**2 - t)**2 + (M**2 - u)**2 &
                             + 2D0*M**2*s ) / s**2
!       f(t) = 12/s^2 * (M^2-t) * (M^2-u)
!              + 8/3 * [ (M^2-t)*(M^2-u) - 2*M^2*(M^2+t) ] / (M^2-t)^2
!              + 8/3 * [ (M^2-t)*(M^2-u) - 2*M^2*(M^2+u) ] / (M^2-u)^2
!              - 2/3*M^2 * (s-4*M^2) / [ (M^2-t)*(M^2-u) ]
!              - 6 * [ (M^2-t)*(M^2-u) + M^2*(u-t) ] / [ s*(M^2-t)]
!              - 6 * [ (M^2-t)*(M^2-u) + M^2*(t-u) ] / [ s*(M^2-u)]
        case(13)
            T1 = 12D0/s**2 * ( M**2 - t ) * ( M**2 - u )
            T2 = 8D0/3D0 * ( (M**2 - t)*(M**2 - u) &
                  - 2D0*M**2*(M**2 + t) ) / (M**2 - t)**2
            T3 = 8D0/3D0 * ( (M**2 - t)*(M**2 - u) &
                  - 2D0*M**2*(M**2 + u) ) / (M**2 - u)**2
            T4 = 2D0/3D0 * M**2 * (s - 4D0*M**2) &
                 / (M**2 - t) / (M**2 - u)
            T5 = 6D0 * ( (M**2 - t)*(M**2 - u) + M**2*(u - t) ) &
                 / s / (M**2 - t)
            T6 = 6D0 * ( (M**2 - t)*(M**2 - u) + M**2*(t - u) ) &
                 / s / (M**2 - u)
            dsigma_dt_func = T1 + T2 + T3 - T4 - T5 - T6
        case default
            dsigma_dt_func = 0D0
        end select


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine collis(ik1,ik2,kf3,kf4,tcp,jjj,kkk,iway,icnew,jcnew,lmn,time)
!!      Performs parton-parton collision & updates particle list
!       If lmn = 4, 6, & 7, updates 'PYJETS',
!        'sbe',diquark list, string list, & lc(1-2,m), too
!       ik1,ik2: line number of the colliding pair in parton list
!       kf3 and kf4: kf code of the colliding pair after collision
!       tcp: collision time
!       jjj: jjj-th loop within a event
!       if kkk=1 throw away current collision
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MPLIS=80000,MSCA=20000)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/sa33/smadel,ecce,secce,parecc,iparres
        common/work7/reac(20),crose(20)
        common/ctllist_p/nreac(20),nrel
        common/show/vip(mplis),xap(mplis)
        common/collist/lc(2,MCLIS),tc(2,MCLIS),icol
        common/e_cons1/ pss(kszj,5), ff(kszj)
!       ifcom(i): line number (in 'sbe') of first component of i-th diquark
!       nreac(i): statistics of the # of successful i-th collision
!       nrel: statistics of the # of collision blocked
!       nsca: number of particles after collision
!       pip(1-3,*): momentum of particle after collision
!       pip(4,*): energy of particle after collision
!       pip(5,*): virtuality of particle after collision
!       pip(6,*): x value of particle after collision
!       kpip(*): flavor code of particle after collision
!       vip: virtuality of parton
!       xap: momentum fraction of parton
!       ppsa(5): cumulate the charge losed in collision processes in an event
!       definition of ppsa(1), ..., and ppsa(4) are given in 'eloss'
        dimension pi(4), pj(4), pii(4), pjj(4), pi_o(4), pj_o(4), ps(4)
        dimension px(4), py(4), pij(4), b(3)


!       The hadronization model.
        i_had_model = INT( adj1(12) )
        kf1 = K(ik1,2)
        kf2 = K(ik2,2)
!       kf1 and kf2: kf code of the colliding pair before interaction
        stili = 0D0
!       stili: statistics of the times calling "ti_li1"
        tl0 = adj1(24)
!       tl0: cut off virtuality of time-like branching, i. e. Mu0**2

!       statistics of the total charge of partons
        ichabe=0
        do i=1,N,1
            ik = K(i,2)
            ichabe = ichabe + ichge(ik)
        end do

        nsca = 0
        pi(4) = P(ik1,4)
        pj(4) = P(ik2,4)
        pi_o(4) = pi(4)
        pj_o(4) = pj(4)
        pij4 = pi(4) + pj(4)
        do i=1,3,1
            pi(i) = P(ik1,i)
            pj(i) = P(ik2,i)
            pi_o(i) = pi(i)
            pj_o(i) = pj(i)
            b(i) = ( pi(i) + pj(i) ) / pij4
        end do
!       Boosts to the current CMS of colliding pair.
        call lorntz(0,b,pi,pj)
        pij(1) = pi(1) + pj(1)
        pij(2) = pi(2) + pj(2)
        pij(3) = pi(3) + pj(3)
        pij(4) = pi(4) + pj(4)
        eiej2 = dot( pij, pij )
!       squared invariant mass of colliding pair
!       call lorntz(1,b,pi,pj)
        if( eiej2 < 0D0 )then
!       iway=1 means collision does not happen
            iway = 1
            return
        endif
!       calculate the total cross section, decide the type of reaction
!        (when ilo=0) and sample the t value as well (when ilo=1)
        ilo=1
        ikf1 = ABS(kf1)
        ikf2 = ABS(kf2)
        if( (ikf1 <= 6 .OR. ikf1 == 21) .AND. (ikf2 <= 6 .OR. ikf2 == 21) )then
!       d, u, s, c, b quarks, their anti quarks, and gluon only
            call fsig(ilo,kf1,kf2,kf3,kf4,eiej2,sig,tsmp,lmn,jjj)
        else
            sig = 0D0
        end if
!       lmn: order number of the process happened
        if( sig <= 0D0 )then
            ! throw away that collision
            jjj = jjj - 1
            kkk = 1
            return
        end if
!       am1 = amass(kf1)
!       am2 = amass(kf2)
!       am3 = amass(kf3)
!       am4 = amass(kf4)
        am1 = 0D0
        am2 = 0D0
        am3 = 0D0
        am4 = 0D0
!       in consistent with sm in "fsig"
!       kinematics of 2->2 process,no matter what is the final state, is
!        treated in zero mass approximation from momentum point of view
        paa = SQRT( pi(1)**2 + pi(2)**2 + pi(3)**2 )
        xs = eiej2
        xt = tsmp

!       Takes heavy masses into account ("fsig_new") now.
        ikf3 = ABS( kf3 )
        ikf4 = ABS( kf4 )
        if( ikf1 == 4 .OR. ikf1 == 5 ) am1 = amass( ikf1 )
        if( ikf2 == 4 .OR. ikf2 == 5 ) am2 = amass( ikf2 )
        if( ikf3 == 4 .OR. ikf3 == 5 ) am3 = amass( ikf3 )
        if( ikf4 == 4 .OR. ikf4 == 5 ) am4 = amass( ikf4 )
        xu = am1**2 + am2**2 + am3**2 + am4**2 - xs - xt
!       zero rest mass approximation
!       xu = - xs - xt

!       square of momentum transfer (Q**2)
        ! pT**2, as the same in PYTHIA (Q^2_{hard})
        qc2 = xt * xu / xs
        ! as the same in PYTHIA
!       qc2 = 4D0 * qc2
!       qc2 = min1( -xt, -xu )

!       calculate the direction cosines of momentum of one colliding
!        particle after scattering relative to the momentum
!        of corresponding particle before scattering
        fi = 2D0 * pio * PYR(1)
!       relation between t and cos(ceta) in zero rest mass approximation
!        (i. e. treated as ela. scattering)
        if( paa < 1D-10 ) paa = 1D-10

!       cctas = 1D0 + xt / 2D0 / paa**2
!       Uses the full form.
        cctas = 1D0 + 2D0 * xs * xt / ( xs - am1 - am2 ) / ( xs - am1 + am2 )

        if( ABS(cctas) > 1D0 )then
            cctas = cctas / ABS( cctas )
        end if
        sctas = SQRT( 1D0 - cctas**2 )
        cfis  = COS(fi)
        sfis  = SIN(fi)
        pc1   = paa
        call rotate(cctas,sctas,cfis,sfis,pc1,pi,pj)
!       pi and pj, as input in 'rotate', are four momentum of colliding
!        pair (before scattering and rotation), as output are four
!        momentum of that pair after scattering and rotation
!       pc1: momentum modulus of pi or pj, both are equal in their cms
!        before and after scattering
!       pi(4) = SQRT( pi(1)*pi(1) + pi(2)*pi(2) + pi(3)*pi(3) + am3*am3 )
!       pj(4) = SQRT( pj(1)*pj(1) + pj(2)*pj(2) + pj(3)*pj(3) + am4*am4 )

!       initial four momentum of parton-parton scattering.
        do i=1,4,1
            pii(i) = pi(i)
            pjj(i) = pj(i)
        end do
!       boost back
        call lorntz(1,b,pi,pj)

!       proceeds for w/o time-like branching

!       final state of two parton scattering (nsca=2)
        nsca = 2
        kpip(1) = kf3
        kpip(2) = kf4
        do i=1,4,1
            pip(i,1) = pi(i)
            pip(i,2) = pj(i)
        end do

!       charge conservation
        ichbe = ichge(kf1) + ichge(kf2)
        ichaf = ichge(kf3) + ichge(kf4)
        if( ichbe /= ichaf )then
            write(19,*) "w/o time-like iii, jjj, kf1, kf2=", iii, jjj, kf1, kf2
            write(19,*) "nsca, kf3, kf4, ichbe, ichaf=",nsca,kf3,kf4,ichbe,ichaf
            ppsa(5) = ppsa(5) - ichbe + ichaf
        end if

!       updates particle list  (note: line numbers of the colliding
!        pair after interaction are the same as ones before interaction)

!       updates 'PYJETS', i.e. feedback final state of parton scattering
!        to ik1 & ik2
        l = ik1
        if( pip(4,1) < 1D-10 ) pip(4,1) = 1D-10
        do k1=1,3,1
            P(l,k1) = pip(k1,1)
        end do
        P(l,4) = pip(4,1)
        K(l,2) = kpip(1)
!       Only for the inelastic processes.
        if( kf3 /= kf1 ) P(l,5) = amass( kpip(1) )

        l = ik2
        if( pip(4,2) < 1D-10 ) pip(4,2) = 1D-10
        do k1=1,3,1
            P(l,k1) = pip(k1,2)
        end do
        P(l,4) = pip(4,2)
        K(l,2) = kpip(2)
!       Only for the inelastic processes.
        if( kf4 /= kf2 ) P(l,5) = amass( kpip(2) )

        reac(lmn)  = reac(lmn)  + 1D0
        crose(lmn) = crose(lmn) + sig
        nreac(lmn) = nreac(lmn) + 1
        ! Only elastic processes.
        if( iparres == 0 .OR. iparres == 2 ) return

        ! Not 4, 6 & 7 processes.
        if( ( iparres == 1 .OR. iparres == 3 ) .AND. &
            ( lmn /= 4 .AND. lmn /= 6 .AND. lmn /= 7) ) return
!#TODO(Lei20241018): need to be extended for heavy quarks.

!       proceeds for lmn=4 or =6 or =7 and hadronization by string
!        fragmentation

        if(lmn == 4 .or. lmn == 6 .or. lmn == 7)then 
!       update collision list
            call update_ctlm(time)
            ! w/ time-like branching
            if( i_time_shower == 1 ) goto 1004

1008        continue
!       copies new string composed of final state parton pair (ik1 & ik2)
!        to the end of 'PYJETS'
            if(ik2 > ik1)then
                K(ik1,1) = 2
                K(ik2,1) = 1
!       copy ik1 to N+1
                call coend(ik1)
                icnew = N
!       copy ik2 to n+1
                call coend(ik2)
                jcnew = N
            end if
            if(ik2 < ik1)then
                K(ik1,1) = 1
                K(ik2,1) = 2
!       copy ik2 to n+1
                call coend(ik2)
                jcnew = N
!       copy ik1 to n+1
                call coend(ik1)
                icnew = N
            end if
            nstr1 = nstr1 + 1
            nstr1a(nstr1) = N - 1
            nstr1v(nstr1) = N

!       handles scattering parton pair (scattered parton pair has been
!        treated as a new string above)

!       7: gg->qqbar
            if(lmn == 7)then
!       throws away scattering parton pair (ik1 & ik2)
                if(ik1 > ik2)then
                    call parmov(ik1,ik1,ik2,lmn)
                    call parmov(ik2,ik1,ik2,lmn)
                end if
                if(ik1 < ik2)then
                    call parmov(ik2,ik1,ik2,lmn)
                    call parmov(ik1,ik1,ik2,lmn)
                endif
                ! New partons are above n-2 w/ or w/o branching.
                n00 = N - 2
                call updpli_p(n00,time)
                return
            end if

!       4: q1q1bar->q2q2bar, 6: qqbar->gg
            if(lmn == 4 .OR. lmn == 6)then
!       hadronization with coalescence
                if( i_had_model /= 0 .AND. i_had_model /= 3 )then
!       throws away scattering parton pair (ik1 & ik2)
                    if(ik1 > ik2)then
                        call parmov(ik1,ik1,ik2,lmn)
                        call parmov(ik2,ik1,ik2,lmn)
                    endif
                    if(ik1 < ik2)then
                        call parmov(ik2,ik1,ik2,lmn)
                        call parmov(ik1,ik1,ik2,lmn)
                    endif
                    ! New partons are above n-2 w/ or w/o branching.
                    n00 = N - 2
                    call updpli_p(n00,time)
                    return
!       hadronization with string fragmentation.
                elseif( i_had_model == 0 .OR. i_had_model == 3 )then
!       is ik1 (ik2) a component of string ?
                    call adjst(ik1,ik1str,ik1sa,ik1sv)
!       ik1str: oredr number of string to which ik1 belongs, equal 0 otherwise
!       ik1sa: line number of first component of above string, equal 0 otherwise
!       ik1sv: line number of last component of above string, equal 0 otherwise
                    call adjst(ik2,ik2str,ik2sa,ik2sv)
!       is ik1 a component of diquark ?
                    call adjdi(ik1,idi1,iway1)
                    ifcom1 = 0
                    npt1   = 0
                    if(iway1 == 2) npt1 = npt(idi1)
!       npt1: line # (in 'sbe') of partner of idi1-th diquark
                    if(iway1 == 1) ifcom1 = ifcom(idi1)
!       ifcom1: line # (in 'sbe') of first component of idi1-th diquark
!       does ik2 is a component of diquark ?
                    call adjdi(ik2,idi2,iway2)
                    ifcom2 = 0
                    npt2 = 0
                    if(iway2 == 2) npt2   = npt(idi2)
                    if(iway2 == 1) ifcom2 = ifcom(idi2)
!       ifcom1: line # of first component of idi1-th diquark
!       removes diquarks (to which ik1 & ik2 belong) from diquark list.
!        It doesn't affect other lists.
                    if( (idi1 /= 0 .and. idi2 /= 0) .and. idi1 == idi2 ) &
                        call diqmov(idi1)
                    if( (idi1 /= 0 .and. idi2 /= 0) .and. idi1 > idi2 )then
                        call diqmov(idi1)
                        call diqmov(idi2)
                    end if
                    if( (idi1 /= 0 .and. idi2 /= 0) .and. idi1 < idi2 )then
                        call diqmov(idi2)
                        call diqmov(idi1)
                    end if
                    if(idi1 /= 0 .and. idi2 == 0) call diqmov(idi1)
                    if(idi1 == 0 .and. idi2 /= 0) call diqmov(idi2)
!       moves conponents of string out & updates lists
                    if(ik1str == 0 .and. ik2str == 0)then
!       moves ik1 & ik2 out as well as updates lists
                        if(ik1 > ik2)then
                            call parmov(ik1,ik1,ik2,lmn)
                            call parmov(ik2,ik1,ik2,lmn)
                        end if
                        if(ik1 < ik2)then
                            call parmov(ik2,ik1,ik2,lmn)
                            call parmov(ik1,ik1,ik2,lmn)
                        end if
                        n00 = N - 2
                        call updpli_p(n00,time)
                        return
                    end if
                    if(ik1str == 0 .and. ik2str /= 0)then
                        if(ik1 < ik2)then
                            call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)
                            call parmov(ik1,ik1,ik2,lmn)
                        end if
                        if(ik1 > ik2)then
                            call parmov(ik1,ik1,ik2,lmn)
                            call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)
                        end if
                        n00 = N - ( ik2sv - ik2sa + 1 ) - 1
                        call updpli_p(n00,time)
                        return
                    end if
                    if(ik1str /= 0 .and. ik2str == 0)then
                        if(ik1 < ik2)then
                            call parmov(ik2,ik1,ik2,lmn)
                            call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)
                        else
                            call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)
                            call parmov(ik2,ik1,ik2,lmn)
                        endif
                        n00 = N - ( ik1sv - ik1sa + 1 ) - 1
                        call updpli_p(n00,time)
                        return
                    end if
!       proceeds for both of ik1str & ik2str are not equal to zero
                    if(ik1str == ik2str)then
                        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)
                        n00 = N - ( ik1sv - ik1sa + 1 )
                        call updpli_p(n00,time)
                    end if
                    if(ik1str > ik2str)then
                        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)
                        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)
                        n00 = N - ( ik1sv - ik1sa + 1 ) - ( ik2sv - ik2sa + 1 )
                        call updpli_p(n00,time)
                    end if
                    if(ik1str < ik2str)then
                        call strmov(ik2str,ik2sa,ik2sv,ik1,ik2,lmn)
                        call strmov(ik1str,ik1sa,ik1sv,ik1,ik2,lmn)
                        n00=n-(ik1sv-ik1sa+1)-(ik2sv-ik2sa+1)
                        call updpli_p(n00,time)
                    end if
                    return
                end if
            end if

1004        continue
!       initial state of time-like branching
            nsca = 2
            kpip(1) = kf3
            kpip(2) = kf4
            do i=1,4
                pip(i,1) = pi(i)
                pip(i,2) = pj(i)
            enddo
            pip(5,1) = qc2
            pip(5,2) = qc2
            pip(6,1) = 1D0
            pip(6,2) = 1D0

!0700   performs time_like branching -----------------------------------
            ichtb = 0
            do i1=1,nsca,1
                kff = kpip(i1)
                ichtb = ichtb + ichge(kff)
            end do
            nsca0 = 0
! 300       continue
            kmsca = nsca
            do i = nsca0+1, kmsca, 1
                if( pip(5,i) > 4*tl0 )then
!       phase space for branching vanishes at a virtuality of 4*tl0
!        cf. Eq.25 in B. R. Webber, Ann. Rev. Nucl. Part. Sci., 36(86)253
                    ea = pip(4,i)
                    qa = pip(5,i)
                    do j=1,3
                        pa(j) = pip(j,i)
                    enddo
                    kfk = kpip(i)
                    xa = pip(6,i)
!       time-like branching for the final state partons & the induced
!        partons formed in time like branchings
                    call ti_li1( stili, i )
                    stili = stili + 1D0
!       stili: statistics of the times calling 'ti_li1'
!-------update of i------------
!                   kpip(i) = kfk
!                   do j=1,3
!                       pip(j,i) = pa(j)
!                   enddo
!                   pip(4,i) = ea
!                   pip(5,i) = qa
!                   pip(6,i) = xa
!-------finished----------------
                end if
            end do
            nsca0 = kmsca
!           if( msca0 /= nsca ) goto 300
            icht = 0
            do i1=1,nsca
                kff  = kpip(i1)
                icht = icht + ichge(kff)
            end do
            if( ichtb /= icht )then
                write(19,*) "af. ti_li1 iii,nsca,ichtb,icht", &
                                       iii,nsca,ichtb,icht
                ppsa(5) = ppsa(5) - ichtb + icht
            end if
!0700   performs time_like branching finished---------------------------

!0800   'nsca' to 'PYJETS' ---------------------------------------------
!       two-body scattering: kf1 + kf2 -> kf3 + kf4
            if(nsca == 2)then
!       updates 'PYJETS', i.e. feedback final state of two parton scattering
!        to ik1 & ik2
                l = ik1
                if( pip(4,1) < 1D-10 ) pip(4,1) = 1D-10
                do k1=1,3
                    P(l,k1) = pip(k1,1)
                end do
                P(l,4)  = pip(4,1)
                K(l,2)  = kpip(1)
                P(l,5)  = amass( kpip(1) )
                taup(l) = 0D0
                l = ik2
                if( pip(4,2) < 1D-10 ) pip(4,2) = 1D-10
                do k1=1,3
                    P(l,k1) = pip(k1,2)
                end do
                p(l,4)  = pip(4,2)
                k(l,2)  = kpip(2)
                p(l,5)  = amass( kpip(2) )
                taup(l) = 0D0
                goto 1008
!       branching: g->gg, g->qqbar, & q-> qg (qbar->qbarg)
!       if considering only g->gg & no g->qqbar,then it is possible to
!        construct chains of kf3->kf3+kf5, kf4->kf4+kf6
!       if nsca=3, makes 'nsca' in order of kf3 (or kf4),kf5,kf6
!       if nsca=4, makes 'nsca' in order of kf3,kf4,kf5,kf6
            else if(nsca == 3)then
                do i1=1,nsca,1
                    if( kpip(i1) == kf3 .or. kpip(i1) == kf4 )then
!       moves particle i1 to the first position 'nsca'
                        kpip1   = kpip(1)
                        pip1    = pip(1,1)
                        pip2    = pip(2,1)
                        pip3    = pip(3,1)
                        pip4    = pip(4,1)
                        pip5    = pip(5,1)
                        pip6    = pip(6,1)
                        kpip(1) = kpip(i1)
                        do i2=1,6,1
                            pip(i2,1) = pip(i2,i1)
                        end do
                        kpip(i1)  = kpip1
                        pip(1,i1) = pip1
                        pip(2,i1) = pip2
                        pip(3,i1) = pip3
                        pip(4,i1) = pip4
                        pip(5,i1) = pip5
                        pip(6,i1) = pip6
                    end if
                end do
            else if(nsca == 4)then
                do i1=1,nsca,1
                    if( kpip(i1) == kf3 )then
!       moves particle i1 to the first position 'nsca'
                        kpip1   = kpip(1)
                        pip1    = pip(1,1)
                        pip2    = pip(2,1)
                        pip3    = pip(3,1)
                        pip4    = pip(4,1)
                        pip5    = pip(5,1)
                        pip6    = pip(6,1)
                        kpip(1) = kpip(i1)
                        do i2=1,6,1
                            pip(i2,1) = pip(i2,i1)
                        end do
                        kpip(i1)  = kpip1
                        pip(1,i1) = pip1
                        pip(2,i1) = pip2
                        pip(3,i1) = pip3
                        pip(4,i1) = pip4
                        pip(5,i1) = pip5
                        pip(6,i1) = pip6
                    end if
                end do
                do i3=2,nsca,1
                    if(kpip(i3) == kf4)then
!       moves particle i3 to the second position in 'nsca'
                        kpip1   = kpip(2)
                        pip1    = pip(1,2)
                        pip2    = pip(2,2)
                        pip3    = pip(3,2)
                        pip4    = pip(4,2)
                        pip5    = pip(5,2)
                        pip6    = pip(6,2)
                        kpip(2) = kpip(i3)
                        do i2=1,6,1
                            pip(i2,2) = pip(i2,i3)
                        end do
                        kpip(i3)  = kpip1
                        pip(1,i3) = pip1
                        pip(2,i3) = pip2
                        pip(3,i3) = pip3
                        pip(4,i3) = pip4
                        pip(5,i3) = pip5
                        pip(6,i3) = pip6
                    end if
                end do
            end if
            do i=1,3,1
                ! three position of ik1
                px(i) = V(ik1,i)
                ! three position of ik2
                py(i) = V(ik2,i)
            end do
!       induced parton is distributed randomly in between colliding pair
            n0 = N
            do i=3,nsca,1
                rl = PYR(1)
                n0 = n0 + 1
                do k1=1,3
                    V(n0,k1) = px(k1) * rl + py(k1) * ( 1D0 - rl )
                end do
                V(n0,4)  = tcp
                taup(n0) = 0D0
                ishp(n0) = 1
            end do

!--------- keeps four momentum conservation ----------------------
            do i2=1,4
                ! pi_o (pj_o) four momentum of colliding pair (ik1 & ik2)
                ps(i2) = pi_o(i2) + pj_o(i2)
            end do
            do i1=1,nsca,1
                do i2=1,4,1
                    pss(i1,i2) = pip(i2,i1)
                enddo
                pss(i1,5) = P(i1,5)
            end do
!       adjust four momentum conservation by iteration,no more than
!        5000 iterations
            call cconse( ps, 1, nsca, 1 )
!       total momentum of partons after scattering (including induced partons
!        in time-like branche) is normalized to the total momentum of
!        scattering pair
            do i1=1,nsca,1
                do i2=1,4,1
                    pip(i2,i1) = pss(i1,i2)
                end do
            end do
!---------- keep four momentum conservation finished -------------

!       a part (first & second) of 'nsca' to 'PYJETS', i.e.
!        feedback first & second partons in 'nsca' to ik1 & ik2
!       updates 'PYJETS'
            l = ik1
            if( pip(4,1) < 1D-10 ) pip(4,1) = 1D-10
            do k1=1,3,1
                P(l,k1) = pip(k1,1)
            end do
            p(l,4)  = pip(4,1)
            k(l,2)  = kpip(1)
            p(l,5)  = amass( kpip(1) )
            taup(l) = 0D0
!       updates 'PYJETS'
            l = ik2
            if( pip(4,2) < 1D-10 ) pip(4,2) = 1D-10
            do k1=1,3,1
                P(l,k1) = pip(k1,2)
            end do
            P(l,4)  = pip(4,2)
            K(l,2)  = kpip(2)
            P(l,5)  = amass( kpip(2) )
            taup(l) = 0D0
!0800   finished--------------------------------------------------------

!       deals with third & fourth partons in 'nsca'
            if(nsca == 3)then
                if( i_had_model == 0 .OR. i_had_model == 3 )then
!       moves third induced partons to 'trs' for sfm
                    ntrs = ntrs + 1
                    ktrs(ntrs,1) = 3
                    ktrs(ntrs,2) = kpip(3)
                    do j=1,4,1
                        ptrs(ntrs,j) = pip(j,3)
                        vtrs(ntrs,j) = V(n0,j)
                    end do
                    ptrs(ntrs,5) = 0D0
                    vtrs(ntrs,5) = 0D0
                elseif( i_had_model /= 0 .AND. i_had_model /= 3 )then
!       Adds this parton (gluon) into collision list for coal.
                    N = N + 1
                    K(N,1) = 2
                    K(N,2) = kpip(3)
                    K(N,3) = 0
                    K(N,4) = 0
                    K(N,5) = 0
                    do j=1,4,1
                        P(N,j) = pip(j,3)
                    end do
                    P(N,5) = 0D0
                end if
                goto 1008
            end if
            if(nsca == 4)then
!       composes third & fourth partons in 'nsca' as a string & moves it
!        to the end of 'PYJETS'
!       updates 'PYJETS
                N = N + 1
                l = N
                if( pip(4,3) < 1D-10 ) pip(4,3) = 1D-10
                do k1=1,3
                    P(l,k1) = pip(k1,3)
                end do
                p(l,4) = pip(4,3)
                k(l,1) = 2
                k(l,2) = kpip(3)
                p(l,5) = amass( kpip(3) )
!       updates 'PYJETS'
                N = N + 1
                l = N
                if( pip(4,4) < 1D-10 ) pip(4,4) = 1D-10
                do k1=1,3,1
                    P(l,k1) = pip(k1,4)
                end do
                P(l,4) = pip(4,4)
                K(l,1) = 1
                K(l,2) = kpip(4)
                P(l,5) = amass( kpip(4) )
                goto 1008
            end if
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ti_li1( stili, ij )
!!      Performs time-like branching (along main chain until tl0=Mn0**2
!!       for branching parton).
!       stili: statistics the times calling 'ti_li1'
!       ij: order # of spliting parton in parton list
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        PARAMETER (pio=3.141592653589793D0)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension pa1(4),pa2(4),p00(4)


        il = 1
!       il=1 means the branching will go on, il=0 means 'stop'.
        do i=1,4
            pa1(i) = 0D0
            pa2(i) = 0D0
            p00(i) = 0D0
        enddo
        sampq = 0D0
10      continue
        sampq = sampq + 1D0

!       varieties in time like branchs are all the squared mass (m**2) and
!        m**2 = E**2 - |p|**2 > 0
!------------------------------------------------------------------
!       branching:   A  -> A1   +   A2
!       Vituality:  qqb    qqa     qqc
!       ID No.   :  kfk    kf1     kf2
!        x       :  xa x'  xa1 x   xa2
!       momentum :  pa     pa1     pa2
!
!       equation : qqb = qqa/z + qqc/(1-z) + pt2 / ( z*(1 - z) ) (*)
!-------------------------------------------------------------------
!       above Eq. is deduced from relation of mass and energy (Einstein
!        relation) and four momentum conservation, cf. note made in
!        Webber's paper

!-------flavor of A1,A2 will be decided in the ensueing module.
        q2max = qa
        kf0   = kfk
        call suda_ti( ilo, kf0, ii, q2max, q2, z )
!       time_like Sudakov factor, input: kf0 & q2max, output: others
!       ilo: =0 fail in sampling q2 and z
!            =1 success in sampling q2 and z
        if( ilo == 0 ) then
            ! forced
            qa = tl0
            return
        end if

        qqb = q2
        if( stili > 2D0 ) qqb = qa
        zab = z
!       'zab' is the sampled ratio of x/x' in this branching.
        if( ABS(zab) <= 1D10 ) return
        xa1 = xa*zab
        xa2 = xa - xa1

        kf1 = 0
        kf2 = 0
!       initial state is quark
        if( kfk /= 21 ) then
!       branch q->qg
            kf1 = kfk
            kf2 = 21
            if( PYR(1) >= 0.5D0 ) then
!       branch q->gq
                kf1 = 21
                kf2 = kfk
            endif
!       initial state is gluon
        else
            if(ii == 2) then
!       branch: g->gg
                kf1 = 21
                kf2 = 21
            else
!       considering g->gg only
                goto 200
!       branch: g->qq(-)
!       ea = pip(4,ij): energy of spliting particle ij
!       pip(1-3,ij): three momentum of spliting particle ij
                eaa = ea*ea - pip(1,ij)*pip(1,ij) - pip(2,ij)*pip(2,ij) &
                    - pip(3,ij)*pip(3,ij)
                ! invariant mass of initial state gluon
                eaa = sqrt(eaa)

                if( eaa < ( 2D0*amass(2) ) ) goto 200

                call break_f(eaa,kff,amq)
                kf1 =  kff
                kf2 = -kf1
                if( PYR(1) >= 0.5D0 ) then
                    kf1 =  kf2
                    kf2 = -kf2
                end if
200         end if
        end if
!----------finish

        q2max = zab*( qqb - tl0 / (1D0 - zab) )
!       the max. value of A1's virtuality decided by the eqution(*) with
!        qqc=tl0 & pt2=0 .
        call suda_ti( ilo, kf1, ii, q2max, q2, z )
        if( ilo == 0 ) then
!       qqa=tl0 means the forward evolution has finished, thus set il=0.
            qqa = tl0
            il  = 0
        else
            qqa = q2
        end if
        q2max = ( qqb - qqa/zab ) * ( 1D0 - zab )
!       the max. value of A2's virtuality decided by the eqution(*) with pt2=0
        call suda_ti( ilo, kf2, ii, q2max, q2, z )
        if( ilo == 0 ) then
            qqc = tl0
        else
            qqc = q2
        endif
        pt2 = ( qqb - qqa/zab - qqc/(1D0 - zab) ) * ( 1D0 - zab ) * zab
!       pt2: the square of transverse momentum decided by the equation(*)

!0130   the following block gives the momentum of A1, A2 i.e. pa1, pa2
!        according to the definition of z & pt.
        pa02 = pa(1)**2 + pa(2)**2 + pa(3)**2
        paz  = dsqrt(pa02)   ! third momentum of A
!       note Eq. (*) with assumption that z axis is set on momentum of A,
!        i. e. z axis of the frame, in branching, is set on the pa
!        direction
        ! third momentum of A1
        pz1 = zab*paz
        ! third momentum of A2
        pz2 = ( 1D0 - zab ) * paz
        ! momentum modulus of A1
        ppa1 = dsqrt( pz1**2 + pt2 )
        ! momentum modulus of A2
        ppa2 = dsqrt( pz2**2 + pt2 )
        if( ppa1 < 1D-10 ) ppa1 = 1D-10
        if( ppa2 < 1D-10 ) ppa2 = 1D-10
        ppt = dsqrt(pt2)
        sctas1 = ppt / ppa1
        cctas1 = pz1 / ppa1
        fi = 2D0*pio*PYR(1)
        cfis = dcos(fi)
        sfis = dsin(fi)
!       direction cosines of A1 relative to momentum of A
!       note : transverse momentum of A1 & A2 must be in same plan in
!        order to keep transverse momentum conservation
        sctas2 = -ppt/ppa2
        cctas2 =  pz2/ppa2
        do i=1,3,1
            pa1(i) = pa(i)
            pa2(i) = pa(i)
        end do
        call rotate( cctas1, sctas1, cfis, sfis, ppa1, pa1, p00 )
!       originally, momentum of A1 is relative to the
!        direction of momentum of A
!       after rotation, momentum of A1 is relative to the cms where A
!        belongs to
        call rotate( cctas2, sctas2, cfis, sfis, ppa2, pa2, p00 )
!0130   finished---------------------------------------

!0140   new induced parton (A2) should be added to particle list.
        nsca = nsca + 1
        kpip(nsca) = kf2
        pip(4,nsca) = amass(kf2)**2
        do i=1,3
            pip(4,nsca) = pip(4,nsca) + pa2(i)**2
            pip(i,nsca) = pa2(i)
        end do
        pip(4,nsca) = dsqrt( pip(4,nsca) )
!       pip(5,nsca) = qa2
        pip(5,nsca) = qqc
        pip(6,nsca) = xa2
!0140   finished--------------------------------------------

!------- update five very important values of KFK,QA,XA,PA,EA (time like)
!        & then forward evolution along main chain further
        kfk = kf1
        qa  = qqa
        do i=1,3
            pa(i) = pa1(i)
        enddo
        ea = dsqrt( pa(1)**2 + pa(2)**2 + pa(3)**2 + amass(kfk)**2 )
        xa = xa1

!       if( sampq > 10D0 ) return
        ! once branching only
        if( sampq >= 1D0 ) return
        if( il == 1 ) goto 10
!       il=0 means the evolution 'stop' while il=1 'continue'
!       the time like branching of A2 parton will consider later, i.e.
!        in 'collis' after statement numbered 300


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine suda_ti( ilo, kf0, ii, q2max, q2, z )
!!      Time-like Sudakov factor.
!!      Performs time-like forward branching, one step only.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        PARAMETER (pio=3.141592653589793D0)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xa,sl0,tl0,qa, &
        ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/sa24/adj1(40),nnstop,non24,zstop
!----------------------------------------------------------------------
!       process    : A  ->   A1  +  A2
!       virtuality : qqb     qqa    qqc
!
!       equation: qqb = qqa / z + qqc / (1 - z) + pt2 / (1 - z) / z
!----------------------------------------------------------------------


        ! q->qg cf. Sa' note p. 44
        yint1(x) = -dlog( 1D0 - x )
        ! g->gg
        yint2(x) =  dlog( 1D0 / (1D0 - x) - 1D0 )
!       yint1, 2: the integral functions of part of the splitting functions
        zsp = 0D0
        adj25 = adj1(25)
!       adj25: lamda_QCD
        adj25 = adj25*adj25
        ilo = 0
        ! cf. W' paper p.264 Eq.25
        if( q2max <= 4*tl0 ) return
        ikf = iabs(kf0)
        ii  = 0
        tmax = dlog( q2max / adj25 )
        ! cf. W' paper p.264
        tmin = dlog( 4.*tl0 / adj25 )
!       tmax, tmin: bounds of t = dlog( Q**2/lamda ) in the Sudakov form factor
        ! cf. W' paper p.263
        zmin = tl0/q2max
        ! cf. W' paper p.263
        zmax = 1D0 - tl0/q2max
!       approximated values of the bounds of z
        if(zmax <= zmin) return
        ! q->qg
        if(ikf /= 21) then
            ! cf. Sa' note p. 43
            ciqg = 4D0/3D0 * ( ( 2D0*yint1(zmax) - zmax - 0.5D0*zmax**2 ) &
                             - ( 2D0*yint1(zmin) - zmin - 0.5D0*zmin**2 ) )
!       ciqg: the integration of splitting function for q->qg
            cisum = ciqg
            ii = 1
        ! g->
        else
            cigg = ( 6D0*yint2(zmax) - 12D0*zmax + 3D0*zmax**2 - 2D0*zmax**3 ) &
                 - ( 6D0*yint2(zmin) - 12D0*zmin + 3D0*zmin**2 - 2D0*zmin**3 )
!       cigg: the integration of splitting function for g->gg
            ciqq = 1D0/2D0 * ( ( zmax - zmax**2 + 2D0/3D0*zmax**3 ) &
                             - ( zmin - zmin**2 + 2D0/3D0*zmin**3 ) )
!       ciqq: the integration of splitting function for g->qqbar
            cisum = cigg + ciqq
            aaa = PYR(1)
            ! g->gg
            if( aaa <= (cigg/cisum) ) then
                ii = 2
                cisum = cigg
            ! g->qqbar
            else
                ii = 3
                cisum = ciqq
            end if
        end if
        ce = 9D0/2D0/cisum
        sampz = 0D0
100     continue
        aaa = PYR(1)
        ! variable alpha_s
!       tt = tmax * aaa**ce
        ! constant alpha_s
        tt = tmax + 2D0*pio / as / cisum * LOG(aaa)
!       tt: the t value sampled from time-like Sudakov factor
        if( tt <= tmin ) return
        q2 = adj25*dexp(tt)
!       adj25: Lambda_s**2
        rzmax = 0.5d0 + 0.5D0 * SQRT( 1D0 - 4D0*tl0/q2 )
!       cf. Eq.24 in B. R. Webber, Ann. Rev. Nucl. Part. Sci. 36(1986)253
!        that journal is simplified as W' elsewhere
        rzmin = 1D0 - rzmax
!       rzmax, rzmin: exact value of zmax & zmin corresponding to the t
!        sampled.
200     continue
!0170   sample z.
        aaa = PYR(1)
        bbb = PYR(1)
        ! q->qg
        if(ii == 1) then
            zsp = aaa*yint1(zmax) + (1D0 - aaa)*yint1(zmin)
            zsp = 1D0 - dexp(-zsp)
            ! accepting probability
            ratio = (1D0 + zsp*zsp) / 2D0
!           if(ratio <= bbb) goto 200
            ! the same as previous statement
            if( bbb > ratio ) goto 200
        end if
        ! g->gg
        if(ii == 2) then
            zsp = aaa*yint2(zmax) + (1. - aaa)*yint2(zmin)
            zsp = 1D0 - 1D0 / (1D0 + EXP(zsp) )
            ratio = ( 1D0 - zsp*(1D0 - zsp) )**2
!           if(ratio <= bbb) goto 200
            ! the same as previous statement
            if(bbb > ratio) goto 200
        end if
        ! g->qqbar
        if(ii == 3) then
            zsp = aaa*zmax + (1. - aaa)*zmin
            ratio = 1D0 - 2D0*zsp + 2D0*zsp**2
!           if(ratio <= bbb) goto 200
            ! the same as previous statement
            if(bbb > ratio) goto 200
        end if
!0170   sample z finished--------------------------------------------------
        if( zsp > rzmax .or. zsp < rzmin ) then
!       if the z sampled fall out of the exact region, reject the sampled t,
!       let tmax equal to the t sampled, and go back to resample.
            sampz = sampz + 1D0
!       forced stopping time-like branching 211004
            if( sampz > 1000D0 )then
                ilo = 0
                return
            end if
            tmax = tt
            goto 100
        end if
        z = zsp
        ilo = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cconse(ps,npl,np,nstep)
!!      Adjusts four momentum conservation by iteration,no more than
!!       5000 iterations
!       pp : four momentum of particles
!       ps : above four momenta should be conserved to ps
!       npl : line # of the first particle
!       np : line # of last particle
!       nstep : interval of the step
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/e_cons1/ pp(kszj,5), ff(kszj)
        dimension ps(4),pxyz(3),arp(3)


        dep  = 1D-5
        ps4  = ps(4)
        pxyz = 0D0
        jj = 0
100     es = 0D0
        do i=npl,np,nstep
            es = es + pp(i,4)
        end do
        fr = es / ps4
!       f(dabs(1.-fr)  <=  dep) return
        do i=npl,np,nstep
            amas  = pp(i,5)
            amas2 = amas*amas
            ppm   = pp(i,4)
            ppf   = ppm / fr
            ff(i) = SQRT( ABS( ppf*ppf - amas2 ) / ( ppm*ppm - amas2 ) )
            do j=1,3,1
                ppp = ff(i) * pp(i,j)
                pp(i,j) = ppp
                pxyz(j) = pxyz(j) + ppp
            end do
        end do
        do i=1,3,1
            arp(i) = ABS( 1D0 - pxyz(i) / ps(i) )
        end do
        if( ABS(1D0 - fr) <= dep .AND. arp(1) <= dep .AND. arp(2) <= dep &
            .AND. arp(3) <= dep ) return
        do i=1,3,1
            pxyz(i) = pxyz(i) - ps(i)
            pxyz(i) = pxyz(i) / ( DBLE(np-npl) / DBLE(nstep) + 1D0 )
        end do
        do i=npl,np,nstep
            do j=1,3,1
                pp(i,j) = pp(i,j) - pxyz(j)
            end do
            pp5  = pp(i,5)
            pp52 = pp5 * pp5
            pp(i,4) = SQRT( pp52 + pp(i,1)**2 + pp(i,2)**2 + pp(i,3)**2 )
        end do
        jj = jj + 1
        if(jj < 5000) goto 100


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_qq(zz)
!!      Splitting function for g-->q*q(-)
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        fap_qq = ( zz**2 + (1D0 - zz)**2 ) / 2D0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_gg(zz)
!!      Splitting function for g-->g*g.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        fap_gg = 6D0 * ( 1D0 - zz * ( 1D0 - zz ) )**2 / ( zz * ( 1D0 - zz ) )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fap_qg(zz)
!!      Splitting function for q-->q*g.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        fap_qg = 4D0 * ( 1D0 + zz**2 ) / 3D0 / ( 1D0 - zz )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_ud(xa)
!!      The struction function of valence u, d quark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        external fun_beta


        ata1   = 0.419D0 + 0.004D0 * dta - 0.007D0 * dta * dta
        ata2   = 3.46D0 + 0.724D0 * dta - 0.066D0 * dta * dta
        rud    = 4.4D0 - 4.86D0 * dta + 1.33D0 * dta * dta
        anud   = 3D0 / ( fun_beta( ata1, ata2 + 1D0 ) &
               * ( 1D0 + rud * ata1 / ( ata1 + ata2 + 1D0 ) ) )
        aaa    = anud * xa**ata1 * ( 1D0 - xa )**ata2 * ( 1D0 + rud * xa ) / xa
        fun_ud = aaa


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_d(xa)
!!      The struction function of d quark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        external fun_beta


        ata3  = 0.763D0 - 0.237D0 * dta + 0.026D0 * dta * dta
        ata4  = 4D0 + 0.627D0 * dta - 0.019D0 * dta * dta
        rd    = -0.421D0 * dta + 0.033D0 * dta * dta
        andd  = 1D0 / ( fun_beta( ata3, ata4 + 1D0 ) &
              * ( 1D0 + rd * ata3 / ( ata3 + ata4 + 1D0 ) ) )
        aaa   = andd * xa**ata3 * ( 1D0 - xa )**ata4 * ( 1D0 + rd * xa ) / xa
        fun_d = aaa


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_u(xa)
!!      The struction function of u quark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        fun_u = fun_ud(xa) - fun_d(xa)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_g(xa)
!!      The struction function of gluon.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)


        aag   = 1.56D0 - 1.71D0 * dta + 0.638D0 * dta * dta
        ag    = -0.949D0 * dta + 0.325D0 * dta * dta
        bg    = 6D0 + 1.44D0 * dta - 1.05D0 * dta * dta
        a1g   = 9D0 - 7.19D0 * dta + 0.255D0 * dta * dta
        a2g   = -16.5D0 * dta + 10.9D0 * dta * dta
        a3g   = 15.3D0 * dta - 10.1D0 * dta * dta
        fun_g = aag * xa**ag * ( 1D0 - xa )**bg &
              * ( 1D0 + a1g*xa + a2g*xa*xa + a3g*xa*xa*xa ) / xa


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_s(xa)
!!      The struction function of sea quark.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xt,sm,as,dta,xxa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)


        aas   = 1.265D0 - 1.132D0 * dta + 0.293D0 * dta * dta
        ass   = -0.372D0 * dta - 0.029D0 * dta * dta
        bs    = 8.05D0 + 1.59D0 * dta - 0.153D0 * dta * dta
        a1s   = 6.31D0 * dta - 0.273D0 * dta * dta
        a2s   = -10.5D0 * dta - 3.17D0 * dta * dta
        a3s   = 14.7D0 * dta + 9.8D0 * dta * dta
        fun_s = aas * xa**ass * ( 1D0 - xa )**bs &
              * ( 1D0 + a1s*xa + a2s*xa*xa + a3s*xa*xa*xa ) / xa / 6D0


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_gal(xa)
!!      Logarithm Gammma function.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension coef(0:6)
        data coef/76.18009173,-86.50532033,24.01409822,-1.231739516, &
         0.00120858003,-5.36382e-6,2.50662827565/


        xx   = xa - 1D0
        temp = xx + 5.5D0
        temp = ( xx + 0.5D0 ) * LOG(temp) - temp
        ser  = 1D0
        do i=0,5,1
            xx  = xx + 1D0
            ser = ser + coef(i) / xx
        end do
        co = coef(6) * ser
        fun_gal = temp + LOG(co)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fun_beta(xa,ya)
!!      Beta function.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)


        b = fun_gal(xa) + fun_gal(ya) - fun_gal(xa+ya)
        fun_beta=dexp(b)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs12_0(xt)
!!      Differential cross section of q1*q2 -> q1*q2, process 1.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11 = adj1(1)
        adj20 = adj1(20)
!       xt = xt - tcut
!       xu = sm - xt - xs - 2D0 * tcut
        xu = sm - (xt - tcut) - xs
!       xs: s
!       sm: ma**2 + mb**2 + mc**2 + md**2
!       xu: u
!       xt: t
!       ccc = LOG( xt * xu / xs / adj1(25)*adj1(25) )
!       ccc = ccc**2
!       fsabcd = pio * as**2 / xs**2 / ccc
        pioa = pio * as**2
        if( INT(adj20) == 1 .or.INT(adj20) == 3 ) goto 100
        fsabcd = pioa / xs**2
        fs12_0 = fsabcd * (xs**2 + xu**2)*4D0 / xt**2 / 9D0 * adj11
!       LO pQCD

        return
100     fs12_0 = pioa * 8D0 / 9D0 / xt**2 * adj11
!       limited and regularized (B. Zhang)


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_0(xt)
!!      Differential cross section of q1*q1 -> q1*q1, process 2.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11=adj1(1)
        adj20=adj1(20)
!210803
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fs11_0=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xs**2+xt**2)/xu**2)/9.- &
         8.*xs**2/27./xu/xt)*adj11   ! 210803
!240803 LO pQCD
        return

100     fs11_0=pioa*8./9./xt**2*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs12_1(xt)
!!      Differential cross section of q1*q2(-) -> q1*q2(-), process 3.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11=adj1(1)
        adj20=adj1(20)
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fs12_1=fsabcd*(4.*(xs**2+xu**2)/xt**2/9.)*adj11   ! 210803
!240803 LO pQCD
        return

100     fs12_1=pioa*8./9./xt**2*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_1(xt)
!!      Differential cross section of q1*q1(-) -> q2*q2(-), process 4.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa33/smadel,ecce,secce,parecc,iparres


        if( iparres == 0 .OR. iparres == 2 &
            .OR. ( (iparres == 1 .OR. iparres == 3) &
            .AND. i_inel_proc == 7 ) )then
            fs11_1 = 0D0   ! 160110
            return   ! 160110
        endif
!240412
        adj11=adj1(1)
        adj20=adj1(20)
!210803
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fs11_1=fsabcd*(4.*(xt**2+xu**2)/xs**2/9.)*adj11   ! 210803
!240803 LO pQCD
        return

100     fs11_1=0.   ! pioa*4./9.*(1.+xt*xt/xs*xs)*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fs11_2(xt)
!!      Differential cross section of q1*q1(-) -> q1*q1(-), process 5.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11=adj1(1)
        adj20=adj1(20)
!210803
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fs11_2=fsabcd*(4.*((xs**2+xu**2)/xt**2+(xu**2+xt**2)/xs**2)/9.- &
         8.*xu**2/27./xs/xt)*adj11   ! 210803
!240803 LO pQCD
        return

100     fs11_2=pioa*8./9./xt**2*adj11
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsqq(xt)
!!      Differential cross section of q*q(-) -> g*g, process 6.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/i_inel_proc,i_time_shower,para1_1,para1_2
        common/sa33/smadel,ecce,secce,parecc,iparres


        if( iparres == 0 .OR. iparres == 2 &
            .OR. ( (iparres == 1 .OR. iparres == 3) &
            .AND. i_inel_proc == 7 ) )then
            fsqq = 0D0   ! 160110
            return   ! 160110
        end if
!240412
        adj11=adj1(1)
        adj20=adj1(20)
!210803
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fsqq=fsabcd*(32.*(xu**2+xt**2)/xu/xt/27.- &
         8.*(xu**2+xt**2)/xs**2/3.)*adj11   ! 210803
!240803 LO pQCD
        return

100     fsqq=-pioa*32./27./xs/xt*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsgg_1(xt)
!!      Differential cross section of g*g ->q*q(-), process 7.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa33/smadel,ecce,secce,parecc,iparres


        if( iparres == 0 .OR. iparres == 2 )then
            fsgg_1 = 0D0   ! 160110
            return   ! 160110
        end if
!240412
        adj11=adj1(1)
        adj20=adj1(20)
!210803
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fsgg_1=fsabcd*((xu**2+xt**2)/xu/xt/6.-3.*(xu**2+xt**2)/xs**2/8.) &
         *adj11   ! 210803
!240803 LO pQCD
        return

100     fsgg_1=-pioa*1./3./xs/xt*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsqg(xt)
!!      Differential cross section of q*g -> q*g, process 8.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11=adj1(1)
        adj20=adj1(20)
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fsqg=fsabcd*((xu**2+xs**2)/xt**2-4.*(xu**2+xs**2)/xu/xs/9.) &
         *adj11   ! 210803
!240803 LO pQCD
        return

100     fsqg=pioa*2./xt**2*adj11   ! 240803 200407
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function fsgg_0(xt)
!!      Differential cross section of g*g -> g*g, process 9.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (MSCA=20000)
        common/papr_p/core,xs,xu,xxt,sm,as,dta,xa,sl0,tl0,qa, &
         ea,sqq,sqg,sgg,pa(3),pip(6,msca),mtime,kfk,nsca,kpip(msca)
        common/syspar_p/rsig1,pio,tcut
        common/sa24/adj1(40),nnstop,non24,zstop


        adj11=adj1(1)
        adj20=adj1(20)
!020903 xt=xt-tcut
!071018 xu=sm-xt-xs-2.*tcut   ! 200407
        xu=sm-(xt-tcut)-xs   ! 071018
!120699 ccc=dlog(xt*xu/xs/adj1(25)/adj1(25))
!120699 ccc=ccc**2
!120699 fsabcd=pio*as**2/xs**2/ccc
        pioa=pio*as**2   ! 240803
        if(INT(adj20) == 1 .or. INT(adj20) == 3) goto 100
        fsabcd=pioa/xs**2   ! 120699
        fsgg_0=fsabcd*(9.*(3.-xu*xt/xs**2-xu*xs/xt**2-xs*xt/xu**2)/2.) &
         *adj11   ! 210803
!240803 LO pQCD
        return

100     fsgg_0=pioa*9./2./xt**2*adj11   ! 240803
!       limited and regularized


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function amass(kf)
!!      Mass of parton.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
        common/sa24/adj1(40),nnstop,non24,zstop


!       Hadronization model.
        i_had_model = INT( adj1(12) )
        i_mass = 0
!       For the debug, uses the current algebra mass.
        if( i_had_model == 1 .AND. i_deex > 99 ) i_mass = 3
        amass = 0D0
        r     = 0D0
        k1    = ABS(kf)
        select case( k1 )
        case( 21 )
        case( 2 )
            r = 0.33D0
            if( i_mass == 3 ) r = 0.0056D0
        case( 1 )
            r = 0.33D0
            if( i_mass == 3 ) r = 0.0099D0
        case( 3 )
            r = 0.5D0
            if( i_mass == 3 ) r = 0.199D0
        case( 4 )
            r = 1.5D0
            if( i_mass == 3 ) r = 1.23D0
        case( 5 )
            r = 4.8D0
            if( i_mass == 3 ) r = 4.17D0
        case( 6 )
            r = 175D0
            if( i_mass == 3 ) r = 165D0
        end select
        amass = r


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine integ(fun,ux,dx,idw,eee,sum)
!       Calculates distribution (eee) and integration (sum) of the 'fun',which
!        is a y=t-tcut distribution
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension eee(0:idw)
        external fun


        dw  = idw
        del = ( ux - dx ) / dw
!       calculate distribution function of external function 'fun'
        sum = 0D0
        do i=0,idw,1
            eee(i) = fun( dx + del*i )
            sum = sum + eee(i)
        end do
        sum = sum - 0.5D0 * ( eee(0) + eee(idw) )
        sum = sum * del


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eee_dd(idw,eee,dd)
!!      Calculate the relative integral distribution function
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension dd(0:idw),eee(0:idw)


        dd(0) = 0D0
        do i=1,idw,1
            dd(i) = dd(i-1) + eee(i) + eee(i-1)
!           dd(i) = dd(i-1) + ABS( eee(i) ) + ABS( eee(i-1) )
        end do
        do i=1,idw,1
            dd(i) = dd(i) / dd(idw)
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine samp_integ(idw,ux,dx,xf,dd)
!!      Samples a value from relative integral distribution function 'dd',
!!       which is steaming from a corresponding differential y=t-tcut
!!       distribution function.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        dimension dd(0:idw)


        dw  = idw*1D0
        del = ( ux - dx ) / dw
        ran = PYR(1)
        ir  = 0
        do i=1,idw,1
            if( dd(i) < ran ) cycle
            ir = i
            exit
        end do
        xf = PYR(1)
        xf = dx + del * ( ir - 1 ) + xf * del


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eloss(dell)
!!      Considers parton energy loss phenomenologically.
!       dell: energy loss per unit distance
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop
        common/e_loss1/ rpo(kszj,5),ppo(kszj,5)
        dimension rr(3),pi(3),pia(3),deltp(3)
!       rpo,ppo: parton four coordinate before Newton motion, four momentum
!        before energy loss
!       ppsa(1): cumulate the px losed in an event
!       ppsa(2): cumulate the py losed in an event
!       ppsa(3): cumulate the pz losed in an event
!       ppsa(4): cumulate the e losed in an event


        adj138 = adj1(38)
        no = N
        do i=1,no,1   ! 1
            pei=p(i,4)   ! energy before loss
            ppo(i,4)=pei
            do j=1,3
                rr(j)=v(i,j)-rpo(i,j)
                pi(j)=p(i,j)   ! three momentum before loss
                ppo(i,j)=pi(j)
            enddo
            pt=pi(1)*pi(1)+pi(2)*pi(2)
            if(pt > 1.e20)  pt=1.e20
            if(pt < 1.e-20) pt=1.e-20
            pt=dsqrt(pt)
            if(pt <= adj138) goto 100
!       below ptmin='adj138' Gev/c jet can not loss energy
            rm=rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3)
            if(rm > 1.e20)  rm=1.e20
            if(rm < 1.e-20) rm=1.e-20
            rm=dsqrt(rm)
            delte=rm*dell
            if(delte >= pei) goto 100
!       'delte' should be less than energy of particle
            if(pei < 1.e-20) pei=1.e-20
            srtf=1.-delte/pei
            p(i,4)=pei-delte   ! energy after loss
            do j=1,3
                pia(j)=srtf*pi(j)   ! three momentum after loss
!       massless and collinear approximations
                deltp(j)=pi(j)-pia(j)   ! three momentum loss
                p(i,j)=pia(j)   ! three momentum after loss
                ppsa(j)=ppsa(j)+pi(j)-pia(j)   ! cumulate three momentum loss
            enddo
            ppsa(4)=ppsa(4)+delte   ! cumulate energy loss
!       tansfer losed energy to a gluon (added to the end of particle list)
            call eto1g(i,delte,deltp)
!       tansfer losed energy to two gluons (added to the end of particle list)
!           call eto2g(i,delte,deltp)
100         continue
        enddo   ! 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ctlcre_para(time)
!!      Creates the collision time list after energy loss.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,cspipiKK,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
!       adj136 = adj1(36)
        icol = 1
        ! upper diagonal
        do i=2,N,1
            kfi = ABS( K(i,2) )
!       consider d, u, s, c, b, their antiquarks, and gluon only
            if(kfi > 6 .and. kfi /= 21) cycle
            ! upper diagonal
            do j=1,i-1,1
                kfj = ABS( K(j,2) )
!       consider d,u,s,c,b, their antiquarks, and gluon only
                if(kfj > 6 .and. kfj /= 21) cycle
                call rsfilt_p(j,i,iflag)
                if(iflag == 0) cycle
                tc(1,icol) = 0D0
                tc(2,icol) = 0D0
                call tcolij_par(i,j,time,icol)
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if(tc1 > tc2) tcicol = tc2
                if( tcicol > 0D0 )then
                    tci = tcicol
                    do j1 = 1, icol-1, 1
                        tcj = tc(1,j1)
                        tck = tc(2,j1)
                        if(tcj > tck) tcj = tck
                        if( ABS( tcj - tci ) < dddt )then
                            icol = icol - 1
                            exit
                        end if
                    end do
                    icol = icol + 1
                end if
            end do
        end do
        icol = icol - 1
        n0 = N


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eto1g(ii,delte,deltp)
!!      Transfer the losed energy to a gluon (added to the end of
!!       particle list).
!       ii: the line number (in 'PYJETS') of parton which lossing energy
!       delte: energy losed
!       deltp(3): three momentum losed
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        dimension deltp(3)


        nl  = N + 1
        pp4 = delte
        p(N+1,4) = pp4
        do i=1,3,1
            ppi = deltp(i)
            P(N+1,i) = ppi
        end do
        P(N+1,5) = 0D0
        K(N+1,2) = 21
        V(N+1,4) = V(ii,4)
        taup(N+1) = 0D0
        ishp(N+1) = 1
        do i=1,3,1
            V(N+1,i) = PYR(1) * V(ii,i)
        end do
        N = N + 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine eto2g(ii,delte,deltp)
!!      Transfers the losed energy to two gluons (added to the end of
!!       particle list).
!       ii: the line number (in 'PYJETS') of parton which losed energy
!       delte: energy losed
!       deltp(3): three momentum losed
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        dimension deltp(3)


        nl = N + 1
        nl = N + 2
        pp4 = PYR(1) * delte
        if(pp4 < 1D-20) pp4 = 1D-20
        P(N+1,4) = pp4
        pp42 = delte - pp4
        if(pp42 < 1D-20) pp42 = 1D-20
        P(N+2,4) = pp42
        do i=1,3,1
            ppi = PYR(1) * deltp(i)
            P(N+1,i) = ppi
            ppi2 = deltp(i) - ppi
            P(N+2,i) = ppi2
        end do
        P(N+1,5)  = 0D0
        P(N+2,5)  = 0D0
        K(N+1,2)  = 21
        K(N+2,2)  = 21
        K(N+1,1)  = 1
        K(N+2,1)  = 1
        V(N+1,4)  = V(ii,4)
        V(N+2,4)  = V(ii,4)
        taup(N+1) = 0D0
        taup(N+2) = 0D0
        ishp(N+1) = 1
        ishp(N+2) = 1
        do i=1,3,1
            rpi = PYR(1) * V(ii,i)
            V(N+1,i) = rpi
            V(N+2,i) = V(ii,i) - rpi
        end do
        N = N + 2


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coend(ii)
!!      Copies parton ii to the end of 'PYJETS'.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio


        N = N + 1
        do i3=1,5,1
            K(N,i3) = K(ii,i3)
            P(N,i3) = P(ii,i3)
            V(N,i3) = V(ii,i3)
        end do
        taup(N) = taup(ii)
        ishp(N) = ishp(ii)
!       ndiq(N) = ndiq(ii)
!       do j1=1,idi,1
!           kk = ifcom(j1)
!           jj = npt(j1)
!           if( kk == ii ) ifcom(j1) = N
!           if( jj == ii ) npt(j1)   = N
!       end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjst(ik1,ik1str,ik1sa,ik1sv)
!!      Finds order number etc. of ik1 string.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)


        ik1str = 0   ! if ik1 not belong to string
        ik1sa  = 0   ! if ik1 not belong to string
        ik1sv  = 0   ! if ik1 not belong to string
        if(ik1 <= nbe)then   ! nbe is given before call break
            do i1 = 1,nstr0,1   ! nstr0 is given after call break
                i1a = nstr1a(i1)   ! i1a: line number of 'A' of i1-th string
                i1v = nstr1v(i1)   ! i1v: line number of 'V' of i1-th string
                if(ik1 >= i1a .and. ik1 < i1v+1)then
                    ik1str = i1    ! order number of string to which ik1 belongs
                    ik1sa  = i1a   ! line number of 'A' of above string
                    ik1sv  = i1v   ! line number of 'V' of above string
                    return
                end if
            end do
        end if
        if(ik1 > nbe .and. nstr1 > nstr0)then
            do i1 = nstr0+1, nstr1, 1
                i1a = nstr1a(i1)
                i1v = nstr1v(i1)
                if( ik1 >= i1a .and. ik1 < i1v+1 )then
                    ik1str = i1
                    ik1sa  = i1a
                    ik1sv  = i1v
                    return
                end if
            end do
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine adjdi(ii,idi1,idway)
!!      Is ii component of a diquark.
!       idway=0: if ii is not first component or it's partner of a diquark
!       idway=1: if ii is first component of a diquark
!       idway=2: if ii is partner of a diquark
!       idi1: order # of diquark to which ii belong
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
!       ifcom(i): line number (in 'sbe') of first component of i-th diquark


        idway = 0
        idi1  = 0
        do i1=1,idi,1
            i2 = ifcom(i1)
            if(ii == i2)then
!       ii is first component of idi1-th diquark
                idi1 = i1
                idway = 1
                return
            end if
        end do

        do i1=1,idi,1
            ii1=npt(i1)   ! line number of partner of i1-th diquark
            if(ii == ii1)then   ! ii is partner of i1-th diquark
!       ii is partner of idi1-th diquark
                idi1=i1   ! ii is partner of idi1-th diquark
                idway=2
                return
            end if
        end do


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine diqmov(idi1)
!!      Removes diquark idi1 out of list.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio


        if(idi > 0)then
            do i1=1,nbe
                if(ndiq(i1) == idi1) ndiq(i1) = 0
!       ndiq(j): = 0 if j is quark (antiquark)
!                = order # of diquark if j is diquark (anti-diquark)
!             j: line number in 'sbe'
            end do
            if(idi1 == idi) idi = idi - 1
            if(idi1 < idi)then
!       updates diquark list
                do i1 = idi1+1, idi, 1
                    i2 = i1 - 1
                    ifcom(i2) = ifcom(i1)
!       ifcom(idi): line number of first component of idi-th diquark
                    npt(i2)=npt(i1)
                end do
            idi = idi - 1
            end if
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine parmov(ii,ik1,ik2,lmn)
!!      If 'ii' is not equal to ik1 (ik2) and lmn=4 or 6, moves parton 'ii'
!!       to 'trs' as well as updates lists.
!!      otherwise, throws away parton ii (ii is equal to ik1 (ik2)),
!!       i.e. updates lists
!       updates 'PYJETS' one step downward from ii+1 to n
!       updates 'sbe' one step downward from ii+1 to nbe if ii <= nbe
!       updates diquark list
!       updates string list
!       updates collision time list
!       ik1 & ik2 are the line number (in parton list) of colliding pair
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
        common/trs/ntrs,nontrs,ktrs(kszj,5),ptrs(kszj,5),vtrs(kszj,5)


!       The hadronization model.
        i_had_model = INT( adj1(12) )
!       moves parton ii (not equal to ik1 or ik2) to 'trs'
        if( ( i_had_model == 0 .or. i_had_model == 3 ) .and. &
            (lmn == 4 .or. lmn == 6) .and. &
            (ii /= ik1 .or. ii /= ik2))then
            ntrs = ntrs + 1
            ktrs(ntrs,1) = 3
            ktrs(ntrs,2) = k(ii,2)
            do j=1,4,1
                ptrs(ntrs,j) = P(ii,j)
                vtrs(ntrs,j) = V(ii,j)
            end do
            ptrs(ntrs,5) = P(ii,5)
            vtrs(ntrs,5) = 0D0
        end if

!       throws away parton ii, i.e. updates lists
!       updates diquark list
        do j1=1,idi,1
            jj = ifcom(j1)
            kk = npt(j1)
            if(jj >= ii) ifcom(j1) = jj - 1
            if(kk >= ii) npt(j1)   = kk - 1
        end do
        do i2 = ii+1, N, 1
            i3 = i2 - 1
            ndiq(i3) = ndiq(i2)
        end do

!       updates particle list 'PYJETS'
        if(ii == N) N = N - 1
        if(ii < N)then
            do i1 = ii+1, N, 1
                i2 = i1 - 1
                do i3=1,5,1
                    K(i2,i3) = K(i1,i3)
                    P(i2,i3) = P(i1,i3)
                    V(i2,i3) = V(i1,i3)
                end do
                ishp(i2) = ishp(i1)
                taup(i2) = taup(i1)
            enddo
            N = N - 1
        end if
!       if(ik1 > ii) ik1 = ik1 - 1
!       if(ik2 > ii) ik2 = ik2 - 1

!       updates particle list 'sbe'
        if(ii == nbe) nbe = nbe - 1
        if(ii < nbe)then
            do i1 = ii+1, nbe, 1
                i2 = i1 - 1
                do i3=1,5,1
                    kbe(i2,i3) = kbe(i1,i3)
                    pbe(i2,i3) = pbe(i1,i3)
                    vbe(i2,i3) = vbe(i1,i3)
                end do
            end do
            nbe = nbe-1
        end if

!       updates string list
        do i1=1,nstr1,1
            nstraa = nstr1a(i1)
            if(nstraa >= ii) nstr1a(i1) = nstraa - 1
            nstrvv = nstr1v(i1)
            if(nstrvv >= ii) nstr1v(i1) = nstrvv - 1
        end do

!       updates the values of lc(1-2,m) if which is  >=  ii
        do m=1,icol,1
            lc1 = lc(1,m)
            if(lc1 >= ii) lc(1,m) = lc1 - 1
            lc2 = lc(2,m)
            if(lc2 >= ii) lc(2,m) = lc2 - 1
        end do

!       Removing colliding pairs with lc=0. This case happens when ii
!        is the 1-st one in PYJETS, especially when removing string
!        after inelastic collisions.
        jcol = 0
!       Loops over colliding pairs.
        do i=1,icol,1
            iic = lc(1,i)
            jjc = lc(2,i)
!           Throw away the pairs with lc=0.
            if(iic /= 0 .AND. jjc /= 0)then
                jcol = jcol + 1
                tc(1,jcol) = tc(1,i)
                tc(2,jcol) = tc(2,i)
                lc(1,jcol) = lc(1,i)
                lc(2,jcol) = lc(2,i)
            end if
        end do
        icol = jcol


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine strmov(istr1,istr1a,istr1v,ik1,ik2,lmn)
!!      Removes all conponents of 'istr1'-th string out.
!       istr1: order number of string (in string list) to be moved
!       istr1a: line number (in parton list) of first component of above string
!       istr1v: line number (in parton list) of last component of above string
!       ik1 & ik2 are line # of colliding pair in parton list 'PYJETS'
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        PARAMETER (MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/collist/lc(2,mclis),tc(2,mclis),icol
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0, &
         nstr1,nstr1a(kszj),nstr1v(kszj)
!       moves components of istr1-th string out of parton list


        do ii = istr1v, istr1a, -1
            call parmov(ii,ik1,ik2,lmn)
        end do
!       moves string istr1-th out of string list
        if(istr1 == nstr1) nstr1 = nstr1 - 1
        if(istr1 < nstr1)then
            jj = istr1v - istr1a + 1
            do ii = istr1+1, nstr1, 1
!               if(ii >= istr1)then
!               nstr1v(ii-jj) = nstr1v(ii)
!               nstr1a(ii-jj) = nstr1a(ii)
                nstr1v(ii-1) = nstr1v(ii)
                nstr1a(ii-1) = nstr1a(ii)
!               end if
            end do
            nstr1 = nstr1 - 1
        end if


        return
        end




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine copl_p(time)
!!      Calculates position of center of mass of the non-freeze-out system
!!      distance of a particle from this cms is used to checke whether
!!       it freezes out or not.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NON1,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)


        tt = time
        coor = 0D0
        samass = 0D0
        do ii=1,N,1
!           if( taup(ii) > tt ) cycle
            if( ishp(ii) == 0 ) cycle
            KF = K(ii,2)
            ! Exludes junctions.
            if( KF == 88 ) cycle
            dmass  = amass( KF )
            samass = samass + dmass
            do jj=1,3,1
                coor(jj) = coor(jj) + dmass * V(ii,jj)
            end do
        end do
        ! 0.0099 : d current algebra mass.
        coor = coor / MAX( 0.0099D0, samass )


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine his_p(t1,rnt,rnp,istop)
!!      Classical Newton motion in AA/NN CMS/Lab. system.
!       istop=1 means all particles have get out of considered volume.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/collist/lc(2,mclis),tc(2,mclis),icol


        istop = 0
        i_n   = 0
        r0    = PARAM(10) * MAX( rnt, rnp )
        do i=1,N,1
!       do move particles which have not been produced
!           if( t1 <= taup(i) ) cycle
            if( ishp(i) == 0 .OR.  K(i,2) == 88 )then
                i_n = i_n + 1
                cycle
            end if
            aa  = 0D0
            pp4 = P(i,4)
            do j=1,3,1
                vp = P(i,j) / pp4
                V(i,j) = V(i,j) + vp * ( t1 - V(i,4) )
                aa = aa + ( V(i,j) - coor(j) )**2
            end do
            aa = SQRT(aa)
            if(aa < r0)then
                V(i,4) = t1
            else
!       if freeze-out already, deduct the distance between the last and
!        current collisions
                do j=1,3,1
                    vp = P(i,j) / pp4
                    V(i,j) = V(i,j) - vp * ( t1 - V(i,4) )
                end do
                ishp(i) = 0
                do il=1,icol,1
                    if( lc(1,il) == i ) tc(1,il) = 0D0
                    if( lc(2,il) == i ) tc(2,il) = 0D0
                end do
            end if
        end do
        if( i_n == N ) istop = 1


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine rsfilt_p( i, j, iflag )
!!      Subroutine rsfilt_p plays the role of the filter.
!       Only quaks and gluons.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000)
        COMMON/PYJETS/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa8_p/taup(kszj),coor(3),ishp(kszj)


        iflag = 0
        ki    = ABS( K(i,2) )
        kj    = ABS( K(j,2) )
        if( i /= j .AND. ishp(i) == 1 .AND. ishp(j) == 1 )then
            if( ( ( ki >= 1 .AND. ki <= 6 ) .OR. ki  == 21 ) .AND. &
                ( ( kj >= 1 .AND. kj <= 6 ) .OR. kj == 21 ) ) &
            iflag = 1
        end if


        return
        end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpli_p(n00,time)
!!      Updates collision time list for new partons when removing string
!!       and diqaruk after inelastic parton collisions.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER (KSZJ=80000,MCLIS=280000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/papr/t0,siig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm, &
         rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen
        common/scatt/pi(4),pj(4),ic,jc,n0
        common/collist/lc(2,mclis),tc(2,mclis),icol


        dddt = adj1(19)
!       Loops over new ii (> n00) and old jj partons (< n00).
        icol = icol + 1
        do ii = n00+1, N, 1
            i1  = ii
            kfi = ABS( K(ii,2) )
            if( kfi > 6 .and. kfi /= 21 ) cycle
!           Consider d,u,s,c,b, their antiquarks, and gluon.
            do jj=1,n00,1
                j1 = jj
                if(i1 == j1) cycle
                call rsfilt_p(j1,i1,iflag)
                if(iflag == 0) cycle
                tc(1,icol) = 0D0
                tc(2,icol) = 0D0
                call tcolij_par(i1,j1,time,icol)
                tc1 = tc(1,icol)
                tc2 = tc(2,icol)
                tcicol = tc1
                if(tc1 > tc2) tcicol = tc2
                if( tcicol > 0D0 )then
                    tci = tcicol
                    do j2 = 1, icol-1, 1
                        tcj = tc(1,j2)
                        tck = tc(2,j2)
                        if(tcj > tck) tcj = tck
                        if( ABS(tcj-tci) < dddt )then
                            icol = icol - 1
                            exit
                        end if
                    end do
                    icol = icol + 1
                end if
            end do
        end do
        icol = icol - 1
        n0 = N


        return
        end



!ccccccccccccccccccccccccccccccccccccc end ccccccccccccccccccccccccccccccccccccc