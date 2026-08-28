#pragma once

#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>

#define NEWTON_MAXITER 3
#define NEWTON_TOL 1E-9
#define FN_THRESH 0.05

/**
 * YEHU SIMULATION CODE IN CPP
 * REQUIRES EIGEN LIBRARY
 * Champ C. Darabundit and Zhen Zhang, CAML 2026
 */
namespace yehu{
    
    template<class T>
    class yehu{
        public:
        yehu() {}
        ~yehu() = default;

        /**
         * @brief prepare for synthesis after parameters are initialized
         * 
         * @param _SR sample rate in Hz
         */
        void prepareToPlay(T _SR) {
            SR = _SR;
            k = 1.0 / SR;
            ksqr = k * k;

            hmin1 = std::sqrt( 0.5 * (c1*c1*k*k+4*sig11*k+std::sqrt( std::pow(c1*c1*k*k+4*sig11*k, 2)+16*Kp1*Kp1*k*k)) ); 
            N1 = static_cast<int>( std::floor ( L1 / hmin1 ) );
            h1 = L1 / N1;

            hmin2 = std::sqrt( 0.5 * (c2*c2*k*k+4*sig11*k+std::sqrt( std::pow(c2*c2*k*k+4*sig12*k, 2)+16*Kp2*Kp2*k*k)) ); 
            N2 = static_cast<int>( std::floor ( L2 / hmin2 ) );
            h2 = L2 / N2;

            // Call discretization to compute matrices
            this->discretize();

            // Allocate states and reset the DSP
            uprev1.resize(N1-1);
            u1.resize(N1-1);
            unext1.resize(N1-1);

            uprev2.resize(N2-1);
            u2.resize(N2-1);
            unext2.resize(N2-1);

            q.resize(Nm);
            qprev.resize(Nm);
            qnext.resize(Nm);

            this->resetDSP();
        }
        
        /**
         * @brief Discretize our system and calculate computation matrices
         * 
         */
        void discretize() {
            // *****************
            // STRINGS
            // *****************

            // STRING 1
            Eigen::SparseMatrix<T> Dxx = this->makeDxx(N1-1, h1);
            Eigen::SparseMatrix<T> Dxxxx = Dxx * Dxx;
            Eigen::SparseMatrix<T> Id = Eigen::SparseMatrix<T>(N1-1, N1-1);
            Id.setIdentity();
            
            os1 = 1+sig01*k;
            As1 = (2*Id +(2*sig11*k+c1*c1*ksqr)*Dxx - Kp1*Kp1*ksqr*Dxxxx)/os1;
            Bs1 = ((sig01*k-1)*Id -(2*sig11*k)*Dxx)/os1;
            Cs1 = k*k/os1;

            oos1 = 2 * SR + 2*sig01;
            AAs1 = ((c1*c1+2*sig11 * SR)*Dxx-Kp1*Kp1*Dxxxx+2 * SR * SR*Id)/oos1;
            BBs1 = ((-2*sig11 * SR)*Dxx-2 * SR * SR*Id)/oos1;

            // STRING 2
            Dxx = this->makeDxx(N2-1, h2);
            Dxxxx = Dxx * Dxx;
            Id.resize(N2-1, N2-1);
            Id.setIdentity();

            os2 = 1+sig02*k;
            As2 = (2*Id +(2*sig12*k+c2*c2*ksqr)*Dxx - Kp2*Kp2*ksqr*Dxxxx)/os2;
            Bs2 = ((sig02*k-1)*Id -(2*sig12*k)*Dxx)/os2;
            Cs2 = k*k/os2;

            oos2 = 2 * SR + 2*sig02;
            AAs2 = ((c2*c2+2*sig12 * SR)*Dxx-Kp2*Kp2*Dxxxx+2 * SR * SR*Id)/oos2;
            BBs2 = ((-2*sig12 * SR)*Dxx-2 * SR * SR*Id)/oos2;
            
            // *****************
            // BOW HAIR
            // *****************
            oh = mh + kh*k*k/4 + rh*k*0.5;
            Ah =  (2*mh - k*k*kh*0.5)/oh;
            Bh = (-mh - k*k*kh/4 + k*rh*0.5)/oh;
            Ch = k*k/oh;

            ooh = 2 * SR + k*kh*0.5/mh +rh/mh;
            AAh = (2 * SR * SR + kh*0.5/mh - kh/mh)/ooh;
            BBh = (-2 * SR * SR - kh*0.5/mh)/ooh;

            // *****************
            // BODY 
            // *****************
            const Eigen::Matrix<T, Eigen::Dynamic, 1> denom = Mq + (ksqr * 0.25) * Kq + (k * 0.5) * Rq;

            Aq = (2 * Mq - (ksqr * 0.5) * Kq).cwiseQuotient(denom);
            Bq = (-Mq - (ksqr * 0.25) * Kq + (k * 0.5) * Rq).cwiseQuotient(denom);
            Cq = denom.cwiseInverse() * ksqr;


            // *****************
            // INTERPOLATORS
            // *****************
            // Icoup1 = itpl(L1, xcoup, h1, N1-1, 1);
            // Ibow1 = itpl(L1, xb, h1, N1-1, 2);
            Icoup1.resize(N1-1);
            Icoup1.reserve(1);
            Ibow1.resize(N1-1);
            Ibow1.reserve(2);
            If1.resize(N1-1);
            If1.reserve(4);

            itpl(Icoup1, L1, xcoup, h1, 1);
            itpl(Ibow1, L1, xb, h1, 2);
            Jcoup1 = Icoup1.transpose() / h1;
            Jbow1 = Ibow1.transpose() / h1;
            IJbow1 = Ibow1.dot(Jbow1)/rhoA1/oos1 +  1.0/mh/ooh;

            // Icoup2 = itpl(L2, xcoup, h2, N2-1, 1);
            // Ibow2 = itpl(L2, xb, h2, N2 - 1, 2);
            Icoup2.resize(N2-1);
            Icoup2.reserve(1);
            Ibow2.resize(N2-1);
            Ibow2.reserve(2);
            If2.resize(N2-1);
            If2.reserve(4);

            itpl(Icoup2, L2, xcoup, h2, 1);
            itpl(Ibow2, L2, xb, h2, 2);
            Jcoup2 = Icoup2.transpose() / h2;
            Jbow2 = Ibow2.transpose() / h2;
            IJbow2 = Ibow2.dot(Jbow2)/rhoA2/oos2 +  1.0/mh/ooh;

            
            // *****************
            // PRECOMPUTE CONSTANTS
            // *****************
            // COUPLING
            IcoupAs1 = Icoup1 * As1;
            IcoupAs2 = Icoup2 * As2;
            IcoupBs1 = Icoup1 * Bs1;
            IcoupBs2 = Icoup2 * Bs2;

            T sumCq = Cq.sum();
            T D1 = Icoup1.dot(Cs1*Jcoup1)/(rhoA1) + sumCq;
            T D2 = Icoup2.dot(Cs2*Jcoup2)/(rhoA2) + sumCq;
            T coupDenom = 1./(sumCq * sumCq - D1*D2);

            coupCs11 = -D2 * coupDenom; 
            coupCs12 = sumCq * coupDenom;
            coupCs21 = sumCq * coupDenom; 
            coupCs22 = -D1 * coupDenom;

            // BOWING/BRIDGE
            CsJcoup1 = Cs1 * Jcoup1 / rhoA1;
            CsJcoup2 = Cs2 * Jcoup2 / rhoA2;
            CsJbow1 = Cs1 * Jbow1 / rhoA1;
            CsJbow2 = Cs2 * Jbow2 / rhoA2;

            // FINGER
            oneOverMf = 1./mf;
            oneOverRhoAOs1 = 1/(rhoA1 * os1);
            oneOverRhoAOs2 = 1/(rhoA2 * os2);

            // BOWING
            oneOverS0 = 1./ s0;
        }
        /**
         * @brief Reset all system states
         * 
         */
        void resetDSP() {
            // *****************
            // STRING 1 
            // *****************
            std::fill(uprev1.data(), uprev1.data() + uprev1.size(), 0.0);
            std::fill(u1.data(), u1.data() + u1.size(), 0.0);
            std::fill(unext1.data(), unext1.data() + unext1.size(), 0.0);
            v1 = vr1 = fn1 = 0.0;
            haprev1 = ha1 = hanext1 = 0.0;
            z1 = zprevh1 = z_av1 = znexth1 = 0.0;
            chiprev1 = 0.0;
            
            // *****************
            // STRING 2
            // *****************
            std::fill(uprev2.data(), uprev2.data() + uprev2.size(), 0.0);
            std::fill(u2.data(), u2.data() + u2.size(), 0.0);
            std::fill(unext2.data(), unext2.data() + unext2.size(), 0.0);
            v2 = vr2 = fn2 = 0.0;
            haprev2 = ha2 = hanext2 = 0.0;
            z2 = zprevh2 = z_av2 = znexth2 = 0.0;
            chiprev2 = 0.0;

            // *****************
            // FINGER
            // *****************
            ufprev = uf = u0f;
            etaprev1 = 0;
            etaprev2 = 0;

            // *****************
            // BODY
            // *****************
            std::fill(q.data(), q.data() + q.size(), 0.0);
            std::fill(qprev.data(), qprev.data() + qprev.size(), 0.0);
            std::fill(qnext.data(), qnext.data() + qnext.size(), 0.0);

        }

        /**
         * @brief Process a single sample
         * 
         * @param out output variable
         * @param fbri output bridge force
         * @param vb bow velocity
         * @param fn bow force
         * @param ff finger force
         */
        void processSample(T & out, T & fbri, T & vb, T & fn, T & ff) {
            // *****************
            // BRIDGE/BODY
            // *****************

            auto b0vec = Aq.cwiseProduct(q) + Bq.cwiseProduct(qprev);
            T b0 = b0vec.sum();
            T b1 = b0 - (IcoupAs1.dot(u1) + IcoupBs1.dot(uprev1));
            T b2 = b0 - (IcoupAs2.dot(u2) + IcoupBs2.dot(uprev2));
            
            T Fs1 = coupCs11 * b1 + coupCs12 * b2;
            T Fs2 = coupCs21 * b1 + coupCs22 * b2;
            fbri = Fs1 + Fs2;

            qnext = b0vec - Cq*(Fs1+Fs2);

            // *****************
            // BOW
            // *****************
            T sn1 = Ibow1.dot(AAs1*u1 + BBs1*uprev1) + AAh*ha1 + BBh*haprev1;
            T sn2 =  Ibow2.dot(AAs2*u2 + BBs2*uprev2) + AAh*ha2 + BBh*haprev2 ; 

            // APPLY BOW HERE
            Eigen::Matrix<T, Eigen::Dynamic, 1> es1, es2;
            if (fn > FN_THRESH) { // Bow the F string
                T fC = muC*fn;
                T fS = muS*fn;
                T z_ba = 0.7*fC * oneOverS0;      
                T epsilon = fC/s1bar;

                // NR CALL
                this->nrBow(vb, vr1, z_av1, zprevh1, fC, fS, z_ba, epsilon, IJbow1, sn1);
               
                T vrEpsPp = 1 / ( vr1 * vr1 + epsilon * epsilon );
                s1 = fC * std::sqrt(vrEpsPp);
                znexth1 = 2.0 * z_av1 - zprevh1;
                T gn = (znexth1 - zprevh1) * SR;

                T F_fr1 =  s0 * z_av1 + s1 * gn;
                es1 = As1 * u1 + Bs1 * uprev1 - CsJbow1*F_fr1 + CsJcoup1*Fs1;
                hanext1 = Ah*ha1 + Bh*haprev1 - Ch*F_fr1;

                es2 = As2 * u2 + Bs2 * uprev2+ CsJcoup2 * Fs2;
                hanext2 = Ah*ha2 + Bh*haprev2;
            }
            else if (fn < -FN_THRESH){ // Bow the C string
                T fC = -muC * fn;
                T fS = -muS * fn;
                T z_ba = 0.4 * fC * oneOverS0;
                T epsilon = fC / s1bar;

                // NR CALL
                this->nrBow(vb, vr2, z_av2, zprevh2, fC, fS, z_ba, epsilon, IJbow2, sn2);
                
                T vrEpsPp = 1.0 / (vr2 * vr2 + epsilon * epsilon);
                s1 = fC * std::sqrt(vrEpsPp);

                znexth2 = 2.0 * z_av2 - zprevh2;
                T gn = (znexth2 - zprevh2) * SR;

                T F_fr2 = s0 * z_av2 + s1 * gn;

                es2 = As2 * u2 + Bs2 * uprev2 - CsJbow2 * F_fr2 + CsJcoup2 * Fs2;
                hanext2 = Ah * ha2 + Bh * haprev2 - Ch * F_fr2;

                es1 = As1 * u1 + Bs1 * uprev1 + CsJcoup1 * Fs1;
                hanext1 = Ah * ha1 + Bh * haprev1;

            } 
            else { // without bowing
                es1 = As1 * u1 + Bs1 * uprev1 + CsJcoup1 * Fs1;
                hanext1 = Ah * ha1 + Bh * haprev1;
                es2 = As2 * u2 + Bs2 * uprev2 + CsJcoup2 * Fs2;
                hanext2 = Ah * ha2 + Bh * haprev2;
            }

            // *****************
            // FINGER COLLISION
            // *****************
            T f1 = 2.0 * uf - ufprev - ff * ksqr * oneOverMf;
            T eta1 = If1.dot( u1 ) - uf;
            T etast1 = If1.dot( es1 ) - f1;
            
            eta1 = (eta1 < 0) ? 0 : eta1;
            T V1 = kf/(alf + 1) * std::pow(eta1, alf+1);
            T dV1deta1 = kf * std::pow(eta1, alf);

            T Gtilde1 = dV1deta1 / std::sqrt( 2 * V1 + eps );
            T rfree1 = etast1 - ( If1.dot(uprev1) - ufprev );

            T g1 = gradientScale( Gtilde1, Gtilde1_prev, chiprev1, rfree1 );

            T eta2 = If2.dot( u2 ) - uf;
            T etast2 = If2.dot( es2 ) - f1;

            eta2 = (eta2 < 0) ? 0 : eta2;
            T V2 = kf/(alf + 1) * std::pow(eta2, alf + 1);
            T dV2deta2 = kf * std::pow( eta2, alf );

            T Gtilde2 = dV2deta2 / std::sqrt( 2 * V2 + eps );
            T rfree2 = etast2 - ( If2.dot(uprev2) - ufprev );
            T g2 = gradientScale( Gtilde2, Gtilde2_prev, chiprev2, rfree2 );


            // SHERMAN-MORRISON 
            T g1sqr = g1 * g1;
            T g2sqr = g2 * g2;

            T pf1 = -0.25 * ksqr * g1sqr * (If1.dot(uprev1) - ufprev) + chiprev1 * g1 * ksqr
                - g1 * rf * k * 0.5 * (If1.dot(uprev1) - ufprev) * chiprev1;

            T pf2 = -0.25 * ksqr * g2sqr * (If2.dot(uprev2) - ufprev) + chiprev2 * g2 * ksqr
                - g2 * rf * k * 0.5 * (If2.dot(uprev2) - ufprev) * chiprev2;

            Eigen::Matrix<T, Eigen::Dynamic, 1> bf1 = es1 - oneOverRhoAOs1 * Jf1 * pf1;
            Eigen::Matrix<T, Eigen::Dynamic, 1> bf2 = es2 - oneOverRhoAOs2 * Jf2 * pf2;
            T bf3 = f1 + oneOverMf * (pf1 + pf2);

            T cf1 = oneOverRhoAOs1 * (0.25 * ksqr * g1sqr + 0.5 * rf * g1 * k * chiprev1);
            T cf2 = oneOverRhoAOs2 * (0.25 * ksqr * g2sqr + 0.5 * rf * g2 * k * chiprev2);
            T denf1 = 1.0 / (1.0 + cf1 * IJf1);
            T denf2 = 1.0 / (1.0 + cf2 * IJf2);

            Eigen::Matrix<T, Eigen::Dynamic, 1> Af1inv_bf1 = bf1 - (cf1 * denf1) * Jf1 * (If1 * bf1);
            Eigen::Matrix<T, Eigen::Dynamic, 1> Af2 = oneOverRhoAOs1 * (-0.25 * ksqr * g1sqr * Jf1 - 0.5 * rf * g1 * k * chiprev1 * Jf1);
            Eigen::Matrix<T, Eigen::Dynamic, 1> Af1inv_Af2 = Af2 - (cf1 * denf1) * Jf1 * (If1 * Af2);
            Eigen::Matrix<T, Eigen::Dynamic, 1> Af3inv_bf2 = bf2 - (cf2 * denf2) * Jf2 * (If2 * bf2);
            Eigen::Matrix<T, Eigen::Dynamic, 1> Af4 = oneOverRhoAOs2 * (-0.25 * ksqr * g2sqr * Jf2 - 0.5 * rf * g2 * k * chiprev2 * Jf2);
            Eigen::Matrix<T, Eigen::Dynamic, 1> Af3inv_Af4 = Af4 - (cf2 * denf2) * Jf2 * (If2 * Af4);

            Eigen::Matrix<T,1, Eigen::Dynamic> Af5 = oneOverMf * (-0.25 * ksqr * g1sqr * If1 - 0.5 * rf * g1 * k * chiprev1 * If1);
            Eigen::Matrix<T,1, Eigen::Dynamic> Af6 = oneOverMf * (-0.25 * ksqr * g2sqr * If2 - 0.5 * rf * g2 * k * chiprev2 * If2);

            T Af7 = 1.0 + oneOverMf * ( 0.25 * ksqr * g1sqr + 0.5 * rf * g1 * chiprev1 * k
                + 0.25 * ksqr * g2sqr + 0.5 * rf * g2 * chiprev2 * k );

            T S = Af7 - Af5 * Af1inv_Af2 - Af6 * Af3inv_Af4;
            T rhsU = bf3 - Af5 * Af1inv_bf1 - Af6 * Af3inv_bf2;

            T ufnext = rhsU / S;

            // UPDATE VARIABLES
            unext1 = Af1inv_bf1 - Af1inv_Af2 * ufnext;
            unext2 = Af3inv_bf2 - Af3inv_Af4 * ufnext;

            chinext1 = 0.5 * g1 * (If1.dot(unext1 - uprev1) - (ufnext - ufprev)) + chiprev1;
            chinext2 = 0.5 * g2 * (If2.dot(unext2 - uprev2) - (ufnext - ufprev)) + chiprev2;

            // OUTPUT
            out = Icoup1.dot(unext1 - uprev1) * SR * 0.5;

            // STATE UPDATE
            zprevh1 = znexth1;
            uprev1 = u1;
            u1 = unext1;
            haprev1 = ha1;
            ha1 = hanext1;

            ufprev = uf;
            uf = ufnext;

            etaprev1 = eta1;
            chiprev1 = chinext1;

            etaprev2 = eta2;
            chiprev2 = chinext2;

            qprev = q;
            q = qnext;

            Gtilde1_prev = Gtilde1;
            Gtilde2_prev = Gtilde2;

            zprevh2 = znexth2;
            uprev2 = u2;
            u2 = unext2;
            haprev2 = ha2;
            ha2 = hanext2;
        }

        /**
         * @brief process a vector of audio
         * 
         * @param out output vector
         * @param fbri output bridge force
         * @param vb vector of bow velocities
         * @param fn vector of bow forces
         * @param xf vector of finger positions
         * @param ff vector of finger forces
         */
        void process(std::vector<T> & out, std::vector<T> & fbri, std::vector<T> & vb, std::vector<T> & fn, std::vector<T> & xf, std::vector<T> & ff)
        {
            double xf_prev = 0.0;
            int Nf = out.size();
            assert( fbri.size() == Nf && vb.size() == Nf && fn.size() == Nf && xf.size() == Nf && ff.size() == Nf);
            for (int samp = 0; samp < Nf; ++samp){
                this->setFingerPos( xf[samp] );
                this->processSample( out[samp], fbri[samp], vb[samp], fn[samp], ff[samp]);
            }
        }

        /**
         * @brief Newton-Raphson call for bowing
         * 
         * @param vb input bow velocity
         * @param vr difference btw bow and mass velocity 
         * @param z_av average bristle deflection
         * @param zprevh previous bristle deflection
         * @param fC dynamic friction force
         * @param fS static friction force
         * @param z_ba breakaway displacement
         * @param epsilon dynamic friction over s1bar
         * @param IJbow coupling coefficient associated with string
         * @param sn auxiliary variable
         */
        void nrBow(T & vb, T & vr, T & z_av, T & zprevh, T & fC, T & fS, T & z_ba, T & epsilon, T& IJbow, T & sn) {
            // Init the NR loop
            bool nr = true;
            T err = 1.0;
            int iter = 0;

            while (nr) {
                T vrEpsPp = 1.0 / (vr*vr + epsilon*epsilon);
                T s1 = fC * std::sqrt(vrEpsPp);
                T d_s1 = -fC * vr * vrEpsPp;

                T dz = 2.0 * SR * (z_av - zprevh);

                T signVr = (vr >= 0.0) ? 1.0 : -1.0;
                signVr = (vr == 0.0) ? 0.0 : signVr;
                T absVr = std::abs(vr);

                T z_ss, dz_ss;
                if (vr == 0.0) {
                    z_ss = fS * oneOverS0;
                    dz_ss = 2.0 * fS * oneOverS0;
                } 
                else {
                    T expVrVs = std::exp(-std::pow(std::abs(vr / vS), Sexp));
                    z_ss = signVr * (fC + (fS - fC) * expVrVs) * oneOverS0;
                    dz_ss = (-Sexp * signVr * std::pow(vr, Sexp - 1.0)) / (std::pow(vS, Sexp) * s0)
                    * ((fS - fC) * expVrVs);
                }

                T signZav = (z_av >= 0.0) ? 1.0 : -1.0;
                signZav = (z_av == 0.0) ? 0.0 : signZav;
                T absZav = std::abs(z_av);
                T absZss = std::abs(z_ss);
                T signZss = (z_ss >= 0.0) ? 1.0 : -1.0;
                signZss = (z_ss == 0.0) ? 0.0 : signZss;

                T alpha, d_alpha_z, d_alpha_v;
                if ((signVr == signZav && absZav <= z_ba) || signVr != signZav) {
                    alpha = 0.0;
                    d_alpha_z = 0.0;
                    d_alpha_v = 0.0;
                } 
                else if (signVr == signZav && absZav >= absZss) {
                    alpha = 1.0;
                    d_alpha_z = 0.0;
                    d_alpha_v = 0.0;
                } 
                else {
                    T oneOverZssZba = 1.0 / (absZss - z_ba);
                    alpha = 0.5 * (1.0 + signZav * std::sin(M_PI * (z_av - 0.5 * signZav * (absZss + z_ba)) * oneOverZssZba));

                    T dz_ss_vAbs = signZss * dz_ss;

                    T theta = M_PI * (absZav - (absZss + z_ba) * 0.5) * oneOverZssZba;
                    d_alpha_z = signZav * (M_PI * 0.5 * std::cos(signZav * theta) * oneOverZssZba);
                    d_alpha_v = dz_ss_vAbs * ((z_ba - absZav) * oneOverZssZba * oneOverZssZba * M_PI * 0.5 * std::cos(signZav * theta));
                }

                    T oneOverZss = 1.0 / z_ss;
                    T gn = vr * (1.0 - alpha * z_av * oneOverZss);
                    T fvz = s0 * z_av + s1 * gn;

                    T Nt1 = vr + IJbow * fvz - sn + vb;
                    T Nt2 = gn - dz;

                    T dgn_z = -vr * oneOverZss * (d_alpha_z * z_av + alpha);
                    T dgn_v = 1.0 - z_av * ((alpha + d_alpha_v * vr) * z_ss - dz_ss * alpha * vr) * oneOverZss * oneOverZss; // diff v

                    T zv_denom = 1.0 / (2.0 * IJbow * d_s1 * gn - dgn_z * k + 2.0 * IJbow * dgn_v * s1
                    + IJbow * dgn_v * k * s0 - IJbow * d_s1 * dgn_z * gn * k + 2.0);
                    T zv11 = -(dgn_z * k - 2.0) * zv_denom;
                    T zv12 = IJbow * k * (s0 + dgn_z * s1) * zv_denom;
                    T zv21 = dgn_v * k * zv_denom;
                    T zv22 = -k * (IJbow * d_s1 * gn + IJbow * dgn_v * s1 + 1.0) * zv_denom;
                    
                    T diff_vr = ( zv11 * Nt1 + zv12 * Nt2 );
                    T diff_zav = (zv21 * Nt1 + zv22 * Nt2 );

                    err = std::sqrt( diff_vr*diff_vr + diff_zav*diff_zav);

                    vr = vr - diff_vr;
                    z_av = z_av - diff_zav;

                    iter++;
                    nr = (err > NEWTON_TOL && iter < NEWTON_MAXITER);
            }
        }

        T gradientScale(T Gtilde, T Gtilde_prev, T chiprev, T rfree) {
            T G = Gtilde;
            T gamma = 1;

            T val = 4 * chiprev + G * rfree;
            T xi, xi_prev;
            if ( val < 0 ) { // Otherwise G = Gtilde
                xi = Gtilde * rfree;
                xi_prev = Gtilde_prev * rfree;
                if (std::abs(xi) > 0 && xi < -4 * chiprev ) {
                    gamma = - 4 * chiprev / xi;
                    G = Gtilde * gamma;
                }
                else if (std::abs(xi) == 0 && std::abs(xi_prev) > 0) {
                    gamma = - 4 * chiprev / xi_prev;
                    G = Gtilde_prev * gamma;
                }
            }
            return G;
        }

        /**
         * @brief Set the finger position
         * 
         * @param _xf new finger position [m]
         */
        void setFingerPos(T _xf) {
            xf = _xf;
            itpl(If1, L1, xf, h1, 3);
            itpl(If2, L2, xf, h2, 3);
            Jf1 = If1.transpose() / h1;
            Jf2 = If2.transpose() / h2;
            IJf1 = If1.dot(Jf1);
            IJf2 = If2.dot(Jf2);
        }

        /**
         * @brief Set the finger parameters
         * 
         * @param _mf finger mass [kg]
         * @param _kf finger stiffness [N/m^alpha]
         * @param _rf finger damping [s/m]
         * @param _alf nonlinear exponent [-]
         * @param _u0f initial finger position [m]
         */
        void setFingerParams(T _mf, T _kf, T _rf, T _alf, T _u0f){
            mf = _mf;
            kf = _kf;
            rf = _rf;
            alf = _alf;
            u0f = _u0f;
        }

        /**
         * @brief Set the bow parameters
         * 
         * @param _mh bow hair mass [kg]
         * @param _kh bow hair stiffness [N/m]
         * @param _rh bow hair damping constant [kg/s]
         * @param _xb bowing position [m]
         */
        void setBowParams(T _mh, T _kh, T _rh, T _xb) {
            mh = _mh;
            kh = _kh;
            rh = _rh;
            xb = _xb;
        }

        /**
         * @brief Set the friction parameters
         * 
         * @param _s0 bristle stiffness [N/m]
         * @param _s1bar damping of bristle [kg/s]
         * @param _vS Stribeck velocity
         * @param _muC Columb friction coefficient
         * @param _muS Stribeck friction coefficient
         * @param _Sexp Stribeck exponent
         */
        void setFrictionParams(T _s0, T _s1bar, T _vS, T _muC, T _muS, T _Sexp){
            s0 = _s0;
            s1bar = _s1bar;
            vS = _vS;
            muC = _muC;
            muS = _muS;
            Sexp = _Sexp;
            oneOverS0 = 1 / s0;
        }

        /**
         * @brief Set first string parameters
         * 
         * @param f0 fundamental string frequency [Hz]
         * @param L string length [m]
         * @param r string radius [m]
         * @param rho string density [kg/m^3]
         * @param E Young's modulus [Pa]
         * @param sig0 freq. independant damping [1/s]
         * @param sig1 freq. dependant damping [m^2/s]
         * @param _xcoup bridge coupling position [m]
         */
        void setString1Params(T f0, T L, T r, T rho, T E, T sig0, T sig1, T _xcoup)
        {
            f01 = f0;
            L1 = L;
            E1 = E;
            rho1 = rho;
            r1 = r;
            sig01 = sig0;
            sig11 = sig1;
            xcoup = _xcoup;
            
            // Derived parameters
            c1 = 2 * f01 * L1 * xcoup;
            rhoA1 = rho1 * ( M_PI * r1 * r1);
            EI1 = E1 * ( 0.25 * M_PI * std::pow( r1, 4 ) );
            Kp1 = std::sqrt( EI1 / rhoA1 );
        }

        /**
         * @brief Set second string parameters
         * 
         * @param f0 fundamental string frequency [Hz]
         * @param L string length [m]
         * @param r string radius [m]
         * @param rho string density [kg/m^3]
         * @param E Young's modulus [Pa]
         * @param sig0 freq. independant damping [1/s]
         * @param sig1 freq. dependant damping [m^2/s]
         * @param _xcoup bridge coupling position [m]
         */
        void setString2Params(T f0, T L, T r, T rho, T E, T sig0, T sig1, T _xcoup)
        {
            f02 = f0;
            L2 = L;
            E2 = E;
            rho2 = rho;
            r2 = r;
            sig02 = sig0;
            sig12 = sig1;
            xcoup = _xcoup;
            
            // Derived parameters
            c2 = 2 * f02 * L2 * xcoup;
            rhoA2 = rho2 * ( M_PI * r2 * r2);
            EI2 = E2 * ( 0.25 * M_PI * std::pow( r2, 4 ) );
            Kp2 = std::sqrt( EI2 / rhoA2 );
        }

        /**
         * @brief Set the body modal parameters
         * 
         * @param _M Modal mass vector [kg]
         * @param _K Modal stiff vector [N/m]
         * @param _R Modal dampings [1/s]
         */
        void setBodyParameters(std::vector<T> & _M, std::vector<T> & _K, std::vector<T> & _R)
        {
            assert(_M.size() == _K.size() && _K.size() == _R.size());
            Nm = _M.size();
            Mq = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(_M.data(), _M.size());
            Kq = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(_K.data(), _K.size());
            Rq = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(_R.data(), _R.size());
        }

        private: 
        // Simulation parameters
        T SR, k, ksqr;

        // String 1 parameters
        T f01, L1, r1, rho1, rhoA1, T1, E1, c1, EI1, Kp1, sig01, sig11;
        T hmin1, h1;
        size_t N1;

        // String 2 parameters
        T f02, L2, r2, rho2, rhoA2, T2, E2, c2, EI2, Kp2, sig02, sig12;
        T hmin2, h2;
        size_t N2;

        // Finger parameters
        // mass, stiffness, nonlinear exponent, position, force
        T mf, kf, rf, alf, xf, u0f;;

        // Friction parameters
        T s0, s1, s1bar, vS, muC, muS, Sexp;

        // Bow hair parameters
        // Mass, stiffness, damping
        T mh, kh, rh;
        T xb;

        // Body parameters
        size_t Nm;
        Eigen::Matrix<T, Eigen::Dynamic, 1> Mq, Kq, Rq; // Mass stiffness damping
        T xcoup;     // coupling point 

        // Computation states
        // For string 1
        Eigen::Matrix<T, Eigen::Dynamic, 1> uprev1, u1, unext1;
        T v1, vr1, fn1;
        T haprev1, ha1, hanext1;
        T z1, zprevh1, z_av1, znexth1, chinext1, chiprev1;
        T etaprev1;

        // For string 2
        Eigen::Matrix<T, Eigen::Dynamic, 1> uprev2, u2, unext2;
        T v2, vr2, fn2;
        T haprev2, ha2, hanext2;
        T z2, zprevh2, z_av2, znexth2, chinext2, chiprev2;
        T etaprev2;

        bool nr;

        // finger position
        T ufprev, uf;
        static constexpr T eps = std::numeric_limits<T>::epsilon();
        T Gtilde1_prev, Gtilde2_prev;
        
        // Body states
        Eigen::Matrix<T, Eigen::Dynamic, 1> q, qprev, qnext;

        // Computation matrices
        // First string
        Eigen::SparseMatrix<T, Eigen::RowMajor> As1, Bs1, AAs1, BBs1;
        T os1, oos1, Cs1;

        // Second string
        Eigen::SparseMatrix<T, Eigen::RowMajor> As2, Bs2, AAs2, BBs2;
        T os2, oos2, Cs2;

        // For bow hair
        T oh, Ah, Bh, Ch, ooh, AAh, BBh; 

        // For body
        Eigen::Matrix<T, Eigen::Dynamic, 1> Aq, Bq, Cq;
        
        // Coupling vectors
        Eigen::SparseVector<T, Eigen::RowMajor> Icoup1, Ibow1, If1, Icoup2, Ibow2, If2;
        Eigen::SparseVector<T> Jcoup1, Jbow1, Jf1, Jcoup2, Jbow2, Jf2;
        T IJbow1, IJbow2, IJcoup1, IJcoup2, IJf1, IJf2;

        // Computation Matrices
        Eigen::SparseVector<T, Eigen::RowMajor> IcoupAs1, IcoupAs2, IcoupBs1, IcoupBs2;

        T coupCs11, coupCs12, coupCs21, coupCs22;

        Eigen::SparseVector<T> CsJcoup1, CsJcoup2, CsJbow1, CsJbow2;

        T oneOverMf, oneOverRhoAOs1, oneOverRhoAOs2;

        T oneOverS0;


        // Private functions
        Eigen::SparseMatrix<T> makeDxx(int N, T h){
            Eigen::SparseMatrix<T> Dxx (N, N);
            std::vector<Eigen::Triplet<T>> triplets;
            triplets.reserve(3*N - 2);
            for (int i = 0; i < N; ++i){
                triplets.emplace_back(i, i, -2.0);
                
                if (i + 1 < N) {
                    triplets.emplace_back(i, i + 1, 1.0); // superdiagonal
                    triplets.emplace_back(i + 1, i, 1.0); // subdiagonal
                }
            }
            Dxx.setFromTriplets(triplets.begin(), triplets.end());
            return Dxx / (h * h);
        }

        void itpl(Eigen::SparseVector<T, Eigen::RowMajor>& I,
            T L, T xin, T h, int order)
        {
            I.data().clear();

            const int N = static_cast<int>(I.size());
            const T in_meter = L * xin;
            const int l0 = static_cast<int>(std::floor(in_meter / h));
            const T al0 = in_meter / h - l0;

                if (order == 1) {
                    if (l0 - 1 < N && l0 - 1 >= 0) {
                    I.insert(l0 - 1) = 1.0;
                    }
                    else {
                        I.insert( std::clamp( l0 - 1, 0, N - 1) ) = 1.0;
                    }
                }
                if (order == 2) {
                    if (l0 - 1 < N && l0 - 1 >= 0) {
                    I.insert(l0 - 1) = 1.0 - al0;
                    I.insert(l0    ) = al0;
                    }
                    else {
                         I.insert( std::clamp( l0 - 1, 0, N - 1) ) = 1.0;
                    }
                }
                if (order == 3) {
                    if (l0 - 2 >= 0 && l0 + 1 < N ) {
                    I.insert(l0 - 2) = (al0 * (al0 - 1) * (al0 - 2)) / -6.0 ;
                    I.insert(l0 - 1) = ((al0 - 1) * (al0 + 1) * (al0 - 2)) / 2.0 ;
                    I.insert(l0    ) = (al0 * (al0 + 1) * (al0 - 2)) / -2.0 ;
                    I.insert(l0 + 1) = (al0 * (al0 + 1) * (al0 - 1)) / 6.0 ;
                    }
                    else {
                         I.insert( std::clamp( l0 - 1, 0, N - 1) ) = 1.0;
                    }

                }
            }

    };
}