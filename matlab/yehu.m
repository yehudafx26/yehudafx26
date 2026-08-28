function [out, Fbri] = yehu(vb, fn, xf, ff, ...
  fStringParams, cStringParams, bodyParams, fingerParams, bowParams, frictionParams, SR, ges_num, aniPlot)
% YEHU.M
% Simulate a yehu by applying bow velocity and force VB and FN and
% controlling the finger position XF and force FF. Returns the bridge
% velocity OUT.
%
% Physical parameters supplied by the FSTRINGPARAMS, CSTRINGPARAMS,
% BODYPARAMS, FINGERPARAMS, BOWPARAMS, and FRICTIONPARAMS. See YEHUPARAMS.M
% for details.
%
% Simulation sampling rate is SR.
%
% Champ C. Darabundit and Zhen Zhang, CAML 2026

%% TIME CONSTANTS %%%%%
k = 1 / SR;           % time step
ksqr = k * k;
Nf = length(vb);  % number of discrete time steps

maxIter = 3;          % max number iterations in NR
tol = 1e-9;           % NR threshold
fNthresh = 0.05;       % Minimum bow force

%% STRING PARAMETERS %%%%%
% [f0, L, r, rho, E, sig0, sig1]
% Load and compute derived parameters
fStringParams = num2cell(fStringParams);
[f01, L1, r1, rho1, E1, sig01, sig11, xcoup] = fStringParams{:};
cf1 = 2*f01*L1*xcoup;
rhoA1 = rho1 * ( pi() * r1 * r1);
EI1 = E1 * ( 0.25 * pi() * r1^4 );
Kp1 = sqrt( EI1 / rhoA1 );
T1 = cf1*cf1*rhoA1;

% Minimum spatial sampling for string 1
hmin1 = sqrt((cf1*cf1*k*k+4*sig11*k+sqrt((cf1*cf1*k*k+4*sig11*k)^2+16*Kp1*Kp1*k*k))*0.5);
N1 = floor(L1/hmin1) ;
h1 = L1 / N1 ;

cStringParams = num2cell(cStringParams);
[f02, L2, r2, rho2, E2, sig02, sig12, xcoup] = cStringParams{:};
cf2 = 2*f02*L2*xcoup;
rhoA2 = rho2 * ( pi() * r2 * r2);
EI2 = E2 * ( 0.25 * pi() * r2^4 );
Kp2 = sqrt( EI2 / rhoA2 );
T2 = cf2*cf2*rhoA2;

% Minimum spatial sampling for string 2
hmin2 = sqrt((cf2*cf2*k*k+4*sig12*k+sqrt((cf2*cf2*k*k+4*sig12*k)^2+16*Kp2*Kp2*k*k))*0.5);
N2 = floor(L2/hmin2) ;
h2 = L2 / N2 ;

%% FINGER PARAMETERS %%%%%
% [mf, kf, rf, alf]
fingerParams = num2cell(fingerParams);
[mf, kf, rf, alf, x0] = fingerParams{:};

%% BOW/FRICTION PARAMETERS %%%%%
% [mh, kh, rh, s0, s1_bar, vS, muC, muS, Sexp]
bowParams = num2cell( bowParams );
[mh, kh, rh, xb] = bowParams{:};
frictionParams = num2cell( frictionParams );
[s0, s1_bar, vS, muC, muS, Sexp] = frictionParams{:};

%% BODY PARAMETERS %%%%%
% [Mq, Kq, Rq]
Mq = bodyParams(:, 1);
Kq = bodyParams(:, 2);
Rq = bodyParams(:, 3);
Nm = length(Mq);

%% COMPUTATION CONSTANTS %%%%%

% STRING 1 %%%%%
Dxx1 = sparse(1/h1^2*toeplitz([-2 1 zeros(1,N1-3)]));
Dxxxx1 = Dxx1*Dxx1;

Id1 = sparse(eye(N1-1));

os1 = 1+sig01*k;
As1 = (2*Id1 +(2*sig11*k+cf1^2*ksqr)*Dxx1 - Kp1^2*ksqr*Dxxxx1)/os1;
Bs1 = ((sig01*k-1)*Id1 -(2*sig11*k)*Dxx1)/os1;
Cs1 = k*k/os1;

oos1 = 2 * SR + 2*sig01;
AAs1 = ((cf1*cf1+2*sig11 * SR)*Dxx1-Kp1*Kp1*Dxxxx1+2 * SR * SR*speye(N1-1))/oos1;
BBs1 = ((-2*sig11 * SR)*Dxx1-2 * SR * SR*speye(N1-1))/oos1;

% STRING 2 %%%%%
Dxx2 = sparse(1/h2^2*toeplitz([-2 1 zeros(1,N2-3)]));
Dxxxx2 = Dxx2*Dxx2;

Id2 = sparse(eye(N2-1));

os2 = 1+sig02*k;
As2 = (2*Id2 +(2*sig12*k+cf2^2*ksqr)*Dxx2 - Kp2^2*ksqr*Dxxxx2)/os2;
Bs2 = ((sig02*k-1)*Id2 -(2*sig12*k)*Dxx2)/os2;
Cs2 = k*k/os2;

oos2 = 2 * SR + 2*sig02;
AAs2 = ((cf2*cf2+2*sig12 * SR)*Dxx2-Kp2*Kp2*Dxxxx2+2 * SR * SR*speye(N2-1))/oos2;
BBs2 = ((-2*sig12 * SR)*Dxx2-2 * SR * SR*speye(N2-1))/oos2;

% BOW HAIR %%%%%
oh = mh + kh*k*k/4 + rh*k*0.5;
Ah =  (2*mh - k*k*kh*0.5)/oh;
Bh = (-mh - k*k*kh/4 + k*rh*0.5)/oh;
Ch = k*k/oh;

ooh = 2 * SR+k*kh*0.5/mh+rh/mh;
AAh = (2 * SR * SR+kh*0.5/mh-kh/mh)/ooh;
BBh = (-2 * SR * SR-kh*0.5/mh)/ooh;

% BODY %%%%%
Aq = (2*Mq-k*k*Kq*0.5)./(Mq+k*k*Kq/4+k*Rq*0.5);
Bq = (-Mq-k*k*Kq/4+k*Rq*0.5)./(Mq+k*k*Kq/4+k*Rq*0.5);
Cq = (k*k)./(Mq+k*k*Kq/4+k*Rq*0.5);

% INTERPOLATORS %%%%%
[Ibow1,Jbow1]=itpl(L1, xb , h1 , N1-1 ,1);
[Icoup1,Jcoup1] = itpl(L1, xcoup, h1, N1-1, 1);
[If1, Jf1] = itpl(L1, xf(1), h1, N1 - 1, 3);

IJbow1 = Ibow1*Jbow1/rhoA1/oos1 +  1/mh/ooh;
IJcoup1 = (Icoup1*Cs1*Jcoup1/(rhoA1)+sum(Cq));

[Ibow2,Jbow2]=itpl(L2, xb ,h2 , N2-1 ,1);
[Icoup2,Jcoup2] = itpl(L2, xcoup, h2, N2-1, 1);
[If2, Jf2] = itpl(L2, xf(1), h2, N2-1, 3);

IJbow2 = Ibow2*Jbow2/rhoA2/oos2 +  1/mh/ooh;
IJcoup2 = (Icoup2*Cs2*Jcoup2/(rhoA2)+sum(Cq));

% Fixed matrices for computing coupling with bridge
IcoupAs1 = Icoup1 * As1;
IcoupAs2 = Icoup2 * As2;
IcoupBs1 = Icoup1 * Bs1;
IcoupBs2 = Icoup2 * Bs2;
sumCq = sum(Cq);
D1 = Icoup1*Cs1*Jcoup1/(rhoA1) + sumCq;
D2 = Icoup2*Cs2*Jcoup2/(rhoA2) + sumCq;
coupDenom = 1./(sum(Cq)^2 - D1*D2);
coupCs11 = -D2 * coupDenom; coupCs12 = sumCq * coupDenom;
coupCs21 = sumCq * coupDenom; coupCs22 = -D1 * coupDenom;

% Fixed matrices for bowing and bridge projection
CsJcoup1 = Cs1 * Jcoup1 / rhoA1;
CsJcoup2 = Cs2 * Jcoup2 / rhoA2;
CsJbow1 = Cs1 * Jbow1 / rhoA1;
CsJbow2 = Cs2 * Jbow2 / rhoA2;

% Finger
oneOverMf = 1./mf;
oneOverRhoAOs1 = 1/(rhoA1 * os1);
oneOverRhoAOs2 = 1/(rhoA2 * os2);
IJf1 = If1 * Jf1;
IJf2 = If2 * Jf2;

% Bowing
oneOverS0 = 1 / s0;

% Gradient
Gtilde1_prev = 0;
Gtilde2_prev = 0;


%% STATE INITIALIZATION %%%%%
uprev1 = zeros(N1-1,1); u1 = uprev1;   unext1 = uprev1;
haprev1 = 0; ha1 = haprev1; hanext1 = 0;
v1 = zeros(1,Nf); v1(1) = 0; vr1 = v1(1);
z1 = 0; zprevh1 = 0; z_av1 = 0.5 * (z1 + zprevh1); znexth1 = 0;
fN1 = 0;

uprev2 = zeros(N2-1,1); u2 = uprev2;   unext2 = uprev2;
haprev2 = 0; ha2 = haprev2; hanext2 = 0;
v2 = zeros(1,Nf); v2(1) = 0; vr2 = v2(1);
z2 = 0; zprevh2 = 0; z_av2 = 0.5 * (z2 + zprevh2); znexth2 = 0;
fN2 = 0;

chiprev1    = 0 ; chiprev2 = 0;

u1l = length(u1); u2l = length(u2);

qprev = zeros(length(Nm),1);
q = qprev; qnext = qprev;

ufprev      = x0 ;
uf      = ufprev ;
etaprev1    =  If1*uprev1 - ufprev;
etaprev2 = If2*uprev2- ufprev;

ub = 0;

%% VIDEO SET

if aniPlot == 1
  if ~exist('../videos', 'dir')
    mkdir('../videos');
  end
  fps = 60;
  frameStep = round(1/(fps*k));
  filename = sprintf('gess%d.mp4', ges_num);
  filepath = fullfile('../videos', filename);
  vid = VideoWriter(filepath, 'MPEG-4');
  vid.FrameRate = fps;
  open(vid);
end


%% PROCESS LOOP
out = zeros(Nf, 1);
Fbri = zeros(Nf, 1);
tic
for n = 2:Nf

  % ---------- BRIDGE ----------
  bf1 = sum(Aq.*q + Bq.*qprev) - (IcoupAs1*u1 + IcoupBs1*uprev1);
  bf2 = sum(Aq.*q + Bq.*qprev) - (IcoupAs2*u2 + IcoupBs2*uprev2);

  Fs1 = coupCs11 * bf1 + coupCs12 * bf2;
  Fs2 = coupCs21 * bf1 + coupCs22 * bf2;

  qnext = Aq.*q+Bq.*qprev - Cq.*(Fs1+Fs2);

  Fbri(n) = Fs1+Fs2; % Store bridge force

  % ---------- BOW ----------
  vb_in = vb(n);
  fn_in = fn(n);


  sn1 = Ibow1 * (AAs1*u1 + BBs1*uprev1) + AAh*ha1 + BBh*haprev1 ;
  sn2 = Ibow2 * (AAs2*u2 + BBs2*uprev2) + AAh*ha2 + BBh*haprev2 ;

  if fn_in > fNthresh    % bow F string

    fC = muC*fn_in;
    fS = muS*fn_in;
    z_ba = 0.7*fC * oneOverS0;
    epsilon = fC/s1_bar;
    [vr1, z_av1] = nr_bow(vb_in, vr1, z_av1, zprevh1, fC, fS, z_ba, epsilon, IJbow1, sn1, s0, oneOverS0, vS, Sexp, SR, k, tol, maxIter);

    vrEpsPp = 1./(vr1*vr1 + epsilon*epsilon);
    s1 = fC .* sqrt(vrEpsPp);

    znexth1 = 2*z_av1 - zprevh1;
    gn = (znexth1 - zprevh1) * SR;

    F_fr1 = s0*z_av1 + s1.*gn;

    es1 = As1 * u1 + Bs1 * uprev1 - CsJbow1*F_fr1+ CsJcoup1*Fs1;

    hanext1 = Ah*ha1 + Bh*haprev1 - Ch*F_fr1;

    es2 = As2 * u2 + Bs2 * uprev2+ CsJcoup2 * Fs2;
    hanext2 = Ah*ha2 + Bh*haprev2;

  elseif fn_in < -fNthresh % bow C string
    fC = -muC*fn_in;
    fS = -muS*fn_in;
    z_ba = 0.4 * fC * oneOverS0;
    epsilon = fC/s1_bar;
    [vr2, z_av2] = nr_bow(vb_in, vr2, z_av2, zprevh2, fC, fS, z_ba, epsilon, IJbow2, sn2, s0, oneOverS0, vS, Sexp, SR, k, tol, maxIter);

    vrEpsPp = 1./(vr2*vr2 + epsilon*epsilon);
    s1 = fC .* sqrt(vrEpsPp);

    znexth2 = 2*z_av2 - zprevh2;
    gn = (znexth2 - zprevh2) * SR;

    F_fr2 = s0 * z_av2 + s1 * gn;


    es2 = As2 * u2 + Bs2 * uprev2 - CsJbow2*F_fr2 + CsJcoup2*Fs2;

    hanext2 = Ah*ha2 + Bh*haprev2 - Ch*F_fr2;

    es1 = As1 * u1 + Bs1 * uprev1+ CsJcoup1*Fs1;
    hanext1 = Ah*ha1 + Bh*haprev1;

  else  % without bowing


    es1 = As1 * u1 + Bs1 * uprev1+ CsJcoup1 * Fs1;
    hanext1 = Ah*ha1 + Bh*haprev1;
    es2 = As2 * u2 + Bs2 * uprev2+ CsJcoup2 * Fs2;
    hanext2 = Ah*ha2 + Bh*haprev2;

  end

  % ---------- FINGER COLLISION ----------
  if xf(n) ~= xf(n-1) % Recompute interpolators if finger moves
    [If1, Jf1] = itpl(L1, xf(n), h1, N1-1, 3);
    [If2, Jf2] = itpl(L2, xf(n), h2, N2-1, 3);
    IJf1 = If1 * Jf1;
    IJf2 = If2 * Jf2;
  end

  f1      = 2*uf - ufprev - ff(n)* ksqr * oneOverMf;
  

  % Constraint mu_t chi1
  eta1    = If1*u1 - uf;
  etast1   = If1*es1 - f1 ;
  eta1 = max(eta1, 0);
  V1 = kf/(alf+1) * eta1^(alf+1);
  dV1deta1 = kf * eta1^alf;

  Gtilde1 = dV1deta1/sqrt(2*V1+eps);

  rfree1 = etast1 - (If1*uprev1-ufprev);
  g1 = gradient_scale(Gtilde1, Gtilde1_prev, chiprev1, rfree1);


  % Constraint mu_t chi2
  eta2    = If2*u2 - uf;
  etast2   =  If2*es2 - f1 ;

  eta2 = max(eta2, 0);
  V2 = kf/(alf+1) * eta2^(alf+1);
  dV2deta2 = kf * eta2^alf;

  Gtilde2 = dV2deta2/sqrt(2*V2+eps);

  rfree2 = etast2 - (If2*uprev2-ufprev);
  g2 = gradient_scale(Gtilde2, Gtilde2_prev, chiprev2, rfree2);

  % SHERMAN-MORRISON COMPUTATION
  g1sqr = g1 * g1;
  g2sqr = g2 * g2;

  pf1 = -0.25*ksqr*g1sqr*(If1*uprev1-ufprev) + chiprev1*g1*ksqr ...
    - g1*rf*k*0.5*(If1*uprev1-ufprev)*chiprev1;

  pf2 = -0.25*ksqr*g2sqr*(If2*uprev2-ufprev) + chiprev2*g2*ksqr ...
    - g2*rf*k*0.5*(If2*uprev2-ufprev)*chiprev2;

  bf1 = es1 - oneOverRhoAOs1  * Jf1 * pf1;
  bf2 = es2 - oneOverRhoAOs2 * Jf2 * pf2;
  bf3 = f1 + oneOverMf * (pf1 + pf2);

  cf1 = oneOverRhoAOs1 * (0.25*ksqr*g1sqr + 0.5*rf*g1*k*chiprev1);
  cf2 = oneOverRhoAOs2 * (0.25*ksqr*g2sqr + 0.5*rf*g2*k*chiprev2);
  denf1 = 1/(1 + cf1*IJf1);
  denf2 = 1/(1 + cf2*IJf2);

  Af1inv_bf1 = bf1 - (cf1 * denf1)*Jf1*(If1*bf1);
  Af2 = oneOverRhoAOs1 * ( -0.25*ksqr*g1sqr*Jf1 - 0.5*rf*g1*k*chiprev1*Jf1);
  Af1inv_Af2 = Af2 - (cf1 * denf1)*Jf1*(If1*Af2);
  Af3inv_bf2 = bf2 - (cf2 * denf2)*Jf2*(If2*bf2);
  Af4 = oneOverRhoAOs2  * ( -0.25*ksqr*g2sqr*Jf2 - 0.5*rf*g2*k*chiprev2*Jf2);
  Af3inv_Af4 = Af4 - (cf2 * denf2)*Jf2*(If2*Af4);

  Af5 = oneOverMf * (-0.25*ksqr*g1sqr*If1 - 0.5*rf*g1*k*chiprev1*If1);
  Af6 = oneOverMf * (-0.25*ksqr*g2sqr*If2 - 0.5*rf*g2*k*chiprev2*If2);

  Af7 = 1 + oneOverMf * ( 0.25*ksqr*g1sqr + 0.5*rf*g1*chiprev1*k ...
    + 0.25*ksqr*g2sqr + 0.5*rf*g2*chiprev2*k );

  S = Af7 - Af5*Af1inv_Af2 - Af6*Af3inv_Af4;

  rhsU = bf3 - Af5*Af1inv_bf1 - Af6*Af3inv_bf2;

  ufnext = rhsU / S;

  unext1 = Af1inv_bf1 - Af1inv_Af2*ufnext;
  unext2 = Af3inv_bf2 - Af3inv_Af4*ufnext;

  chinext1 = 0.5 * g1 * ( If1*(unext1-uprev1) - (ufnext-ufprev) ) + chiprev1;
  chinext2 = 0.5 * g2 * ( If2*(unext2-uprev2) - (ufnext-ufprev) ) + chiprev2;

  muchi1 = 0.5 * (chinext1 + chiprev1);
  muchi2 = 0.5 * (chinext2 + chiprev2);

  % Output
  out(n)   = (Icoup1*unext1-Icoup1*uprev1) * SR*0.5;

  % Check if computation errors
  if muchi1 < -eps
    warning('mu_t chi1 < 0, %.2e', muchi1);
  end

  if muchi2 < -eps
    warning('mu_t chi2 < 0, %.2e', muchi2);
  end

  ub = ub + vb_in*k;

  if aniPlot == 1
    animation2(aniPlot, vid, frameStep, Nf, k, n, N1, N2, u1, u2, h1, h2, xf, L1, ...
      uf, ub, xb, fn, xcoup, Icoup1 );
  end

  % Update the states
  zprevh1 = znexth1;
  uprev1 = u1;
  u1 = unext1;
  haprev1 = ha1;
  ha1 = hanext1;

  ufprev      = uf ;
  uf      = ufnext ;

  etaprev1    = eta1 ;
  chiprev1    = chinext1 ;

  etaprev2    = eta2 ;
  chiprev2    = chinext2 ;

  qprev = q;
  q = qnext;

  Gtilde1_prev = Gtilde1;
  Gtilde2_prev = Gtilde2;

  % update the variables
  zprevh2 = znexth2;
  uprev2 = u2;
  u2 = unext2;
  haprev2 = ha2;
  ha2 = hanext2;
end
toc

  function [vr, z_av] = nr_bow(vb, vr, z_av, zprevh, fC, fS, z_ba, epsilon, IJbow, sn, s0, oneOverS0, vS, Sexp, SR, k, tol, maxIter)

    nr = true;
    err = 1;
    iter = 0;

    while nr
      vrEpsPp = 1./(vr*vr + epsilon*epsilon);
      s1 = fC * sqrt(vrEpsPp);
      d_s1 = - fC * vr * vrEpsPp ;

      dz = 2 * SR * (z_av - zprevh);

      signVr = sign(vr);
      absVr = abs(vr);

      if vr == 0
        z_ss = fS * oneOverS0;
        dz_ss = 2*fS * oneOverS0;
      else
        expVrVs = exp(-abs(vr/vS).^Sexp);
        % z_ss = signVr * (fC + (fS - fC)*expVrVs + s2*absVr) * oneOverS0;
        z_ss = signVr * (fC + (fS - fC)*expVrVs) * oneOverS0;
        % dz_ss = (-Sexp * signVr * vr^(Sexp - 1))./(vS^Sexp*s0) * ((fS-fC)*expVrVs) + signVr*s2 * oneOverS0;
        dz_ss = (-Sexp * signVr * vr^(Sexp - 1))./(vS^Sexp*s0) * ((fS-fC)*expVrVs);
      end

      signZav = sign(z_av);
      absZav = abs(z_av);
      absZss = abs(z_ss);
      signZss = sign(z_ss);


      if (signVr == signZav && absZav <= z_ba) || signVr ~= signZav
        alpha = 0;
        d_alpha_z = 0;
        d_alpha_v = 0;
      elseif signVr == signZav && absZav >= absZss
        alpha = 1;
        d_alpha_z = 0;
        d_alpha_v = 0;
      else
        oneOverZssZba = 1 / (absZss-z_ba);
        alpha = 0.5*(1+signZav.*sin(pi*(z_av - 0.5*signZav.*(absZss+z_ba)) * oneOverZssZba ));

        dz_ss_vAbs = signZss.*dz_ss;

        theta = pi*(absZav-(absZss+z_ba)*0.5) * oneOverZssZba;
        d_alpha_z = signZav * (pi*0.5*cos(signZav*theta) * oneOverZssZba);
        d_alpha_v = dz_ss_vAbs.*( (z_ba-absZav) * oneOverZssZba * oneOverZssZba * pi*0.5.*cos(signZav*theta));
      end

      oneOverZss = 1 / z_ss;
      gn = vr.*(1 - alpha.*z_av * oneOverZss );

      fvz = s0*z_av + s1.*gn;

      Nt1 = vr + IJbow*fvz - sn + vb;
      Nt2 = gn - dz;

      dgn_z = -vr * oneOverZss *(d_alpha_z.*z_av + alpha);
      dgn_v = 1 - z_av.*((alpha + d_alpha_v.*vr).*z_ss - dz_ss.*alpha.*vr) * oneOverZss * oneOverZss;      % diff v

      zv_denom = 1 ./ (2*IJbow*d_s1*gn - dgn_z*k + 2*IJbow*dgn_v*s1 + IJbow*dgn_v*k*s0 - IJbow*d_s1*dgn_z*gn*k + 2);
      zv11 = -(dgn_z*k - 2) * zv_denom;
      zv12 = IJbow*k*(s0 + dgn_z*s1) * zv_denom;
      zv21 = dgn_v*k *zv_denom;
      zv22 = -k*(IJbow*d_s1*gn + IJbow*dgn_v*s1 + 1) * zv_denom;

      diff_vr = ( zv11 * Nt1 + zv12 * Nt2 );
      diff_zav = (zv21 * Nt1 + zv22 * Nt2 );
      % vr1new = vr - ( zv11 * Nt1 + zv12 * Nt2 );
      % z_av1new = z_av - (zv21 * Nt1 + zv22 * Nt2 );

      err = norm([diff_vr, diff_zav]);

      vr = vr - diff_vr;
      z_av = z_av - diff_zav;

      iter = iter + 1;
      nr = (err > tol && iter < maxIter);
    end
  end

  function G = gradient_scale(Gtilde, Gtilde_prev, chiprev, rfree)
    G = Gtilde;
    gamma = 1;

    val = 4*chiprev + G*rfree;

    if val < 0

      xi = Gtilde.' * rfree;
      xiPrev = Gtilde_prev.' * rfree;
      if abs(xi)>0 && xi<-4*chiprev
        gamma = -4*chiprev / xi;
        G = Gtilde*gamma;
      elseif abs(xi)==0 && abs(xiPrev)>0
        gamma = -4*chiprev / xiPrev;

        lambda = 1;

        G = Gtilde_prev*gamma*lambda;

      else
        G = Gtilde*gamma;
      end
    end

  end
end