clear
close all
clc

%% Geometry and fluid properties

r1 = [0.0198, 0.0380, 0.0198, 0.0198, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.15, 0.1, 0.0198];
r2 = [0.0380, 0.0760, 0.0380, 0.0380, 0.3, 0.5, 0.3, 0.5, 0.3, 0.4, 0.6, 0.3, 0.3, 0.3, 0.0380];
lp = [0.325, 0.315, 0.325, 0.325, 1, 1 , 1, 1, 1, 1, 1, 1, 1, 1, 0.325];
lavalues = [0.00, 0.00, 0.00, 0.02, 0.00, 0.00, 0.05, 0.05, 0, 0, 0, 0, 0, 0, 0.02];
lbvalues = [0.00, 0.00, 0.00, 0.03, 0.00, 0.00, 0.10, 0.10, 0, 0, 0, 0, 0, 0, 0.03];

rho0 = 1.18;
c0 = 346;

sigma=0.263;
dh=0.0035; % Given in Kirby [12]
tw=0.001; % Given in Kirby [12]

Rfvalues = [5000, 30700, 5000, 5000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 5000]; % Flow resistivity

Mvalues = [0, 0.163, 0.15];

freq = 1:1:2000;

TLthird = zeros(length(Mvalues),length(r1),length(freq));

TM = zeros(2,2,length(freq),length(Mvalues),length(r1));

TLfirst = zeros(length(Mvalues),length(r1),length(freq));

TLsecond = zeros(length(Mvalues),length(r1),length(freq));

%% Main loop
for i = 1:length(Mvalues)

    M = Mvalues(i);

    for r=1:length(r1)
        if r >= 5 && M ~= 0 
           continue
        end

        r1current = r1(r);
        r2current = r2(r);
        lpcurrent = lp(r); 
        Rf = Rfvalues(r);
        la = lavalues(r);
        lb = lbvalues(r);

        S1 = pi*r1current^2;
        Z01 = rho0*c0/S1;

        %% Build the symbolic eigenequations 

        syms Skz Sk0 Skw Srhow SM SZp
        
        % Forward wave
        krplus  = sqrt(Sk0^2-2*SM*Sk0*Skz-(1-SM^2)*Skz^2);
        krwplus = sqrt(Skw^2-Skz^2);                     
        
        Cbplus = -besselj(1,krwplus*r2current)/bessely(1,krwplus*r2current);

        Fplus = (krplus/(rho0*(1-SM*Skz/Sk0)) * besselj(1,krplus*r1current)) / (besselj(0,krplus*r1current) + 1i*SZp*krplus*besselj(1,krplus*r1current)/(Sk0*(1-SM*Skz/Sk0))) - (krwplus/Srhow) * (besselj(1,krwplus*r1current)+Cbplus*bessely(1,krwplus*r1current)) / (besselj(0,krwplus*r1current)+Cbplus*bessely(0,krwplus*r1current));

        % Backward wave
        krminus  = sqrt(Sk0^2+2*SM*Sk0*Skz-(1-SM^2)*Skz^2);
        krwminus = sqrt(Skw^2-Skz^2);                    
        
        Cbminus = -besselj(1,krwminus*r2current)/bessely(1,krwminus*r2current);
        
        Fminus = (krminus/(rho0*(1+SM*Skz/Sk0)) * besselj(1,krminus*r1current)) / (besselj(0,krminus*r1current) + 1i*SZp*krminus*besselj(1,krminus*r1current)/(Sk0*(1+SM*Skz/Sk0))) - (krwminus/Srhow) * (besselj(1,krwminus*r1current)+Cbminus*bessely(1,krwminus*r1current)) / (besselj(0,krwminus*r1current)+Cbminus*bessely(0,krwminus*r1current));

        dFplus  = diff(Fplus,Skz);
        dFminus = diff(Fminus,Skz);
        
        Fplusfun   = matlabFunction(Fplus,'Vars',[Skz Sk0 Skw Srhow SM SZp]);
        dFplusfun  = matlabFunction(dFplus,'Vars',[Skz Sk0 Skw Srhow SM SZp]);
        Fminusfun  = matlabFunction(Fminus,'Vars',[Skz Sk0 Skw Srhow SM SZp]);
        dFminusfun = matlabFunction(dFminus,'Vars',[Skz Sk0 Skw Srhow SM SZp]);
    
    
        for n = 1:length(freq)
        
            f = freq(n);
        
            omega = 2*pi*f;
            k0 = omega/c0;
        
        
            %% Absorbent material properties
        
            X = rho0*f/Rf;
        
            Zcratio = 1+0.0571*X^(-0.754)-1i*0.087*X^(-0.732); % From Cummings and Chang [9]
        
            kwratio = 1+0.0978*X^(-0.700)-1i*0.189*X^(-0.595); % From Cummings and Chang [9]
        
            Zcw = Zcratio*rho0*c0;
        
            kw = kwratio*k0;
        
            cw = omega/kw;
        
            rhow = Zcw/cw;
           
            Y0 = 1/(rho0*c0);

            Yw= 1/Zcw;
            
            if M == 0
                Zp = (0.006 + 1i*k0*(tw + 0.375*dh*(1 + (Yw/Y0)*(kw/k0))))/sigma;
            else
                Zp = (1/sigma) * ( 7.337e-3*(1 + 72.23*M) + 1i*2.2245e-5*(1 + 51*tw) * (1 + 102*dh*(1 + (Yw/Y0)*(kw/k0))) * f );
            end

            if r==1 || r==2 
                Zp=0;
            end

            %% Initial guess
            
            if n == 1
                    
                a = -1i*Zp*rho0/k0;

                b = 2*rho0/r1current - 2*rhow*r1current/(r1current^2-r2current^2)- 1i*Zp*rho0/k0*(kw^2-k0^2);
                
                c = 2*rho0/r1current*(kw^2-k0^2);

                kr2roots = roots([a b c]);
                
                kzcandidates = sqrt(k0^2 - kr2roots);
                
                if length(kzcandidates())<2 || imag(kzcandidates(1))<=0
                    kz0 = kzcandidates(1);
                else
                    kz0 = kzcandidates(2);
                end
                
                %% Converting the M = 0 radial wavenumber to kz+ and kz-
                        
                kr2 = k0^2-kz0^2;
                        
                kzplus = (-M*k0+sqrt(k0^2-(1-M^2)*kr2))/(1-M^2);
                        
                kzminus = (M*k0+sqrt(k0^2-(1-M^2)*kr2))/(1-M^2);
             
            else
                
                 kzplus = kzguessplus;
                 kzminus = kzguessminus;
                
            end
        
        
             %% Newton-Raphson
        
             tol = 1e-10;
             maxIter = 200;
             iter=0;
        
             while true
                   iter = iter +1;
        
                   Fp  = Fplusfun(kzplus,k0,kw,rhow,M,Zp);
                   dFp = dFplusfun(kzplus,k0,kw,rhow,M,Zp);
        
                   Fm  = Fminusfun(kzminus,k0,kw,rhow,M,Zp);
                   dFm = dFminusfun(kzminus,k0,kw,rhow,M,Zp);
        
                   kznewplus = kzplus - Fp/dFp;
        
                   kznewminus = kzminus - Fm/dFm;
        
        
                   if abs(kznewplus-kzplus)<tol*abs(kzplus) && abs(kznewminus-kzminus)<tol*abs(kzminus)
        
                      kzplus=kznewplus;
                      kzminus=kznewminus;
        
                      break
        
                   end
        
                   kzplus = kznewplus;
                   kzminus = kznewminus;

                  if iter>=maxIter
                
                     fprintf('No convergence: r=%d, f=%d Hz, M=%g\n',r,f,M);
                     break
        
                  end
        
            end
        
        
            %% Saving roots for next frequency
        
            kzguessplus = kzplus;
            kzguessminus = kzminus;
        
        
            %% Four-pole parameters
        
            alpha2plus = kzplus/(rho0*c0*(k0-M*kzplus));
        
            alpha2minus = kzminus/(rho0*c0*(k0+M*kzminus));

            Yplus = alpha2plus;

            Yminus = alpha2minus;

            expplus = exp(1i*kzplus*lpcurrent);
        
            expminus = exp(-1i*kzminus*lpcurrent);
        

            %% Transfer matrix
            % Third Approach
            if la==0 && lb==0
               A11 = (alpha2minus*expplus+alpha2plus*expminus)/(alpha2minus+alpha2plus);
          
               A12 = (expplus-expminus)/(S1*(alpha2minus+alpha2plus));
          
               A21 = (alpha2minus*alpha2plus*S1*(expplus-expminus))/(alpha2minus+alpha2plus);
          
               A22 = (alpha2plus*expplus+alpha2minus*expminus)/(alpha2minus+alpha2plus);
               
               TM(:,:,n,i,r) = [A11 A12; A21 A22];
            else
               krplus  = sqrt(k0^2 - 2*M*k0*kzplus  - (1-M^2)*kzplus^2);
               krminus = sqrt(k0^2 + 2*M*k0*kzminus - (1-M^2)*kzminus^2);
   
               krwplus  = sqrt(kw^2 - kzplus^2);
               krwminus = sqrt(kw^2 - kzminus^2);
   
               Cbplus  = -besselj(1,krwplus*r2current)/bessely(1,krwplus*r2current);
               Cbminus = -besselj(1,krwminus*r2current)/bessely(1,krwminus*r2current);

               Za = -1i*Zcw*cot(kw*la);
               Zb =  1i*Zcw*cot(kw*lb);

               Bexpplus = exp(-1i*kzplus*lpcurrent);
               Bexpminus = exp(1i*kzminus*lpcurrent);

               psi2pplus = @(r)((besselj(0,krplus*r)).*(besselj(0,krwplus*r1current) + Cbplus*bessely(0,krwplus*r1current)) ./ (besselj(0,krwplus*r1current) + 1i*Zp*krplus*k0*besselj(1,krplus*r1current)/(k0-M*kzplus)^2)) ;
              
               psi2pminus = @(r)((besselj(0,krminus*r)).*(besselj(0,krwminus*r1current) + Cbminus*bessely(0,krwminus*r1current)) ./ (besselj(0,krminus*r1current) + 1i*Zp*krminus*k0*besselj(1,krminus*r1current)/(k0+M*kzminus)^2 )) ;
              
               psi3pplus = @(r)(besselj(0,krwplus*r1current) + Cbplus*bessely(0,krwplus*r1current));
              
               psi3pminus = @(r)(besselj(0,krwminus*r1current) + Cbminus*bessely(0,krwminus*r1current));
            
               S2plus  = 2*pi * integral(@(r) psi2pplus(r).*r, 0, r1current);
               S2minus = 2*pi * integral(@(r) psi2pminus(r).*r, 0, r1current);
              
               S3plus  = 2*pi * integral(@(r) psi3pplus(r).*r, r1current, r2current);
               S3minus = 2*pi * integral(@(r) psi3pminus(r).*r, r1current, r2current);

               A = [ S2plus/S1, S2minus/S1;alpha2plus*S2plus + S3plus/Za, -(alpha2minus*S2minus + S3minus/Za) ];

               B = [ (S2plus/S1)*Bexpplus, (S2minus/S1)*Bexpminus;(alpha2plus*S2plus - S3plus/Zb)*Bexpplus, -(alpha2minus*S2minus - S3minus/Zb)*Bexpminus ];

               TM(:,:,n,i,r) = A/B;
               A11 = TM(1,1,n,i,r);
               A12 = TM(1,2,n,i,r);
               A21 = TM(2,1,n,i,r);
               A22 = TM(2,2,n,i,r);

            end

   
            % First approach 

            T11 = (Yminus*expplus+Yplus*expminus)/(Yminus+Yplus);

            T12 = (expplus-expminus)/(S1*(Yminus+Yplus));

            T21 = (Yminus*Yplus*S1*(expplus-expminus))/(Yminus+Yplus);

            T22 = (Yplus*expplus+Yminus*expminus)/(Yminus+Yplus);
            
            TMfirst = [T11 T12;T21 T22];

            S2=pi*(r2current^2-r1current^2);
            Z02 = rho0*c0/S2;

            TMla = [cos(k0*la), 1i*Z02*sin(k0*la);1i*sin(k0*la)/Z02, cos(k0*la)];

            TMlb = [cos(k0*lb), 1i*Z02*sin(k0*lb);1i*sin(k0*lb)/Z02, cos(k0*lb)];

            TMfirsttotal =  TMla * TMfirst * TMlb ;

            A11fullthree = TMfirsttotal(1,1);
            A12fullthree = TMfirsttotal(1,2);
            A21fullthree = TMfirsttotal(2,1);
            A22fullthree = TMfirsttotal(2,2);

            if r==15 && M==0 %% Procedure on page 380 of Textbook (Chapter 3) 
                d1 = 2*r1current; 
                d2 = 2*r2current; 
                zeta = Zp/(rho0*c0);
                ka2 = k0^2 - 1i*4*k0/(d1*zeta); 
                kb2 = k0^2 - 1i*4*k0*d1/((d2^2-d1^2)*zeta); 
                alpha1 = 0; 
                alpha2 = ka2; 
                alpha3 = 0; 
                alpha4 = -(ka2-k0^2); 
                alpha5 = 0; 
                alpha6 = -(kb2-k0^2); 
                alpha7 = 0; 
                alpha8 = kb2; 
                
                p = [1, 0, alpha2 + alpha8, 0, alpha2*alpha8 - alpha4*alpha6]; 
                beta = roots(p); 
                
                A0 = zeros(4,4); 
                Al = zeros(4,4); 

                for m = 1:4 
                    bm = beta(m); 
                    psi1 = 1; 
                    psi2 = -(bm^2 + alpha1*bm + alpha2) / (alpha3*bm + alpha4); 
                    psi3 = 1/bm; 
                    psi4 = psi2/bm; 
                    
                    % z = 0 
                    A0(1,m) = psi3; 
                    A0(2,m) = psi4; 
                    A0(3,m) = -1/(1i*k0 + M*bm); 
                    A0(4,m) = -psi2/(1i*k0 + M*bm);

                    % z = lp 
                    e = exp(bm*lpcurrent); 
                    Al(1,m) = psi3*e; 
                    Al(2,m) = psi4*e; 
                    Al(3,m) = -e/(1i*k0 + M*bm); 
                    Al(4,m) = -psi2*e/(1i*k0 + M*bm); 
                end 

                T4 = A0/Al; 

                T11 = T4(1,1); 
                T12 = T4(1,2); 
                T13 = T4(1,3); 
                T14 = T4(1,4); 
                T21 = T4(2,1); 
                T22 = T4(2,2); 
                T23 = T4(2,3); 
                T24 = T4(2,4); 
                T31 = T4(3,1); 
                T32 = T4(3,2); 
                T33 = T4(3,3); 
                T34 = T4(3,4); 
                T41 = T4(4,1); 
                T42 = T4(4,2); 
                T43 = T4(4,3); 
                T44 = T4(4,4); 

                % 4x4 to 2x2 
                X1 = -1i*tan(k0*la); 
                X2 = 1i*tan(k0*lb); 
                F1 = T42 + X2*T44 - X1*(T22 + X2*T24); 
                A1 = (X1*T21 - T41)/F1; 
                B1 = (X1*T23 - T43)/F1; 
                A2 = T12 + X2*T14; 
                B2 = T32 + X2*T34; 
                Ta = T11 + A1*A2; 
                Tb = T13 + B1*A2; 
                Tc = T31 + A1*B2; 
                Td = T33 + B1*B2; 

                TMsecondrhoc0 = [Ta Tb; Tc Td]; 
                TMsecond = [TMsecondrhoc0(1,1), TMsecondrhoc0(1,2)*rho0*c0; TMsecondrhoc0(2,1)/(rho0*c0), TMsecondrhoc0(2,2)]; 
                
                A = TMsecond(1,1); 
                B = TMsecond(1,2); 
                C = TMsecond(2,1); 
                D = TMsecond(2,2); 
                
                TLsecond(i,r,n) = 20*log10(abs((A + B/Z01 + C*Z01 + D)/2)); 
            end
            %% Transmission loss
            
            TLthird(i,r,n) = 20*log10(abs((A11+A12/Z01+A21*Z01+A22)/2));
            TLfirst(i,r,n) = 20*log10(abs((A11fullthree + A12fullthree/Z01 + A21fullthree*Z01 + A22fullthree)/2));

         end
        
     end
end

%% Plot

figure(3)
hold on
plot(freq,squeeze(TLthird(1,1,:))) %If an array has dimensions like m × 1 × n × 1, squeeze will remove the 1-sized dimensions, resulting in m × n (was not letting me input directly into plot)
plot(freq,squeeze(TLthird(2,1,:)))
xlabel('Frequency (Hz)')
ylabel('Transmission loss (dB)')
legend('M=0','M=0.163','Location','best')
title('Comparison of the Third Approach M=0, 0.163')
hold off

figure(4)
hold on
plot(freq,squeeze(TLthird(1,2,:)))
plot(freq,squeeze(TLthird(3,2,:)))
xlabel('Frequency (Hz)')
ylabel('Transmission loss (dB)')
legend('M=0','M=0.15','Location','best')
title('Comparison of the Third Approach M=0, 0.15')
hold off

figure(5)
hold on
plot(freq, squeeze(TLthird(1,3,:)), '-')
plot(freq, squeeze(TLfirst(1,3,:)), '--')
plot(freq, squeeze(TLthird(2,3,:)), '-')
plot(freq, squeeze(TLfirst(2,3,:)), '--')
xlabel('Frequency (Hz)')
ylabel('TLa (dB)')
legend('Third approach, M=0.0','First approach, M=0.0', 'Third approach, M=0.163', 'First approach, M=0.163', 'Location','best')
title('Comparison of the first and third approaches')
hold off

figure(6)
hold on
plot(freq, squeeze(TLthird(1,4,:)), '-')
plot(freq, squeeze(TLfirst(1,4,:)), '--')
plot(freq, squeeze(TLthird(2,4,:)), '-')
plot(freq, squeeze(TLfirst(2,4,:)), '--')
xlabel('Frequency (Hz)')
ylabel('Transmission loss (dB)')
legend('Third approach, M=0.0', 'First approach, M=0.0', 'Third approach, M=0.163', 'First approach, M=0.163', 'Location','best')
title('Comparison of the first and third approaches, la=2 cm, lb=3 cm')
hold off

figure(7)
hold on
plot(freq,squeeze(TLsecond(1,15,:)),'-')
plot(freq,squeeze(TLthird(1,15,:)),'--')
xlabel('Frequency (Hz)')
ylabel('Transmission loss (dB)')
legend('Second approach','Third approach','Location','best')
title('Comparison of the second and third approaches, M=0, la=2 cm, lb=3 cm')
hold off

figure(8)
hold on
plot(freq,squeeze(TLthird(1,5,:)))
plot(freq,squeeze(TLthird(1,6,:)))
xlabel('Frequency (Hz)')
ylabel('TLa (dB)')
legend('OD=60 cm, ID=40 cm (d/h=1, l/h=10)','OD=100 cm, ID=80 cm (d/h=0.5, l/h=5)','Location','best')
title('Transmission loss for two large mufflers')
hold off

figure(9)
hold on
plot(freq,squeeze(TLthird(1,5,:)), '-')
plot(freq,squeeze(TLthird(1,6,:)), '--')
plot(freq,squeeze(TLthird(1,7,:)), '-')
plot(freq,squeeze(TLthird(1,8,:)), '--')
xlabel('Frequency (Hz)')
ylabel('TLa (dB)')
legend('(a) r2=30 cm, la=lb=0', '(b) r2=50 cm, la=lb=0', '(c) r2=30 cm, la=5 cm, lb=10 cm', '(d) r2=50 cm, la=5 cm, lb=10 cm', 'Location','best')
title('Effect of annular cavities')
hold off

figure(10)
hold on
plot(freq,squeeze(TLthird(1,9,:)), '-')
plot(freq,squeeze(TLthird(1,10,:)), '--')
plot(freq,squeeze(TLthird(1,11,:)), '.-')
xlabel('Frequency (Hz)')
ylabel('TLa (dB)')
legend('d/h = 1.0', 'd/h = 2.0','d/h = 4.0','Location','best')
title('Transmission loss curves for different d/h values')
hold off

figure(11)
hold on
plot(freq,squeeze(TLthird(1,12,:)), '-')
plot(freq,squeeze(TLthird(1,13,:)), '--')
plot(freq,squeeze(TLthird(1,14,:)), '-')
xlabel('Frequency (Hz)')
ylabel('TLa (dB)')
legend('d/h = 1.0', 'd/h = 2.0', 'd/h = 4.0', 'Location','best')
title('Transmission loss curves for different d/h values')
hold off