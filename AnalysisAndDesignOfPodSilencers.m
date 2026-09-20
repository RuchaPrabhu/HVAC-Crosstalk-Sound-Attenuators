clear % clear all clears all variables and functions but makes performance slower
clc % clears command window
close all % closes figures
%% Defining Constants
% Values given in Analysis and design of pod silencers
E=15000;
ru=0.05;
rd=0.05;
r3=0.15;
lp=0.50;
lec=0.15;
lcc=0.15;
M=0.05;
rho0=1.18;
c0=346;
r1values=[0.025 0.050 0.075];
r2values5 = [0.075 0.100 0.125];
r2values6 = [0.050 0.075 0.100];

% Perturbation offsets used to hunt for nearby roots of the
% characteristic equation at each frequency, so we don't accidentally
% lock onto the wrong (higher-order) root as frequency increases.
offsets = [0, 0.3, -0.3, 1i*0.3, -1i*0.3, 0.3+0.3i, -0.3-0.3i, 0.6, -0.6, 1i*0.6, -1i*0.6];

syms kzS k0S kwS r1S r2S r3S Y0S YwS

kr1 = sqrt(kwS^2 - kzS^2);
kr2 = sqrt(k0S^2 - kzS^2);
kr3 = kr1;
Asym = sym(zeros(5,5));
Asym(1,1) = -kr1*besselj(1,kr1*r1S)/(kwS*YwS);
Asym(1,2) =  kr2*besselj(1,kr2*r1S)/(k0S*Y0S);
Asym(1,3) =  kr2*bessely(1,kr2*r1S)/(k0S*Y0S);
Asym(2,1) =  besselj(0,kr1*r1S);
Asym(2,2) = -besselj(0,kr2*r1S);
Asym(2,3) = -bessely(0,kr2*r1S);
Asym(3,2) = -kr2*besselj(1,kr2*r2S)/(k0S*Y0S);
Asym(3,3) = -kr2*bessely(1,kr2*r2S)/(k0S*Y0S);
Asym(3,4) =  kr3*besselj(1,kr3*r2S)/(kwS*YwS);
Asym(3,5) =  kr3*bessely(1,kr3*r2S)/(kwS*YwS);
Asym(4,2) =  besselj(0,kr2*r2S);
Asym(4,3) =  bessely(0,kr2*r2S);
Asym(4,4) = -besselj(0,kr3*r2S);
Asym(4,5) = -bessely(0,kr3*r2S);
Asym(5,4) =  kr3*besselj(1,kr3*r3S)/(kwS*YwS);
Asym(5,5) =  kr3*bessely(1,kr3*r3S)/(kwS*YwS);
dAsym = diff(Asym, kzS);

Build_A_num  = matlabFunction(Asym, 'Vars', [kzS k0S kwS r1S r2S r3S Y0S YwS]);
dBuild_A_num = matlabFunction(dAsym, 'Vars', [kzS k0S kwS r1S r2S r3S Y0S YwS]);

kr1np = sqrt(k0S^2 - kzS^2);
kr2np = sqrt(kwS^2 - kzS^2);
AsymNP = sym(zeros(3,3));
AsymNP(1,1) =  besselj(0,kr1np*r2S);
AsymNP(1,2) = -besselj(0,kr2np*r2S);
AsymNP(1,3) = -bessely(0,kr2np*r2S);
AsymNP(2,1) =  kr1np*besselj(1,kr1np*r2S)/(k0S*Y0S);
AsymNP(2,2) = -kr2np*besselj(1,kr2np*r2S)/(kwS*YwS);
AsymNP(2,3) = -kr2np*bessely(1,kr2np*r2S)/(kwS*YwS);
AsymNP(3,2) =  kr2np*besselj(1,kr2np*r3S)/(kwS*YwS);
AsymNP(3,3) =  kr2np*bessely(1,kr2np*r3S)/(kwS*YwS);
dAsymNP = diff(AsymNP, kzS);

Build_A_No_Pod_num  = matlabFunction(AsymNP, 'Vars',[kzS k0S kwS r2S r3S Y0S YwS]);
dBuild_A_No_Pod_num = matlabFunction(dAsymNP,'Vars', [kzS k0S kwS r2S r3S Y0S YwS]);

%% Initialising the Matrices
freq = 1:1:5000;
TLpod = zeros(3,length(freq));
TLfull = zeros(3,length(freq));
TLpodr2 = zeros(3,length(freq));
TLnoPod = zeros(3,length(freq));

for k = 1:3
    %% Initial Guess
    kzguess = 0.98*(2*pi*1/c0) - 0.02i*(2*pi*1/c0);
    kzguess5 = kzguess;
    kzguess6 = kzguess;
    kzguessOld = kzguess;
    kzguess5Old = kzguess;
    kzguess6Old = kzguess;

    for n = 1:length(freq)

        %% For Fig. 3 & 4  (r1 = r1_values(k), r2 = 0.100)

        r1 = r1values(k);
        r2 = 0.100;

        %% Defining Properties

        k0 = 2*pi*freq(n)/c0;
        Y0=rho0*c0;
        lambda = c0/freq(n);

        A = E*lambda/(rho0*c0); % Page 238 of ACOUSTICS OF DUCTS AND MUFFLERS

        Yw = rho0*c0*(1 + 0.0485*A^(0.754) - 1i*0.087*A^(0.732));

        kw = k0*(1 + 0.0978*A^(0.6929) - 1i*0.189*A^(0.6185));

        %% Newton-Raphson Scheme 
        if n <= 2
            pred = kzguess;
        else
            pred = 2*kzguess - kzguessOld;
        end

        candidates = [];
        for oi = 1:length(offsets)

            kz = pred + offsets(oi)*max(abs(pred),1)*0.1;
            tol = 1e-11;
            maxIter = 150;
            iter = 0;
            converged = false;

            while true

                iter = iter + 1;

                Amat = Build_A_num(kz, k0, kw, r1, r2, r3, Y0, Yw);
                dAmat = dBuild_A_num(kz, k0, kw, r1, r2, r3, Y0, Yw);

                F = det(Amat);

                dF = 0;
                for row = 1:5
                    Atemp = Amat;
                    Atemp(row,:) = dAmat(row,:);
                    dF = dF + det(Atemp);
                end

                kznew = kz - (F/dF);

                if abs(kznew - kz) < tol * max(abs(kz),1)
                    kz = kznew;
                    converged = true;
                    break;
                end

                kz = kznew;

                if iter >= maxIter
                    break;
                end
            end

            if converged
                resid = abs(det(Build_A_num(kz, k0, kw, r1, r2, r3, Y0, Yw)));
                if resid < 1e-6 && real(kz) > 0
                    candidates(end+1) = kz; 
                end
            end
        end

        if isempty(candidates)
            fprintf('Frequency = %d (Fig.3&4, k=%d): no clean root found, using predicted continuation\n', freq(n), k);
            kz = pred;
        else
            [~, idx] = min(abs(imag(candidates)));
            kz = candidates(idx);
        end

        kzguessOld = kzguess;
        kzguess = kz;

        %% Defining Transfer Matrices

        %% Transfer Matrix for Pod
        % Defined in Analysis and design of pod silencers
        Y0 = rho0*c0/(pi*(r2^2-r1^2));

        Y = Y0*k0/kz;

        Tpod = [cos(kz*lp) 1j*Y*sin(kz*lp) ; 1j/Y*sin(kz*lp) cos(kz*lp)];

        %% Transfer Matrix for Cone Expansion
        % Page 79 of ACOUSTICS OF DUCTS AND MUFFLERS
        d0=2*ru;
        d1=2*r3;
        L=lec;
        Y0 = rho0*c0/(pi*d0^2/4);
        C = cos(k0*L);
        S = sin(k0*L);
        m = (d1-d0)/L;
        T11 = (d1/d0)*C - (m/(k0*d0))*S;
        T12 = 1i*Y0*(d0/d1)*S;
        T21 = 1i/Y0*(1 + m^2/(k0^2*d0*d1))*S -1i*(m/(k0*d0*Y0))*(1-d0/d1)*C;
        T22 = (m/(k0*d1))*S + (d0/d1)*C;
        Texp = [T11 T12;T21 T22];

        %% Transfer Matrix for Cone Contraction
        % Page 79 of ACOUSTICS OF DUCTS AND MUFFLERS
        d0=2*r3;
        d1=2*rd;
        L=lcc;
        C = cos(k0*L);
        S = sin(k0*L);
        Y0 = rho0*c0/(pi*d0^2/4);
        m = (d1-d0)/L;
        T11 = (d1/d0)*C - (m/(k0*d0))*S;
        T12 = 1i*Y0*(d0/d1)*S;
        T21 = 1i/Y0*(1 + m^2/(k0^2*d0*d1))*S -1i*(m/(k0*d0*Y0))*(1-d0/d1)*C;
        T22 = (m/(k0*d1))*S + (d0/d1)*C;
        Tcon = [T11 T12;T21 T22];

        %% Overall Transfer Matrix
        T = Texp*Tpod*Tcon;

        %% Calculation of Transmission Loss
        Am = T(1,1);
        Bm = T(1,2);
        Cm = T(2,1);
        Dm = T(2,2);
        Fm = Tpod(1,1);
        Gm = Tpod(1,2);
        Hm = Tpod(2,1);
        Im = Tpod(2,2);

        Zf= rho0*c0/(pi*ru^2);
        Zp= rho0*c0/(pi*(r2^2-r1^2));

        TLfull(k,n) = 20*log10(abs((Am + Bm/Zf + Cm*Zf + Dm)/2));
        TLpod(k,n) = 20*log10(abs((Fm + Gm/Zp + Hm*Zp + Im)/2));


        %% For Fig. 5  (r1 = 0.050, r2 = r2_values(k))
        r2 = r2values5(k);
        r1 = 0.050;

        %% Newton-Raphson 
        if n <= 2
            pred = kzguess5;
        else
            pred = 2*kzguess5 - kzguess5Old;
        end

        candidates = [];
        for oi = 1:length(offsets)

            kz = pred + offsets(oi)*max(abs(pred),1)*0.1;
            tol = 1e-11;
            maxIter = 150;
            iter = 0;
            converged = false;

            while true

                iter = iter + 1;

                Amat  = Build_A_num(kz,k0,kw,r1,r2,r3,rho0*c0,Yw);
                dAmat = dBuild_A_num(kz,k0,kw,r1,r2,r3,rho0*c0,Yw);

                F = det(Amat);

                dF = 0;
                for row = 1:5
                    Atemp = Amat;
                    Atemp(row,:) = dAmat(row,:);
                    dF = dF + det(Atemp);
                end

                kznew = kz - F/dF;

                if abs(kznew-kz) < tol*max(abs(kz),1)
                    kz = kznew;
                    converged = true;
                    break
                end

                kz = kznew;

                if iter>=maxIter
                    break
                end
            end

            if converged
                resid = abs(det(Build_A_num(kz,k0,kw,r1,r2,r3,rho0*c0,Yw)));
                if resid < 1e-6 && real(kz) > 0
                    candidates(end+1) = kz; 
                end
            end
        end

        if isempty(candidates)
            fprintf('Frequency = %d (Fig.5, k=%d): no clean root found, using predicted continuation\n', freq(n), k);
            kz = pred;
        else
            [~, idx] = min(abs(imag(candidates)));
            kz = candidates(idx);
        end

        kzguess5Old = kzguess5;
        kzguess5 = kz;

        %% Pod transfer matrix

        Zp = rho0*c0/(pi*(r2^2-r1^2));
        Y = Zp*k0/kz;

        Tpod = [cos(kz*lp) 1i*Y*sin(kz*lp); 1i/Y*sin(kz*lp) cos(kz*lp)];

        %% Transmission Loss

        Am = Tpod(1,1);
        Bm = Tpod(1,2);
        Cm = Tpod(2,1);
        Dm = Tpod(2,2);

        TLpodr2(k,n)=20*log10(abs((Am+Bm/Zp+Cm*Zp+Dm)/2));

        %% For Fig. 6  (r1 = 0, r2 = r2_values(k))
        r1=0;
        r2 = r2values6(k);

        %% Newton-Raphson 
        if n <= 2
            pred = kzguess6;
        else
            pred = 2*kzguess6 - kzguess6Old;
        end

        candidates = [];
        for oi = 1:length(offsets)

            kz = pred + offsets(oi)*max(abs(pred),1)*0.1;
            tol = 1e-11;
            maxIter = 150;
            iter = 0;
            converged = false;

            while true

                iter = iter + 1;

                Amat  = Build_A_No_Pod_num(kz,k0,kw,r2,r3,rho0*c0,Yw);
                dAmat = dBuild_A_No_Pod_num(kz,k0,kw,r2,r3,rho0*c0,Yw);

                F = det(Amat);

                dF = 0;
                for row = 1:3
                    Atemp = Amat;
                    Atemp(row,:) = dAmat(row,:);
                    dF = dF + det(Atemp);
                end

                kznew = kz - F/dF;

                if abs(kznew-kz) < tol*max(abs(kz),1)
                    kz = kznew;
                    converged = true;
                    break
                end

                kz = kznew;

                if iter >= maxIter
                    break
                end
            end

            if converged
                resid = abs(det(Build_A_No_Pod_num(kz,k0,kw,r2,r3,rho0*c0,Yw)));
                if resid < 1e-6 && real(kz) > 0
                    candidates(end+1) = kz;
                end
            end
        end

        if isempty(candidates)
            fprintf('Frequency = %d (Fig.6, k=%d): no clean root found, using predicted continuation\n', freq(n), k);
            kz = pred;
        else
            [~, idx] = min(abs(imag(candidates)));
            kz = candidates(idx);
        end

        kzguess6Old = kzguess6;
        kzguess6 = kz;

        %% Transfer matrix of absorptive section

        Zp = rho0*c0/(pi*r2^2);      % r1 = 0

        Y = Zp*k0/kz;

        Tpod = [ cos(kz*lp)  1i*Y*sin(kz*lp); 1i/Y*sin(kz*lp)  cos(kz*lp)];

        %% Transmission loss

        A11 = Tpod(1,1);
        A12 = Tpod(1,2);
        A21 = Tpod(2,1);
        A22 = Tpod(2,2);

        TLnoPod(k,n) = 20*log10(abs((A11 + A12/Zp + A21*Zp + A22)/2));

    end
end

%% Plots

figure(3)
hold on
plot(freq,TLfull(1,:),'--')
plot(freq,TLpod(1,:),'-')
title('Comparison of TL of the (full) pod silencer with that of the absorptive section alone')
axis([0 5000 0 50])
xlabel('Frequency(Hz)')
ylabel('Transmission Loss (dB)')
legend('TL of the (full) pod silencer', 'TL of the absorptive section alone', 'Location', 'southeast')
hold off

styles = {'--','-','.-'};

figure(4)
hold on
for k = 1:3
    plot(freq, TLpod(k,:), styles{k});
end
xlabel('Frequency (Hz)');
ylabel('Transmission loss (dB)');
legend('r1 = 25 mm', 'r1 = 50 mm', 'r1 = 75 mm', 'Location', 'best');
title('Effect of the pod silencer on TL of the absorptive section of the silencer');
axis([0 5000 0 100])
hold off

figure(5)
hold on
for k = 1:3
    plot(freq, TLpodr2(k,:), styles{k});
end
xlabel('Frequency (Hz)');
ylabel('Transmission loss (dB)');
legend('r2 = 75 mm', 'r2 = 100 mm', 'r2 = 125 mm', 'Location', 'best');
title('Effect of radius r2 on TL of the absorptive section of the silencer');
axis([0 5000 0 100])
hold off

figure(6)
hold on
for k = 1:3
    plot(freq, TLnoPod(k,:), styles{k});
end
xlabel('Frequency (Hz)');
ylabel('Transmission loss (dB)');
legend('r2 = 50 mm', 'r2 = 75 mm', 'r2 = 100 mm', 'Location', 'best');
title('Effect of radius r2 on TL of the absorptive section of the silencer without pod (r1=0)');
axis([0 5000 0 100])
hold off