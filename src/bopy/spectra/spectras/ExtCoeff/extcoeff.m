% extcoeff.m : compare extinction coefficient of HbO2, Hb, and H2O
% March 4, 2003 

clear all;

ini = 650;
fin = 1000;

%sfactor=input('scaling factor for hemoglobin(1000 for mM, 1000000 for uM):');
sfactor = 1000000; % scale hemoglobin extinction coefficient to be cm-1/uM

% From Steve Jacques' web site : Gratzer, Kollias
% Extinction Coefficients at 750 nm, 786 nm, 830 nm
% in [cm-1/M]
load hemoglobin.dat;
lamdaH = hemoglobin(:,1);

ind1 = find(lamdaH == ini);
ind2 = find(lamdaH == fin);

eHBO2 = hemoglobin(ind1:ind2,2) ./sfactor;   
eHB   = hemoglobin(ind1:ind2,3) ./sfactor;  

clear hemoglobin;
clear ind1; clear ind2;

% Mua from water (Segelstein) : interpolated data [cm-1]
load newwater.dat;
lamdaW = newwater(:,1);

ind1 = find(lamdaW == ini);
ind2 = find(lamdaW == fin);

Muawater = newwater(ind1:ind2,2);

eH2O = Muawater ./log(10)./55.4;    % e = mua / ln10 / C where C = 55.4 M for pure water
clear newwater;
clear ind1; clear ind2;

wave = [ini:2:fin];

figure,plot(wave,eHBO2,'r.-',wave,eHB,'m.-',wave,eH2O,'b.-');
title('Extinction coefficient: H2O is in cm-1/M whereas hemoglobins are in cm-1/\muM');
xlabel('wavelength (nm)');
ylabel('extinction coefficient');
legend('eHBO2 (cm-1/\muM)','eHB (cm-1/\muM)','eH2O (cm-1/M)');
