%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Weak quadratic chirp detection through the adaptive Liènard oscillator
%    array in chaotic intermittence regime, and time-frequency
%                              representation
%
%                          January, 10th 2023
%       Author: Maribel Tello-Bello / Pedro Pancóatl-Bortolotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc;
clear;
format long;
%close all;
load Data/LinearChirp.mat signal
signal = interp1([-1, 1], [-0.004, 0.004], signal);
%% Runge-Kutta Parameters
h=1.5e-4;  %Step size
a=0;      %Initial Time
b=1; %Final Time for Integration

xl=0;
yl=0;

%% Chirp Parameters

fini=1000;
Wini=fini*(2*pi);
ffin=7000;
Wfin=ffin*(2*pi);
K=(ffin-fini)/(b-a);
A = 0.01; % Amplitud
W=Wini;
tic

%% Window Size 
tv= sqrt(0.14/K);
TV=a:tv:b;  % Array for overlaping
global ts
ts=h;        %Step size
NVent=length(TV);

%% Liènard-type System Parameters
F=0.7;
alfa=2;
beta=1.8654;
epsilon=2.5263;

%% Delta Omega distribution
DW1=0.92;
DW2=0.96;
DW3=1.04;
DW4=1.08;

%% SNR and chirp amplitude 
dB=-15;
Ampli=75;
F1=Ampli*abs(max(signal)); %Set amplitude F1=0.3
sigma=sqrt( (F1^2)/(2*10^(dB/10)) ); %Standard deviation for adding noise
Umbral=1.9; %Upper threshold

MonteCarlo=50; %Number of Monte Carlo iterations
tic
for MC=1:1:MonteCarlo
    
    ruido = normrnd(0,sigma,[1,length(signal)]); %Noise
    CHIRP= Ampli*signal + ruido; %Chirp with additive noise
    
    freq(1)=fini; %set initial frequency
    Tfreq(1)=a;
    W=Wini;
    
    for Z=2:1:NVent-1
            
            ini=Z-1;           
            fin=Z+1;
            SampleIni=round(TV(ini)/h)+1;
            SampleFin=round(TV(fin)/h);
            Signal=CHIRP(SampleIni:SampleFin); %Signal Windowing 
            tc=TV(1):h:TV(2)/20; 
            
            %%Oscillator 1    
            W1=DW1*W;
            ODEl=@(t1,x1,y1) (W1)*y1;
            ODE2=@(t1,x1,y1)  (W1)*((-epsilon*(-1+alfa*tan(x1)^2)*y1)-beta*sinh(x1) + F*cos((W1)*t1) );           
            [t1,x1,y1]=RK4B( ODEl,ODE2,xl,yl,h,TV(ini),TV(fin),Signal,W1 );          
            fc1=cos(W1*tc);
            xc01=xcorr(fc1,x1); %cross correlation for improving the intermittent periods
            xc1=xc01(1:length(xc01)/2).^2;
            xc1=(2*xc1/max(xc1));
            
            %%Oscillator 2  
            W2=DW2*W;
            ODE3=@(t2,x2,y2) (W2)*y2;
            ODE4=@(t2,x2,y2)  (W2)*((-epsilon*(-1+alfa*tan(x2)^2)*y2)-beta*sinh(x2) + F*cos((W2)*t2) ); 
            [t2,x2,y2]=RK4B( ODE3,ODE4,xl,yl,h,TV(ini),TV(fin),Signal,W2 );   
            fc2=cos(W2*tc);
            xc02=xcorr(fc2,x2);
            xc2=xc02(1:length(xc02)/2).^2;
            xc2=(2*xc2/max(xc2));
          
            
            %%Oscillator 3
            W3=W;
            ODE5=@(t3,x3,y3) (W3)*y3;
            ODE6=@(t3,x3,y3)  (W3)*((-epsilon*(-1+alfa*tan(x3)^2)*y3)-beta*sinh(x3) + F*cos((W3)*t3) ); %Chirp version 1 !!!BUENO PARA DUFFING
            [t3,x3,y3]=RK4B( ODE5,ODE6,xl,yl,h,TV(ini),TV(fin),+Signal,W3);  
            fc3=cos(W3*tc);
            xc03=xcorr(fc3,x3);
            xc3=xc03(1:length(xc03)/2).^2;  
            xc3=(2*xc3/max(xc3));
            
            
            %%Oscillator 4
            W4=DW3*W;
            ODE7=@(t4,x4,y4) (W4)*y4;
            ODE8=@(t4,x4,y4)  (W4)*((-epsilon*(-1+alfa*tan(x4)^2)*y4)-beta*sinh(x4) + F*cos((W4)*t4) ); %Chirp version 1 !!!BUENO PARA DUFFING
            [t4,x4,y4]=RK4B( ODE7,ODE8,xl,yl,h,TV(ini),TV(fin),Signal,(W4)  ); 
            fc4=cos(W4*tc);
            xc04=xcorr(fc4,x4);
            xc4=xc04(1:length(xc04)/2).^2;
            xc4=(2*xc4/max(xc4));
            
            %%Oscillator 5
            W5=DW4*W;
            ODE9=@(t5,x5,y5) (W5)*y5;
            ODE10=@(t5,x5,y5)  (W5)*((-epsilon*(-1+alfa*tan(x5)^2)*y5)-beta*sinh(x5) + F*cos((W5)*t5) ); %Chirp version 1 !!!BUENO PARA DUFFING                
            [t5,x5,y5]=RK4B( ODE9,ODE10,xl,yl,h,TV(ini),TV(fin),Signal,(W5) ); 
            fc5=cos(W5*tc);
            xc05=xcorr(fc5,x5);
            xc5=xc05(1:length(xc05)/2).^2;  
            xc5=(2*xc5/max(xc5));
            
            txc=t1(1:length(xc1));
            
            %Detect intermittent periods

            [Pic1,delT1]=findpeaks(xc1,'MinPeakDistance', 0.6*(2*pi)/(h*(W-W1)),'MinPeakHeight',Umbral );
            [Pic2,delT2]=findpeaks(xc2,'MinPeakDistance', 0.6*(2*pi)/(h*(W-W2)),'MinPeakHeight',Umbral );
            [Pic3,delT3]=findpeaks(xc3,'MinPeakDistance', length(xc3)/1.02 ,'MinPeakHeight',Umbral);
            [Pic4,delT4]=findpeaks(xc4,'MinPeakDistance', 0.6*(2*pi)/(h*(W4-W)),'MinPeakHeight',Umbral);
            [Pic5,delT5]=findpeaks(xc5,'MinPeakDistance', 0.6*(2*pi)/(h*(W5-W)),'MinPeakHeight',Umbral );


            [ deltaT1 ] = Delta4( delT1 );
            [ deltaT2 ] = Delta4( delT2 );
            [ deltaT3 ] = Delta4( delT3 );
            [ deltaT4 ] = Delta4( delT4 );
            [ deltaT5 ] = Delta4( delT5 );
           
            mediana1=mean(deltaT1);
            mediana2=mean(deltaT2);
            mediana3=mean(deltaT3);
            mediana4=mean(deltaT4);
            mediana5=mean(deltaT5);
                        
     VecDeltasP=[mediana1 mediana2 mediana3 mediana4 mediana5];
     VecDeltasP(VecDeltasP==Inf)=0;
     [elem,pos] = find (VecDeltasP==0);

    if pos == 1
            VecDeltasP(pos)=[W*DW1];
            VecDeltasP(pos+1)=[(W*DW2) - (2*pi/VecDeltasP(2))];
            VecDeltasP(pos+2)=[(W) - (2*pi/VecDeltasP(3))];
            VecDeltasP(pos+3)=[(W*DW3) - (2*pi/VecDeltasP(4))];
            VecDeltasP(pos+4)=[(W*DW4) - (2*pi/VecDeltasP(5))];
            FrecWT1=mean(VecDeltasP);
    elseif  pos == 2
            VecDeltasP(pos-1)=[W*DW1 + (2*pi/VecDeltasP(1))];
            VecDeltasP(pos)=[(W*DW2)];
            VecDeltasP(pos+1)=[(W) - (2*pi/VecDeltasP(3))];
            VecDeltasP(pos+2)=[(W*DW3) - (2*pi/VecDeltasP(4))];
            VecDeltasP(pos+3)=[(W*DW4) - (2*pi/VecDeltasP(5))];
            FrecWT1=mean(VecDeltasP);
    elseif  pos == 3
            VecDeltasP(pos-2)=[(W*DW1) + (2*pi/VecDeltasP(1))];
            VecDeltasP(pos-1)=[(W*DW2) + (2*pi/VecDeltasP(2))];
            VecDeltasP(pos)=[W];
            VecDeltasP(pos+1)=[(W*DW3) - (2*pi/VecDeltasP(4))];
            VecDeltasP(pos+2)=[(W*DW4) - (2*pi/VecDeltasP(5))];
            FrecWT1=mean(VecDeltasP);
    elseif  pos == 4
            VecDeltasP(pos-3)=[(W*DW1) + (2*pi/VecDeltasP(1))];
            VecDeltasP(pos-2)=[(W*DW2) + (2*pi/VecDeltasP(2))];
            VecDeltasP(pos-1)=[W + (2*pi/VecDeltasP(3))];
            VecDeltasP(pos)=[(W*DW3)];
            VecDeltasP(pos+1)=[(W*DW4) - (2*pi/VecDeltasP(5))];
            FrecWT1=mean(VecDeltasP);
    elseif  pos == 5
            VecDeltasP(pos-4)=[(W*DW1) + (2*pi/VecDeltasP(1))];
            VecDeltasP(pos-3)=[(W*DW2) + (2*pi/VecDeltasP(2))];
            VecDeltasP(pos-2)=[W + (2*pi/VecDeltasP(3))];
            VecDeltasP(pos-1)=[(W*DW3) + (2*pi/VecDeltasP(4))];
            VecDeltasP(pos)=[(W*DW4)];
            FrecWT1=mean(VecDeltasP);
    elseif length(pos)>1 
            disp('Hay mas de un cero');
            %FrecWT1=wact.*(1.08);
            VecDeltasP(VecDeltasP==0)=[];
            VecDeltasP=W + 2.*pi./VecDeltasP;
            FrecWT1=mean(VecDeltasP);
    else
            VecDeltasP(1)=[(W*DW1) + (2*pi/VecDeltasP(1))];
            VecDeltasP(2)=[(W*DW2) + (2*pi/VecDeltasP(2))];
            VecDeltasP(3)=[W + (2*pi/VecDeltasP(3))];
            VecDeltasP(4)=[(W*DW3) - (2*pi/VecDeltasP(4))];
            VecDeltasP(5)=[(W*DW4)- (2*pi/VecDeltasP(5))];
            FrecWT1=mean(VecDeltasP);
    end                        
     %}     
            
            freq(MC,Z)=FrecWT1/(2*pi);
            W=FrecWT1; 
                       
            Tfreq(Z)=(TV(fin)+TV(ini))/2;              
                                        
    end
end

frec2=transpose(freq);

for c=1:1:(Z)    
    frecMonteCarlo(c)=mean(frec2(c,:));   
end

frecMonteCarlo(1)=fini;
Tfreq(1)=a;

toc    
figure();
hold on;
for B=1:1:MC
    stairs(Tfreq,freq(B,:),'g');
end
stairs(Tfreq,frecMonteCarlo,'r',LineWidth=2)
ChirpIdeal=(((ffin-fini)/b).*Tfreq+fini);
plot(Tfreq,ChirpIdeal,'-k');
xlabel('t/s');
ylabel('Frequency/Hz');
title(['TF Representation: linear chirp at ',num2str(dB), ' dB' ] )
axis([0 1 0 1e4]);
grid();
savefig('Plots/ALI/linear_-15db.fig')      

Error=100*(abs(frecMonteCarlo-ChirpIdeal)./ChirpIdeal); 
mse_media=(mean(Error));
disp('Percentage of mean absolute error %: ' ); disp(mse_media);

% run('QuadraticChirpDetection.m')
