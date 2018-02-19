%****************************************************************************************
%  This program is a demostration of digital filter, including 
%  Butterworth IIR, FIR, median filter, Wiener filter, adaptive filter and Wavelet filter
%
%  Created on: Dec 13, 2015
%  Author: Adam
%
%***************************************************************************************

% Generate two signals: Mix_Signal_1 and Mix_Signal_2
Fs = 1000;    % sample rates
N  = 1000;    % samples
n  = 0:N-1;
t  = 0:1/Fs:1-1/Fs;  
Signal_Original_1 =sin(2*pi*10*t)+sin(2*pi*20*t)+sin(2*pi*30*t); 
Noise_White_1    = [0.3*randn(1,500), rand(1,500)];        % Gauss white noise and uniform distribution white noise
Mix_Signal_1   = Signal_Original_1 + Noise_White_1;        % mixed signal 

Signal_Original_2  =  [zeros(1,100), 20*ones(1,20), -2*ones(1,30), 5*ones(1,80), -5*ones(1,30), 9*ones(1,140), -4*ones(1,40), 3*ones(1,220), 12*ones(1,100), 5*ones(1,20), 25*ones(1,30), 7 *ones(1,190)]; 
Noise_White_2     =  0.5*randn(1,1000);                    % Gauss white noise
Mix_Signal_2        =  Signal_Original_2 + Noise_White_2;  % mixed signal 

%****************************************************************************************
%  
%                Butterworth lowpass filter
%
%***************************************************************************************

% for Mix_Signal_1
figure(1);
Wc=2*50/Fs;         %cut off frequence: 50Hz
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal_1);

subplot(4,1,1);     % before filter                
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('Original Signal  1');

subplot(4,1,2);           
plot(Signal_Filter);
axis([0,1000,-4,4]);
title('after Butterworth LPF');

% for Mix_Signal_2
Wc=2*100/Fs;        %cut off frequence: 100Hz
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal_2);

subplot(4,1,3);             
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('Original Signal  2');

subplot(4,1,4);          
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('after Butterworth LPF');

%****************************************************************************************
%  
%                FIR LPF
%
%***************************************************************************************
figure(2);
F   =  [0:0.05:0.95]; 
A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
b  =  firls(20,F,A);
Signal_Filter = filter(b,1,Mix_Signal_1);

subplot(4,1,1);     % before filter                
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('Original Signal  1');

subplot(4,1,2);       
plot(Signal_Filter);
axis([0,1000,-5,5]);
title('after FIR LPF');

F   =  [0:0.05:0.95]; 
A  =  [1    1      1     1     1    0      0    0     0    0     0     0     0     0     0     0    0   0   0   0] ;
b  =  firls(20,F,A);
Signal_Filter = filter(b,1,Mix_Signal_2);
subplot(4,1,3);                                                    
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('Original Signal  2');

subplot(4,1,4);                     
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('after FIR LPF');

%****************************************************************************************
%  
%                moving average filter
%
%***************************************************************************************

figure(3);
b  =  [1 1 1 1 1 1]/6;
Signal_Filter = filter(b,1,Mix_Signal_1);

subplot(4,1,1);             
plot(Mix_Signal_1);
axis([0,1000,-4,4]);
title('Original Signal 1');

subplot(4,1,2);         
plot(Signal_Filter);
axis([0,1000,-4,4]);
title('moving average filter');

%混合信号 Mix_Signal_2  移动平均滤波
b  =  [1 1 1 1 1 1]/6;
Signal_Filter = filter(b,1,Mix_Signal_2);
subplot(4,1,3);                           
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('Original Signal 2');

subplot(4,1,4);                        
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('moving average filter');

%****************************************************************************************
%  
%                median filter
%
%***************************************************************************************
figure(4);
Signal_Filter=medfilt1(Mix_Signal_1,10);

subplot(4,1,1);                  
plot(Mix_Signal_1);
axis([0,1000,-5,5]);
title('Original Signal 1');

subplot(4,1,2);   
plot(Signal_Filter);
axis([0,1000,-5,5]);
title(' median filter');

%混合信号 Mix_Signal_2  中值滤波
Signal_Filter=medfilt1(Mix_Signal_2,10);
subplot(4,1,3);                     
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('Original Signal 2');

subplot(4,1,4);         
plot(Signal_Filter);
axis([0,1000,-10,30]);
title(' median filter');

%****************************************************************************************
%  
%               Wiener filter
%
%***************************************************************************************
figure(5);
Rxx=xcorr(Mix_Signal_1,Mix_Signal_1);              % autocorrelation function
M=100;                                             % order
for i=1:M                                          % autocorrelation matrix
    for j=1:M
        rxx(i,j)=Rxx(abs(j-i)+N);
    end
end
Rxy=xcorr(Mix_Signal_1,Signal_Original_1);         % Cross correlation function
for i=1:M
    rxy(i)=Rxy(i+N-1);
end                                                % Cross correlation vector
h = inv(rxx)*rxy';                                 % coefficient of Wiener filter
Signal_Filter=filter(h,1, Mix_Signal_1);           

subplot(4,1,1);                     
plot(Mix_Signal_1);
axis([0,1000,-5,5]);
title('Original Signal 1');

subplot(4,1,2);              
plot(Signal_Filter);
axis([0,1000,-5,5]);
title('Wiener filter');

Rxx=xcorr(Mix_Signal_2,Mix_Signal_2);  
M=500;                             
for i=1:M                            
    for j=1:M
        rxx(i,j)=Rxx(abs(j-i)+N);
    end
end
Rxy=xcorr(Mix_Signal_2,Signal_Original_2);  
for i=1:M
    rxy(i)=Rxy(i+N-1);
end                                          
h=inv(rxx)*rxy';                                      
Signal_Filter=filter(h,1, Mix_Signal_2);            

subplot(4,1,3);                                                    
plot(Mix_Signal_2);
axis([0,1000,-10,30]);
title('Original Signal 2');

subplot(4,1,4);                                       
plot(Signal_Filter);
axis([0,1000,-10,30]);
title('Wiener filter');

%****************************************************************************************
%  
%                adaptive filter
%
%***************************************************************************************

figure(6);
N=1000;                                             % samples
k=100;                                              % order of LMS filter
u=0.001;                                            % Step Length Factor

% initial
yn_1=zeros(1,N);                                    % output signal
yn_1(1:k)=Mix_Signal_1(1:k);                        % 
w=zeros(1,k);                                       % initial valve
e=zeros(1,N);                                       % error

% LMS
for i=(k+1):N
        XN=Mix_Signal_1((i-k+1):(i));
        yn_1(i)=w*XN';
        e(i)=Signal_Original_1(i)-yn_1(i);
        w=w+2*u*e(i)*XN;
end

subplot(4,1,1);
plot(Mix_Signal_1);                 
axis([k+1,1000,-4,4]);
title('Original Signal 1');

subplot(4,1,2);
plot(yn_1);                                        
axis([k+1,1000,-4,4]);
title('adaptive filter');

N=1000;                                           
k=500;                                          
u=0.000011;                                       

yn_1=zeros(1,N);                                 
yn_1(1:k)=Mix_Signal_2(1:k);                  
w=zeros(1,k);                                 
e=zeros(1,N);                                   

for i=(k+1):N
        XN=Mix_Signal_2((i-k+1):(i));
        yn_1(i)=w*XN';
        e(i)=Signal_Original_2(i)-yn_1(i);
        w=w+2*u*e(i)*XN;
end

subplot(4,1,3);
plot(Mix_Signal_2);                         
axis([k+1,1000,-10,30]);
title('Original Signal ');

subplot(4,1,4);
plot(yn_1);                                        
axis([k+1,1000,-10,30]);
title('adaptive filter');

%****************************************************************************************
%  
%                wavelet filter
%
%***************************************************************************************

figure(7);
subplot(4,1,1);
plot(Mix_Signal_1);                 
axis([0,1000,-5,5]);
title('Original Signal ');

subplot(4,1,2);
[xd,cxd,lxd] = wden(Mix_Signal_1,'sqtwolog','s','one',2,'db3');
plot(xd);                    
axis([0,1000,-5,5]);
title('wavelet filter');

subplot(4,1,3);
plot(Mix_Signal_2);                          
axis([0,1000,-10,30]);
title('Original Signal ');

subplot(4,1,4);
[xd,cxd,lxd] = wden(Mix_Signal_2,'sqtwolog','h','sln',3,'db3');
plot(xd);                          
axis([0,1000,-10,30]);
title('wavelet filter');