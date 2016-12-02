function [ avg, max_coeff ,sigma,min1,dc,firstharmonic ] = fourier_desc( binaryImage )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Elliptical description of boundary curve is used( points are represented
%as complex coordinates)
% zero component is the center of gravity
% first harmonic is the radius/area
% higher harmonics contain finer detail but are more susceptible to noise
% calculating mean, maximum and variance - should be close for similar
% images?

[x y] = find(binaryImage == 1); % White hot?
 
% No of samples
 m = length([x y]);

 %Complex Contour Representation
 c =[x + i*y];
 
 %Discrete Fourier Tranform
 F = fft(c);
 
 %Normalize with zero component for scale invariance
 F_norm = F/abs(F(1));
%  figure;
% plot(x,y,'o');
% xlabel('Real Axis')
% ylabel('Imaginary Axis')
% figure;
% plot(x,log(real(F_norm+1)),'r');
% figure;
% plot(log(imag(F_norm+1)),y,'b');
 %Descriptors
 avg= mean(F_norm);
 max_coeff = F_norm(end);
 %max1 = max(F_norm)
 min1 = min(F_norm);
 dc= F_norm(1);
 firstharmonic = F_norm(2);
 sigma = sqrt(sum((F_norm-avg).^2));

end

