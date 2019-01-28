function y = rest_filter(TR,filter,x) 

%Usage: "TR" is the repeition time, "filter" is the band to keep i.e.
%filter = [0.01, 0.1], "x" is the timeseries that will be filtered, and "y" is
%the filtered output

%Find the Fourier transform of each column of the timeseries
y = fft(x, [ ], 1);

%Find 1 less than the number of rows that are in the timeseries data
f = (0:size(x,1)-1); 

%Create an array of the smallest values from f
f = min(f, size(x,1)-f);

%Convolution: Find all the values from f that are less than the smaller
%value in the filter range and multiplies it by the time series.
idx = find(f<filter(1)*(TR*size(x,1))); 

%Find every value greater than 1 from the filtered data
idx = idx(idx>1); 

%Set a portion of the Fourier transformed data equal to 0
y(idx,:)=0; 

%Takes the inverse Fourier transform of the filtered data and extracts only
%the real values and now is in the form that can be used for the rest of
%MATLAB
y=real(ifft(y, [ ], 1)); 
end