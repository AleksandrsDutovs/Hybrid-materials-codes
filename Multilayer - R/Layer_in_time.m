pkg load optim; #The package needed for the least square fitting to work

clear all #Clears all variables
close all #Closes all windows
#clc #Clears command window

format none; #The number of digits after the decimal separator does not change when the number goes from 3 to 4 digits before decimal point

lambda_min=480; #nanometers, Starting value of wavelength range
lambda_max=800; #nanometers, Last value of wavelength range
lambda_range=[lambda_min:lambda_max]; #Full wavelength range with the step of 1 nm written in a row
lambda_size=size(lambda_range);

ref_thickness=320; #nanometers, Approximate thickness of oxide layer
ref_multiplier=1; #Multiplier for the fitting data to match experimental data
#guessed_values=[ref_thickness ref_multiplier]; #variable used in dielectricfilm function

file_names=flipud(importdata('filelist.txt')); #Creates a variable with all the file names (contained in text file) to be used for data analysis
number_files=size(file_names,1); #Total number of files to be investigated

global layers=1; #Number of layers in the structure between air and substrate (e.g. 1 if there is only oxide and no barrier layers). Optical parameters must be defined for all these layers

time_interval=0.8; #seconds, Time interval between each measurement
data_time=0:time_interval:time_interval*(number_files-1); #The time of each measurement file

#Setting up reference spectrum
reference=dlmread('ref_spektrs.txt','\t',[17,0,3664,1]); #Reads the specified file from row 17 to 3664 and from column 0 to 1.
lambda=reference(:,1); #Wavelengths from the reference file.
intensity_reference=reference(:,2); #Intensity data of the reference spectrum.
intensity_reference_interp=interp1(lambda,intensity_reference,lambda_range,'spline'); #Light intensity of the reference is recalculated for desired wavelength range

#Defining the optical parameters of the layers stack
water_data=dlmread('Water.txt','\t'); #Loading the data for water from M. Daimon and A. Masumura.
                                      #Measurement of the refractive index of distilled water from the near-infrared region to the ultraviolet region,
                                      #Appl. Opt. 46, 3811-3820 (2007)
                                      #High performance liquid chromatography (HPLC) distilled water at 19.0 °C.
lambda_water=water_data(:,1)*1000; #Wavelength values from Daimon data
n_water=water_data(:,2); # Refractive index of water from Daimon data
n_air=interp1(lambda_water,n_water,lambda_range,'spline'); #Refractive index data is recalculated for desired wavelength range
k_air=0*ones(lambda_size); #Extinction coefficient of water (=0)
nk_air=n_air+j*k_air; #Complex refractive index of water

n_oxide=1.6*ones(lambda_size); #Refractive coefficient of oxide (=1.6), considered to be constant
k_oxide=zeros(lambda_size); #Extinction coefficient of oxide (=0)
nk_oxide=n_oxide+j*k_oxide; #Complex refractive index of oxide

Al=load('AlRakic.dat'); #Loading the data for aluminum from A. D. Raki?, A. B. Djurišic, J. M. Elazar, and M. L. Majewski.
                        #Optical properties of metallic films for vertical-cavity optoelectronic devices,
                        #http://dx.doi.org/10.1364/AO.37.005271\, Appl. Opt. 37, 5271-5283 (1998)
                        #Fit of the experimental data from several sources to Brendel-Bormann (BB) model
                        #First column is wavelength in micrometers, second - refractive index, third - extinction coefficient
lambda_Al=Al(:,1)*1000; #Wavelength values from Rakic data
n_AlRakic=Al(:,2); #Refractive index of aluminum from Rakic data
k_AlRakic=Al(:,3); #Extinction coefficient of aluminum from Rakic data
n_Al=interp1(lambda_Al,n_AlRakic,lambda_range,'spline'); #Refractive index data is recalculated for desired wavelength range
k_Al=interp1(lambda_Al,k_AlRakic,lambda_range,'spline'); #Extinction coefficient data is recalculated for desired wavelength range
nk_Al=n_Al+j*k_Al; #Complex refractive index of aluminum

global nk_all=[nk_air;nk_oxide;nk_Al]; #All refractive indices and extinction coefficients in one place: 1st row - air, 2nd row - oxide, 3rd row - aluminum
                                       #It is a global variable so that the function can use it without the need to define it within the function

function R=dielectricfilm(wavelength,variables) #Function to calculate reflected light ratio from the oxide
                                                #variables=[thickness multiplier]
                                                #Wavelength could be a single wavelength or a range for full spectra
  thickness=variables(1);                       #Thickness is oxide layer thickness in nm
  multiplier=variables(2);                      #Multiplier is a constant number for normalization
  R=zeros(size(wavelength));
  lambda_0=wavelength*1e-9; #Wavelength is converted to meters
  L=length(lambda_0); #The amount of wavelengths of interest
  l=thickness*1e-9; #Thickness is converted to meters
    
  for a=1:L #Function goes through all of wavelengths of interest
    global nk_all #Calls globally defined variable
    nk=nk_all(:,a); #Picks complex refractive index values for a certain wavelength
    global layers #Calls globally defined variable
    #Main body of function
    for b=1:layers #Function goes through all defined layers in the structure
      
      
      tau=2*nk(b)/(nk(b)+nk(b+1)); #Elementary transmission coefficient (from Orfanidis eq. (5.2.7))
      rho=(nk(b)-nk(b+1))/(nk(b)+nk(b+1)); #Elementary reflection coefficient (from Orfanidis eq. (5.2.7))
      T=[1 rho;rho 1]/tau; #Matching matrix constant (from Orfanidis eq. (5.2.3))
      
      k=nk(b+1)*2*pi/lambda_0(a); #Propagation wavenumber in oxide at specific wavelength in air (from Orfanidis Chapter 5.4 (bottom of page 163))
      P=[exp(i*k*l(b)) 0;0 exp(-i*k*l(b))]; #Propagation matrix constant (from Orfanidis eq. (5.1.11))
    endfor
    Ef=[1;(nk(2)-nk(3))/(nk(2)+nk(3))]; #Reflection from aluminum into oxide
    
    E0=T*P*Ef; #Backwards calculation of initial intensity, taking into account transmission through air-oxide interface, propagation in oxide, and reflection at oxide-aluminum interface
    R(a)=abs(E0(2)/E0(1))^2*multiplier; #Reflection coefficient at air-oxide interface (ratio of reflected light E0(2) with incident light E0(1))
  endfor
endfunction

R0=(dielectricfilm(lambda_range,[0.000000001 1]))'; #Some arbitrary reflection from the sample with a very very small oxide layer, which represents sample before anodisation

spectrum=zeros(lambda_size(2),number_files+1);
spectrum(:,1)=intensity_reference_interp;
R=zeros(lambda_size(2),number_files);
fittedR=zeros(lambda_size(2),number_files);
fittedNM=zeros(number_files+1,2);
fittedNM(1,1)=ref_thickness;
fittedNM(:,2)=ref_multiplier;

#Data analysis
for n=1:number_files
  filename=char(file_names(n,1)); #Determining filename to be read
  spectrum0=dlmread(filename,'\t',[17,1,3664,1]); #Reading spectrum file
  spectrum(:,n+1)=interp1(lambda,spectrum0,lambda_range,'spline'); #Interpolating data from the file
  R(:,n)=(spectrum(:,n+1)./intensity_reference_interp').*R0; #Calculating reflection of data in the file
  
  guessed_values=[fittedNM(n,1) fittedNM(n,2)]; #variable used in dielectricfilm function
  
  [R_fit,var_fit]=leasqr(lambda_range,(R(:,n))',guessed_values,@dielectricfilm); #Least square fitting of reflection data with guessed thickness and multiplier values
  fittedR(:,n)=R_fit; #Fitted reflection values
  fittedNM(n,:)=var_fit; #Obtained actual thickness (1st column) and multiplier (2nd column)
  fittedNM(n+1,1)=var_fit(1);
  
  #Plotting measured and fitted reflection and saving plot as jpg file
  figure(n);
  plot(lambda_range,R(:,n),lambda_range,fittedR(:,n));
  title(['Measured reflection and fitting ' filename]);
  xlabel('Wavelength (nm)');
  ylabel('Reflection (a.u.)');
  axis([lambda_min lambda_max 0.75 0.95]);
  saveas(figure(n),['ReflectionWithFit' mat2str(n) ', d ' mat2str(fittedNM(n,1)) ', m ' mat2str(fittedNM(n,2)) '.jpg']);
  close(figure(n));
end

fittedNM=flipud(fittedNM(1:end-1,1:2));
save reflection_experiment_all.txt R;
save reflection_fitted_all.txt fittedR;
save fitted_parameters(thickness_multiplier).txt fittedNM;

figure(6000)
plot(data_time,fittedNM(:,1),'.');
title('PAAO thickness dependence on anodization time');
xlabel('Time (s)');
ylabel('Thickness (nm)');
saveas(figure(6000),'Thickness_on_Time','jpg');