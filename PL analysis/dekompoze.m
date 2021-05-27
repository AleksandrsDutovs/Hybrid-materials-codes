function dekompoze(datfilename,usec)
  pkg load optim;
  
  w=50;
  datfilename;                                        #here we define name of file (for php automatization script)
				#usec
  sec=usec/1000000;                                   #here we define PL intensity units - photons per second
  
  
  
  indata=load('-ascii',datfilename);                  
  
  imin=21;
  imax=1678;
  indata_trim=indata(imin:imax,:);
  i400=247;
  dark=indata(imin:i400,2);                           #here we define (and in line 32 remove) dark spectra (if spectrometer software is not able to do it)
 mean_dark=mean(dark);


  
  nm_original=indata_trim(:,1);                       #from line 22 to 32 we remove signal noises, pikes and other processing artifacts
  intens_original=(indata_trim(:,2)-mean_dark)/sec;
  
  trim_w=5;
  rempleaks=removepeaks(indata_trim,trim_w);
  for q=1:20
    rempleaks=removepeaks(rempleaks,trim_w);
  endfor
  
  nm=rempleaks(:,1);
  intens=(rempleaks(:,2)-mean_dark)/sec;
  

  eV=1240./nm;

  eV_interp=1.3:0.01:3.2;
  intens_interp = interp1 (eV, intens, eV_interp, "spline");
  
  fitfilename = ['peaks_' strrep(datfilename, '.dat', '.fit')];   
  p=[ 2.07      0.27    477 ...                                   # here we define fitting gaussian parameters - central energy, HWHM, intensity
      2.45      0.20    1047 ...                                  # It is possible to fit with 5 gaussians
      #2.70      0.14    883 ...
      2.80      0.11    839 ... 
      #];
  
   #p = load ('-ascii', fitfilename);                             #to achieve better fitting, after first run comment lines with p=[] and uncommet this line

  [F, P2]=leasqr(eV_interp,intens_interp,p,@gauss3);              #F is sum of gaussians and we fit it to experimental spectrum with least square method
  P2'
  
                                                                   #then we plot everytning and save data
  p1 = [P2(1) P2(2) P2(3)];
  G1 = gauss1(eV_interp,p1);
 # plot(eV,G1);
  p2 = [P2(4) P2(5) P2(6)];
  G2 = gauss1(eV_interp,p2);
#  plot(eV,G2);
  p3 = [P2(7) P2(8) P2(9)];
  G3 = gauss1(eV_interp,p3);
#  plot(eV,G3);
  #p4 = [P2(10) P2(11) P2(12)];
  #G4 = gauss1(eV_interp,p4);
 # plot(eV,G4);
  #p5 = [P2(13) P2(14) P2(15)];
  #G5 = gauss1(eV_interp,p5);
  
    save ("-ascii",  fitfilename, "P2");

#  plot(nm_original,intens_original,nm,intens);
  plot(eV,intens, "linewidth", 3,...
       eV_interp,G1, "linewidth", 2, "linestyle", "-.", "r",...
       eV_interp,G2, "linewidth", 2, "linestyle", "-.", "g",...
       eV_interp,G3, "linewidth", 2, "linestyle", "-.", "b",...
       #eV_interp,G4, "linewidth", 2, "linestyle", "-.",...
       #eV_interp,G5, "linewidth", 1, "linestyle", "-.",...
        eV_interp,F, "linewidth", 2, "k");


  grid on;  
				  axis ([1.2 3.3 0 600 ]);
 # xlabel ("wavelength (nm)");
  set(gca, "linewidth", 3, "fontsize", 20)
  xlabel ("E (eV)", "fontsize", 18);
  ylabel ("Intensity (counts/s)", "fontsize", 18);
  legend ({"experimental", "1. Gaussian", "2. Gaussian", "3. Gaussian","fitted line"},"fontsize", 14, "location", "northwest");
  #sltitle=strrep (datfilename, '_', '\_');
  #title (sltitle)
  saveas(1,['glud_' strrep(datfilename, '.dat', '.png')]);
  
				# pause();

endfunction;