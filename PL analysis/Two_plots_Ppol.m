clear all;
close all;
pkg load optim;
lambda=[420:1:800]';   
 sec=5;
sec_1 = 1000000/10000;
#set reference
refdata=dlmread("Ppol_REF.txt", '\t', [18,0,3664,13]);
ref_lambda = refdata(:,1);
ref_i = refdata(:,2);
ref_i_sec =  ref_i/sec_1;
i_ref_i = interp1(ref_lambda, ref_i_sec, lambda, "spline");
#figure(10)
#plot (lambda, i_ref_i);
#/set reference

 rezmape='Prezulc/';


faili=[ 'L1_Ppol.txt';
        'L2_Ppol.txt';
        'L3_Ppol.txt';
        'L4_Ppol.txt';
        'L5_Ppol.txt';
        'L6_Ppol.txt';
        'L7_Ppol.txt';
        'L8_Ppol.txt';
        'L9_Ppol.txt';
        'L10_Ppol.txt'];
       

mapes=[
       'TPpol';
       'GPpol'
       ]


[N Pvienalga]=size(faili);
[M Pvienalga]=size(mapes);


for n=1:N
  h=figure(n);
  clear rez;
  for m=1:M
    nosaukums=[mapes(m,:) '/' faili(n,:)];
    indata=dlmread(nosaukums, '\t', [18,0,3664,13]);
    #indata=load('-ascii',nosaukums);
    imin=21;
    imax=1678;
    indata_trim=indata(imin:imax,:);
    #i400=247;
    #dark=indata(imin:i400,2);
    #mean_dark=mean(dark);
    #x=indata(:,1);
    #y=indata(:,2);
    
      #nm_original=indata_trim(:,1);
     
      #intens_original=(indata_trim(:,2)-mean_dark)/sec;
      trim_w=5;
      rempleaks=removepeaks(indata_trim,trim_w);
        for q=1:20
         rempleaks=removepeaks(rempleaks,trim_w);
        endfor
      nm=rempleaks(:,1);
      intens=(rempleaks(:,2))/sec;
      intens_interp = interp1 (nm, intens, lambda, "spline");
      y = intens_interp./i_ref_i;
    if m==1
       rez=[lambda y];
    else
       rez=[rez y];
    endif
    
  endfor
  
  plot (rez (:,1), rez(:,2), "linewidth", 3, "r", rez(:,1), rez(:,3), "linewidth", 3, "b");
  hold on; 
  axis ([500 700 0 15]);
  set(gca,'XTick',500:50:700);
  set(gca,'YTick',0:5:15);
  #set(h, "papersize", [1000, 1000]);
  set(gca, "linewidth", 2, "fontsize", 20)
  
  xlabel ('\lambda, nm', "fontsize", 20);
  ylabel ("Izkliede, a.u.", "fontsize", 20);
  legend ({"Bez grafena", "Ar grafenu"},"fontsize", 13, "location", "southeast");
  grid on;
  izejas_bilde=[rezmape strrep(faili(n,:),".txt",".png")];
  saveas(h, izejas_bilde);
  
  izejas_dati=[rezmape strrep(faili(n,:),".txt",".dat")];
  save ("-ascii", izejas_dati, "rez");
  close all;
endfor