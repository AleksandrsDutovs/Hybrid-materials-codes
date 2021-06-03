# 1. USE .TXT FILES WITH HEADER ONLY!
# 2. REFERENCE SPECTRUM FILE NAME MUST BE: ref_spektrs.txt
# 3. SPECTRUM FILES MUST BE STARTING FROM: R00001.txt
close all;
pkg load optim;

# mainigo definesana
lambda_s = 420;                     # gaismas diapazona sakums
lambda_b = 800;                     # gaismas diapazona beigas
biezums_teoretiski = 260;           # aptuvena biezuma vertiba, kadai vajadzetu but, nanometros
spektru_skaits = 598;                # ieguto spektru skaits
intervals = 0.4;                      # spektru uznemsanas intervals

lambda=[lambda_s:lambda_b];
lambda_prim = lambda';

n0=ones(size(lambda));              # atstarosanas koeficients gaisam
n_Al2O3=1.60*n0;                    # atstarosanas koeficients Al2O3
epsair=n0.^2;
#[nmAlTab nkAlTab] = AlRakic1998();  # complex refractive index of Al (Aluminium) - Rakic
#nkAl=interp1 (nmAlTab,              # interpolated nk for aluminum
#              nkAlTab,
#              lambda,
#              "spline");
[nmAl2O3Tab nkAl2O3Tab]=Malitson1972();
nkAl2O3 = interp1 (nmAl2O3Tab, nkAl2O3Tab, lambda, "spline");
epsAl2O3=nkAl2O3.^2;

[nmAlTab nkAlTab] = AlRakic1998(); # complex refractive index of Al (Aluminium) - Rakic
nkAl=interp1 (nmAlTab, nkAlTab, lambda, "spline"); # interpolated nk for aluminum
epsAl=nkAl.^2;
AlConcInBarrLayer=0.000001;

n_barjer = bruggeman(AlConcInBarrLayer,epsAl,epsAl2O3);


global ntilde=[                     # nk values for multilayer
	       n0;
	       n_Al2O3;
         n_barjer;
	       nkAl
	       ];              

function R = multilay (nm,p)        # p(1) thickness of Al2O3
				                            # p(2) constant multiplier for normalization   
  global ntilde;
  R=zeros(size(nm));
  l=[                               # multilayer thickness in meters
     p(1);
     40;
     p(2)
     ]*1E-9;
  
  lambda_0 = nm * 1E-9;             # wavelength in meters
  L=length(lambda_0);  
  
  for a=1:L
    n=ntilde(:,a);
    M=[1 0 ; 0 1];
    
    for m=1:2
      tau=2*n(m)/(n(m)+n(m+1));
      rho=(n(m)-n(m+1))/(n(m)+n(m+1));
      T= [1 rho ; rho 1 ]/tau;
      M=M*T;
      
      k=n(m+1)*2*pi/lambda_0(a);
      P=[
	 exp(i * k *l(m))  0;
	 0 exp(-i * k *l(m))
	 ];
      M=M*P;
    endfor
    Ef=[                            # Reflection from the last surface
	1;
	(n(m+1)-n(m+2))/(n(m+1)+n(m+2))
	];
    E_0=M*Ef;
    R(a)=abs( E_0(2)/ E_0(1))^2*p(2);
  endfor
endfunction

biezumi_anod_laika = NaN(spektru_skaits:1);

#pedeja merijuma iegusana vispirms
ref_spektrs = dlmread("ref_spektrs.txt",'\t',[18,0,3664,13]);
lambda_ref_spektrs = ref_spektrs(:,1);
intens_ref_spektrs = ref_spektrs(:,2);
intens_ref_spektrs_interp = interp1(lambda_ref_spektrs,intens_ref_spektrs,lambda_prim,"spline");

beigu_spektrs = dlmread(sprintf("R%05d.txt", spektru_skaits),'\t',[18,0,3664,13]);
lambda_spektrs = beigu_spektrs(:,1);
intens_spektrs = beigu_spektrs(:,2);  
intens_spektrs_interp = interp1(lambda_spektrs,intens_spektrs,lambda_prim,"spline");

R = multilay(lambda,[0 1]);                                
R_spektrs = ((intens_spektrs_interp ./ intens_ref_spektrs_interp ).*R')';
minetas_vertibas = [biezums_teoretiski 1];
parbaude = multilay(lambda,minetas_vertibas);

[appr,p]= leasqr(lambda,R_spektrs,minetas_vertibas,@multilay);
#normesana
norm1= R_spektrs./p(2);
norm2 = appr./p(2);                 
                               
biezumi_anod_laika(spektru_skaits,1) = p(1)+40;   

figure(spektru_skaits);
plot(lambda,norm1,'*',lambda,norm2,'k',"linewidth",5);
set(gca,"linewidth", 4,"fontsize", 24);                # smukakam grafikam. labak piemerots attelu liksanai darbos 
xlabel ("\\lambda ,nm");
ylabel ("R(\\lambda)");
legend (sprintf ("h(PAAO) = %d nm", p(1)));
saveas (spektru_skaits,sprintf("R_anod_laika_%05d.png",spektru_skaits));

i = 1;
while i <= (spektru_skaits)                                    
   dati = dlmread(sprintf("R%05d.txt", (i)),'\t',[18,0,3664,13]);
   lambda_spektrs = dati(:,1);
   intens_spektrs = dati(:,2);  
   intens_spektrs_interp = interp1(lambda_spektrs,intens_spektrs,lambda_prim,"spline");
   R = multilay(lambda,[0 1]);                                
   R_spektrs = ((intens_spektrs_interp ./ intens_ref_spektrs_interp ).*R')';
   minetas_vertibas = [(0.41*i) 0.9];
   #if i < 64 
    # minetas_vertibas = [200 0.9];
  # endif
   #if i < 50
    # minetas_vertibas = [150 0.9];
   #endif
   parbaude = multilay(lambda,minetas_vertibas);
   [appr,p]= leasqr(lambda,R_spektrs,minetas_vertibas,@multilay);   
   #normesana
   norm1= R_spektrs ./ p(2);
   norm2 = appr ./ p(2);
   
   biezumi_anod_laika((i),1) = p(1)+40;
   figure((i));
   plot(lambda,norm1,'*',lambda,norm2,'k',"linewidth",5);
   set(gca,"linewidth", 4,"fontsize", 24);            # smukakam grafikam. labak piemerots attelu liksanai darbos 
   title (i);
   xlabel ("\\lambda ,nm");
   ylabel ("R(\\lambda)");
   legend (sprintf ("h(PAAO) = %d nm, %d", p(1)+40, p(2)));
   saveas ((i),sprintf("R_anod_laika_%05d.png",(i)));
   
   i = i+1;
  close all;
endwhile
save biezumi_anod_laika.txt biezumi_anod_laika;   

j=1;
while j <= (spektru_skaits)
  laiks (1, j) = j*intervals;  
  f = laiks';
  j = j+1;
endwhile
#for j = 1:spektru_skaits
  #laiks(1,j+1) = j*intervals;
  #f = laiks'
#endfor

figure(1)
plot(f,biezumi_anod_laika,"g","linewidth",5);
set(gca,"linewidth", 4,"fontsize", 24);
xlabel ("t ,s");
ylabel ("h(AAO) ,nm");
legend({'h(PAAO)'},"location", "northwest");
saveas(1,"parauga_augsanas_dinamika.png");
