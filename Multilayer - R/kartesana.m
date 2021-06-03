# 1. USE .TXT FILES WITH HEADER ONLY!
# 2. REFERENCE SPECTRUM FILE NAME MUST BE: ref_spektrs.txt
# 3. SPECTRUM FILES MUST BE STARTING FROM: R00001.txt

pkg load optim; 		#Interpretatation of spectral data obtained from measurments in PAAO synthesis process

# setting up variables
steps = 19;                         # steps for 1 line
lambda_s = 450;                     # wavelenght start
lambda_b = 800;                     # wavelenght end
biezums_teoretiski = 270;           # theoretical layer, nm
mala = 0.9;                         # sample lenght


lambda=[lambda_s:lambda_b];
lambda_prim = lambda';

n0=ones(size(lambda));              # reflectance coeff for air
n_Al2O3=1.60*n0;                    # reflectance coeff Al2O3
epsair=n0.^2;

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


#Ref spektra ielase
ref_spektrs = dlmread("MC_ref_10xBF_cam_on.txt",'\t',[18,0,3664,13]); #reference spectum file name
lambda_ref_spektrs = ref_spektrs(:,1);
intens_ref_spektrs = ref_spektrs(:,2);
intens_ref_spektrs_interp = interp1(
                                    lambda_ref_spektrs,
                                    intens_ref_spektrs,
                                    lambda_prim,
                                    "spline");

all_steps = steps^2;                              
biezumi = NaN(all_steps:1);   #biezumi = layers

for i = 1:all_steps
   dati = dlmread(
        sprintf("R%05d.txt", i),
        '\t',[18,0,3664,13]);
   lambda_spektrs = dati(:,1);
   intens_spektrs = dati(:,2);  
   intens_spektrs_interp = interp1(
                                   lambda_spektrs,
                                   intens_spektrs,
                                   lambda_prim,
                                   "spline");
   R = multilay(lambda,[0 1]);                                
   R_spektrs = ((intens_spektrs_interp ./ intens_ref_spektrs_interp ).*R')';
   minetas_vertibas = [biezums_teoretiski 1];           #minetas vertibas = guessed values
   parbaude = multilay(lambda,minetas_vertibas);

   [appr,p]= leasqr(
                    lambda,
                    R_spektrs,
                    minetas_vertibas,
                    @multilay);
    
    biezumi(i,1) = p(1)+40;
endfor

save biezumi.txt biezumi;

#STATISTIKA
minimum = min(biezumi);
maximum = max(biezumi);
vid = mean(biezumi);
sd = std(biezumi);

statistika = [minimum;
              maximum;
              vid;
              sd];

save statistika.txt statistika;
              
# matrix of thickness for map
for z = 1:steps
  for n =1:steps
    biezumu_matrix(n,z)=biezumi(n*steps-(steps-z),1);
  endfor 
endfor

# map plotting
grid = linspace(0,mala,steps);
xx=grid';
yy=xx;

figure(1);
[x,y]=meshgrid(xx,yy);
surf(x,
     y,
     biezumu_matrix);
colorbar ("eastoutside");
set(gca,                             
    "linewidth",1,
    "fontsize", 14);
view (-40, 60);
xlabel ("X ass, cm");
ylabel ("Y ass, cm");
zlabel("h(PAAO), nm");

saveas (1, "parauga_biezuma_karte.png");
save chart_data.txt x y biezumu_matrix;
