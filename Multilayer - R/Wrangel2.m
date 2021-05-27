function Wrangel2(datfilename) 	#this function is used for porous anodized alumina layer thickness determination using light electromagnetic wave transfer model in multilayer environment
				#name of function is dedicated to baron Pyotr Wrangel

lambda=[450:800];              # wavelength range in nm

n0=ones(size(lambda));          # refractive index in air 
n_ZnO=1.95*n0;                  # refractive index of ZnO
n_Al2O3=1.60*n0;                # refractive index of Al2O3
epsair=n0.^2;
#

[nmAl2O3Tab nkAl2O3Tab]=Malitson1972();
nkAl2O3 = interp1 (nmAl2O3Tab, nkAl2O3Tab, lambda, "spline");
epsAl2O3=nkAl2O3.^2;

[nmAlTab nkAlTab] = AlRakic1998();                  # complex refractive index of Al 
nkAl=interp1 (nmAlTab, nkAlTab, lambda, "spline");  # interpolated nk for aluminum
epsAl=nkAl.^2;
AlConcInBarrLayer=0.000001;

n_barjer = bruggeman(AlConcInBarrLayer,epsAl,epsAl2O3); #here we define refractive index of alumina layer (40 nm) using Bruggeman approximation

        
        # material properties 
global ntilde=[
	       n0;
	       n_ZnO;
	       n_Al2O3;
         n_barjer;
         nkAl
	       ]; 

pkg load optim;
w=50;
  datfilename;
  indata=load('-ascii',datfilename);
nm_exp=indata(:,1);
nhy_exp=indata(:,2)*0.95; 
nhy_exp_interp = interp1(nm_exp, nhy_exp, lambda, "linear");


fitfilename = ['wrangel2_' strrep(datfilename, '.dat', '.fit')];  #here we are saving fitted data (PAAO layer thickness) in separate file
  p= [426 0.98]'                                                  #here we guess layer thickness for fitting
  #p = load ('-ascii', fitfilename);                              #for better fitting after first run, comment p=[] line and uncomment this line
  [F, P2]=leasqr(lambda, nhy_exp_interp, p, @fixlay);             #here we approximate light tranmition and propagation model (function named fixlay) described in sections 1.6 and 2.5 of Master Thesys
  P2                                                              #in actual model layer thickness of ZnO and alumina barrier layer (35 and 40 nm resp.) were set as constant values

R1=fixlay(lambda,p); 
save ("-ascii",  fitfilename, "P2");

h = figure();
plot( ...
    #lambda,R1, "linewidth", 3,"g", ...
     lambda, F, "linewidth", 3,"r", ...
     lambda,nhy_exp_interp);
set(h, "papersize", [100, 100]);
set (gca, "fontsize", 28)
#title (strrep (datfilename, '_', '\_' ));
xlabel ("Wavelenght (nm)", "fontsize", 24);
ylabel ("Reflectance", "fontsize", 24);
legend ({'Simulated', 'Experimental'}, "fontsize", 18, "location", "northwest");
grid on;
hold on;

#plot(nm_exp,nhy_exp)
axis ([450 750 0.6 1.14]);
set(gca,'XTick',450:100:750);
set(gca,'YTick',0.6:0.2:1.14);
saveas(1,['wrangel2_' strrep(datfilename, '.dat', '.png')]);
savedat=[lambda; 
          R1]';
save ("-ascii", ['wrangel2_' strrep(datfilename, '.dat', '.txt')], "savedat");

endfunction

