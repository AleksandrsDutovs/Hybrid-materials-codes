function R = fixlay (nm,p) 	# this is supportive function for PAAO layer thickness determination
				# p(1) thickness of ZnO
				# p(2) thickness of Al2O3
				# p(3) barrier layer
        
  global ntilde;
  R=zeros(size(nm)); 
  l=[                           # multilayer thickness in meters
     35;
     p(1);
     40;
     0.94
               ]*1E-9;
  
  lambda_0 = nm * 1E-9;         # wavelength in meters
  L=length(lambda_0);
  
  for a=1:L
    n=ntilde(:,a);
    M=[1 0 ; 0 1];
    
    for m=1:3                    # m1 = Gaiiz/ZnO, m2 = ZnO/Al2O3
      tau=2*n(m)/(n(m)+n(m+1));   # tau, rho = Frenel coefficients, 90 degree falling light
      rho=(n(m)-n(m+1))/(n(m)+n(m+1)); 
      T= [1 rho ; rho 1 ]/tau;    # matching matrix, transmitatnce and reflection are one-time processes
      M=M*T;
      
      k=n(m+1)*2*pi/lambda_0(a);
      P=[
	        exp(i * k *l(m))  0;
	        0 exp(-i * k *l(m))
	      ];                        # propagation matrix 
      M=M*P;
    endfor
    Ef=[     # Reflection from the last surface (Al)
	      1;   
	      (n(m+1)-n(m+2))/(n(m+1)+n(m+2)) #intensity of reflected wave
	      ];
    E_0=M*Ef;    # = kopejie lauki m1 m2 reizinati ar Al reflection laukiem     
    R(a)=abs( E_0(2)/ E_0(1))^2*1; 
        # R = reflectance coefficient as square of absolute value of two e-fields relation
  endfor
  
endfunction
