function R = multilay (nm,p) #считает коэф отражения гибрида в завис от длины волны
				# p(1) thickness of ZnO
				# p(2) thickness of Al2O3
				# p(3) barrier 
        # p(4) коэф отражения полированного алюминия, принять как 1, 
        # но Ал не ровный (https://en.wikipedia.org/wiki/Reflectance) 
  global ntilde;
  R=zeros(size(nm)); #пустой массив с нулями
  l=[                           # multilayer thickness in meters
     p(1);
     p(2);
     p(3);
     p(4)
               ]*1E-9;
  
  lambda_0 = nm * 1E-9;         # wavelength in meters
  L=length(lambda_0);
  
  for a=1:L
    n=ntilde(:,a);
    M=[1 0 ; 0 1];
    
    for m=1:3                     # m1 = Gaiiz/ZnO, m2 = ZnO/Al2O3
      tau=2*n(m)/(n(m)+n(m+1));   # tau, rho = Frenel coeff, 90 grad falling light
      rho=(n(m)-n(m+1))/(n(m)+n(m+1)); 
      T= [1 rho ; rho 1 ]/tau;    # sophocles orfanidis e-mag waves
                                  # matching matrix, transm and reflection are one-time processes
      M=M*T;
      
      k=n(m+1)*2*pi/lambda_0(a);
      P=[
	        exp(i * k *l(m))  0;
	        0 exp(-i * k *l(m))
	      ];                        # propagation matrix vilna izplat vid?
      M=M*P;
    endfor
    Ef=[     # Reflection from the last surface (Al)
	      1;   # krito?a el lauka intens
	      (n(m+1)-n(m+2))/(n(m+1)+n(m+2)) #atstarota el lauka intens (att. 5.2.1 b orfanidis)
	      ];
    E_0=M*Ef;    # = kopejie lauki m1 m2 reizinati ar Al reflection laukiem     
    R(a)=abs( E_0(2)/ E_0(1))^2*p(4); 
        # R = atstarosanas koef ka el lauka attiecibas modula kvadrats
  endfor
  
endfunction