function nt = bruggeman(fa,ea,eb)   #Simulation of PAAO barrier layer complex refractive index by Bruggeman approximation 

e1 = - ((sqrt((9*eb.^2-18*ea.*eb+9*ea.^2)*fa.^2 + ((-12*eb.^2)+18*ea.*eb-6*ea.^2)*fa+4*eb.^2+4*ea.*eb+ea.^2))
+((3*eb-3*ea)*fa-2*eb+ea))/4 ;
e2 =   ((sqrt((9*eb.^2-18*ea.*eb+9*ea.^2)*fa.^2 + ((-12*eb.^2)+18*ea.*eb-6*ea.^2)*fa+4*eb.^2+4*ea.*eb+ea.^2))
+((3*ea-3*eb)*fa+2*eb-ea))/4 ;

  em = e1;
  

  nt=sqrt(em);
  
