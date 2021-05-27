function rez=gauss1(eV,p)
  rez=p(3)*exp(-((eV-p(1)).^2)/(2*p(2).^2));
endfunction