function rez=gauss1(eV,p)
  rez=p(3)*exp(-((eV-p(1)).^2)/(2*p(2).^2)) ...
    + p(6)*exp(-((eV-p(4)).^2)/(2*p(5).^2));
endfunction