# Optical constants of Al2O3 (Aluminium oxide, Sapphire)
# I. H. Malitson, M. J. Dodge, Refractive index and bire-
# fringence of synthetic sapphire, J. Opt. Soc. Am. 62 (11)
# (1972) 1405.
# http://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o
function [nm nk]=Malitson1972()
  data =load("-ascii",'Malitson-o.txt');
  nm=data(:,1)*1000;
  nk=data(:,2);
endfunction