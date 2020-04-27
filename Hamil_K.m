function result = Hamil_K(kx,ky)

 hbar = 1.054572e-34; % reduced Planck constant ( J . s )
 vf=1e6;  % Fermi velocity m/sec
 
 fk=kx+i*ky;
 result=hbar.*vf.*[1e-9 fk; fk' -1e-9]; % Very small band gap is introduced to avoid NAN resluts
 %result=hbar.*vf.*[0 fk; fk' 0];
 
end