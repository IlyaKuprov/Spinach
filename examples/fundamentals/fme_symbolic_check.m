% Exact integral-versus-Fourier checks of low-order Floquet coefficients.
% Requires Symbolic Math Toolbox. Related to Appendix A and supplement
% Section S5 of Bock, Ernst, and Uhrig, DOI:10.5194/mr-2026-10.
% A two-by-two noncommuting model is used, not their MAS calculation.
% The mean Hamiltonian and the first inverse-frequency AHT correction
% are independently obtained from the first two Magnus integrals.
% Leading AHT and Van Vleck kicks are checked by exact integration;
% all these symbolic coefficients are then compared with production
% floquet_magnus coefficients after a numerical substitution.
%
% talos@spindynamics.org

function fme_symbolic_check()

% Dimensionless phase x=omega*t and exactly represented spin matrices
syms x y a b c real
p=sym(pi);
S=pauli(2);
Sx=sym(full(S.x)); Sy=sym(full(S.y)); Sz=sym(full(S.z));
H0=a*Sz; Hp=(b*Sx-1i*c*Sy)/2; Hm=Hp';
Hx=H0+Hp*exp(1i*x)+Hm*exp(-1i*x);
Hy=subs(Hx,x,y);

% First two Magnus integrals over one period; omega^-1 factored out
F0=simplify(int(Hx,x,0,2*p)/(2*p));
F1=simplify(-1i*int(int(Hx*Hy-Hy*Hx,y,0,x),x,0,2*p)/(4*p));
F1_fourier=Hp*Hm-Hm*Hp+H0*(Hp-Hm)-(Hp-Hm)*H0;
assert(isequal(simplify(F0-H0),sym(zeros(2))));
assert(isequal(simplify(F1-F1_fourier),sym(zeros(2))));

% Integrating the oscillatory Hamiltonian fixes the AHT endpoint gauge
Ka=simplify(int(Hy-H0,y,0,x));
Kv=simplify(Ka-int(Ka,x,0,2*p)/(2*p));
Kv_fourier=-1i*Hp*exp(1i*x)+1i*Hm*exp(-1i*x);
assert(isequal(simplify(Kv-Kv_fourier),sym(zeros(2))));
assert(isequal(simplify(subs(Ka,x,0)),sym(zeros(2))));

% Exact expressions meet the production routine at nonzero coefficients
values=[3 2 5]; omega=17;
H={double(subs(Hm,[a b c],values)),...
   double(subs(H0,[a b c],values)),...
   double(subs(Hp,[a b c],values))};
[Fa,~]=floquet_magnus(H,omega,1,0,'aht');
[Fzero,Ka_num]=floquet_magnus(H,omega,0,0,'aht');
[~,Kv_num]=floquet_magnus(H,omega,0,0,'van_vleck');
Fa_exact=double(subs(F0+F1/omega,[a b c],values));
errors=zeros(1,8);
errors(1:2)=[norm(Fa-Fa_exact,'fro'),norm(Fzero-H{2},'fro')];

% Project exact kicks onto every retained harmonic, not selected times
for n=-1:1
    ka_coeff=simplify(int(Ka*exp(-1i*n*x),x,0,2*p)/(2*p));
    kv_coeff=simplify(int(Kv*exp(-1i*n*x),x,0,2*p)/(2*p));
    ka_exact=double(subs(ka_coeff,[a b c],values))/omega;
    kv_exact=double(subs(kv_coeff,[a b c],values))/omega;
    errors(2*n+5:2*n+6)=[norm(Ka_num{n+2}-ka_exact,'fro'),...
                        norm(Kv_num{n+2}-kv_exact,'fro')];
end

% Allow 100 ulps for the finite matrix sums and symbolic-to-double conversion
bound=100*eps(max(1,norm(Fa_exact,'fro')));
assert(max(errors)<bound,'Symbolic and production coefficients disagree.');
fprintf('Exact Magnus/Fourier identities passed; numerical error %.3g.\n',...
        max(errors));

end


