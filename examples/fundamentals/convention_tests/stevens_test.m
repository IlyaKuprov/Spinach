% Tests of Spinach Stevens operator function against
% explicit expressions from the literature.
%
% i.kuprov@soton.ac.uk

function stevens_test()

% Spin and multiplicity
sqn=12; mult=2*sqn+1;

% Spinach results, rank 6
O_6_p6_Spinach=stevens(mult,6,+6);
O_6_p5_Spinach=stevens(mult,6,+5);
O_6_p4_Spinach=stevens(mult,6,+4);
O_6_p3_Spinach=stevens(mult,6,+3);
O_6_p2_Spinach=stevens(mult,6,+2);
O_6_p1_Spinach=stevens(mult,6,+1);
O_6_0_Spinach =stevens(mult,6, 0);
O_6_m1_Spinach=stevens(mult,6,-1);
O_6_m2_Spinach=stevens(mult,6,-2);
O_6_m3_Spinach=stevens(mult,6,-3);
O_6_m4_Spinach=stevens(mult,6,-4);
O_6_m5_Spinach=stevens(mult,6,-5);
O_6_m6_Spinach=stevens(mult,6,-6);

% Explicit formulae, rank 6
S=pauli(mult); s=sqn*(sqn+1);
O_6_p6_Stevens=(1/2)*(S.p^6+S.m^6);
O_6_p5_Stevens=(1/2)*acomm(S.z,S.p^5+S.m^5)/2;
O_6_p4_Stevens=(1/2)*acomm(11*S.z^2-(s+38)*S.u,S.p^4+S.m^4)/2;
O_6_p3_Stevens=(1/2)*acomm(11*S.z^3-(3*s+59)*S.z,S.p^3+S.m^3)/2;
O_6_p2_Stevens=(1/2)*acomm(33*S.z^4-(18*s+123)*S.z^2+(s^2+10*s+102)*S.u,S.p^2+S.m^2)/2;
O_6_p1_Stevens=(1/2)*acomm(33*S.z^5-(30*s-15)*S.z^3+(5*s^2-10*s+12)*S.z,S.p+S.m)/2;
O_6_0_Stevens=231*S.z^6-(315*s-735)*S.z^4+(105*s^2-525*s+294)*S.z^2-(5*s^3-40*s^2+60*s)*S.u;
O_6_m1_Stevens=(1/2i)*acomm(33*S.z^5-(30*s-15)*S.z^3+(5*s^2-10*s+12)*S.z,S.p-S.m)/2;
O_6_m2_Stevens=(1/2i)*acomm(33*S.z^4-(18*s+123)*S.z^2+(s^2+10*s+102)*S.u,S.p^2-S.m^2)/2;
O_6_m3_Stevens=(1/2i)*acomm(11*S.z^3-(3*s+59)*S.z,S.p^3-S.m^3)/2;
O_6_m4_Stevens=(1/2i)*acomm(11*S.z^2-(s+38)*S.u,S.p^4-S.m^4)/2;
O_6_m5_Stevens=(1/2i)*acomm(S.z,S.p^5-S.m^5)/2;
O_6_m6_Stevens=(1/2i)*(S.p^6-S.m^6);

% Differences, rank 6
disp([norm(O_6_p6_Spinach-O_6_p6_Stevens,1) norm(O_6_p5_Spinach-O_6_p5_Stevens,1) ...
      norm(O_6_p4_Spinach-O_6_p4_Stevens,1) norm(O_6_p3_Spinach-O_6_p3_Stevens,1) ...
      norm(O_6_p2_Spinach-O_6_p2_Stevens,1) norm(O_6_p1_Spinach-O_6_p1_Stevens,1) ...
      norm(O_6_0_Spinach-O_6_0_Stevens,1) ...
      norm(O_6_m1_Spinach-O_6_m1_Stevens,1) norm(O_6_m2_Spinach-O_6_m2_Stevens,1) ...
      norm(O_6_m3_Spinach-O_6_m3_Stevens,1) norm(O_6_m4_Spinach-O_6_m4_Stevens,1) ...
      norm(O_6_m5_Spinach-O_6_m5_Stevens,1) norm(O_6_m6_Spinach-O_6_m6_Stevens,1)]);

% Spinach results, rank 4
O_4_p4_Spinach=stevens(mult,4,+4);
O_4_p3_Spinach=stevens(mult,4,+3);
O_4_p2_Spinach=stevens(mult,4,+2);
O_4_p1_Spinach=stevens(mult,4,+1);
O_4_0_Spinach =stevens(mult,4, 0);
O_4_m1_Spinach=stevens(mult,4,-1);
O_4_m2_Spinach=stevens(mult,4,-2);
O_4_m3_Spinach=stevens(mult,4,-3);
O_4_m4_Spinach=stevens(mult,4,-4);

% Explicit formulae, rank 4
S=pauli(mult); s=sqn*(sqn+1);
O_4_p4_Stevens=(1/2)*(S.p^4+S.m^4);
O_4_p3_Stevens=(1/2)*acomm(S.z,S.p^3+S.m^3)/2;
O_4_p2_Stevens=(1/2)*acomm(7*S.z^2-(s+5)*S.u,S.p^2+S.m^2)/2;
O_4_p1_Stevens=(1/2)*acomm(7*S.z^3-(3*s+1)*S.z,S.p+S.m)/2;
O_4_0_Stevens=35*S.z^4-(30*s-25)*S.z^2+(3*s^2-6*s)*S.u;
O_4_m1_Stevens=(1/2i)*acomm(7*S.z^3-(3*s+1)*S.z,S.p-S.m)/2;
O_4_m2_Stevens=(1/2i)*acomm(7*S.z^2-(s+5)*S.u,S.p^2-S.m^2)/2;
O_4_m3_Stevens=(1/2i)*acomm(S.z,S.p^3-S.m^3)/2;
O_4_m4_Stevens=(1/2i)*(S.p^4-S.m^4);

% Differences, rank 4
disp([norm(O_4_p4_Spinach-O_4_p4_Stevens,1) norm(O_4_p3_Spinach-O_4_p3_Stevens,1) ...
      norm(O_4_p2_Spinach-O_4_p2_Stevens,1) norm(O_4_p1_Spinach-O_4_p1_Stevens,1) ...
      norm(O_4_0_Spinach-O_4_0_Stevens,1) ...
      norm(O_4_m1_Spinach-O_4_m1_Stevens,1) norm(O_4_m2_Spinach-O_4_m2_Stevens,1) ...
      norm(O_4_m3_Spinach-O_4_m3_Stevens,1) norm(O_4_m4_Spinach-O_4_m4_Stevens,1)]);

% Spinach results, rank 2
O_2_p2_Spinach=stevens(mult,2,+2);
O_2_p1_Spinach=stevens(mult,2,+1);
O_2_0_Spinach =stevens(mult,2, 0);
O_2_m1_Spinach=stevens(mult,2,-1);
O_2_m2_Spinach=stevens(mult,2,-2);

% Explicit formulae, rank 2
S=pauli(mult); s=sqn*(sqn+1);
O_2_p2_Stevens=(1/2)*(S.p^2+S.m^2);
O_2_p1_Stevens=(1/2)*acomm(S.z,S.p+S.m)/2;
O_2_0_Stevens=3*S.z^2-s*S.u;
O_2_m1_Stevens=(1/2i)*acomm(S.z,S.p-S.m)/2;
O_2_m2_Stevens=(1/2i)*(S.p^2-S.m^2);

% Differences, rank 2
disp([norm(O_2_p2_Spinach-O_2_p2_Stevens,1) norm(O_2_p1_Spinach-O_2_p1_Stevens,1) ...
      norm(O_2_0_Spinach-O_2_0_Stevens,1) ...
      norm(O_2_m1_Spinach-O_2_m1_Stevens,1) norm(O_2_m2_Spinach-O_2_m2_Stevens,1)]);

% IST formulae, rank 2
T=irr_sph_ten(mult,2);
T2p2=T{1}; T2p1=T{2}; T20=T{3}; T2m1=T{4}; T2m2=T{5};
O_2_p2_IST=T2m2+T2p2;
O_2_p1_IST=(T2m1-T2p1)/2;
O_2_0_IST=sqrt(6)*T20;
O_2_m1_IST=1i*(T2m1+T2p1)/2;
O_2_m2_IST=1i*(T2m2-T2p2);

% Differences, rank 2
disp([norm(O_2_p2_Spinach-O_2_p2_IST,1) norm(O_2_p1_Spinach-O_2_p1_IST,1) ...
      norm(O_2_0_Spinach-O_2_0_IST,1) ...
      norm(O_2_m1_Spinach-O_2_m1_IST,1) norm(O_2_m2_Spinach-O_2_m2_IST,1)]);

end

