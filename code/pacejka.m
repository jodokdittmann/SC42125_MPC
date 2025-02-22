function y = pacejka(B, C, D, E, x)

	y = D*sin(C*atan(B*x - E*(B*x - atan(B*x))));

end
