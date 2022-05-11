syms E I K L w(x) k L

Dw      = diff(w, x);
DDw     = diff(w, x, 2);
DDDw    = diff(w, x, 3);

eqn     = E*I * diff(w, x, 4) == k * x;

cond_1  = w(0) == 0;
cond_2  = Dw(0) == 0;
cond_3  = DDw(L) == 0;
cond_4  = DDDw(L) == 0;

solution = dsolve(eqn, cond_1, cond_2, cond_3, cond_4);

disp(solution)