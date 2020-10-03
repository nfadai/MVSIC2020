function [Si_a, Sti_a] = reorder(Si, Sti)

Si_a = Si;
Sti_a = Sti;

Si_a(1) = Si(12);
Sti_a(1) = Sti(12);

Si_a(2) = Si(13);
Sti_a(2) = Sti(13);

Si_a(3) = Si(1);
Sti_a(3) = Sti(1);

Si_a(4) = Si(2);
Sti_a(4) = Sti(2);

Si_a(5) = Si(4);
Sti_a(5) = Sti(4);

Si_a(6) = Si(5);
Sti_a(6) = Sti(5);

Si_a(7) = Si(8);
Sti_a(7) = Sti(8);

Si_a(8) = Si(3);
Sti_a(8) = Sti(3);

Si_a(9) = Si(6);
Sti_a(9) = Sti(6);

Si_a(10) = Si(7);
Sti_a(10) = Sti(7);

Si_a(11) = Si(11);
Sti_a(11) = Sti(11);

Si_a(12) = Si(9);
Sti_a(12) = Sti(9);

Si_a(13) = Si(10);
Sti_a(13) = Sti(10);

Si_a(14) = Si(14);
Sti_a(14) = Sti(14);

Si_a(15) = Si(15);
Sti_a(15) = Sti(15);

Si_a(16) = Si(16);
Sti_a(16) = Sti(16);

end
