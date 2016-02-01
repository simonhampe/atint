my $m1 = new matroid::Matroid(N_ELEMENTS=>4,BASES=>[[0,1],[0,2],[1,3],[2,3]]);
my $m2 = matroid::uniform_matroid(2,4);
my $m3 = new matroid::Matroid(N_ELEMENTS=>4,BASES=>[[0,1],[0,2],[0,3],[1,3],[2,3]]);
my $m4 = new matroid::Matroid(N_ELEMENTS=>4,BASES=>[[0,1],[0,2],[1,2],[1,3],[2,3]]);
my @r = map { matroid_ring_cycle<Min>($_)} ($m1,$m2,$m3,$m4);

my $n1 = matroid_ring_cycle<Max>(matroid::uniform_matroid(2,3));
my $n2 = zero_in_matroid_ring<Max>(3);

compare_data("1", matroid_ring_linear_space(@r));
compare_data("2", matroid_ring_linear_space($n1,$n2));
