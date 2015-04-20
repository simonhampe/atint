for my $i (1 .. 8) {
	my $f = load($i."f");
	my $g = load($i."g");
	if(!compare_object($i."", morphism_composition($f,$g))) {
		return 0;
	}
}
return 1;
