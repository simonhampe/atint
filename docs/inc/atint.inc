<div xmlns:html="http://www.w3.org/1999/xhtml" xmlns:pm="http://www.polymake.de/ns/docs#3" id="content"><h1>application: atint</h1><div class="undoc">UNDOCUMENTED</div><br/><b>imports from: </b><a href="common.html">common</a>, <a href="graph.html">graph</a><br/><b>uses: </b><a href="fan.html">fan</a>, <a href="group.html">group</a>, <a href="matroid.html">matroid</a>, <a href="polytope.html">polytope</a>, <a href="topaz.html">topaz</a>, <a href="tropical.html">tropical</a><br/><div class="level2"><h2><b>Objects</b></h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i21" onclick="swap_content('span:i21'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i21" onclick="swap_content('span:i21'); return false;"><h3><div class="object">MinMaxFunction</div></h3></a><span style="display: inline" id="span:i21"><div class="undoc"><p>UNDOCUMENTED</p></div><b>derived from: </b><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a><br/><div class="level3"/><h4>Properties of MinMaxFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i101" onclick="swap_content('span:i101'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i101" onclick="swap_content('span:i101'); return false;"><h3><div class="category">Basic properties</div></h3></a><span style="display: inline" id="span:i101"><li><div class="li"><div class="icon" id="icon:i144" onclick="swap_content('span:i144'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i144" onclick="swap_content('span:i144'); return false;"><b>CONSTANT_COEFFICIENTS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i144"><div class="descr"><p> A vector of the constant coefficients of the function. The i-th entry is the constant term of the function 
 represented by the i-th row of <a href="#i143" onclick="unfold('span:i143');">LINEAR_COEFFICIENTS</a>.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i138" onclick="swap_content('span:i138'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i138" onclick="swap_content('span:i138'); return false;"><b>DOMAIN_DIMENSION</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i138"><div class="descr"><p> The dimension n of the domain Rn of the function. Is equal to the columns of <a href="#i143" onclick="unfold('span:i143');">LINEAR_COEFFICIENTS</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i142" onclick="swap_content('span:i142'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i142" onclick="swap_content('span:i142'); return false;"><b>FUNCTION_MATRIX</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i142"><div class="descr"><p> This is simply the concatenation of <a href="#i143" onclick="unfold('span:i143');">LINEAR_COEFFICIENTS</a> and <a href="#i144" onclick="unfold('span:i144');">CONSTANT_COEFFICIENTS</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i140" onclick="swap_content('span:i140'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i140" onclick="swap_content('span:i140'); return false;"><b>INPUT_STRING</b>:
		</a>common::String<br/><span style="display: inline" id="span:i140"><div class="descr"><p> The syntax for STRING in pseudo-regex is the following (after deleting any whitespace): 
 "(min|max)\((SUM)+(,SUM)*\)"
 where SUM = "(TERM)+([\+\-](TERM))*"
 where TERM = (RAT)|(RAT[:alpha:][^\+,\)]*)
 where RAT is any string that can be parsed to a polymake Rational and does not contain any letters.
 To make this more comprehensible: STRING takes anything of the form min(...) or max(...), 
 where ... is a comma-separated list of strings of the form a + bx + cy + ... , where a,b,c are rational 
 numbers and x,y are variables. Any string that starts with a letter and does not contain "+" or "," can be 
 a variable. Such a sum can contain several such terms for the same variable and they need not be in any order.
 A term can be of the form "3x", "3*x", but "x3" will be interpreted as 1 * "x3". Coefficients should not contain letters,
 i.e. 3E10*x would be interpreted as 3 * "E10*x". Empty terms are allowed.
 If no INPUT_VARS are given, MinMaxFunction will recognize the dimension of the domain from the number of variables and will
 assign the variables to the standard basis vectors by sorting them alphabetically. 
 E.g. INPUT_STRING=&gt;"min(x+2z,3+z)",INPUT_VARS=&gt;["x","y","z"] is a different function than INPUT_STRING=&gt;"min(x+2z,3+z)",
 since the second is the function f:R² -&gt; R, (a,b) -&gt; min(a+2b,3+b) and the first is f: R³-&gt;R, (a,b,c)-&gt;min(a+2c,3+c)
 If however, <a href="#i141" onclick="unfold('span:i141');">INPUT_VARS</a> is given, the parse will *ignore* any terms containing variables that are not 
 defined in <a href="#i141" onclick="unfold('span:i141');">INPUT_VARS</a>. The variables will be assigned to the standard basis vectors in order of their
 appearance in <a href="#i141" onclick="unfold('span:i141');">INPUT_VARS</a>.
 E.g. INPUT_STRING=&gt;"min(x+2c,3+x)",INPUT_VARS=&gt;["x","y"] will give the function f:R²-&gt;R, (a,b) -&gt; min(a,3+a)
 Note that ANY whitespace in <a href="#i140" onclick="unfold('span:i140');">INPUT_STRING</a> is ignore, as well as any brackets '(',')' besides the ones enclosing
 the function list.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i141" onclick="swap_content('span:i141'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i141" onclick="swap_content('span:i141'); return false;"><b>INPUT_VARS</b>:
		</a>common::Array&lt;String&gt;<br/><span style="display: inline" id="span:i141"><div class="descr"><p> A list of strings that are variables. Any string that starts with a letter (i.e. [a-zA-Z]) and does not contain
 '+' or ',' can be a variable. Illegal variable names will not cause an error, but the parser will not detect them
 in <a href="#i140" onclick="unfold('span:i140');">INPUT_STRING</a>, so their coefficients will be 0. If you assign a value to this, only terms with variables in 
 this array will be parsed. Note that ANY whitespace in a variable name will be ignored. The i-th element of this array
 will be assigned to the i-th standard basis vector.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i143" onclick="swap_content('span:i143'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i143" onclick="swap_content('span:i143'); return false;"><b>LINEAR_COEFFICIENTS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i143"><div class="descr"><p> A matrix of the linear coefficients of the functions. Each row corresponds to a linear function within
 the min or max and the entry in column j is the coefficient of variable x_j</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i137" onclick="swap_content('span:i137'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i137" onclick="swap_content('span:i137'); return false;"><b>NORMAL_FAN</b>:
		</a><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a><br/><span style="display: inline" id="span:i137"><div class="descr"><p> This computes the normal fan of the polytope conv( (c_i,a_i) ) in Q^{n+1}, where the function is 
 max/min(a_i*x + c_i,i=1,..,k) with a_i in Q^n and c_i in Q. (For min it actually computes the negative
 of the normal fan</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i145" onclick="swap_content('span:i145'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i145" onclick="swap_content('span:i145'); return false;"><b>STANDARD_FORM</b>:
		</a>common::String<br/><span style="display: inline" id="span:i145"><div class="descr"><p> This will give a nice string representation of the function, where all variables are of the form 'x_i',
 i = 1,..,<a href="#i138" onclick="unfold('span:i138');">DOMAIN_DIMENSION</a>. There are no redundancies and terms are sorted in ascending alphabetic order of their variable,
 with the constant term at the end</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i139" onclick="swap_content('span:i139'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i139" onclick="swap_content('span:i139'); return false;"><b>USES_MIN</b>:
		</a>common::Bool<br/><span style="display: inline" id="span:i139"><div class="descr"><p> Indicates whether this function uses min (TRUE) or max (FALSE). The default (if you directly give a coefficient matrix,
 but not this value) is FALSE</p></div></span></div></div></li></span></div></div></li></ul><h4>User Methods of MinMaxFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i99" onclick="swap_content('span:i99'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i99" onclick="swap_content('span:i99'); return false;"><b>pow</b></a>
			()  

			<span style="display: inline" id="span:i99"><br/><div class="descr_func"><div class="undoc"><p>UNDOCUMENTED</p></div></div><div class="level3"/></span></div></div></li><li><div class="li"><div class="icon" id="icon:i100" onclick="swap_content('span:i100'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i100" onclick="swap_content('span:i100'); return false;"><h3><div class="category">Tropical geometry</div></h3></a><span style="display: inline" id="span:i100"><li><div class="li"><div class="icon" id="icon:i136" onclick="swap_content('span:i136'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i136" onclick="swap_content('span:i136'); return false;"><b>linearityDomains</b></a>
			()  

			
				→ WeightedComplex<span style="display: inline" id="span:i136"><br/><div class="descr_func"><div class="descr"><p> This computes the (homogenized version of the) domains of linearity of the function. More precisely, it computes the 
 normal fan of the polytope defined by this function in the following way:
 If the standard form reads max(a_i*x + c_i,i=1,..,k) with a_i in Q^n and c_i in Q, then the polytope is
 conv( (c_i,a_i) ) in Q^{n+1} (in homog. coordinates its actually in n+2 coordinates).
 In fact, this returns the property <a href="#i160" onclick="unfold('span:i160');">DOMAIN</a>.</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i135" onclick="swap_content('span:i135'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i135" onclick="swap_content('span:i135'); return false;"><b>valueAt</b></a>
			()  

			
				→ Rational<span style="display: inline" id="span:i135"><br/><div class="descr_func"><div class="descr"><p> Computes the value of this function at a given point. If the point has dimension DOMAIN_DIMENSION, it simply computes the value. 
 If it has dimension DOMAIN_DIMENSION +1, it assumes the point is given in homogeneous coordinates (i.e. the first coordinate should be one or zero, 
 this is not checked, however) and computes the function value of the last DOMAIN_DIMENSION coordinates.</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td><div class="descr"><div>The function value at the given point.</div></div></td></tr></table></div></span></div></div></li></span></div></div></li></ul></span></div></div></li><li><div class="li"><div class="icon" id="icon:i22" onclick="swap_content('span:i22'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i22" onclick="swap_content('span:i22'); return false;"><h3><div class="object">RationalCurve</div></h3></a><span style="display: inline" id="span:i22"><div class="undoc"><p>UNDOCUMENTED</p></div><div class="level3"/><h4>Properties of RationalCurve</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i104" onclick="swap_content('span:i104'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i104" onclick="swap_content('span:i104'); return false;"><h3><div class="category">Basic properties</div></h3></a><span style="display: inline" id="span:i104"><li><div class="li"><div class="icon" id="icon:i151" onclick="swap_content('span:i151'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i151" onclick="swap_content('span:i151'); return false;"><b>COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i151"><div class="descr"><p> A list of positive rational coefficients. The list should have the same length as <a href="#i155" onclick="unfold('span:i155');">SETS</a> and
 contain only entries &gt; 0. The i-th entry then gives the length of the bounded edge defined by
 the i-th partition. If you're not sure if all your coefficients are &gt; 0, use <a href="#i154" onclick="unfold('span:i154');">INPUT_SETS</a> and
 <a href="#i149" onclick="unfold('span:i149');">INPUT_COEFFS</a> instead.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i153" onclick="swap_content('span:i153'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i153" onclick="swap_content('span:i153'); return false;"><b>GRAPH</b>:
		</a>graph::Graph&lt;Undirected&gt;<br/><span style="display: inline" id="span:i153"><div class="descr"><p> Contains the abstract graph (non-metric) corresponding to the curve. All unbounded leaves are modelled
 as bounded edges.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i152" onclick="swap_content('span:i152'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i152" onclick="swap_content('span:i152'); return false;"><b>GRAPH_EDGE_LENGTHS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i152"><div class="descr"><p> Contains the lengths of the edges of <a href="#i153" onclick="unfold('span:i153');">GRAPH</a> that represent bounded edges of the curve.
 The coefficients appear in the order that the corr. edges appear in <a href="#" onclick="unfold('span:');">EDGES</a>.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i149" onclick="swap_content('span:i149'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i149" onclick="swap_content('span:i149'); return false;"><b>INPUT_COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i149"><div class="descr"><p> Same as <a href="#i151" onclick="unfold('span:i151');">COEFFS</a>, except that entries may be &lt;=0. This should have the same length as <a href="#i154" onclick="unfold('span:i154');">INPUT_SETS</a>.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i154" onclick="swap_content('span:i154'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i154" onclick="swap_content('span:i154'); return false;"><b>INPUT_SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i154"><div class="descr"><p> Same as <a href="#i155" onclick="unfold('span:i155');">SETS</a>, except that sets may appear several times.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i148" onclick="swap_content('span:i148'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i148" onclick="swap_content('span:i148'); return false;"><b>INPUT_STRING</b>:
		</a>common::String<br/><span style="display: inline" id="span:i148"><div class="descr"><p> This property can also be used to define a rational curve: A linear combination of partitions is
 given as a string, using the following syntax:
 A partition is given as a subset of {1,..,n} and written as a comma-separated list of leaf
 indices in round brackets, e.g. "(1,2,5)"
 A linear combination can be created using rational numbers, "+","+" and "-" in the obvious way, 
 e.g. "2*(1,2,5) + 1*(3,4,7) - 2(1,2) (The "*" is optional)
 Of course, each set should contain at least two elements. If you don't specify N_LEAVES, it is
 set to be the largest leaf index occuring in the sets.
 Partitions needn't be irredundant and coefficients can be any rational number. If the resulting 
 element is not in the moduli space, an error is thrown.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i150" onclick="swap_content('span:i150'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i150" onclick="swap_content('span:i150'); return false;"><b>N_LEAVES</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i150"><div class="descr"><p> The number of leaves of the rational curve.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i155" onclick="swap_content('span:i155'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i155" onclick="swap_content('span:i155'); return false;"><b>SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i155"><div class="descr"><p> A list of partitions of [n] that define the tree of the curve: For each bounded edge we have 
 the corresponding partition of the n leaves. These should be irredundant. If you want to input
 a possibly redundant list, use <a href="#i154" onclick="unfold('span:i154');">INPUT_SETS</a> and <a href="#i149" onclick="unfold('span:i149');">INPUT_COEFFS</a> instead. The number of marked leaves
 should always be given by <a href="#i150" onclick="unfold('span:i150');">N_LEAVES</a>. The sets are subsets of {1,...,n} (NOT {0,..,n-1}!)</p></div></span></div></div></li></span></div></div></li></ul><h4>User Methods of RationalCurve</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i102" onclick="swap_content('span:i102'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i102" onclick="swap_content('span:i102'); return false;"><b>VISUAL</b></a>
			()  

			<span style="display: inline" id="span:i102"><br/><div class="descr_func"><div class="undoc"><p>UNDOCUMENTED</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td>option list:</td><td>Visual::RationalCurve::decorations</td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i103" onclick="swap_content('span:i103'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i103" onclick="swap_content('span:i103'); return false;"><h3><div class="category">Tropical geometry</div></h3></a><span style="display: inline" id="span:i103"><li><div class="li"><div class="icon" id="icon:i146" onclick="swap_content('span:i146'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i146" onclick="swap_content('span:i146'); return false;"><b>matroid_vector</b></a>
			()  

			
				→ Vector&lt;Rational&gt;<span style="display: inline" id="span:i146"><br/><div class="descr_func"><div class="descr"><p> Computes the coordinates of the curve in the moduli space given in matroid coordinates.</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i147" onclick="swap_content('span:i147'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i147" onclick="swap_content('span:i147'); return false;"><b>metric_vector</b></a>
			()  

			<span style="display: inline" id="span:i147"><br/><div class="descr_func"><div class="descr"><p> Returns the (n over 2) metric vector of the rational n-marked curve</p></div></div><div class="level3"/></span></div></div></li></span></div></div></li></ul></span></div></div></li><li><div class="li"><div class="icon" id="icon:i23" onclick="swap_content('span:i23'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i23" onclick="swap_content('span:i23'); return false;"><h3><div class="object">RationalFunction</div></h3></a><span style="display: inline" id="span:i23"><div class="undoc"><p>UNDOCUMENTED</p></div><div class="level3"/><h4>Properties of RationalFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i107" onclick="swap_content('span:i107'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i107" onclick="swap_content('span:i107'); return false;"><h3><div class="category">Tropical geometry</div></h3></a><span style="display: inline" id="span:i107"><li><div class="li"><div class="icon" id="icon:i160" onclick="swap_content('span:i160'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i160" onclick="swap_content('span:i160'); return false;"><b>DOMAIN</b>:
		</a><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a><br/><span style="display: inline" id="span:i160"><div class="descr"><p> This property describes the linearity domains of the function. I.e. the function is
 affine integral linear on each maximal cone of <a href="#i160" onclick="unfold('span:i160');">DOMAIN</a>. It is of type WeightedComplex, but
 actually only the basic geometric properties <a href="#" onclick="unfold('span:');">RAYS</a>/<a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>,<a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>/<a href="#i123" onclick="unfold('span:i123');">CMPLX_MAXIMAL_CONES</a>,<a href="#" onclick="unfold('span:');">LINEALITY_SPACE</a> and <a href="#i120" onclick="unfold('span:i120');">USES_HOMOGENEOUS_C</a> are of any interest. Further properties are ignored and can be omitted upon generation.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i161" onclick="swap_content('span:i161'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i161" onclick="swap_content('span:i161'); return false;"><b>LIN_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i161"><div class="descr"><p> The value at index i describes the function value of <a href="#i160" onclick="unfold('span:i160');">DOMAIN</a>-&gt;LINEALITY_SPACE-&gt;row(i)</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i159" onclick="swap_content('span:i159'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i159" onclick="swap_content('span:i159'); return false;"><b>RAY_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><span style="display: inline" id="span:i159"><div class="descr"><p> The value at index i describes the function value of ray <a href="#i160" onclick="unfold('span:i160');">DOMAIN</a>-&gt;CMPLX_RAYS-&gt;row(i). More precisely, if the
 corresponding ray is a vertex, it describes its function value. If it is a directional ray, it describes 
 the slope on that ray.</p></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i106" onclick="swap_content('span:i106'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i106" onclick="swap_content('span:i106'); return false;"><h3><div class="category">Tropical geometry / Rational functions</div></h3></a><span style="display: inline" id="span:i106"><li><div class="li"><div class="icon" id="icon:i158" onclick="swap_content('span:i158'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i158" onclick="swap_content('span:i158'); return false;"><b>POWER</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i158"><div class="descr"><p> This property describes how many times a function should be applied, if it occurs in a divisor term
 E.g. if you create a MinMaxFunction f and you set its POWER property to 2, a call to
 divisor(X,f) actually computes f*f*X</p></div></span></div></div></li></span></div></div></li></ul><h4>User Methods of RationalFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i105" onclick="swap_content('span:i105'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i105" onclick="swap_content('span:i105'); return false;"><h3><div class="category">Rational functions</div></h3></a><span style="display: inline" id="span:i105"><li><div class="li"><div class="icon" id="icon:i157" onclick="swap_content('span:i157'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i157" onclick="swap_content('span:i157'); return false;"><b>homogenize</b></a>
			()  

			<span style="display: inline" id="span:i157"><br/><div class="descr_func"><div class="descr"><p> If the domain of the function is in homogeneous coordinates, it returns the function. Otherwise it homogenizes
 the domain and assigns the value 0 to the additional vertex at the origin</p></div></div><div class="level3"/></span></div></div></li><li><div class="li"><div class="icon" id="icon:i156" onclick="swap_content('span:i156'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i156" onclick="swap_content('span:i156'); return false;"><b>pow</b></a>
			(k)  

			
				→ Returns<span style="display: inline" id="span:i156"><br/><div class="descr_func"><div class="descr"><p> This creates a RationalFunction from this one, simply by changing the power</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td>Returns</td><td><div class="descr"><div>the function f with a power set to k. Note that the original power of f is ignored.</div></div></td></tr></table></div></span></div></div></li></span></div></div></li></ul></span></div></div></li><li><div class="li"><div class="icon" id="icon:i20" onclick="swap_content('span:i20'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i20" onclick="swap_content('span:i20'); return false;"><h3><div class="object">WeightedComplex</div></h3></a><span style="display: inline" id="span:i20"><div class="undoc"><p>UNDOCUMENTED</p></div><b>derived from: </b>fan::PolyhedralFan&lt;Rational&gt;<br/><div class="level3"/><h4>Properties of WeightedComplex</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i96" onclick="swap_content('span:i96'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i96" onclick="swap_content('span:i96'); return false;"><h3><div class="category">Basic properties</div></h3></a><span style="display: inline" id="span:i96"><li><div class="li"><div class="icon" id="icon:i128" onclick="swap_content('span:i128'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i128" onclick="swap_content('span:i128'); return false;"><b>CODIM_1_FACES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i128"><div class="descr"><p> An incidence matrix of all codimension one faces. Each row is a face and column indices refer to <a href="#" onclick="unfold('span:');">RAYS</a>.
 If the fan uses homogeneous coordinates, only those faces that intersect {x0 = 1} are listed here.
 If there is a <a href="#i134" onclick="unfold('span:i134');">LOCAL_RESTRICTION</a>, only compatible codimension one faces are kept.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i130" onclick="swap_content('span:i130'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i130" onclick="swap_content('span:i130'); return false;"><b>CODIM_1_IN_MAXIMAL_CONES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i130"><div class="descr"><p> An incidence matrix indicating which codim 1 faces are contained in which maximal cone. Rows refer
 to (rows of) <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a>, columns to (rows of) <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i129" onclick="swap_content('span:i129'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i129" onclick="swap_content('span:i129'); return false;"><b>DESCRIPTION</b>:
		</a>common::String<br/><span style="display: inline" id="span:i129"><div class="descr"><p> An arbitrary textual description of the complex, to be specified by the user. Some
 operations that produce a WeightedComplex automatically set this property to something
 informative.</p></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i94" onclick="swap_content('span:i94'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i94" onclick="swap_content('span:i94'); return false;"><h3><div class="category">Polyhedral complex</div></h3></a><span style="display: inline" id="span:i94"><li><div class="li"><div class="icon" id="icon:i122" onclick="swap_content('span:i122'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i122" onclick="swap_content('span:i122'); return false;"><b>CMPLX_AMBIENT_DIM</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i122"><div class="descr"><p> Returns the ambient dimension of the polyhedral complex. In the case of non-homogeneous coordinates, this is 
 just FAN_AMBIENT_DIM, otherwise FAN_AMBIENT_DIM-1</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i112" onclick="swap_content('span:i112'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i112" onclick="swap_content('span:i112'); return false;"><b>CMPLX_CODIMENSION</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i112"><div class="descr"><p> Returns the codimension of the polyhedral complex in its ambient vector space. More precisely, this
 value is equal to <a href="#i122" onclick="unfold('span:i122');">CMPLX_AMBIENT_DIM</a> - <a href="#i114" onclick="unfold('span:i114');">CMPLX_DIM</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i121" onclick="swap_content('span:i121'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i121" onclick="swap_content('span:i121'); return false;"><b>CMPLX_CODIM_1_FACES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i121"><div class="descr"><p> An incidence matrix describing which codim 1 cone in the complex at x0 = 1 is generated by which rays. 
 Each row corresponds to a codimension one cone (More precisely, the i-th element represents the
 same codim 1 cone as the i-th element of <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a>). The indices in a row refer to rows of <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>, i.e. the cone
 cone described by the i-th element is generated by the rays corresponding to these row indices.
 If the fan does not use homogeneous coordinates, it is just a copy of <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a>.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i113" onclick="swap_content('span:i113'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i113" onclick="swap_content('span:i113'); return false;"><b>CMPLX_CONVERSION_VECTOR</b>:
		</a>common::Vector&lt;Int&gt;<br/><span style="display: inline" id="span:i113"><div class="descr"><p> A vector with an entry for each row in <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>. More precisely, the i-th entry gives the row index
 of the ray in <a href="#" onclick="unfold('span:');">RAYS</a> that is equal to the i-th row of <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i114" onclick="swap_content('span:i114'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i114" onclick="swap_content('span:i114'); return false;"><b>CMPLX_DIM</b>:
		</a>common::Int<br/><span style="display: inline" id="span:i114"><div class="descr"><p> Returns the dimension of the polyhedral complex. In the case of non-homogeneous coordinates, this is 
 just DIM, otherwise DIM-1. When computing this property, it is assumed that the complex is pure.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i123" onclick="swap_content('span:i123'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i123" onclick="swap_content('span:i123'); return false;"><b>CMPLX_MAXIMAL_CONES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i123"><div class="descr"><p> An incidence matrix describing which maximal cone in the complex at x0 = 1 is generated by which rays. 
 Each row corresponds to a maximal cone (More precisely, the i-th element represents the
 same maximal cone as the i-th element of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>). The indices in a row refer to rows of <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>, i.e. the maximal
 cone described by the i-th element is generated by the rays corresponding to these row indices.
 If the fan does not use homogeneous coordinates, it is just a copy of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i115" onclick="swap_content('span:i115'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i115" onclick="swap_content('span:i115'); return false;"><b>CMPLX_RAYS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i115"><div class="descr"><p> If the fan does not use homogeneous coordinates, this is just a copy of RAYS and actually needn't be requested. Otherwise this computes
 a matrix of rays of the complex obtained by intersecting the fan with {x0 = 1}. More precisely, each ray r from <a href="#" onclick="unfold('span:');">RAYS</a> occurs  
 as a row in this matrix...
 - once, if r_0 = 1
 - k times, if r_0 = 0 and k is the number of equivalence classes of maximal cones containing r with respect to the following relation:
 Two maximal cones m, m' containing r are equivalent, if they are equal or there exists a sequence of maximal cones
 m = m_1,...m_r = m', such that r is contained in each m_i and each intersection m_i cap m_i+1 contains at least one ray s with s_0 = 1.
 The reason for this is that, when specifying a piecewise affine linear function on a polyhedral complex, the same directional ray with 
 x0 = 0 might be assigned two different values, if it is contained in two "non-connected" maximal cones (where connectedness is to be 
 understood as described above).
 If there is a <a href="#i134" onclick="unfold('span:i134');">LOCAL_RESTRICTION</a> the above equivalence relation is changed in such a 
 way that the affine ray s with s_0 = 1 that must be contained in the intersection of
 two subsequent cones must be a compatible ray. </p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i118" onclick="swap_content('span:i118'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i118" onclick="swap_content('span:i118'); return false;"><b>DIRECTIONAL_RAYS</b>:
		</a>common::Set&lt;Int&gt;<br/><span style="display: inline" id="span:i118"><div class="descr"><p> Computes the subset of row indices of RAYS, such that the corresponding rows describe a directional ray of the complex.
 This set contains all row indices  if <a href="#i120" onclick="unfold('span:i120');">USES_HOMOGENEOUS_C</a> is false and is the complement of <a href="#i117" onclick="unfold('span:i117');">VERTICES</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i119" onclick="swap_content('span:i119'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i119" onclick="swap_content('span:i119'); return false;"><b>INPUT_HOM_LINEALITY</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i119"><div class="descr"><p> A matrix of generators of the lineality space, but in homogeneous coordinates (i.e. x0 should be 0, otherwise 
 this throws an error). This should be used *instead of* <a href="#" onclick="unfold('span:');">INPUT_LINEALITY</a> / <a href="#" onclick="unfold('span:');">LINEALITY_SPACE</a></p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i116" onclick="swap_content('span:i116'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i116" onclick="swap_content('span:i116'); return false;"><b>INPUT_HOM_RAYS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i116"><div class="descr"><p> A matrix of the rays of the fan. This should be used *instead of* <a href="#" onclick="unfold('span:');">INPUT_RAYS</a>, if the coordinates of the
 fan should be interpreted as homogeneous coordinates and the fan as the complex at x0 = 1.
 The input rays are normalized to x0 = 1 and the fan is then initialized with the normalized rays as
 <a href="#" onclick="unfold('span:');">INPUT_RAYS</a>.
 <a href="#" onclick="unfold('span:');">INPUT_CONES</a> can be combined with this as usual. If your rays are already homogenized and you are sure they
 are also your rays, you can directly put this into RAYS and manually set USES_HOMOGENEOUS_C=&gt;TRUE</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i120" onclick="swap_content('span:i120'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i120" onclick="swap_content('span:i120'); return false;"><b>USES_HOMOGENEOUS_C</b>:
		</a>common::Bool<br/><span style="display: inline" id="span:i120"><div class="descr"><p> Indicates whethere this fan's coordinates should be interpreted as affine coordinates or as homogeneous coordinates.
 In the latter case the fan would actually represent the polyhedral complex that comes from intersecting the fan
 (in affine coordinates) with the affine hyperplane x0 = 1. This property should not be set manually. It will be set
 as true, if <a href="#i116" onclick="unfold('span:i116');">INPUT_HOM_RAYS</a> or <a href="#i119" onclick="unfold('span:i119');">INPUT_HOM_LINEALITY</a> are used and to false, if <a href="#" onclick="unfold('span:');">INPUT_RAYS</a> or <a href="#" onclick="unfold('span:');">INPUT_LINEALITY</a>
 are used.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i117" onclick="swap_content('span:i117'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i117" onclick="swap_content('span:i117'); return false;"><b>VERTICES</b>:
		</a>common::Set&lt;Int&gt;<br/><span style="display: inline" id="span:i117"><div class="descr"><p> Computes the subset of row indices of RAYS, such that the corresponding rows describe a vertex of the complex.
 This set is empty if <a href="#i120" onclick="unfold('span:i120');">USES_HOMOGENEOUS_C</a> is false and is the complement of <a href="#i118" onclick="unfold('span:i118');">DIRECTIONAL_RAYS</a></p></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i95" onclick="swap_content('span:i95'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i95" onclick="swap_content('span:i95'); return false;"><h3><div class="category">Tropical geometry</div></h3></a><span style="display: inline" id="span:i95"><li><div class="li"><div class="icon" id="icon:i125" onclick="swap_content('span:i125'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i125" onclick="swap_content('span:i125'); return false;"><b>BALANCED_FACES</b>:
		</a>common::Vector&lt;Bool&gt;<br/><span style="display: inline" id="span:i125"><div class="descr"><p> A vector whose entries correspond to the rows of CODIM_1_FACES. The i-th entry is true, if and only 
 if the complex is balanced at that face</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i124" onclick="swap_content('span:i124'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i124" onclick="swap_content('span:i124'); return false;"><b>LATTICE_NORMAL_FCT_VECTOR</b>:
		</a>common::Map&lt;Int, Map&lt;Int, Vector&lt;Rational&gt;&gt;&gt;<br/><span style="display: inline" id="span:i124"><div class="descr"><p> For each lattice normal vector, this gives a vector of length (number of rays) + (lineality dim.), such that
 if a rational function is given by values on the rays and lin space generators, the value of the corresponding
 normal <a href="#i131" onclick="unfold('span:i131');">LATTICE_NORMALS</a>-&gt;{i}-&gt;{j} can be computed by multiplying the function value vector with 
 the vector <a href="#i124" onclick="unfold('span:i124');">LATTICE_NORMAL_FCT_VECTOR</a>-&gt;{i}-&gt;{j}. This is done in the following way:  If this fan does not 
 use homogeneous coordinates, it computes a representation of the lattice normal in the generating system 
 consisting of the rays of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>-&gt;row(j) and the generators 
 of the lineality space. It then inserts the coefficients of the representation at the right position. 
 If the fan uses homog. coordinates, however, we use a different generating system (and indices refer to CMPLX_RAYS)
 &lt;(r_i-r_0)_i&gt;0, s_j, l_k&gt;, where r_0 is the ray of the maximal cone 
 with the lowest index in <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>, such that it fulfills x0 = 1, r_i are the remaining rays with x0 = 1, ordered
 according to their index in <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a>, s_j are the rays of the cone with x0 = 0 and l_k are the lineality space 
 generators. We will then store the coefficients a_i of (r_i - r_0) at the index of r_i, then - sum(a_i) at 
 the index of r_0 and the remaining coefficients at the appropriate places.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i127" onclick="swap_content('span:i127'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i127" onclick="swap_content('span:i127'); return false;"><b>LATTICE_NORMAL_SUM</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i127"><div class="descr"><p> Rows of this matrix correspond to <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a>, and each row contains the weighted sum:
 sum_{cone &gt; codim-1-face}( weight(cone) * <a href="#i131" onclick="unfold('span:i131');">LATTICE_NORMALS</a>-&gt;{codim-1-face}-&gt;{cone})</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i126" onclick="swap_content('span:i126'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i126" onclick="swap_content('span:i126'); return false;"><b>LATTICE_NORMAL_SUM_FCT_VECTOR</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><span style="display: inline" id="span:i126"><div class="descr"><p> Rows of this matrix correspond to <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a> (or <a href="#i121" onclick="unfold('span:i121');">CMPLX_CODIM_1_FACES</a> in the homog. case), and each row 
 contains a function vector for the corresponding row of <a href="#i127" onclick="unfold('span:i127');">LATTICE_NORMAL_SUM</a>. This function vector is computed 
 in the same way as described under <a href="#i124" onclick="unfold('span:i124');">LATTICE_NORMAL_FCT_VECTOR</a>.
 Note that for any codim-1-faces at which the complex is not balanced, the corresponding row is a zero row
 If a face is balanced can be checked under <a href="#i125" onclick="unfold('span:i125');">BALANCED_FACES</a>.</p></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i98" onclick="swap_content('span:i98'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i98" onclick="swap_content('span:i98'); return false;"><h3><div class="category">Tropical geometry/Local computations</div></h3></a><span style="display: inline" id="span:i98"><li><div class="li"><div class="icon" id="icon:i134" onclick="swap_content('span:i134'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i134" onclick="swap_content('span:i134'); return false;"><b>LOCAL_RESTRICTION</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><span style="display: inline" id="span:i134"><div class="descr"><p> This contains a list of sets of ray indices (referring to RAYS). All of these sets should 
 describe cones of the polyhedral complex (though not necessarily maximal ones).
 A cone is now called compatible with this property, if contains one of these cones or is 
 contained in one of these cones. If this list is not empty, all computations will be
 done only on (or around) compatible cones. The documentation of each property will explain
 in what way this restriction is enforced. 
 If this list is empty or not defined, there is no restriction.
 Careful: The implementation assumes that ALL maximal cones are compatible. If in doubt,
 you can create a complex with a local restriction from a given complex by using one of
 the "local_..." creation methods</p></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i97" onclick="swap_content('span:i97'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i97" onclick="swap_content('span:i97'); return false;"><h3><div class="category">Tropical properties</div></h3></a><span style="display: inline" id="span:i97"><li><div class="li"><div class="icon" id="icon:i133" onclick="swap_content('span:i133'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i133" onclick="swap_content('span:i133'); return false;"><b>IS_BALANCED</b>:
		</a>common::Bool<br/><span style="display: inline" id="span:i133"><div class="descr"><p> Indicates whether the fan is balanced with the given <a href="#i132" onclick="unfold('span:i132');">TROPICAL_WEIGHTS</a>. If there is a 
 <a href="#i134" onclick="unfold('span:i134');">LOCAL_RESTRICTION</a>, this will only be checked at the compatible codim one faces.</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i131" onclick="swap_content('span:i131'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i131" onclick="swap_content('span:i131'); return false;"><b>LATTICE_NORMALS</b>:
		</a>common::Map&lt;Int, Map&lt;Int, Vector&lt;Integer&gt;&gt;&gt;<br/><span style="display: inline" id="span:i131"><div class="descr"><p> A list of lattice normal vectors. They are associated to the faces and maximal cones in the following 
 way: For codimension one face i and adjacent maximal cone j, the element <a href="#i131" onclick="unfold('span:i131');">LATTICE_NORMALS</a>-&gt;{i}-&gt;{j} is the corresp.
 lattice normal vector. There is a convenience method lattice_normal(i,j)
 If <a href="#i120" onclick="unfold('span:i120');">USES_HOMOGENEOUS_C</a> is TRUE, this actually computes a lattice normal of the cone (face i)\cap {x0 =1} wrt the
 cone (maximal cone j) \cap {x0=1} to ensure that the normal is of the form (0,...)
 If there is a <a href="#i134" onclick="unfold('span:i134');">LOCAL_RESTRICTION</a>, only lattice normals for compatible codim one faces
 are computed (naturally, since <a href="#i128" onclick="unfold('span:i128');">CODIM_1_FACES</a> only contains those)</p></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i132" onclick="swap_content('span:i132'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i132" onclick="swap_content('span:i132'); return false;"><b>TROPICAL_WEIGHTS</b>:
		</a>common::Vector&lt;Integer&gt;<br/><span style="display: inline" id="span:i132"><div class="descr"><p> Vector of (integer) weights for maximal cones. Indices refer to (the rows of) <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>
 Note that, if the fan only consists of a lineality space, it specifies the weight of the single empty maximal
 cone </p></div></span></div></div></li></span></div></div></li></ul><h4>User Methods of WeightedComplex</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i92" onclick="swap_content('span:i92'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i92" onclick="swap_content('span:i92'); return false;"><b>BB_VISUAL</b></a>
			()  

			<span style="display: inline" id="span:i92"><br/><div class="descr_func"><div class="descr"><p> ---------------------------------------------------------------------------------------------------------</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td>option list:</td><td>Visual::WeightedComplex::BoundingDecorations</td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i93" onclick="swap_content('span:i93'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i93" onclick="swap_content('span:i93'); return false;"><b>characteristic</b></a>
			()  

			<span style="display: inline" id="span:i93"><br/><div class="descr_func"><div class="descr"><blockquote><div>user_method adjHomology {</div><div>  my $this = shift;</div><div>  #Create set graph</div><div>  my $facets = new Vector&lt;Set&lt;Int&gt; &gt;();</div><div>  my $maximalInCodim = $this-&gt;CODIM_1_IN_MAXIMAL_CONES;</div><div>  my $codimInMaximal = transpose($this-&gt;CODIM_1_IN_MAXIMAL_CONES);</div><div>  for my $mc (0 .. $this-&gt;MAXIMAL_CONES-&gt;rows() -1) {</div><div>      print "Finding adjacencies of cone $mc...\n";</div><div>      my $cdset = $codimInMaximal-&gt;row($mc);</div><div>      for my $cd (@{$cdset}) {</div><div>	my $mcset = ($maximalInCodim-&gt;row($cd) - $mc);</div><div>	for my $othermc (@{$mcset}) {</div><div>	  if($othermc &gt; $mc) {</div><div>	    my $s = new Set&lt;Int&gt;($mc,$othermc);</div><div>	    $facets = $facets | $s;</div><div>	  }</div><div>	}</div><div>      }</div><div>  }</div><div>  my $cplx = new topaz::SimplicialComplex(FACETS=&gt;$facets);</div><div>  $cplx-&gt;VISUAL;</div><div>  return $cplx-&gt;HOMOLOGY;</div><div>}</div></blockquote></div></div><div class="level3"/></span></div></div></li><li><div class="li"><div class="icon" id="icon:i90" onclick="swap_content('span:i90'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i90" onclick="swap_content('span:i90'); return false;"><h3><div class="category">Polyhedral geometry</div></h3></a><span style="display: inline" id="span:i90"><li><div class="li"><div class="icon" id="icon:i110" onclick="swap_content('span:i110'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i110" onclick="swap_content('span:i110'); return false;"><b>homogenize</b></a>
			()  

			
				→ If<span style="display: inline" id="span:i110"><br/><div class="descr_func"><div class="descr"><p> Transforms the fan into homogeneous coordinates</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td>If</td><td><div class="descr"><div>this fan uses homog. coordinates, returns the fan itself, Otherwise it will replace the ray matrix
 by the following matrix: [0,ray] for any ray of the original fan and one affine ray [1,0,...0] appended. In <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>
 all get the additional affine ray and <a href="#" onclick="unfold('span:');" class="invalid">TROPICAL WEIGHTS</a> is simply copied</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i91" onclick="swap_content('span:i91'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i91" onclick="swap_content('span:i91'); return false;"><h3><div class="category">Tropical methods</div></h3></a><span style="display: inline" id="span:i91"><li><div class="li"><div class="icon" id="icon:i111" onclick="swap_content('span:i111'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i111" onclick="swap_content('span:i111'); return false;"><b>lattice_normal</b></a>
			(i, j)  

			
				→ The<span style="display: inline" id="span:i111"><br/><div class="descr_func"><div class="descr"><p> Convenience method, returns LATTICE_NORMALS-&gt;{i}-&gt;{j}</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">i</td><td><div class="descr"><div>The index of the facet</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">j</td><td><div class="descr"><div>The index of the maximal cone</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>lattice normal vector of cone j wrt. facet i</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i89" onclick="swap_content('span:i89'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i89" onclick="swap_content('span:i89'); return false;"><h3><div class="category">Visualization</div></h3></a><span style="display: inline" id="span:i89"><li><div class="li"><div class="icon" id="icon:i109" onclick="swap_content('span:i109'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i109" onclick="swap_content('span:i109'); return false;"><b>boundingBox</b></a>
			(distance)  

			
				→ Matrix&lt;Rational&gt;<span style="display: inline" id="span:i109"><br/><div class="descr_func"><div class="descr"><p> Takes a tropical variety and a (positive) Rational as input and computes the relative bounding box of the visualization with the Rational as BoundingDistance (see also the description of <a href="#i92" onclick="unfold('span:i92');">BB_VISUAL</a></p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">distance</td><td><div class="descr"><div>The distance of the bounding box from the affine part of the variety (Can also be a standard perl number type)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td><div class="descr"><div>A 2x<a href="#i122" onclick="unfold('span:i122');">CMPLX_AMBIENT_DIM</a>-matrix specifying the bounding box (it can be used as input for BoundingBox in <a href="#i92" onclick="unfold('span:i92');">BB_VISUAL</a></div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i108" onclick="swap_content('span:i108'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i108" onclick="swap_content('span:i108'); return false;"><b>CMPLX_VISUAL</b></a>
			()  

			<span style="display: inline" id="span:i108"><br/><div class="descr_func"><div class="descr"><p> Displays a weighted polyhedral complex in ambient dimension up to 3 (4 in homog. coordinates) by visualizing all its
 cells in the following way: For each cell, every directional ray is added to every affine ray (which is only the origin
 in the non-homog. case) and the resulting polytope is displayed. The visualization has the following options:
 * DirScale: A rational number with which the directional rays are multiplied before added to the affine rays
 * WeightLabels:  if set to "hidden", the labels indicating the weights are hidden 
 * CoordLabels: If set to "show", the labels indicating the vertex coordinates are displayed, otherwise they are not. Note that this is expensive and significantly increases computation time.</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td>option list:</td><td>Visual::WeightedComplex::decorations</td></tr></table></div></span></div></div></li></span></div></div></li></ul></span></div></div></li></ul></div><div class="level2"><h2>User Functions</h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i13" onclick="swap_content('span:i13'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i13" onclick="swap_content('span:i13'); return false;"><h3><div class="category">Arithmetic</div></h3></a><span style="display: inline" id="span:i13"><li><div class="li"><div class="icon" id="icon:i71" onclick="swap_content('span:i71'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i71" onclick="swap_content('span:i71'); return false;"><b>gcdext</b></a>
			(a, b, s, t)  

			
				→ Integer<span style="display: inline" id="span:i71"><br/><div class="descr_func"><div class="descr"><p> Computes the gcd of <em class="param">a</em> and <em class="param">b</em> and returns it. <em class="param">s</em> and <em class="param">t</em> are set
 such that gcd(<em class="param">a</em>,<em class="param">b</em>) = <em class="param">s</em> * <em class="param">a</em> + <em class="param">t</em> * <em class="param">b</em></p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">a</td><td><div class="descr"><div>first argument of gcd(,)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">b</td><td><div class="descr"><div>second argument of gcd(,)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">s</td><td><div class="descr"><div>coefficient of a (will be set)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">t</td><td><div class="descr"><div>coefficient of b (will be set)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td/></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i10" onclick="swap_content('span:i10'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i10" onclick="swap_content('span:i10'); return false;"><h3><div class="category">Graph theory</div></h3></a><span style="display: inline" id="span:i10"><li><div class="li"><div class="icon" id="icon:i32" onclick="swap_content('span:i32'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i32" onclick="swap_content('span:i32'); return false;"><b>spanning_complete</b></a>
			(n)  

			
				→ perl::ListReturn<span style="display: inline" id="span:i32"><br/><div class="descr_func"><div class="descr"><p> Computes the set of all spanning trees of the complete graph on n nodes</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The number of nodes</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>perl::ListReturn</td><td><div class="descr"><div>An array of Set&lt;Int&gt; The spanning trees, each encoded as
 the set of edges, where the edge indices are determined as follows: For i,j =
 0,..,n-1 and i &lt; j we sort the edges (i,j) in ascending lexicographical ordering
 , i.e, (0,1) &lt; (0,2) &lt; .. &lt; (n-2,n-1).</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i8" onclick="swap_content('span:i8'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i8" onclick="swap_content('span:i8'); return false;"><h3><div class="category">Homology</div></h3></a><span style="display: inline" id="span:i8"><li><div class="li"><div class="icon" id="icon:i30" onclick="swap_content('span:i30'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i30" onclick="swap_content('span:i30'); return false;"><b>adjacencyComplex</b></a>
			(fan)  

			
				→ topaz::SimplicialComplex<span style="display: inline" id="span:i30"><br/><div class="descr_func"><div class="descr"><p> Computes the codimension one adjacency graph of a tropical variety and returns it as a 
 simplicial complex. The nodes of the graph are the maximal cones, which are connected, iff 
 they share a codimension one face</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A tropical variety</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">topaz::SimplicialComplex</a></td><td><div class="descr"><div>The adjacency graph as a simplicial complex</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i29" onclick="swap_content('span:i29'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i29" onclick="swap_content('span:i29'); return false;"><b>equivalenceComplex</b></a>
			(fan)  

			
				→ topaz::SimplicialComplex<span style="display: inline" id="span:i29"><br/><div class="descr_func"><div class="descr"><p> Computes the canonical complex of a tropical variety. The one-skeleton is determined 
 in the following way: Two cones a,b are equivalent, if there exists a sequence of 
 maximal cones a=c0,...,cr=b, such that ci,ci+1 intersect in a codimension one face of 
 which they are the only neighbours
 The vertices are the equivalence classes and two of them are connected if any cones in 
 them are connected.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A tropical variety</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">topaz::SimplicialComplex</a></td><td><div class="descr"><div>The canonical complex of the variety.</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i14" onclick="swap_content('span:i14'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i14" onclick="swap_content('span:i14'); return false;"><h3><div class="category">Linear algebra</div></h3></a><span style="display: inline" id="span:i14"><li><div class="li"><div class="icon" id="icon:i76" onclick="swap_content('span:i76'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i76" onclick="swap_content('span:i76'); return false;"><b>affineTransform</b></a>
			(complex, translate, matrix)  

			
				→ WeightedComplex<span style="display: inline" id="span:i76"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral complex and applies an affine linear transformation, given by a translate 
 vector and a matrix. The method assumes the function is bijective and preserves cones, i.e. it 
 just applies the transformation to the rays and lineality space and leaves the cones and weights
 unchanged.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>The complex to be transformed, supposed to be in homogeneous 
 coordinates (if translate != 0)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">translate</td><td><div class="descr"><div>A vector whose dimension should be equal to the column dimension
 of the transformation matrix</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">matrix</td><td><div class="descr"><div>An integer (square) invertible matrix. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The transformed complex</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i74" onclick="swap_content('span:i74'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i74" onclick="swap_content('span:i74'); return false;"><b>computeMatrixBases</b></a>
			(m)  

			
				→ Set&lt;int&gt;<span style="display: inline" id="span:i74"><br/><div class="descr_func"><div class="descr"><p> Computes a list of sets of column indices of a matrix such that
 the corresponding column sets form a basis of the column space</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Set&lt;int&gt;</a></td><td><div class="descr"><div>An array of Set&lt;int&gt;</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i73" onclick="swap_content('span:i73'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i73" onclick="swap_content('span:i73'); return false;"><b>computeMatrixColoops</b></a>
			(m)  

			
				→ Set&lt;int&gt;<span style="display: inline" id="span:i73"><br/><div class="descr_func"><div class="descr"><p> Computes a list of column indices of a matrix such that 
 the corresponding columns are contained in every column basis of the
 column space. Is equal to the intersection of all sets returned by
 <a href="#i74" onclick="unfold('span:i74');">computeMatrixBases</a></p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Set&lt;int&gt;</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i75" onclick="swap_content('span:i75'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i75" onclick="swap_content('span:i75'); return false;"><b>linearRepresentation</b></a>
			(v, generators)  

			
				→ Vector<span style="display: inline" id="span:i75"><br/><div class="descr_func"><div class="descr"><p> Takes a vector v and a matrix with column dimension equal to the dimension of v. Assuming that 
 v is in the row span of the matrix, it computes one(!) possible representation of v in these 
 generators</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td class="param">v</td><td><div class="descr"><div>The vector supposed to be contained in the row span of the generators</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">generators</td><td><div class="descr"><div>A set of row vectors whose linear span should contain v</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td><div class="descr"><div>A vector (a1,..,an) such that v = (a1,...,an) * generators. It returns a vector of
 dimension 0, if v is not in the span of the generators. An error is thrown if the dimensions of 
 v and the generators mismatch</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i72" onclick="swap_content('span:i72'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i72" onclick="swap_content('span:i72'); return false;"><b>lllHNF</b></a>
			(matrix, tfmatrix, kdim)  

			
				→ Matrix<span style="display: inline" id="span:i72"><br/><div class="descr_func"><div class="descr"><p> Computes a row HNF of an integer matrix A
 It returns the normal form and stores the unimodular transformation matrix and the kernel dimension  of the transposed matrix in
 the last two parameters. The algorithm is the LLL-based HNF alg. by Havas, Majevski, Matthews</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">matrix</td><td><div class="descr"><div>the matrix for which the transformation is computed</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">tfmatrix</td><td><div class="descr"><div>The matrix that will contain the transformation matrix</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">kdim</td><td><div class="descr"><div>This will be set to dim Ker(A)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td/></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i6" onclick="swap_content('span:i6'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i6" onclick="swap_content('span:i6'); return false;"><h3><div class="category">Polyhedral geometry</div></h3></a><span style="display: inline" id="span:i6"><li><div class="li"><div class="icon" id="icon:i25" onclick="swap_content('span:i25'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i25" onclick="swap_content('span:i25'); return false;"><b>complexify</b></a>
			(rays, cones, weights, uses_homog)  

			
				→ perl::Object<span style="display: inline" id="span:i25"><br/><div class="descr_func"><div class="descr"><p> Takes a collection of weighted cones of the same dimension and refines them such that they 
 form a weighted complex. Weights of cones lying one over the other add up.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rays</td><td><div class="descr"><div>A matrix of rays of the cones. Has to be irredundant and 
 normalized.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Set&lt;int&gt; &gt;</a></td><td class="param">cones</td><td><div class="descr"><div>The maximal cones, given in terms of row indices of their 
 rays </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Integer&gt;</a></td><td class="param">weights</td><td><div class="descr"><div>The i-th element is the weight of the cone described by the 
 i-th row of cones</div></div></td></tr><tr><td>bool</td><td class="param">uses_homog</td><td><div class="descr"><div>Whether the rays are given in homogeneous coordinates </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>perl::Object</td><td><div class="descr"><div>A WeightedComplex object, whose support is the union of the cones. 
 The weight of a cone is the sum of the weights of the orginal cones containing it. The 
 result is homogeneous iff uses_homog is true</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i26" onclick="swap_content('span:i26'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i26" onclick="swap_content('span:i26'); return false;"><b>product_complex</b></a>
			(A)  

			
				→ the<span style="display: inline" id="span:i26"><br/><div class="descr_func"><div class="descr"><p> Computes the cartesian product of a set of polyhedral fans. If any of them uses homogeneous coordinates, so will the result</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>complexes</td><td class="param">A</td><td><div class="descr"><div>list of WeightedComplex objects</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>the</td><td><div class="descr"><div>cartesian product of the complexes. If any of the complexes use homogeneous coordinates, so does the result. If any of the complexex carries weights, the result carries the sum of the weights of its factors (Complexes without weights are treated as complexes with weight 1).</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i27" onclick="swap_content('span:i27'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i27" onclick="swap_content('span:i27'); return false;"><b>skeleton_complex</b></a>
			(fan, k, preserveRays)  

			
				→ The<span style="display: inline" id="span:i27"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral complex and computes the k-skeleton. Will return an empty fan, if k is larger then the dimension of the given complex or smaller than 0.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A fan (or polyhedral complex)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the skeleton that should be computed</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">preserveRays</td><td><div class="descr"><div>When true, the function assumes that all rays of the fan remain in the k-skeleton, so it just copies the RAYS, instead of computing an irredundant list. This property can always be set to true, if fan is not in homogeneous coordinates or if the corresponding complex at x0 = 1 only has vertices. By default, this property is false.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>k-skeleton of the fan (or complex, if USES_HOMOGENEOUS_C is true)</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i7" onclick="swap_content('span:i7'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i7" onclick="swap_content('span:i7'); return false;"><h3><div class="category">Random number generators</div></h3></a><span style="display: inline" id="span:i7"><li><div class="li"><div class="icon" id="icon:i28" onclick="swap_content('span:i28'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i28" onclick="swap_content('span:i28'); return false;"><b>randomInteger</b></a>
			(max_arg, n)  

			
				→ Array&lt;Integer&gt;<span style="display: inline" id="span:i28"><br/><div class="descr_func"><div class="descr"><p> Returns n random integers in the range 0.. (max_arg-1),inclusive
 Note that this algorithm is not optimal for real randomness:
 If you change the range parameter and then change it back, you will
 usually get the exact same sequence</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">max_arg</td><td><div class="descr"><div>The upper bound for the random integers</div></div></td></tr><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The number of integers to be created</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Array&lt;Integer&gt;</a></td><td/></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i15" onclick="swap_content('span:i15'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i15" onclick="swap_content('span:i15'); return false;"><h3><div class="category">Rational functions</div></h3></a><span style="display: inline" id="span:i15"><li><div class="li"><div class="icon" id="icon:i77" onclick="swap_content('span:i77'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i77" onclick="swap_content('span:i77'); return false;"><b>function_mult</b></a>
			(f, r)  

			
				→ RationalFunction<span style="display: inline" id="span:i77"><br/><div class="descr_func"><div class="descr"><p> This function takes a RationalFunction f and a Rational value r and returns the RationalFunction r*f</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a></td><td class="param">f</td><td/></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">r</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i78" onclick="swap_content('span:i78'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i78" onclick="swap_content('span:i78'); return false;"><b>function_value</b></a>
			(X, An)  

			
				→ RationalFunction<span style="display: inline" id="span:i78"><br/><div class="descr_func"><div class="descr"><p> This creates a RationalFunction by directly defining it on the rays of a given domain. It is simply shorter then
 using the constructor.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>The DOMAIN of the function</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">An</td><td><div class="descr"><div>arbitrarily long list of function values. It should be as long as 
 X-&gt;CMPXL_RAYS-&gt;rows + X-&gt;LINEALITY_DIM. Function values will be assigned in this order 
 (first rays, then linspace). Superfluous values are ignored, missing ones replaced by zero</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i79" onclick="swap_content('span:i79'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i79" onclick="swap_content('span:i79'); return false;"><b>wfan</b></a>
			(fan, homog)  

			
				→ WeightedComplex<span style="display: inline" id="span:i79"><br/><div class="descr_func"><div class="descr"><p> This function takes a fan::PolyhedralFan and converts it into a WeightedComplex object.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">fan::PolyhedralFan</a></td><td class="param">fan</td><td><div class="descr"><div>The fan to be converted</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">homog</td><td><div class="descr"><div>Whether the fan is to be interpreted in homogeneous coordinates</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>A complex with no weights but with the same geometry as fan.</div></div></td></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i12" onclick="swap_content('span:i12'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i12" onclick="swap_content('span:i12'); return false;"><h3><div class="category">Testing</div></h3></a><span style="display: inline" id="span:i12"><li><div class="li"><div class="icon" id="icon:i70" onclick="swap_content('span:i70'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i70" onclick="swap_content('span:i70'); return false;"><b>atint_testrun</b></a>
			()  

			<span style="display: inline" id="span:i70"><br/><div class="descr_func"><div class="descr"><p> Tests the extension a-tint by calling (more or less) all functions and methods
 Needs no parameters and returns nothing.</p></div></div><div class="level3"/></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i16" onclick="swap_content('span:i16'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i16" onclick="swap_content('span:i16'); return false;"><h3><div class="category">Tropial geometry</div></h3></a><span style="display: inline" id="span:i16"><li><div class="li"><div class="icon" id="icon:i80" onclick="swap_content('span:i80'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i80" onclick="swap_content('span:i80'); return false;"><b>rational_curve_from_metric</b></a>
			(v)  

			
				→ RationalCurve<span style="display: inline" id="span:i80"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^(n over 2) that describes an n-marked rational abstract
 curve as a distance vector between its leaves. It then computes the 
 curve corresponding to this vector.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector of length (n over 2). Its entries are 
 interpreted as the distances d(i,j) ordered lexicographically according to i,j. However, they need not be positive, as long as v is equivalent to a proper 
 metric modulo leaf lengths.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td/></tr></table></div></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i11" onclick="swap_content('span:i11'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i11" onclick="swap_content('span:i11'); return false;"><h3><div class="category">Tropical geometry</div></h3></a><span style="display: inline" id="span:i11"><li><div class="li"><div class="icon" id="icon:i50" onclick="swap_content('span:i50'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i50" onclick="swap_content('span:i50'); return false;"><b>bergman_fan</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<span style="display: inline" id="span:i50"><br/><div class="descr_func"><div class="descr"><p> Creates the bergman fan of a given matroid fan.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">m</td><td><div class="descr"><div>A matroid</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i40" onclick="swap_content('span:i40'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i40" onclick="swap_content('span:i40'); return false;"><b>bergman_fan_linear</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<span style="display: inline" id="span:i40"><br/><div class="descr_func"><div class="descr"><p> Uses an algorithm by Felipe Rincón to compute the bergman fan of a linear matroid</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td><div class="descr"><div>Any rational matrix</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i44" onclick="swap_content('span:i44'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i44" onclick="swap_content('span:i44'); return false;"><b>bergman_fan_matroid</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<span style="display: inline" id="span:i44"><br/><div class="descr_func"><div class="descr"><p> Uses an algorithm by Felipe Rincón to compute the bergman fan of a general matroid</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">m</td><td><div class="descr"><div>Any matroid</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i33" onclick="swap_content('span:i33'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i33" onclick="swap_content('span:i33'); return false;"><b>count_mn_cones</b></a>
			(n)  

			
				→ Integer<span style="display: inline" id="span:i33"><br/><div class="descr_func"><div class="descr"><p> Computes the number of maximal cones of the tropical moduli space M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be &gt;= 3</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The number of maximal cones</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i36" onclick="swap_content('span:i36'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i36" onclick="swap_content('span:i36'); return false;"><b>count_mn_rays</b></a>
			(n)  

			
				→ Integer<span style="display: inline" id="span:i36"><br/><div class="descr_func"><div class="descr"><p> Computes the number of rays of the tropical moduli space M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be &gt;= 3</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The number of rays</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i65" onclick="swap_content('span:i65'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i65" onclick="swap_content('span:i65'); return false;"><b>curve_graph_from_metric</b></a>
			(v)  

			
				→ An<span style="display: inline" id="span:i65"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^(n over 2) that describes an n-marked rational abstract
 curve as a distance vector between its leaves. It then computes the 
 graph of the curve corresponding to this vector.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector of length (n over 2). Its entries are 
 interpreted as the distances d(i,j) ordered lexicographically according to i,j. However, they need not be positive, as long as v is equivalent to a proper 
 metric modulo leaf lengths.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>An</td><td><div class="descr"><div>array containing first the graph::Graph and then a Vector&lt;Rational&gt;, containing 
 the lengths of the bounded edges (in the order they appear in EDGES)</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i57" onclick="swap_content('span:i57'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i57" onclick="swap_content('span:i57'); return false;"><b>cycle_intersection</b></a>
			(X, Y)  

			
				→ WeightedComplex<span style="display: inline" id="span:i57"><br/><div class="descr_func"><div class="descr"><p> Computes the intersection product of two tropical cycles in a common ambient vector space V.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>The first tropical variety</div></div></td></tr><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>The second tropical variety. Should have the same actual ambient 
 dimension (homogeneous coordinates don't count) as X.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>Z The intersection product X*Y in V = R^n, where n is the ambient 
 dimension of X and Y. The result has homogeneous coordinates in any case</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i49" onclick="swap_content('span:i49'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i49" onclick="swap_content('span:i49'); return false;"><b>degree</b></a>
			(complex)  

			
				→ Int<span style="display: inline" id="span:i49"><br/><div class="descr_func"><div class="descr"><p> Computes the degree of a tropical variety as the degree of the 
 0-dimensional complex obtained when intersecting 
 the variety with an appropriate linear space L^n_k</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">complex</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i64" onclick="swap_content('span:i64'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i64" onclick="swap_content('span:i64'); return false;"><b>diagonal_functions</b></a>
			(n)  

			
				→ MinMaxFunction<span style="display: inline" id="span:i64"><br/><div class="descr_func"><div class="descr"><p> Produces the rational functions that cut out the diagonal of the vector space R^n.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>Them dimension of the vector space (and the number of functions)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i21" onclick="unfold('span:i21');">MinMaxFunction</a></td><td><div class="descr"><div>A list of the functions (max(x_n,y_n),...,max(x_1,y_1)), where x_i is the coordinate function in the first part of the cartesian product R^n x R^n and y_i is the coordinate function of the second half.</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i58" onclick="swap_content('span:i58'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i58" onclick="swap_content('span:i58'); return false;"><b>divisor</b></a>
			(X, f, k)  

			
				→ WeightedComplex<span style="display: inline" id="span:i58"><br/><div class="descr_func"><div class="descr"><p> Takes as input a tropical variety X , a rational function f and optionally a positive integer k. It then computes
 f^k * X
Note that this still produces a meaningful result, if the WeightedComplex is not balanced: The "divisor" of a given function is computed by taking all codim-1-faces, at which the complex is balanced and computing weights there.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical variety</div></div></td></tr><tr><td><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a></td><td class="param">f</td><td><div class="descr"><div>A rational function, defined on a DOMAIN that contains X. It can be in non-homog. coordinates, even if X is in homog. coordinates, but then its AMBIENT_DIM has to be equal to X-&gt;AMBIENT_DIM. Otherwise you will have to homogenize the function first.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>Optional. A positive integer. Note that if k is given, this will actually ignore f-&gt;POWER. Is 1 by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The divisor f^k * X. It is given in homogeneous coordinates, if and only if: Either X is already in homogeneous coordinates or f is defined on a homogeneous DOMAIN (which is implicitly the case for ALL MinMaxFunctions)</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i55" onclick="swap_content('span:i55'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i55" onclick="swap_content('span:i55'); return false;"><b>divisorByPLF</b></a>
			(fan, A)  

			
				→ The<span style="display: inline" id="span:i55"><br/><div class="descr_func"><div class="descr"><p> NOTE: Deprecated. Use divisor(..) instead
 Computes the divisor of a MinMaxFunction on a given tropical variety. The result will be 
 in homogeneous coordinates, whether the tropical variety uses them or not. The function 
 should be given on the affine coordinates of the variety, NOT the homogeneous ones.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A tropical variety, on which the divisor is computed</div></div></td></tr><tr><td><a href="#i21" onclick="unfold('span:i21');">MinMaxFunction</a></td><td class="param">A</td><td><div class="descr"><div>function whose DOMAIN should be equal to the affine coordinate 
 space of the variety, i.e. AMBIENT_DIM-1, if the variety uses homogeneous coordinates, 
 AMBIENT_DIM otherwise.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>corresponding divisor as a tropical variety in homogeneous coordinates.</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i39" onclick="swap_content('span:i39'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i39" onclick="swap_content('span:i39'); return false;"><b>divisorByValue</b></a>
			(A, An)  

			
				→ The<span style="display: inline" id="span:i39"><br/><div class="descr_func"><div class="descr"><p> NOTE: Deprecated. Use divisor(...) or divisor_nr(...) in conjunction with function_value instead
 Takes as input a tropical fan / tropical variety and an array of rational values. The array length should coincide 
 with the number of <a href="#i115" onclick="unfold('span:i115');">CMPLX_RAYS</a> of the fan plus the dimension of the lineality space and will be interpreted as a rational
 function, where each value has been assigned to the rays given by $fan-&gt;CMPLX_RAYS and to the generators given by $fan-&gt;LINEALITY_SPACE  
 (in that order). Missing values will be filled up by 0's, superfluous ones will be ignored. 
 The function will then compute the corresponding Weil divisor and return it as a 
 tropical variety given as a fan. The fan uses homogeneous coordinates, if and only the input fan does.
 Note that this still produces a meaningful result, if the WeightedComplex is not balanced: The "divisor" of a given function is computed by taking all codim-1-faces, at which the complex is balanced and computing weights there.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">A</td><td><div class="descr"><div>tropical variety on which the divisor is computed. </div></div></td></tr><tr><td>values</td><td class="param">An</td><td><div class="descr"><div>array of rational values that define an integer affine map on the fan. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>divisor of the function defined by values on the given fan, as a tropical variety.</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i53" onclick="swap_content('span:i53'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i53" onclick="swap_content('span:i53'); return false;"><b>divisor_list</b></a>
			(X, A)  

			
				→ WeightedComplex<span style="display: inline" id="span:i53"><br/><div class="descr_func"><div class="descr"><p> Takes as input a tropical variety and an arbitrarily long list of RationalFunction objects f_1,...,f_k. It will 
# then compute f_k^a_k * ... * f_1^a_1 * X, where a_i = f_i-&gt;POWER</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">X</td><td/></tr><tr><td><a href="#i23" onclick="unfold('span:i23');">RationalFunction</a></td><td class="param">A</td><td><div class="descr"><div>list f_1,...,f_k</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The divisor f_k * ... * f_1 * X. It is in homogeneous coordinates, if and only if: 
 Either X is in homog. coordinates or one of the functions is defined on a homogeneous domain.</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i61" onclick="swap_content('span:i61'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i61" onclick="swap_content('span:i61'); return false;"><b>divisor_list_nr</b></a>
			()  

			<span style="display: inline" id="span:i61"><br/><div class="descr_func"><div class="descr"><p> Works exactly as divisor_list(WeightedComplex;RationalFunction,...). Should be called
 ONLY,when every occuring function f is defined on a DOMAIN equal to X (in the sense that
 all properties like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not
 sufficient). In this case this function will in general be faster.</p></div></div><div class="level3"/></span></div></div></li><li><div class="li"><div class="icon" id="icon:i43" onclick="swap_content('span:i43'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i43" onclick="swap_content('span:i43'); return false;"><b>divisor_nr</b></a>
			()  

			<span style="display: inline" id="span:i43"><br/><div class="descr_func"><div class="descr"><p> Works exactly as divisor(WeightedComplex, RationalFunction;Int). Should be called ONLY,
 when the function f is defined on a DOMAIN equal to X (in the sense that all properties
 like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not sufficient). In this
 case this function will in general be faster.</p></div></div><div class="level3"/></span></div></div></li><li><div class="li"><div class="icon" id="icon:i68" onclick="swap_content('span:i68'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i68" onclick="swap_content('span:i68'); return false;"><b>fan_decomposition</b></a>
			(complex)  

			
				→ perl::ListReturn<span style="display: inline" id="span:i68"><br/><div class="descr_func"><div class="descr"><p> Take a polyhedral complex and returns a list of all the local vertex fans, 
 i.e. for each affine ray r, the list contains the fan Star_complex(r) 
 (in non-homogeneous coordinates)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>A tropical variety</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>perl::ListReturn</td><td><div class="descr"><div>A list of WeightedComplex objects in 
 non-homogeneous coordinates. The i-th complex corresponds to the i-th 
 affine ray ( vertex). If the complex is not in homogeneous coordinates, 
 the list contains just the complex itself </div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i35" onclick="swap_content('span:i35'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i35" onclick="swap_content('span:i35'); return false;"><b>halfspace_complex</b></a>
			(constant, equation)  

			
				→ WeightedComplex<span style="display: inline" id="span:i35"><br/><div class="descr_func"><div class="descr"><p> Creates the halfspace complex defined by an rational vector g and a rational b, i.e. the 
 complex consisting of the two maximal cones g &gt;= a and g &lt;= a</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">constant</td><td><div class="descr"><div>The constant translation a</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">equation</td><td><div class="descr"><div>The defining equation g</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The resultin halfspace complex</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i60" onclick="swap_content('span:i60'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i60" onclick="swap_content('span:i60'); return false;"><b>intersect_complete_fan</b></a>
			(fan, completeFan)  

			
				→ WeightedComplex<span style="display: inline" id="span:i60"><br/><div class="descr_func"><div class="descr"><p> Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete 
 (i.e. its support is the whole ambient space) to compute the intersection correctly.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>An arbitrary weighted polyhedral fan</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">fan::PolyhedralFan</a></td><td class="param">completeFan</td><td><div class="descr"><div>A complete polyhedral fan</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The intersection of both fans (whose support is equal to the support of fan). The 
 resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, 
 the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned.</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i56" onclick="swap_content('span:i56'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i56" onclick="swap_content('span:i56'); return false;"><b>is_zerocycle</b></a>
			(c)  

			
				→ Bool<span style="display: inline" id="span:i56"><br/><div class="descr_func"><div class="descr"><p> Checks whether a cycle is the 0 cycle</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">c</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i63" onclick="swap_content('span:i63'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i63" onclick="swap_content('span:i63'); return false;"><b>latticeBasis</b></a>
			(cone)  

			
				→ Matrix&lt;Integer&gt;<span style="display: inline" id="span:i63"><br/><div class="descr_func"><div class="descr"><p> Takes a cone and computes a Z-basis of the vector space spanned by the cone,
 returned as row vectors of a matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">cone</td><td><div class="descr"><div>A cone for which a Z-basis is to be computed</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i66" onclick="swap_content('span:i66'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i66" onclick="swap_content('span:i66'); return false;"><b>latticeNormal</b></a>
			(taumatrix, sigmamatrix, additionalRay)  

			
				→ Vector<span style="display: inline" id="span:i66"><br/><div class="descr_func"><div class="descr"><p> Takes two matrices whose rows define the dual of the linear span of cone tau and sigma.
 Assuming that tau is a codimension one face of sigma, computes a representative of
 the primitive lattice normal vector of sigma with respect to tau</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">taumatrix</td><td><div class="descr"><div>a codimension one face of sigma, given as a matrix defining its linear span</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">sigmamatrix</td><td><div class="descr"><div>an arbitrary cone, given as a matrix defining its linear span</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td class="param">additionalRay</td><td><div class="descr"><div>A ray that is contained in sigma, but not in tau. Used to calculate proper orientation of the normal vector.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i37" onclick="swap_content('span:i37'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i37" onclick="swap_content('span:i37'); return false;"><b>latticeNormalByCone</b></a>
			(tau, sigma)  

			
				→ Vector<span style="display: inline" id="span:i37"><br/><div class="descr_func"><div class="descr"><p> Assuming that tau is a codimension one face of sigma, computes a representative of
 the primitive lattice normal vector of sigma with respect to tau</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">tau</td><td><div class="descr"><div>a codimension one face of tau, given as a cone</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">sigma</td><td><div class="descr"><div>an arbitrary cone</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i59" onclick="swap_content('span:i59'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i59" onclick="swap_content('span:i59'); return false;"><b>linear_nspace</b></a>
			(n, weight)  

			
				→ WeightedComplex<span style="display: inline" id="span:i59"><br/><div class="descr_func"><div class="descr"><p> Creates the tropical variety R^n as a single lineality space with a given weight </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The dimension of the linear space</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i69" onclick="swap_content('span:i69'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i69" onclick="swap_content('span:i69'); return false;"><b>linear_space_by_basis</b></a>
			(rowspace, weight)  

			
				→ WeightedComplex<span style="display: inline" id="span:i69"><br/><div class="descr_func"><div class="descr"><p> Creates a tropical variety consisting of a single linear space defined by a matrix and assigns it a given weight.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rowspace</td><td><div class="descr"><div>A matrix whose rows span the lineality space of the fan (hence the whole fan) and
 which are assumed to be linearly independent. If you are not sure of this, call linear_space_by_matrix instead.</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i45" onclick="swap_content('span:i45'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i45" onclick="swap_content('span:i45'); return false;"><b>linear_space_by_matrix</b></a>
			(rowspace, weight)  

			
				→ WeightedComplex<span style="display: inline" id="span:i45"><br/><div class="descr_func"><div class="descr"><p> Creates a tropical variety consisting of a single linear space defined by a matrix and assigns it a given weight.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rowspace</td><td><div class="descr"><div>A matrix whose rows span the lineality space of the fan (hence the whole fan).
 If you already know that your rows are a basis, call linear_space_by_basis instead.</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i52" onclick="swap_content('span:i52'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i52" onclick="swap_content('span:i52'); return false;"><b>moduliDimensionFromLength</b></a>
			(k)  

			
				→ Int<span style="display: inline" id="span:i52"><br/><div class="descr_func"><div class="descr"><p> Takes a positive length (of a vector) and assumes it is of the form (n over 2)
 It then computes n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The length = (n over 2)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td><div class="descr"><div>n. If k is not of the form (n over 2), the result is arbitrary</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i51" onclick="swap_content('span:i51'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i51" onclick="swap_content('span:i51'); return false;"><b>random_minmaxdivisor</b></a>
			(ambient_dim, equations, upperCoeffBound, onlyFunction)  

			
				→ MinMaxFunction<span style="display: inline" id="span:i51"><br/><div class="descr_func"><div class="descr"><p> Creates a divisor of a random MinMaxFunction (in n-space) with integer coefficients.
 By default this uses maximum for the function. It returns an array containing the function created and its
 divisor in n-space.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">ambient_dim</td><td><div class="descr"><div>The ambient dimension of the divisor</div></div></td></tr><tr><td>int</td><td class="param">equations</td><td><div class="descr"><div>The number of equations to use</div></div></td></tr><tr><td>int</td><td class="param">upperCoeffBound</td><td><div class="descr"><div>The upper bound  of the coefficients (All coefficients are positive)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">onlyFunction</td><td><div class="descr"><div>If true, only the random function is returned and the divisor is not computed. FALSE by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i21" onclick="unfold('span:i21');">MinMaxFunction</a></td><td><div class="descr"><div>if onlyFunction is TRUE, otherwise an array containing first the function, then the divisor</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i47" onclick="swap_content('span:i47'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i47" onclick="swap_content('span:i47'); return false;"><b>rational_curve_from_moduli</b></a>
			(v)  

			
				→ RationalCurve<span style="display: inline" id="span:i47"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^((n over 2) - n) that lies in M_0,n (in its matroid coordinates 
 and computes the corresponding rational curve.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector in the moduli space</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td/></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i41" onclick="swap_content('span:i41'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i41" onclick="swap_content('span:i41'); return false;"><b>rational_curve_list_from_metric</b></a>
			(m)  

			
				→ RationalCurve<span style="display: inline" id="span:i41"><br/><div class="descr_func"><div class="descr"><p> Takes a matrix whose rows are metrics of rational n-marked curves.
 Returns a list, where the i-th element is the curve corr. to 
 the i-th row in the matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td><div class="descr"><div>: An array of RationalCurves</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i38" onclick="swap_content('span:i38'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i38" onclick="swap_content('span:i38'); return false;"><b>rational_curve_list_from_moduli</b></a>
			(m)  

			
				→ RationalCurve<span style="display: inline" id="span:i38"><br/><div class="descr_func"><div class="descr"><p> Takes a matrix whose rows are elements in the moduli space M_0,n in matroid 
 coordinates. Returns a list, where the i-th element is the curve corr. to 
 the i-th row in the matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td><div class="descr"><div>: An array of RationalCurves</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i48" onclick="swap_content('span:i48'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i48" onclick="swap_content('span:i48'); return false;"><b>recession_fan</b></a>
			(complex)  

			
				→ WeightedComplex<span style="display: inline" id="span:i48"><br/><div class="descr_func"><div class="descr"><p> Computes the recession fan of a tropical variety</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>A tropical variety. If it is a fan, the complex itself is returned</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>A tropical fan, the recession fan of the complex</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i54" onclick="swap_content('span:i54'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i54" onclick="swap_content('span:i54'); return false;"><b>sum_curves</b></a>
			(An, v)  

			
				→ RationalCurve<span style="display: inline" id="span:i54"><br/><div class="descr_func"><div class="descr"><p> This function takes a vector of coefficients a_i and a list of RationalCurves c_i and computes
 sum(a_i * c_i). In particular, it also checks, whether the result lies in M_0,n. If not, it returns undef</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td class="param">An</td><td><div class="descr"><div>arbitrary list of RationalCurve objects</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A list of coefficients. Superfluous coefficients are ignored, missing ones replaced
 by +1(!)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i22" onclick="unfold('span:i22');">RationalCurve</a></td><td><div class="descr"><div>The linear combination of the curves defined by the coefficients or undef, if the result
 is not in M_0,n. The history of the operation is kept in INPUT_SETS and INPUT_COEFFS</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i42" onclick="swap_content('span:i42'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i42" onclick="swap_content('span:i42'); return false;"><b>testFourPointCondition</b></a>
			(v)  

			
				→ Int<span style="display: inline" id="span:i42"><br/><div class="descr_func"><div class="descr"><p> Takes a metric vector in Q^{(n over 2)} and checks whether it fulfills 
 the four-point condition, i.e. whether it lies in M_0,n. More precisely 
 it only needs to be equivalent to such a vector</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>The vector to be checked</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td><div class="descr"><div>A quadruple (array) of indices, where the four-point condition 
 is violated or an empty list, if the vector is indeed in M_0,n</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i34" onclick="swap_content('span:i34'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i34" onclick="swap_content('span:i34'); return false;"><b>tropical_cube</b></a>
			(n, k)  

			
				→ WeightedComplex<span style="display: inline" id="span:i34"><br/><div class="descr_func"><div class="descr"><p> Creates the tropical cube T^n_k, i.e. the tropical variety obtained by
 glueing together the 2^n L^n_k obtained by applying all possible sign 
 changes, such that the vertices form the k-skeleton of 
 the n-dimensional cube</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The ambient dimension</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the cube</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>cube</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i46" onclick="swap_content('span:i46'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i46" onclick="swap_content('span:i46'); return false;"><b>tropical_lnk</b></a>
			(n, k)  

			
				→ WeightedComplex<span style="display: inline" id="span:i46"><br/><div class="descr_func"><div class="descr"><p> Creates the linear tropical space L^n_k. This tropical fan is defined in the following way: 
 As rays we take -e_i,i=1,...,n, where e_i is the i-th standard basis vector of R^n and 
 e_0 = e_1 + ... + e_n. As maximal cones we take the cones generated by rays {e_i, i in S}, where
 S runs over all k-element subsets of {0,..,n}.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The ambient dimension of the fan.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the fan (should be smaller equal n, otherwise an error is thrown).</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>A PolyhedralFan object representing L^n_k</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i67" onclick="swap_content('span:i67'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i67" onclick="swap_content('span:i67'); return false;"><b>tropical_m0n</b></a>
			(n)  

			
				→ WeightedComplex<span style="display: inline" id="span:i67"><br/><div class="descr_func"><div class="descr"><p> Creates the moduli space of abstract rational n-marked curves. Its coordinates are
 given as the coordinates of the bergman fan of the matroid of the complete graph on 
 n-1 nodes (but not computed as such)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be at least 4</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The tropical moduli space M_0,n</div></div></td></tr></table></div></span></div></div></li><li><div class="li"><div class="icon" id="icon:i62" onclick="swap_content('span:i62'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i62" onclick="swap_content('span:i62'); return false;"><b>zero_cycle</b></a>
			()  

			<span style="display: inline" id="span:i62"><br/><div class="descr_func"><div class="descr"><p> Returns the zero cycle, i.e. a weighted complex with no rays and cones</p></div></div><div class="level3"/></span></div></div></li></span></div></div></li><li><div class="li"><div class="icon" id="icon:i9" onclick="swap_content('span:i9'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i9" onclick="swap_content('span:i9'); return false;"><h3><div class="category">Tropical geometry /Local geometry</div></h3></a><span style="display: inline" id="span:i9"><li><div class="li"><div class="icon" id="icon:i31" onclick="swap_content('span:i31'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i31" onclick="swap_content('span:i31'); return false;"><b>local_restrict</b></a>
			(complex, cones)  

			
				→ WeightedComplex<span style="display: inline" id="span:i31"><br/><div class="descr_func"><div class="descr"><p> This takes a tropical variety and an IncidenceMatrix describing a set
 of cones (not necessarily maximal ones) of this variety. It will then
 create a variety that contains all compatible maximal cones and is
 locally restricted to the given cone set.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>An arbitrary weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">IncidenceMatrix</a></td><td class="param">cones</td><td><div class="descr"><div>A set of cones, indices refer to RAYS</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i20" onclick="unfold('span:i20');">WeightedComplex</a></td><td><div class="descr"><div>The same complex, locally restricted to the given
 cones</div></div></td></tr></table></div></span></div></div></li></span></div></div></li></ul></div><div class="level3"><h2>Common Option Lists</h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i17" onclick="swap_content('span:i17'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i17" onclick="swap_content('span:i17'); return false;"><b>Visual::RationalCurve::decorations</b></a><span style="display: inline" id="span:i17"><div class="undoc"><p>UNDOCUMENTED</p></div><b>imports from: </b>Visual::Graph::decorations, Visual::Wire::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">LengthLabels</td><td><div class="descr"><div>if set to "hidden", the labels indicating the lengths are hidden</div></div></td></tr></table></span></div></div></li><li><div class="li"><div class="icon" id="icon:i18" onclick="swap_content('span:i18'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i18" onclick="swap_content('span:i18'); return false;"><b>Visual::WeightedComplex::BoundingDecorations</b></a><span style="display: inline" id="span:i18"><div class="undoc"><p>UNDOCUMENTED</p></div><b>imports from: </b>Visual::WeightedComplex::decorations, Visual::Polygons::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Flexible&lt;Rational&gt;</a></td><td class="param">BoundingDistance</td><td><div class="descr"><div>The distance of the border of the bounding box from the smallest box 
 containing the affine points of the complex. This is only relevant, if BoundingMode is "relative"</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">BoundingBox</td><td><div class="descr"><div>A fixed bounding box, determined by two row vectors that specify two of its vertices (on "on top" and one "at the bottom"). Is only relevant, if BoundingMode is "absolute"</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">BoundingMode</td><td><div class="descr"><div>If set to "relative", the function determines the smallest possible box containing all affine points of the complex and then enlarges the box by BoundingDistance to all sides. If set to "absolute", BoundingBox must be specified and the complex will be intersected with that box. By default this is set to "relative".</div></div></td></tr></table></span></div></div></li><li><div class="li"><div class="icon" id="icon:i19" onclick="swap_content('span:i19'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" name="i19" onclick="swap_content('span:i19'); return false;"><b>Visual::WeightedComplex::decorations</b></a><span style="display: inline" id="span:i19"><div class="descr"><p> require Visual;
 require Visual::PointSet;
 require Visual::Wire;
 require Visual::Polygon;</p></div><b>imports from: </b>Visual::Polygons::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Flexible&lt;Rational&gt;</a></td><td class="param">DirScale</td><td><div class="descr"><div>length of the directional rays</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">WeightLabels</td><td><div class="descr"><div>if set to "hidden", the labels indicating the weights are hidden</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">CoordLabels</td><td><div class="descr"><div>If set to "show", the labels indicating the vertex coordinates are displayed, otherwise they are not. Note that this is expensive and significantly increases computation time.
 String If set to "pinned", the weight labels are computed in such a way that they stick to the cells when the complex
 is exploded. Note that this is expensive and significantly increases computation time.
 WeightLabelStyle =&gt; "",</div></div></td></tr></table></span></div></div></li></ul></div></div>
