<div xmlns:html="http://www.w3.org/1999/xhtml" xmlns:pm="http://www.polymake.de/ns/docs#3" id="content"><h1>application: atint</h1><div class="descr"><p> a-tint is an application for polymake (see <a href="http://polymake.mathematik.tu-darmstadt.de" onclick="unfold('span:');">http://polymake.mathematik.tu-darmstadt.de</a>). 
 It aims to implement algorithmic aspects of tropical intersection theory. 
 Tropical varieties are realized as weighted polyhedral complexes. This project is part of and
 financed by the DFG Priority Project SPP 1489 (see also <a href="http://www.computeralgebra.de/index.php" onclick="unfold('span:');">http://www.computeralgebra.de/index.php</a>)
 For more information see <a href="https://bitbucket.org/hampe/atint/" onclick="unfold('span:');">https://bitbucket.org/hampe/atint/</a></p></div><br/><b>imports from: </b><a href="common.html">common</a>, <a href="graph.html">graph</a><br/><b>uses: </b><a href="fan.html">fan</a>, <a href="group.html">group</a>, <a href="matroid.html">matroid</a>, <a href="polytope.html">polytope</a>, <a href="topaz.html">topaz</a><br/><div class="level2"><h2><b>Objects</b></h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i27" onclick="swap_content('span:i27'); return false;">
		 	
	</div><div class="unfoldable"><div class="object"><h3><a href="#" class="javalink" id="i27" onclick="swap_content('span:i27'); return false;"><span class="object">MinMaxFunction</span></a></h3></div><div style="display: inline" id="span:i27" class="foldit"><div class="descr"><p> A special case of a rational function, a tropical polynomial (using either min or max)</p></div><b>derived from: </b><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a><br/><div class="level3"/><h4>Properties of MinMaxFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i158" onclick="swap_content('span:i158'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i158" onclick="swap_content('span:i158'); return false;"><span class="category">Basic properties</span></a></h3></div><div style="display: inline" id="span:i158" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i229" onclick="swap_content('span:i229'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i229" onclick="swap_content('span:i229'); return false;"><b>CONSTANT_COEFFICIENTS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i229" class="foldit"><div class="descr"><p> A vector of the constant coefficients of the function. The i-th entry is the constant term of the function 
 represented by the i-th row of <a href="#i228" onclick="unfold('span:i228');">LINEAR_COEFFICIENTS</a>.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i223" onclick="swap_content('span:i223'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i223" onclick="swap_content('span:i223'); return false;"><b>DOMAIN_DIMENSION</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i223" class="foldit"><div class="descr"><p> The dimension n of the domain R<sup>n</sup> of the function. Is equal to the columns of <a href="#i228" onclick="unfold('span:i228');">LINEAR_COEFFICIENTS</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i227" onclick="swap_content('span:i227'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i227" onclick="swap_content('span:i227'); return false;"><b>FUNCTION_MATRIX</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i227" class="foldit"><div class="descr"><p> This is simply the concatenation of <a href="#i228" onclick="unfold('span:i228');">LINEAR_COEFFICIENTS</a> and <a href="#i229" onclick="unfold('span:i229');">CONSTANT_COEFFICIENTS</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i225" onclick="swap_content('span:i225'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i225" onclick="swap_content('span:i225'); return false;"><b>INPUT_STRING</b>:
		</a>common::String<br/><div style="display: inline" id="span:i225" class="foldit"><div class="descr"><p> The syntax for STRING in pseudo-regex is the following (after deleting any whitespace): 
 "(min|max)(\(|\{)(SUM)+(,SUM)*(\)|\})"
 where SUM = "(TERM)+([\+\-](TERM))*"
 where TERM = (RAT)|(RAT[:alpha:][^\+,\)]*)
 where RAT is any string that can be parsed to a polymake Rational and does not contain any letters.
 To make this more comprehensible: STRING takes anything of the form min(...) or max(...) (or min{...}, max{...}), 
 where ... is a comma-separated list of strings of the form a + bx + cy + ... , where a,b,c are rational 
 numbers and x,y are variables. Any string that starts with a letter and does not contain "+" or "," can be 
 a variable. Such a sum can contain several such terms for the same variable and they need not be in any order.
 A term can be of the form "3x", "3*x", but "x3" will be interpreted as 1 * "x3". Coefficients should not contain letters,
 i.e. 3E10*x would be interpreted as 3 * "E10*x". Empty terms are allowed.
 If no INPUT_VARS are given, MinMaxFunction will recognize the dimension of the domain from the number of variables and will
 assign the variables to the standard basis vectors by sorting them alphabetically. 
 E.g. INPUT_STRING=&gt;"min(x+2z,3+z)",INPUT_VARS=&gt;["x","y","z"] is a different function than INPUT_STRING=&gt;"min(x+2z,3+z)",
 since the second is the function f:R² -&gt; R, (a,b) -&gt; min(a+2b,3+b) and the first is f: R³-&gt;R, (a,b,c)-&gt;min(a+2c,3+c)
 If however, <a href="#i226" onclick="unfold('span:i226');">INPUT_VARS</a> is given, the parser will *ignore* any terms containing variables that are not 
 defined in <a href="#i226" onclick="unfold('span:i226');">INPUT_VARS</a>. The variables will be assigned to the standard basis vectors in order of their
 appearance in <a href="#i226" onclick="unfold('span:i226');">INPUT_VARS</a>.
 E.g. INPUT_STRING=&gt;"min(x+2c,3+x)",INPUT_VARS=&gt;["x","y"] will give the function f:R²-&gt;R, (a,b) -&gt; min(a,3+a)
 Note that ANY whitespace in <a href="#i225" onclick="unfold('span:i225');">INPUT_STRING</a> is ignored, as well as any brackets '(',')' besides the ones enclosing
 the function list.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i226" onclick="swap_content('span:i226'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i226" onclick="swap_content('span:i226'); return false;"><b>INPUT_VARS</b>:
		</a>common::Array&lt;String&gt;<br/><div style="display: inline" id="span:i226" class="foldit"><div class="descr"><p> A list of strings that are variables. Any string that starts with a letter (i.e. [a-zA-Z]) and does not contain
 '+' or ',' can be a variable. Illegal variable names will not cause an error, but the parser will not detect them
 in <a href="#i225" onclick="unfold('span:i225');">INPUT_STRING</a>, so their coefficients will be 0. If you assign a value to this, only terms with variables in 
 this array will be parsed. Note that ANY whitespace in a variable name will be ignored. The i-th element of this array
 will be assigned to the i-th standard basis vector.
 You might want to use <a href="#i45" onclick="unfold('span:i45');">input_vars</a> for defining this property.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i228" onclick="swap_content('span:i228'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i228" onclick="swap_content('span:i228'); return false;"><b>LINEAR_COEFFICIENTS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i228" class="foldit"><div class="descr"><p> A matrix of the linear coefficients of the functions. Each row corresponds to a linear function within
 the min or max and the entry in column j is the coefficient of variable x_j</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i222" onclick="swap_content('span:i222'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i222" onclick="swap_content('span:i222'); return false;"><b>NORMAL_FAN</b>:
		</a><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a><br/><div style="display: inline" id="span:i222" class="foldit"><div class="descr"><p> This computes the normal fan of the polytope conv( (c_i,a_i) ) in Q^{n+1}, where the function is 
 max/min(a_i*x + c_i,i=1,..,k) with a_i in Q^n and c_i in Q. (For min it actually computes the negative
 of the normal fan)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i230" onclick="swap_content('span:i230'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i230" onclick="swap_content('span:i230'); return false;"><b>STANDARD_FORM</b>:
		</a>common::String<br/><div style="display: inline" id="span:i230" class="foldit"><div class="descr"><p> This will give a nice string representation of the function, where all variables are of the form 'x_i',
 i = 1,..,<a href="#i223" onclick="unfold('span:i223');">DOMAIN_DIMENSION</a>. There are no redundancies and terms are sorted in ascending alphabetic order of their variable,
 with the constant term at the end</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i224" onclick="swap_content('span:i224'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i224" onclick="swap_content('span:i224'); return false;"><b>USES_MIN</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i224" class="foldit"><div class="descr"><p> Indicates whether this function uses min (TRUE) or max (FALSE). The default (if you directly give a coefficient matrix,
 but not this value) is FALSE</p></div></div></div></div></li></ul></div></div></div></li></ul><h4>User Methods of MinMaxFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i157" onclick="swap_content('span:i157'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i157" onclick="swap_content('span:i157'); return false;"><span class="category">Rational functions</span></a></h3></div><div style="display: inline" id="span:i157" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i221" onclick="swap_content('span:i221'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i221" onclick="swap_content('span:i221'); return false;"><b>linearityDomains</b></a>
			()  

			
				→ WeightedComplex<div style="display: inline" id="span:i221" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This computes the (homogenized version of the) domains of linearity of the function. More precisely, it computes the 
 normal fan of the polytope defined by this function in the following way:
 If the standard form reads max(a_i*x + c_i,i=1,..,k) with a_i in Q^n and c_i in Q, then the polytope is
 conv( (c_i,a_i) ) in Q^{n+1} (in homog. coordinates its actually in n+2 coordinates).
 In fact, this simply returns the property <a href="#i250" onclick="unfold('span:i250');">DOMAIN</a>. It is kept for backwards compatibility.</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i219" onclick="swap_content('span:i219'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i219" onclick="swap_content('span:i219'); return false;"><b>pow</b></a>
			(k)  

			
				→ MinMaxFunction<div style="display: inline" id="span:i219" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Does the same as pow for RationalFunction, i.e. returns the same function but with a given internal power</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The new power of this function, i.e. using it in a divisorial expression will be interpretes as 
 applying it k times</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i220" onclick="swap_content('span:i220'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i220" onclick="swap_content('span:i220'); return false;"><b>valueAt</b></a>
			(point)  

			
				→ Rational<div style="display: inline" id="span:i220" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the value of this function at a given point. If the point has dimension DOMAIN_DIMENSION, it simply computes the value. 
 If it has dimension DOMAIN_DIMENSION +1, it assumes the point is given in homogeneous coordinates (i.e. the first coordinate should be one or zero, 
 this is not checked, however) and computes the function value of the last DOMAIN_DIMENSION coordinates.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">point</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td><div class="descr"><div>The function value at the given point.</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i156" onclick="swap_content('span:i156'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i156" onclick="swap_content('span:i156'); return false;"><span class="category">Reduction</span></a></h3></div><div style="display: inline" id="span:i156" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i218" onclick="swap_content('span:i218'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i218" onclick="swap_content('span:i218'); return false;"><b>reduce</b></a>
			()  

			
				→ MinMaxFunction<div style="display: inline" id="span:i218" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a MinMaxFunction and reduces its terms to the necessary ones. That means it computes the newton polytope
 and returns a function that has only terms corresponding to the vertices.</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td><div class="descr"><div>The reduced version of the function</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i25" onclick="swap_content('span:i25'); return false;">
		 	
	</div><div class="unfoldable"><div class="object"><h3><a href="#" class="javalink" id="i25" onclick="swap_content('span:i25'); return false;"><span class="object">Morphism</span></a></h3></div><div style="display: inline" id="span:i25" class="foldit"><div class="descr"><p> This represents a morphism as a "rational function to R<sup>n</sup>", i.e. it is given via its <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a>
 and <a href="#i174" onclick="unfold('span:i174');">RAY_VALUES</a> and <a href="#i177" onclick="unfold('span:i177');">LIN_VALUES</a> now lying in some R<sup>n</sup>  instead of R. The user has to take
 care that the function is locally linear, this is not checked in any way.</p></div><div class="level3"/><h4>Properties of Morphism</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i145" onclick="swap_content('span:i145'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i145" onclick="swap_content('span:i145'); return false;"><span class="category">Basic properties</span></a></h3></div><div style="display: inline" id="span:i145" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i176" onclick="swap_content('span:i176'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i176" onclick="swap_content('span:i176'); return false;"><b>DOMAIN</b>:
		</a><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a><br/><div style="display: inline" id="span:i176" class="foldit"><div class="descr"><p> This property describes the domain of the morphism. I.e. the morphism is
 defined on this complex and is locally affine integral linear. It is of type WeightedComplex, but
 actually only the basic geometric properties <a href="#" onclick="unfold('span:');">RAYS</a>/<a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>,<a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>/<a href="#i201" onclick="unfold('span:i201');">CMPLX_MAXIMAL_CONES</a>,<a href="#" onclick="unfold('span:');">LINEALITY_SPACE</a> and <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> are of any interest. Further properties are ignored and can be omitted during generation.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i173" onclick="swap_content('span:i173'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i173" onclick="swap_content('span:i173'); return false;"><b>IS_GLOBAL</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i173" class="foldit"><div class="descr"><p> Is TRUE, iff the morphism is a global affine linear integral map
 In this case, the morphism's description can be read off from the properties
 <a href="#i175" onclick="unfold('span:i175');">MATRIX</a> and <a href="#i172" onclick="unfold('span:i172');">TRANSLATE</a>. Setting any of these properties during initialization 
 will make this TRUE by default.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i177" onclick="swap_content('span:i177'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i177" onclick="swap_content('span:i177'); return false;"><b>LIN_VALUES</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i177" class="foldit"><div class="descr"><p> The vector in row i describes the function value (slope) of <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a>-&gt;LINEALITY_SPACE-&gt;row(i)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i175" onclick="swap_content('span:i175'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i175" onclick="swap_content('span:i175'); return false;"><b>MATRIX</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i175" class="foldit"><div class="descr"><p> If the morphism is a global affine linear map x |-&gt; Ax+v, i.e. IS_GLOBAL is or should be true, 
 then this contains the matrix A (in non-homog. coordinates). If <a href="#i172" onclick="unfold('span:i172');">TRANSLATE</a> is set, but this property is not set, then it
 is the identity by default.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i174" onclick="swap_content('span:i174'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i174" onclick="swap_content('span:i174'); return false;"><b>RAY_VALUES</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i174" class="foldit"><div class="descr"><p> The vector at row i describes the function value of ray <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a>-&gt;CMPLX_RAYS-&gt;row(i). More precisely, if the
 corresponding ray is a vertex, it describes its function value. If it is a directional ray, it describes 
 the slope on that ray.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i172" onclick="swap_content('span:i172'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i172" onclick="swap_content('span:i172'); return false;"><b>TRANSLATE</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i172" class="foldit"><div class="descr"><p> If the morphism is a global affine linear map x |-&gt; Ax+v, i.e. IS_GLOBAL is or should be true, 
 then this contains the translation vector v (in non-homog. coordinates). If <a href="#i175" onclick="unfold('span:i175');">MATRIX</a> is set, but this property is not set,
 then it is the zero vector by default.</p></div></div></div></div></li></ul></div></div></div></li></ul><h4>User Methods of Morphism</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i144" onclick="swap_content('span:i144'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i144" onclick="swap_content('span:i144'); return false;"><span class="category">Morphisms</span></a></h3></div><div style="display: inline" id="span:i144" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i168" onclick="swap_content('span:i168'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i168" onclick="swap_content('span:i168'); return false;"><b>after</b></a>
			(g)  

			
				→ Morphism<div style="display: inline" id="span:i168" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the composition of this morphism with another morphism. The other morphism should either be global 
 and surjective or its image should be contained in the <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a> of this morphism</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td class="param">g</td><td><div class="descr"><div>Another morphism. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>The composition "this after g" (The converse composition is computed using "before")</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i171" onclick="swap_content('span:i171'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i171" onclick="swap_content('span:i171'); return false;"><b>before</b></a>
			(g)  

			
				→ Morphism<div style="display: inline" id="span:i171" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the composition of another morphism with this morphism. This morphism should either be global 
 and surjective or its image should be contained in the <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a> of the other morphism</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td class="param">g</td><td><div class="descr"><div>Another morphism. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>The composition "g after this" (The converse composition is computed using "after")</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i169" onclick="swap_content('span:i169'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i169" onclick="swap_content('span:i169'); return false;"><b>homogenize</b></a>
			()  

			<div style="display: inline" id="span:i169" class="foldit"><br/><div class="descr_func"><div class="descr"><p>   If the domain of the function is in homogeneous coordinates, it returns the function. Otherwise it homogenizes
   the domain and assigns the value (0,..,0) to the additional vertex at the origin</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i170" onclick="swap_content('span:i170'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i170" onclick="swap_content('span:i170'); return false;"><b>restrict</b></a>
			(f, domain)  

			
				→ Morphism<div style="display: inline" id="span:i170" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Restricts a morphism / rational function to a given domain (which should be contained in the <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a> of the function)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td class="param">f</td><td><div class="descr"><div>(or RationalFunction) The function to be restricted</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">domain</td><td><div class="descr"><div>The domain to which the function should be restricted</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>or RationalFunction The restricted function</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i143" onclick="swap_content('span:i143'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i143" onclick="swap_content('span:i143'); return false;"><span class="category">Visualization</span></a></h3></div><div style="display: inline" id="span:i143" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i167" onclick="swap_content('span:i167'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i167" onclick="swap_content('span:i167'); return false;"><b>BB_VISUAL</b></a>
			()  

			<div style="display: inline" id="span:i167" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Visualizes the domain of the morphism. Works exactly as BB_VISUAL of WeightedComplex, but has additional option
 @options</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">FunctionLabels</td><td><div class="descr"><div>If set to "show", textual function representations are diplayed on cones. False by default</div></div></td></tr><tr><td>option list:</td><td>Visual::WeightedComplex::FunctionDecorations</td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i28" onclick="swap_content('span:i28'); return false;">
		 	
	</div><div class="unfoldable"><div class="object"><h3><a href="#" class="javalink" id="i28" onclick="swap_content('span:i28'); return false;"><span class="object">RationalCurve</span></a></h3></div><div style="display: inline" id="span:i28" class="foldit"><div class="descr"><p> An n-marked rational curve, identified by its <a href="#i236" onclick="unfold('span:i236');">SETS</a>, i.e. its partitions of {1,...,n}
 and its <a href="#i235" onclick="unfold('span:i235');">COEFFS</a>, i.e. the length of the corresponding edge</p></div><div class="level3"/><h4>Properties of RationalCurve</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i161" onclick="swap_content('span:i161'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i161" onclick="swap_content('span:i161'); return false;"><span class="category">Basic properties</span></a></h3></div><div style="display: inline" id="span:i161" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i235" onclick="swap_content('span:i235'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i235" onclick="swap_content('span:i235'); return false;"><b>COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i235" class="foldit"><div class="descr"><p> A list of positive rational coefficients. The list should have the same length as <a href="#i236" onclick="unfold('span:i236');">SETS</a> and
 contain only entries &gt; 0. The i-th entry then gives the length of the bounded edge defined by
 the i-th partition. If you're not sure if all your coefficients are &gt; 0, use <a href="#i243" onclick="unfold('span:i243');">INPUT_SETS</a> and
 <a href="#i244" onclick="unfold('span:i244');">INPUT_COEFFS</a> instead.
 Note that the zero curve (i.e. no bounded edges, only leaves) is represented by one empty set
 with corresponding lenghth 0.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i234" onclick="swap_content('span:i234'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i234" onclick="swap_content('span:i234'); return false;"><b>N_LEAVES</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i234" class="foldit"><div class="descr"><p> The number of leaves of the rational curve.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i236" onclick="swap_content('span:i236'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i236" onclick="swap_content('span:i236'); return false;"><b>SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i236" class="foldit"><div class="descr"><p> A list of partitions of [n] that define the tree of the curve: For each bounded edge we have 
 the corresponding partition of the n leaves. These should be irredundant. If you want to input
 a possibly redundant list, use <a href="#i243" onclick="unfold('span:i243');">INPUT_SETS</a> and <a href="#i244" onclick="unfold('span:i244');">INPUT_COEFFS</a> instead. The number of marked leaves
 should always be given by <a href="#i234" onclick="unfold('span:i234');">N_LEAVES</a>. The sets are subsets of {1,...,n} (NOT {0,..,n-1}!)
 Note that the zero curve (i.e. no bounded edges, only leaves) is represented by one empty set
 with corresponding lenghth 0.</p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i162" onclick="swap_content('span:i162'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i162" onclick="swap_content('span:i162'); return false;"><span class="category">Graph properties</span></a></h3></div><div style="display: inline" id="span:i162" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i239" onclick="swap_content('span:i239'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i239" onclick="swap_content('span:i239'); return false;"><b>GRAPH</b>:
		</a>graph::Graph&lt;Undirected&gt;<br/><div style="display: inline" id="span:i239" class="foldit"><div class="descr"><p> Contains the abstract graph (non-metric) corresponding to the curve. All unbounded leaves are modelled
 as bounded edges.
 The vertices at the ends of the "leaves" are always the first <a href="#i234" onclick="unfold('span:i234');">N_LEAVES</a> vertices.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i237" onclick="swap_content('span:i237'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i237" onclick="swap_content('span:i237'); return false;"><b>GRAPH_EDGE_LENGTHS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i237" class="foldit"><div class="descr"><p> Contains the lengths of the edges of <a href="#i239" onclick="unfold('span:i239');">GRAPH</a> that represent bounded edges of the curve.
 The coefficients appear in the order that the corr. edges appear in <a href="#" onclick="unfold('span:');">EDGES</a>.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i241" onclick="swap_content('span:i241'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i241" onclick="swap_content('span:i241'); return false;"><b>NODES_BY_LEAVES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i241" class="foldit"><div class="descr"><p> This incidence matrix gives a list of the vertices of the curve
 Each row corresponds to a vertex and contains as a set the [[LEAVES] that are
 attached to that vertex (again, counting from 1!)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i238" onclick="swap_content('span:i238'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i238" onclick="swap_content('span:i238'); return false;"><b>NODES_BY_SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i238" class="foldit"><div class="descr"><p> This incidence matrix gives a list of the vertices of the curve
 Each row corresponds to a vertex and contains as a set the row indices of the <a href="#i236" onclick="unfold('span:i236');">SETS</a> that
 correspond to edges attached to that vertex</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i240" onclick="swap_content('span:i240'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i240" onclick="swap_content('span:i240'); return false;"><b>NODE_DEGREES</b>:
		</a>common::Vector&lt;Int&gt;<br/><div style="display: inline" id="span:i240" class="foldit"><div class="descr"><p> This gives a list of the vertices of the curve in terms of their valences
 They appear in the same order as in <a href="#i241" onclick="unfold('span:i241');">NODES_BY_LEAVES</a> or <a href="#i238" onclick="unfold('span:i238');">NODES_BY_SETS</a></p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i163" onclick="swap_content('span:i163'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i163" onclick="swap_content('span:i163'); return false;"><span class="category">Input properties</span></a></h3></div><div style="display: inline" id="span:i163" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i244" onclick="swap_content('span:i244'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i244" onclick="swap_content('span:i244'); return false;"><b>INPUT_COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i244" class="foldit"><div class="descr"><p> Same as <a href="#i235" onclick="unfold('span:i235');">COEFFS</a>, except that entries may be &lt;=0. This should have the same length as <a href="#i243" onclick="unfold('span:i243');">INPUT_SETS</a>.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i243" onclick="swap_content('span:i243'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i243" onclick="swap_content('span:i243'); return false;"><b>INPUT_SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i243" class="foldit"><div class="descr"><p> Same as <a href="#i236" onclick="unfold('span:i236');">SETS</a>, except that sets may appear several times.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i242" onclick="swap_content('span:i242'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i242" onclick="swap_content('span:i242'); return false;"><b>INPUT_STRING</b>:
		</a>common::String<br/><div style="display: inline" id="span:i242" class="foldit"><div class="descr"><p> This property can also be used to define a rational curve: A linear combination of partitions is
 given as a string, using the following syntax:
 A partition is given as a subset of {1,..,n} and written as a comma-separated list of leaf
 indices in round brackets, e.g. "(1,2,5)"
 A linear combination can be created using rational numbers, "+","+" and "-" in the obvious way, 
 e.g. "2*(1,2,5) + 1*(3,4,7) - 2(1,2) (The "*" is optional)
 Of course, each set should contain at least two elements. If you don't specify N_LEAVES, it is
 set to be the largest leaf index occuring in the sets.
 Partitions needn't be irredundant and coefficients can be any rational number. If the resulting 
 element is not in the moduli space, an error is thrown.</p></div></div></div></div></li></ul></div></div></div></li></ul><h4>User Methods of RationalCurve</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i160" onclick="swap_content('span:i160'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i160" onclick="swap_content('span:i160'); return false;"><span class="category">Conversion</span></a></h3></div><div style="display: inline" id="span:i160" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i232" onclick="swap_content('span:i232'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i232" onclick="swap_content('span:i232'); return false;"><b>matroid_vector</b></a>
			()  

			
				→ Vector&lt;Rational&gt;<div style="display: inline" id="span:i232" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the coordinates of the curve in the moduli space given in matroid coordinates.
 In the isomorphism of the metric curve space and the moduli coordinates
 the last leaf is considered as the special leaf</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i233" onclick="swap_content('span:i233'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i233" onclick="swap_content('span:i233'); return false;"><b>metric_vector</b></a>
			()  

			<div style="display: inline" id="span:i233" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Returns the (n over 2) metric vector of the rational n-marked curve</p></div></div><div class="level3"/></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i159" onclick="swap_content('span:i159'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i159" onclick="swap_content('span:i159'); return false;"><span class="category">Visualization</span></a></h3></div><div style="display: inline" id="span:i159" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i231" onclick="swap_content('span:i231'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i231" onclick="swap_content('span:i231'); return false;"><b>VISUAL</b></a>
			()  

			<div style="display: inline" id="span:i231" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Visualizes a RationalCurve object. This visualization uses the VISUAL method of its GRAPH, so 
 it accepts all the options of Visual::Graph::decorations. In addition it has another option
 @options</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">LengthLabels</td><td><div class="descr"><div>If "hidden", the edges are not labelled with their lengths. Any other text is ignored. Not set to "hidden" by default.</div></div></td></tr><tr><td>option list:</td><td>Visual::RationalCurve::decorations</td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i29" onclick="swap_content('span:i29'); return false;">
		 	
	</div><div class="unfoldable"><div class="object"><h3><a href="#" class="javalink" id="i29" onclick="swap_content('span:i29'); return false;"><span class="object">RationalFunction</span></a></h3></div><div style="display: inline" id="span:i29" class="foldit"><div class="descr"><p> A rational function on a polyhedral complex. It is given by its <a href="#i250" onclick="unfold('span:i250');">DOMAIN</a>, a WeightedComplex
 and values on this domain, i.e. <a href="#i248" onclick="unfold('span:i248');">RAY_VALUES</a> and <a href="#i251" onclick="unfold('span:i251');">LIN_VALUES</a></p></div><div class="level3"/><h4>Properties of RationalFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i166" onclick="swap_content('span:i166'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i166" onclick="swap_content('span:i166'); return false;"><span class="category">Basic properties</span></a></h3></div><div style="display: inline" id="span:i166" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i250" onclick="swap_content('span:i250'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i250" onclick="swap_content('span:i250'); return false;"><b>DOMAIN</b>:
		</a><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a><br/><div style="display: inline" id="span:i250" class="foldit"><div class="descr"><p> This property describes the linearity domains of the function. I.e. the function is
 affine integral linear on each maximal cone of <a href="#i250" onclick="unfold('span:i250');">DOMAIN</a>. It is of type WeightedComplex, but
 actually only the basic geometric properties <a href="#" onclick="unfold('span:');">RAYS</a>/<a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>,<a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>/<a href="#i201" onclick="unfold('span:i201');">CMPLX_MAXIMAL_CONES</a>,<a href="#" onclick="unfold('span:');">LINEALITY_SPACE</a> and <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> are of any interest. Further properties are ignored and can be omitted upon generation.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i251" onclick="swap_content('span:i251'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i251" onclick="swap_content('span:i251'); return false;"><b>LIN_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i251" class="foldit"><div class="descr"><p> The value at index i describes the function value of <a href="#i250" onclick="unfold('span:i250');">DOMAIN</a>-&gt;LINEALITY_SPACE-&gt;row(i)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i249" onclick="swap_content('span:i249'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i249" onclick="swap_content('span:i249'); return false;"><b>POWER</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i249" class="foldit"><div class="descr"><p> This property describes how many times a function should be applied, if it occurs in a divisor term
 E.g. if you create a MinMaxFunction f and you set its POWER property to 2, a call to
 divisor(X,f) actually computes f*f*X</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i248" onclick="swap_content('span:i248'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i248" onclick="swap_content('span:i248'); return false;"><b>RAY_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:i248" class="foldit"><div class="descr"><p> The value at index i describes the function value of ray <a href="#i250" onclick="unfold('span:i250');">DOMAIN</a>-&gt;CMPLX_RAYS-&gt;row(i). More precisely, if the
 corresponding ray is a vertex, it describes its function value. If it is a directional ray, it describes 
 the slope on that ray.</p></div></div></div></div></li></ul></div></div></div></li></ul><h4>User Methods of RationalFunction</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i165" onclick="swap_content('span:i165'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i165" onclick="swap_content('span:i165'); return false;"><span class="category">Rational functions</span></a></h3></div><div style="display: inline" id="span:i165" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i247" onclick="swap_content('span:i247'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i247" onclick="swap_content('span:i247'); return false;"><b>homogenize</b></a>
			()  

			<div style="display: inline" id="span:i247" class="foldit"><br/><div class="descr_func"><div class="descr"><p> If the domain of the function is in homogeneous coordinates, it returns the function. Otherwise it homogenizes
 the domain and assigns the value 0 to the additional vertex at the origin</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i246" onclick="swap_content('span:i246'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i246" onclick="swap_content('span:i246'); return false;"><b>pow</b></a>
			(k)  

			
				→ Returns<div style="display: inline" id="span:i246" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This creates a RationalFunction from this one, simply by changing the power</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td>Returns</td><td><div class="descr"><div>the function f with a power set to k. Note that the original power of f is ignored.</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i164" onclick="swap_content('span:i164'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i164" onclick="swap_content('span:i164'); return false;"><span class="category">Visualization</span></a></h3></div><div style="display: inline" id="span:i164" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i245" onclick="swap_content('span:i245'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i245" onclick="swap_content('span:i245'); return false;"><b>BB_VISUAL</b></a>
			()  

			<div style="display: inline" id="span:i245" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Visualizes the domain of the function. Works exactly as BB_VISUAL of WeightedComplex, but has additional option
 @options</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">FunctionLabels</td><td><div class="descr"><div>If set to "show", textual function representations are diplayed on cones. False by default</div></div></td></tr><tr><td>option list:</td><td>Visual::WeightedComplex::FunctionDecorations</td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i26" onclick="swap_content('span:i26'); return false;">
		 	
	</div><div class="unfoldable"><div class="object"><h3><a href="#" class="javalink" id="i26" onclick="swap_content('span:i26'); return false;"><span class="object">WeightedComplex</span></a></h3></div><div style="display: inline" id="span:i26" class="foldit"><div class="descr"><p> A weighted complex, derived from fan::PolyhedralFan. Complexes are realized by
 using homogeneous coordinates (i.e. adding a 1 to the front of all vertices and a 0 to all rays)
 and setting USES_HOMOGENEOUS_C to true</p></div><b>derived from: </b>fan::PolyhedralFan&lt;Rational&gt;<br/><div class="level3"/><h4>Properties of WeightedComplex</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i151" onclick="swap_content('span:i151'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i151" onclick="swap_content('span:i151'); return false;"><span class="category">Input properties</span></a></h3></div><div style="display: inline" id="span:i151" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i186" onclick="swap_content('span:i186'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i186" onclick="swap_content('span:i186'); return false;"><b>INPUT_HOM_LINEALITY</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i186" class="foldit"><div class="descr"><p> USE IS NOT RECOMMENDED, since <a href="#i203" onclick="unfold('span:i203');">TROPICAL_WEIGHTS</a> are always given with respect to <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a> 
 anyway. 
 A matrix of generators of the lineality space, but in homogeneous coordinates (i.e. x0 should be 0, otherwise 
 this throws an error). This should be used *instead of* <a href="#" onclick="unfold('span:');">INPUT_LINEALITY</a> / <a href="#" onclick="unfold('span:');">LINEALITY_SPACE</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i187" onclick="swap_content('span:i187'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i187" onclick="swap_content('span:i187'); return false;"><b>INPUT_HOM_RAYS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i187" class="foldit"><div class="descr"><p> USE IS NOT RECOMMENDED, since <a href="#i203" onclick="unfold('span:i203');">TROPICAL_WEIGHTS</a> are always given with respect to <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a> 
 anyway. 
 A matrix of the rays of the fan. This should be used *instead of* <a href="#" onclick="unfold('span:');">INPUT_RAYS</a>, if the coordinates of the
 fan should be interpreted as homogeneous coordinates and the fan as the complex at x0 = 1.
 The input rays are normalized to x0 = 1 and the fan is then initialized with the normalized rays as
 <a href="#" onclick="unfold('span:');">INPUT_RAYS</a>.
 <a href="#" onclick="unfold('span:');">INPUT_CONES</a> can be combined with this as usual. If your rays are already homogenized and you are sure they
 are also your rays, you can directly put this into RAYS and manually set USES_HOMOGENEOUS_C=&gt;TRUE</p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i155" onclick="swap_content('span:i155'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i155" onclick="swap_content('span:i155'); return false;"><span class="category">Lattice properties</span></a></h3></div><div style="display: inline" id="span:i155" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i217" onclick="swap_content('span:i217'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i217" onclick="swap_content('span:i217'); return false;"><b>IS_UNIMODULAR</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i217" class="foldit"><div class="descr"><p> This property is not computed anywhere but should be set to TRUE by the user upon 
 initialization, if all cones of the fan are unimodular in the sense that 
 ANY linearly independent subset of their rays generates 
 the lattice of the cone. Note that this property is only considered if <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> is false.
 If <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> is true, this property is implicitly assumed to be false.
 If this property is TRUE, <a href="#i215" onclick="unfold('span:i215');">LATTICE_GENERATORS</a> and <a href="#i212" onclick="unfold('span:i212');">LATTICE_BASES</a> can be computed
 much faster.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i212" onclick="swap_content('span:i212'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i212" onclick="swap_content('span:i212'); return false;"><b>LATTICE_BASES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i212" class="foldit"><div class="descr"><p> This incidence matrix gives a lattice basis for each maximal cone. More precisely it 
 gives a lattice basis whose span contains the lattice of the maximal cone. Row i 
 corresponds to cone i and gives lattice generator indices referring to
 <a href="#i215" onclick="unfold('span:i215');">LATTICE_GENERATORS</a>
 If this property is computed via rules, it does indeed give a lattice basis for the
 cone lattice, but when it is computed during an operation like refinement or divisor
 it will in general be larger.
 If this property exists, lattice normals might be computed faster.
 If you create a variety with the property <a href="#i217" onclick="unfold('span:i217');">IS_UNIMODULAR</a> set to TRUE, then
 this will simply be a copy of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i215" onclick="swap_content('span:i215'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i215" onclick="swap_content('span:i215'); return false;"><b>LATTICE_GENERATORS</b>:
		</a>common::Matrix&lt;Integer, NonSymmetric&gt;<br/><div style="display: inline" id="span:i215" class="foldit"><div class="descr"><p> This is an irredundant list of all lattice generators of all maximal cones
 If this property exists, lattice normals might be computed faster
 If you create a variety with the property <a href="#i217" onclick="unfold('span:i217');">IS_UNIMODULAR</a> set to TRUE, then
 this will simply be a copy of <a href="#" onclick="unfold('span:');">RAYS</a> (made primitive integer)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i211" onclick="swap_content('span:i211'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i211" onclick="swap_content('span:i211'); return false;"><b>LATTICE_NORMALS</b>:
		</a>common::Map&lt;Int, Map&lt;Int, Vector&lt;Integer&gt;&gt;&gt;<br/><div style="display: inline" id="span:i211" class="foldit"><div class="descr"><p> A list of lattice normal vectors. They are associated to the faces and maximal cones in the following 
 way: For codimension one face i and adjacent maximal cone j, the element <a href="#i211" onclick="unfold('span:i211');">LATTICE_NORMALS</a>-&gt;{i}-&gt;{j} is the corresp.
 lattice normal vector. There is a convenience method lattice_normal(i,j)
 If <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> is TRUE, this actually computes a lattice normal of the cone (face i)\cap {x0 =1} wrt the
 cone (maximal cone j) \cap {x0=1} to ensure that the normal is of the form (0,...)
 If there is a <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a>, only lattice normals for compatible codim one faces
 are computed (naturally, since <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a> only contains those)</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i213" onclick="swap_content('span:i213'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i213" onclick="swap_content('span:i213'); return false;"><b>LATTICE_NORMAL_FCT_VECTOR</b>:
		</a>common::Map&lt;Int, Map&lt;Int, Vector&lt;Rational&gt;&gt;&gt;<br/><div style="display: inline" id="span:i213" class="foldit"><div class="descr"><p> For each lattice normal vector, this gives a vector of length (number of rays) + (lineality dim.), such that
 if a rational function is given by values on the rays and lin space generators, the value of the corresponding
 normal <a href="#i211" onclick="unfold('span:i211');">LATTICE_NORMALS</a>-&gt;{i}-&gt;{j} can be computed by multiplying the function value vector with 
 the vector <a href="#i213" onclick="unfold('span:i213');">LATTICE_NORMAL_FCT_VECTOR</a>-&gt;{i}-&gt;{j}. This is done in the following way:  If this fan does not 
 use homogeneous coordinates, it computes a representation of the lattice normal in the generating system 
 consisting of the rays of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>-&gt;row(j) and the generators 
 of the lineality space. It then inserts the coefficients of the representation at the right position. 
 If the fan uses homog. coordinates, however, we use a different generating system (and indices refer to CMPLX_RAYS)
 &lt;(r_i-r_0)_i&gt;0, s_j, l_k&gt;, where r_0 is the ray of the maximal cone 
 with the lowest index in <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>, such that it fulfills x0 = 1, r_i are the remaining rays with x0 = 1, ordered
 according to their index in <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>, s_j are the rays of the cone with x0 = 0 and l_k are the lineality space 
 generators. We will then store the coefficients a_i of (r_i - r_0) at the index of r_i, then - sum(a_i) at 
 the index of r_0 and the remaining coefficients at the appropriate places.
 In particular, the value of a lattice normal under a rational function can be computed
 simply by taking the scalar product of RAY_VALUES | LIN_VALUES with ist FCT_VECTOR</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i216" onclick="swap_content('span:i216'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i216" onclick="swap_content('span:i216'); return false;"><b>LATTICE_NORMAL_SUM</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i216" class="foldit"><div class="descr"><p> Rows of this matrix correspond to <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a>, and each row contains the weighted sum:
 sum_{cone &gt; codim-1-face}( weight(cone) * <a href="#i211" onclick="unfold('span:i211');">LATTICE_NORMALS</a>-&gt;{codim-1-face}-&gt;{cone})</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i214" onclick="swap_content('span:i214'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i214" onclick="swap_content('span:i214'); return false;"><b>LATTICE_NORMAL_SUM_FCT_VECTOR</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i214" class="foldit"><div class="descr"><p> Rows of this matrix correspond to <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a> (or <a href="#i198" onclick="unfold('span:i198');">CMPLX_CODIM_1_FACES</a> in the homog. case), and each row 
 contains a function vector for the corresponding row of <a href="#i216" onclick="unfold('span:i216');">LATTICE_NORMAL_SUM</a>. This function vector is computed 
 in the same way as described under <a href="#i213" onclick="unfold('span:i213');">LATTICE_NORMAL_FCT_VECTOR</a>.
 Note that for any codim-1-faces at which the complex is not balanced, the corresponding row is a zero row
 If a face is balanced can be checked under <a href="#i208" onclick="unfold('span:i208');">BALANCED_FACES</a>.</p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i152" onclick="swap_content('span:i152'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i152" onclick="swap_content('span:i152'); return false;"><span class="category">Polyhedral properties</span></a></h3></div><div style="display: inline" id="span:i152" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i200" onclick="swap_content('span:i200'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i200" onclick="swap_content('span:i200'); return false;"><b>CMPLX_AMBIENT_DIM</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i200" class="foldit"><div class="descr"><p> Returns the ambient dimension of the polyhedral complex. In the case of non-homogeneous coordinates, this is 
 just FAN_AMBIENT_DIM, otherwise FAN_AMBIENT_DIM-1</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i189" onclick="swap_content('span:i189'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i189" onclick="swap_content('span:i189'); return false;"><b>CMPLX_CODIMENSION</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i189" class="foldit"><div class="descr"><p> Returns the codimension of the polyhedral complex in its ambient vector space. More precisely, this
 value is equal to <a href="#i200" onclick="unfold('span:i200');">CMPLX_AMBIENT_DIM</a> - <a href="#i191" onclick="unfold('span:i191');">CMPLX_DIM</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i198" onclick="swap_content('span:i198'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i198" onclick="swap_content('span:i198'); return false;"><b>CMPLX_CODIM_1_FACES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i198" class="foldit"><div class="descr"><p> An incidence matrix describing which codim 1 cone in the complex at x0 = 1 is generated by which rays. 
 Each row corresponds to a codimension one cone (More precisely, the i-th element represents the
 same codim 1 cone as the i-th element of <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a>). The indices in a row refer to rows of <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>, i.e. the cone
 cone described by the i-th element is generated by the rays corresponding to these row indices.
 If the fan does not use homogeneous coordinates, it is just a copy of <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a>.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i199" onclick="swap_content('span:i199'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i199" onclick="swap_content('span:i199'); return false;"><b>CMPLX_CONES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i199" class="foldit"><div class="descr"><p> This incidence matrix describes the complete polyhedral complex, i.e. it contains all faces
 of the polyhedral complex</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i190" onclick="swap_content('span:i190'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i190" onclick="swap_content('span:i190'); return false;"><b>CMPLX_CONVERSION_VECTOR</b>:
		</a>common::Vector&lt;Int&gt;<br/><div style="display: inline" id="span:i190" class="foldit"><div class="descr"><p> A vector with an entry for each row in <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>. More precisely, the i-th entry gives the row index
 of the ray in <a href="#" onclick="unfold('span:');">RAYS</a> that is equal to the i-th row of <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i191" onclick="swap_content('span:i191'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i191" onclick="swap_content('span:i191'); return false;"><b>CMPLX_DIM</b>:
		</a>common::Int<br/><div style="display: inline" id="span:i191" class="foldit"><div class="descr"><p> Returns the dimension of the polyhedral complex. In the case of non-homogeneous coordinates, this is 
 just DIM, otherwise DIM-1. When computing this property, it is assumed that the complex is pure.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i201" onclick="swap_content('span:i201'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i201" onclick="swap_content('span:i201'); return false;"><b>CMPLX_MAXIMAL_CONES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i201" class="foldit"><div class="descr"><p> An incidence matrix describing which maximal cone in the complex at x0 = 1 is generated by which rays. 
 Each row corresponds to a maximal cone (More precisely, the i-th element represents the
 same maximal cone as the i-th element of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>). The indices in a row refer to rows of <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a>, i.e. the maximal
 cone described by the i-th element is generated by the rays corresponding to these row indices.
 If the fan does not use homogeneous coordinates, it is just a copy of <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i192" onclick="swap_content('span:i192'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i192" onclick="swap_content('span:i192'); return false;"><b>CMPLX_RAYS</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i192" class="foldit"><div class="descr"><p> If the fan does not use homogeneous coordinates, this is just a copy of RAYS and actually needn't be requested. Otherwise this computes
 a matrix of rays of the complex obtained by intersecting the fan with {x0 = 1}. More precisely, each ray r from <a href="#" onclick="unfold('span:');">RAYS</a> occurs  
 as a row in this matrix...
 - once, if r_0 = 1
 - k times, if r_0 = 0 and k is the number of equivalence classes of maximal cones containing r with respect to the following relation:
 Two maximal cones m, m' containing r are equivalent, if they are equal or there exists a sequence of maximal cones
 m = m_1,...m_r = m', such that r is contained in each m_i and each intersection m_i cap m_i+1 contains at least one ray s with s_0 = 1.
 The reason for this is that, when specifying a piecewise affine linear function on a polyhedral complex, the same directional ray with 
 x0 = 0 might be assigned two different values, if it is contained in two "non-connected" maximal cones (where connectedness is to be 
 understood as described above).
 If there is a <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a> the above equivalence relation is changed in such a 
 way that the affine ray s with s_0 = 1 that must be contained in the intersection of
 two subsequent cones must be a compatible ray</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i196" onclick="swap_content('span:i196'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i196" onclick="swap_content('span:i196'); return false;"><b>CODIM_1_FACES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i196" class="foldit"><div class="descr"><p> An incidence matrix of all codimension one faces. Each row is a face and column indices refer to <a href="#" onclick="unfold('span:');">RAYS</a>.
 If the fan uses homogeneous coordinates, only those faces that intersect {x0 = 1} are listed here.
 If there is a <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a>, only compatible codimension one faces are kept.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i193" onclick="swap_content('span:i193'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i193" onclick="swap_content('span:i193'); return false;"><b>CODIM_1_IN_MAXIMAL_CONES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i193" class="foldit"><div class="descr"><p> An incidence matrix indicating which codim 1 faces are contained in which maximal cone. Rows refer
 to (rows of) <a href="#i196" onclick="unfold('span:i196');">CODIM_1_FACES</a>, columns to (rows of) <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i195" onclick="swap_content('span:i195'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i195" onclick="swap_content('span:i195'); return false;"><b>DIRECTIONAL_RAYS</b>:
		</a>common::Set&lt;Int&gt;<br/><div style="display: inline" id="span:i195" class="foldit"><div class="descr"><p> Computes the subset of row indices of RAYS, such that the corresponding rows describe a directional ray of the complex.
 This set contains all row indices  if <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> is false and is the complement of <a href="#i194" onclick="unfold('span:i194');">VERTICES</a></p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i188" onclick="swap_content('span:i188'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i188" onclick="swap_content('span:i188'); return false;"><b>LOCAL_RESTRICTION</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:i188" class="foldit"><div class="descr"><p> This contains a list of sets of ray indices (referring to RAYS). All of these sets should 
 describe cones of the polyhedral complex (though not necessarily maximal ones).
 A cone is now called compatible with this property, if contains one of these cones
 If this list is not empty, all computations will be
 done only on (or around) compatible cones. The documentation of each property will explain
 in what way this restriction is enforced. 
 If this list is empty or not defined, there is no restriction.
 Careful: The implementation assumes that ALL maximal cones are compatible. If in doubt,
 you can create a complex with a local restriction from a given complex by using one of
 the "local_..." creation methods
 This list is assumed to be irredundant, i.e. there are no doubles (though this should not
 break anything, it is simply less efficient). It is, however, possible that one element
 is a subset of another.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i197" onclick="swap_content('span:i197'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i197" onclick="swap_content('span:i197'); return false;"><b>USES_HOMOGENEOUS_C</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i197" class="foldit"><div class="descr"><p> Indicates whethere this fan's coordinates should be interpreted as affine coordinates or as homogeneous coordinates.
 In the latter case the fan would actually represent the polyhedral complex that comes from intersecting the fan
 (in affine coordinates) with the affine hyperplane x0 = 1. This property should not be set manually. It will be set
 as true, if <a href="#i187" onclick="unfold('span:i187');">INPUT_HOM_RAYS</a> or <a href="#i186" onclick="unfold('span:i186');">INPUT_HOM_LINEALITY</a> are used and to false, if <a href="#" onclick="unfold('span:');">INPUT_RAYS</a> or <a href="#" onclick="unfold('span:');">INPUT_LINEALITY</a>
 are used.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i194" onclick="swap_content('span:i194'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i194" onclick="swap_content('span:i194'); return false;"><b>VERTICES</b>:
		</a>common::Set&lt;Int&gt;<br/><div style="display: inline" id="span:i194" class="foldit"><div class="descr"><p> Computes the subset of row indices of RAYS, such that the corresponding rows describe a vertex of the complex.
 This set is empty if <a href="#i197" onclick="unfold('span:i197');">USES_HOMOGENEOUS_C</a> is false and is the complement of <a href="#i195" onclick="unfold('span:i195');">DIRECTIONAL_RAYS</a></p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i153" onclick="swap_content('span:i153'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i153" onclick="swap_content('span:i153'); return false;"><span class="category">Tropical properties</span></a></h3></div><div style="display: inline" id="span:i153" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i208" onclick="swap_content('span:i208'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i208" onclick="swap_content('span:i208'); return false;"><b>BALANCED_FACES</b>:
		</a>common::Vector&lt;Bool&gt;<br/><div style="display: inline" id="span:i208" class="foldit"><div class="descr"><p> A vector whose entries correspond to the rows of CODIM_1_FACES. The i-th entry is true, if and only 
 if the complex is balanced at that face</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i204" onclick="swap_content('span:i204'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i204" onclick="swap_content('span:i204'); return false;"><b>DEGREE</b>:
		</a>common::Integer<br/><div style="display: inline" id="span:i204" class="foldit"><div class="descr"><p> The degree of the tropical variety, i.e. the weight of the intersection product with an 
 appropriate tropical linear space.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i207" onclick="swap_content('span:i207'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i207" onclick="swap_content('span:i207'); return false;"><b>IS_BALANCED</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i207" class="foldit"><div class="descr"><p> Indicates whether the fan is balanced with the given <a href="#i203" onclick="unfold('span:i203');">TROPICAL_WEIGHTS</a>. If there is a 
 <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a>, this will only be checked at the compatible codim one faces.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i209" onclick="swap_content('span:i209'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i209" onclick="swap_content('span:i209'); return false;"><b>IS_IRREDUCIBLE</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:i209" class="foldit"><div class="descr"><p> Whether this complex is irreducible.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i203" onclick="swap_content('span:i203'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i203" onclick="swap_content('span:i203'); return false;"><b>TROPICAL_WEIGHTS</b>:
		</a>common::Vector&lt;Integer&gt;<br/><div style="display: inline" id="span:i203" class="foldit"><div class="descr"><p> Vector of (integer) weights for maximal cones. Indices refer to (the rows of) <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>
 Note that, if the fan only consists of a lineality space, it specifies the weight of the single empty maximal
 cone </p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i205" onclick="swap_content('span:i205'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i205" onclick="swap_content('span:i205'); return false;"><b>WEIGHT_CONE</b>:
		</a>polytope::Cone&lt;Rational&gt;<br/><div style="display: inline" id="span:i205" class="foldit"><div class="descr"><p> The intersection of <a href="#i206" onclick="unfold('span:i206');">WEIGHT_SPACE</a> with the positive orthant.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i206" onclick="swap_content('span:i206'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i206" onclick="swap_content('span:i206'); return false;"><b>WEIGHT_SPACE</b>:
		</a>common::Matrix&lt;Integer, NonSymmetric&gt;<br/><div style="display: inline" id="span:i206" class="foldit"><div class="descr"><p> A Z-basis (as rows) for the space of weight distributions on this 
 tropical cycle making it balanced (i.e. this cycle is irreducible, if and only if 
 <a href="#i206" onclick="unfold('span:i206');">WEIGHT_SPACE</a> has only one row and the gcd of <a href="#i203" onclick="unfold('span:i203');">TROPICAL_WEIGHTS</a> is 1.</p></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i202" onclick="swap_content('span:i202'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i202" onclick="swap_content('span:i202'); return false;"><b>WEIGHT_SYSTEM</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:i202" class="foldit"><div class="descr"><p> The dual of <a href="#i206" onclick="unfold('span:i206');">WEIGHT_SPACE</a>.</p></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i154" onclick="swap_content('span:i154'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i154" onclick="swap_content('span:i154'); return false;"><span class="category">Unused</span></a></h3></div><div style="display: inline" id="span:i154" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i210" onclick="swap_content('span:i210'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i210" onclick="swap_content('span:i210'); return false;"><b>GLUING_FACES</b>:
		</a>common::Map&lt;Set&lt;Int&gt;, Set&lt;Int&gt;&gt;<br/><div style="display: inline" id="span:i210" class="foldit"><div class="descr"><p> Not used at the moment</p></div></div></div></div></li></ul></div></div></div></li></ul><h4>User Methods of WeightedComplex</h4><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i150" onclick="swap_content('span:i150'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i150" onclick="swap_content('span:i150'); return false;"><span class="category">Lattice properties</span></a></h3></div><div style="display: inline" id="span:i150" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i185" onclick="swap_content('span:i185'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i185" onclick="swap_content('span:i185'); return false;"><b>lattice_normal</b></a>
			(i, j)  

			
				→ The<div style="display: inline" id="span:i185" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Convenience method, returns LATTICE_NORMALS-&gt;{i}-&gt;{j}</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">i</td><td><div class="descr"><div>The index of the facet</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">j</td><td><div class="descr"><div>The index of the maximal cone</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>lattice normal vector of cone j wrt. facet i</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i147" onclick="swap_content('span:i147'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i147" onclick="swap_content('span:i147'); return false;"><span class="category">Polyhedral geometry</span></a></h3></div><div style="display: inline" id="span:i147" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i181" onclick="swap_content('span:i181'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i181" onclick="swap_content('span:i181'); return false;"><b>dehomogenize</b></a>
			()  

			<div style="display: inline" id="span:i181" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Assuming that this weighted complex is actually a fan, given in homogeneous coordinates, it will 
 convert it to non-homogeneous coordinates. This will return an error, if the complex has more than
 one vertex. Local restriction is also adapted.</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i182" onclick="swap_content('span:i182'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i182" onclick="swap_content('span:i182'); return false;"><b>homogenize</b></a>
			()  

			
				→ If<div style="display: inline" id="span:i182" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Transforms the fan into homogeneous coordinates</p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td>If</td><td><div class="descr"><div>this fan uses homog. coordinates, returns the fan itself, Otherwise it will replace the ray matrix
 by the following matrix: [0,ray] for any ray of the original fan and one affine ray [1,0,...0] appended. In <a href="#" onclick="unfold('span:');">MAXIMAL_CONES</a>
 all get the additional affine ray and <a href="#" onclick="unfold('span:');" class="invalid">TROPICAL WEIGHTS</a> is simply copied. Local restriction is adapted as well.</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i149" onclick="swap_content('span:i149'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i149" onclick="swap_content('span:i149'); return false;"><span class="category">Tropical geometry / Basic operations</span></a></h3></div><div style="display: inline" id="span:i149" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i184" onclick="swap_content('span:i184'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i184" onclick="swap_content('span:i184'); return false;"><b>invert</b></a>
			()  

			
				→ WeightedComplex<div style="display: inline" id="span:i184" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a WeightedComplex and multiplies all rays / vertices with -1. </p></div></div><div class="level3"><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>w The inverted complex</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i148" onclick="swap_content('span:i148'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i148" onclick="swap_content('span:i148'); return false;"><span class="category">Tropical geometry / Local computations</span></a></h3></div><div style="display: inline" id="span:i148" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i183" onclick="swap_content('span:i183'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i183" onclick="swap_content('span:i183'); return false;"><b>delocalize</b></a>
			(w)  

			
				→ WeightedComplex<div style="display: inline" id="span:i183" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a WeightedComplex and returns the same complex, but without any 
 <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a></p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">w</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td/></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i146" onclick="swap_content('span:i146'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i146" onclick="swap_content('span:i146'); return false;"><span class="category">Visualization</span></a></h3></div><div style="display: inline" id="span:i146" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i180" onclick="swap_content('span:i180'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i180" onclick="swap_content('span:i180'); return false;"><b>BB_VISUAL</b></a>
			()  

			<div style="display: inline" id="span:i180" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Displays a (possibly weighted) polyhedral complex by intersecting it with a bounding box. This bounding box
 is either defined by the vertices of the complex and the option "BoundingDistance" or explicitly given by
 "BoundingBox" and by setting "BoundingMode" to "absolute"
 @options</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">WeightLabels</td><td><div class="descr"><div>If "hidden", no weight labels are displayed. Not hidden by default.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">CoordLabels</td><td><div class="descr"><div>If "show", coordinate labels are displayed at vertices. Hidden by default.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">BoundingMode</td><td><div class="descr"><div>Can be "relative" (intersects with the bounding box returned by the method  boundingBox(BoundingDistance)) or "absolute" (intersects with the given BoundingBox). "relative" by default.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">BoundingDistance</td><td><div class="descr"><div>The distance parameter for relative bounding mode</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">BoundingBox</td><td><div class="descr"><div>The bounding parameter for absolute bounding mode</div></div></td></tr><tr><td>option list:</td><td>Visual::WeightedComplex::BoundingDecorations</td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i179" onclick="swap_content('span:i179'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i179" onclick="swap_content('span:i179'); return false;"><b>boundingBox</b></a>
			(distance)  

			
				→ Matrix&lt;Rational&gt;<div style="display: inline" id="span:i179" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a tropical variety and a (positive) Rational as input and computes the relative bounding box of the visualization with the Rational as BoundingDistance. I.e. will compute for each coordinate the minimum and
 maximum over all vertices and then add/substract the parameter distance. This is returned as a matrix. (see also the description of <a href="#i180" onclick="unfold('span:i180');">BB_VISUAL</a>),</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">distance</td><td><div class="descr"><div>The distance of the bounding box from the affine part of the variety (Can also be a standard perl number type)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td><div class="descr"><div>A 2x<a href="#i200" onclick="unfold('span:i200');">CMPLX_AMBIENT_DIM</a>-matrix specifying the bounding box (it can be used as input for BoundingBox in <a href="#i180" onclick="unfold('span:i180');">BB_VISUAL</a>)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i178" onclick="swap_content('span:i178'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i178" onclick="swap_content('span:i178'); return false;"><b>CMPLX_VISUAL</b></a>
			()  

			<div style="display: inline" id="span:i178" class="foldit"><br/><div class="descr_func"><div class="descr"><p> DEPRECATED. Use <a href="#i180" onclick="unfold('span:i180');">BB_VISUAL</a> instead.
 Displays a weighted polyhedral complex in ambient dimension up to 3 (4 in homog. coordinates) by visualizing all its
 cells in the following way: For each cell, every directional ray is added to every affine ray (which is only the origin
 in the non-homog. case) and the resulting polytope is displayed. The visualization has the following options:
 @options</p></div></div><div class="level3"><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">DirScale</td><td><div class="descr"><div>A rational number with which the directional rays are multiplied before added to the affine rays</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">WeightLabels</td><td><div class="descr"><div>if set to "hidden", the labels indicating the weights are hidden </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">CoordLabels</td><td><div class="descr"><div>If set to "show", the labels indicating the vertex coordinates are displayed, otherwise they are not. Note that this is expensive and significantly increases computation time.</div></div></td></tr><tr><td>option list:</td><td>Visual::WeightedComplex::decorations</td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div></div></div></li></ul></div><div class="level2"><h2>User Functions</h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i16" onclick="swap_content('span:i16'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i16" onclick="swap_content('span:i16'); return false;"><span class="category">Abstract rational curves</span></a></h3></div><div style="display: inline" id="span:i16" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i101" onclick="swap_content('span:i101'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i101" onclick="swap_content('span:i101'); return false;"><b>curve_graph_from_metric</b></a>
			(v)  

			
				→ An<div style="display: inline" id="span:i101" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^(n over 2) that describes an n-marked rational abstract
 curve as a distance vector between its leaves. It then computes the 
 graph of the curve corresponding to this vector.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector of length (n over 2). Its entries are 
 interpreted as the distances d(i,j) ordered lexicographically according to i,j. However, they need not be positive, as long as v is equivalent to a proper 
 metric modulo leaf lengths.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>An</td><td><div class="descr"><div>array containing first the graph::Graph and then a Vector&lt;Rational&gt;, containing 
 the lengths of the bounded edges (in the order they appear in EDGES)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i94" onclick="swap_content('span:i94'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i94" onclick="swap_content('span:i94'); return false;"><b>insert_leaves</b></a>
			(curve, nodes)  

			<div style="display: inline" id="span:i94" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a RationalCurve and a list of node indices. Then inserts additional 
 leaves (starting from N_LEAVES+1) at these nodes and returns the resulting 
 RationalCurve object</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td class="param">curve</td><td><div class="descr"><div>A RationalCurve object</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Int&gt;</a></td><td class="param">nodes</td><td><div class="descr"><div>A list of node indices of the curve</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i97" onclick="swap_content('span:i97'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i97" onclick="swap_content('span:i97'); return false;"><b>rational_curve_embedding</b></a>
			(delta, type)  

			
				→ WeightedComplex<div style="display: inline" id="span:i97" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function creates an embedding of a rational tropical curve using
 a given abstract curve and degree</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">delta</td><td><div class="descr"><div>The degree of the curve in non-homogeneous 
 coordinates. The number of rows
 should correspond to the number of leaves of type and the number of columns
 is the dimension of the space in which the curve should be realized</div></div></td></tr><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td class="param">type</td><td><div class="descr"><div>An abstract rational curve</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The corresponding embedded complex, weighted
 with 1. The position of the curve is determined by the first node, 
 which is always placed at the origin</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i100" onclick="swap_content('span:i100'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i100" onclick="swap_content('span:i100'); return false;"><b>rational_curve_from_cone</b></a>
			(X, n_leaves, coneIndex)  

			
				→ RationalCurve<div style="display: inline" id="span:i100" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a weighted complex X that is supposed to be of the form
 M_0,n x Y for some Y (It assumes that M_0,n occupies the first 
 coordinates) and an index of a maximal cone of that complex.
 It then computes a rational curve corresponding to an interior point of
 that cone (ignoring the second component Y)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A weighted complex of the form M_0,n x Y</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n_leaves</td><td><div class="descr"><div>The n in M_0,n. Needed to determine the dimension of
 the M_0,n component</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">coneIndex</td><td><div class="descr"><div>The index of the maximal cone</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td><div class="descr"><div>c The curve corresponding to an interior point</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i102" onclick="swap_content('span:i102'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i102" onclick="swap_content('span:i102'); return false;"><b>rational_curve_from_metric</b></a>
			(v)  

			
				→ RationalCurve<div style="display: inline" id="span:i102" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^(n over 2) that describes an n-marked rational abstract
 curve as a distance vector between its leaves. It then computes the 
 curve corresponding to this vector.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector of length (n over 2). Its entries are 
 interpreted as the distances d(i,j) ordered lexicographically according to i,j. However, they need not be positive, as long as v is equivalent to a proper 
 metric modulo leaf lengths.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i92" onclick="swap_content('span:i92'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i92" onclick="swap_content('span:i92'); return false;"><b>rational_curve_from_moduli</b></a>
			(v)  

			
				→ RationalCurve<div style="display: inline" id="span:i92" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a vector from Q^((n over 2) - n) that lies in M_0,n (in its matroid coordinates 
 and computes the corresponding rational curve.
 In the isomorphism of the metric curve space and the moduli coordinates
 the last leaf is considered as the special leaf</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector in the moduli space</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i93" onclick="swap_content('span:i93'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i93" onclick="swap_content('span:i93'); return false;"><b>rational_curve_from_rays</b></a>
			(rays, uses_homog)  

			
				→ RationalCurve<div style="display: inline" id="span:i93" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a matrix of rays of a given cone that is supposed to lie
 in a moduli space M_0,n and computes the rational curve corresponding
 to an interior point. More precisely, if there are k vertices in 
 homogeneous coordinates, it computes 1/k * (sum of these vertices),
 then it adds each directional ray. It then returns the curve corresponding
 to this point</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rays</td><td><div class="descr"><div>The rays of the cone, in homog. or non-homog.
 coordinates.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">uses_homog</td><td><div class="descr"><div>Whether the rays are given in homog. coordinates
 False by default</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td><div class="descr"><div>c The curve corresponding to an interior point</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i95" onclick="swap_content('span:i95'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i95" onclick="swap_content('span:i95'); return false;"><b>rational_curve_list_from_metric</b></a>
			(m)  

			
				→ RationalCurve<div style="display: inline" id="span:i95" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a matrix whose rows are metrics of rational n-marked curves.
 Returns a list, where the i-th element is the curve corr. to 
 the i-th row in the matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td><div class="descr"><div>: An array of RationalCurves</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i98" onclick="swap_content('span:i98'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i98" onclick="swap_content('span:i98'); return false;"><b>rational_curve_list_from_moduli</b></a>
			(m)  

			
				→ RationalCurve<div style="display: inline" id="span:i98" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a matrix whose rows are elements in the moduli space M_0,n in matroid 
 coordinates. Returns a list, where the i-th element is the curve corr. to 
 the i-th row in the matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td><div class="descr"><div>: An array of RationalCurves</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i96" onclick="swap_content('span:i96'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i96" onclick="swap_content('span:i96'); return false;"><b>sum_curves</b></a>
			(An, v)  

			
				→ RationalCurve<div style="display: inline" id="span:i96" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function takes a vector of coefficients a_i and a list of RationalCurves c_i and computes
 sum(a_i * c_i). In particular, it also checks, whether the result lies in M_0,n. If not, it returns undef</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td class="param">An</td><td><div class="descr"><div>arbitrary list of RationalCurve objects</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A list of coefficients. Superfluous coefficients are ignored, missing ones replaced
 by +1(!)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td><div class="descr"><div>The linear combination of the curves defined by the coefficients or undef, if the result
 is not in M_0,n. The history of the operation is kept in INPUT_SETS and INPUT_COEFFS</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i99" onclick="swap_content('span:i99'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i99" onclick="swap_content('span:i99'); return false;"><b>testFourPointCondition</b></a>
			(v)  

			
				→ Int<div style="display: inline" id="span:i99" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a metric vector in Q^{(n over 2)} and checks whether it fulfills 
 the four-point condition, i.e. whether it lies in M_0,n. More precisely 
 it only needs to be equivalent to such a vector</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>The vector to be checked</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td><div class="descr"><div>A quadruple (array) of indices, where the four-point condition 
 is violated or an empty list, if the vector is indeed in M_0,n</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i12" onclick="swap_content('span:i12'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i12" onclick="swap_content('span:i12'); return false;"><span class="category">Basic polyhedral operations</span></a></h3></div><div style="display: inline" id="span:i12" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i66" onclick="swap_content('span:i66'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i66" onclick="swap_content('span:i66'); return false;"><b>affineTransformation</b></a>
			(complex, translate, matrix)  

			
				→ WeightedComplex<div style="display: inline" id="span:i66" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral complex and applies an affine linear transformation, given by a
 translate vector and a matrix. The method assumes the function is bijective (on the 
 complex) and preserves cones, i.e. it just applies the transformation to the rays and 
 lineality space and leaves the cones and weights unchanged.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>The complex to be transformed</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">translate</td><td><div class="descr"><div>A vector whose dimension should be equal to the column 
 dimension of the transformation matrix</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">matrix</td><td><div class="descr"><div>An integer matrix </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The transformed complex, in homogeneous coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i70" onclick="swap_content('span:i70'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i70" onclick="swap_content('span:i70'); return false;"><b>coarsen</b></a>
			(complex)  

			
				→ WeightedComplex<div style="display: inline" id="span:i70" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a tropical variety on which a coarsest polyhedral structure exists
 and computes this structure.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>A tropical variety which has a unique 
 coarsest polyhedral structre </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The corresponding coarse complex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i75" onclick="swap_content('span:i75'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i75" onclick="swap_content('span:i75'); return false;"><b>fan_decomposition</b></a>
			(complex)  

			
				→ perl::ListReturn<div style="display: inline" id="span:i75" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Take a polyhedral complex and returns a list of all the local vertex fans, 
 i.e. for each affine ray r, the list contains the fan Star_complex(r) 
 (in non-homogeneous coordinates)
 If the complex has a non-trivial <a href="#i188" onclick="unfold('span:i188');">LOCAL_RESTRICTION</a>, only the local fans
 at compatible vertices are computed</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>A tropical variety</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>perl::ListReturn</td><td><div class="descr"><div>A list of WeightedComplex objects in 
 non-homogeneous coordinates. The i-th complex corresponds to the i-th 
 affine ray ( vertex). If the complex is not in homogeneous coordinates, 
 the list contains just the complex itself </div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i64" onclick="swap_content('span:i64'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i64" onclick="swap_content('span:i64'); return false;"><b>insert_rays</b></a>
			(F, R)  

			
				→ WeightedComplex<div style="display: inline" id="span:i64" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral fan wihtout lineality and a list of rays and triangulates the fan
 such that it contains these rays</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A polyhedral fan (i.e. USES_HOMOGENEOUS_C = FALSE),
 possibly with weights but with no lineality space</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">R</td><td><div class="descr"><div>A list of normalized rays (as row vectors), which will 
 be added to the fan</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A triangulation of F that contains all the 
 original rays of F plus the ones in R</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i72" onclick="swap_content('span:i72'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i72" onclick="swap_content('span:i72'); return false;"><b>intersect_complete_fan</b></a>
			()  

			<div style="display: inline" id="span:i72" class="foldit"><br/><div class="descr_func"><div class="descr"><p> DEPRECATED. Use intersect_container instead.</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i71" onclick="swap_content('span:i71'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i71" onclick="swap_content('span:i71'); return false;"><b>intersect_container</b></a>
			(fan, container, forceLatticeComputation)  

			
				→ WeightedComplex<div style="display: inline" id="span:i71" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes two fans and computes the intersection of both. The function
 relies on the fact that the latter fan contains the first fan to 
 compute the refinement correctly
 The function copies <a href="#i203" onclick="unfold('span:i203');">TROPICAL_WEIGHTS</a> and <a href="#i212" onclick="unfold('span:i212');">LATTICE_BASES</a>
 if they exist</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>An arbitrary weighted polyhedral fan</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">fan::PolyhedralFan</a></td><td class="param">container</td><td><div class="descr"><div>A polyhedral fan containing the 
 first one (as a set)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">forceLatticeComputation</td><td><div class="descr"><div>Whether the properties
 <a href="#i212" onclick="unfold('span:i212');">LATTICE_BASES</a> and <a href="#i215" onclick="unfold('span:i215');">LATTICE_GENERATORS</a> of fan should be computed
 before refining. False by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection of both fans (whose support is equal to the support of fan). The 
 resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, 
 the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i67" onclick="swap_content('span:i67'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i67" onclick="swap_content('span:i67'); return false;"><b>product_complex</b></a>
			(A)  

			
				→ the<div style="display: inline" id="span:i67" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the cartesian product of a set of polyhedral fans. If any of them uses homogeneous coordinates, so will the result</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>complexes</td><td class="param">A</td><td><div class="descr"><div>list of WeightedComplex objects</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>the</td><td><div class="descr"><div>cartesian product of the complexes. If any of the complexes use homogeneous coordinates, so does the result. If any of the complexex carries weights, the result carries the sum of the weights of its factors (Complexes without weights are treated as complexes with weight 1). If you want to compute the properties <a href="#i212" onclick="unfold('span:i212');">LATTICE_BASES</a> and <a href="#i215" onclick="unfold('span:i215');">LATTICE_GENERATORS</a> as cartesian products of the corr. properties of the complexes, use product_complex_lattice instead.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i74" onclick="swap_content('span:i74'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i74" onclick="swap_content('span:i74'); return false;"><b>product_complex_lattice</b></a>
			(A)  

			
				→ the<div style="display: inline" id="span:i74" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the cartesian product of a set of polyhedral fans. If any of them uses homogeneous coordinates, so will the result</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>complexes</td><td class="param">A</td><td><div class="descr"><div>list of WeightedComplex objects</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>the</td><td><div class="descr"><div>cartesian product of the complexes. If any of the complexes use homogeneous coordinates, so does the result. If any of the complexex carries weights, the result carries the sum of the weights of its factors (Complexes without weights are treated as complexes with weight 1). This function will precompute the properties <a href="#i212" onclick="unfold('span:i212');">LATTICE_BASES</a> and <a href="#i215" onclick="unfold('span:i215');">LATTICE_GENERATORS</a> for  all complexes and the result will have the corr. properties set to the cartesian product of these</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i65" onclick="swap_content('span:i65'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i65" onclick="swap_content('span:i65'); return false;"><b>recession_fan</b></a>
			(complex)  

			
				→ WeightedComplex<div style="display: inline" id="span:i65" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the recession fan of a tropical variety. WARNING: This is a highly experimental
 function. If it works at all, it is likely to take a very long time.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>A tropical variety. If it is a fan, the complex itself is returned</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A tropical fan, the recession fan of the complex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i68" onclick="swap_content('span:i68'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i68" onclick="swap_content('span:i68'); return false;"><b>skeleton_complex</b></a>
			(fan, k, preserveRays)  

			
				→ The<div style="display: inline" id="span:i68" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral complex and computes the k-skeleton. Will return an empty fan, if k is larger then the dimension of the given complex or smaller than 0.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A fan (or polyhedral complex)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the skeleton that should be computed</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">preserveRays</td><td><div class="descr"><div>When true, the function assumes that all rays of the fan remain in the k-skeleton, so it just copies the RAYS, instead of computing an irredundant list. This property can always be set to true, if fan is not in homogeneous coordinates or if the corresponding complex at x0 = 1 only has vertices. By default, this property is false.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>k-skeleton of the fan (or complex, if USES_HOMOGENEOUS_C is true)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i73" onclick="swap_content('span:i73'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i73" onclick="swap_content('span:i73'); return false;"><b>triangulateFan</b></a>
			(F)  

			
				→ WeightedComplex<div style="display: inline" id="span:i73" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a polyhedral fan and computes a triangulation</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A polyhedral fan (i.e. USES_HOMOGENEOUS_C = FALSE),
 possibly with weights</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A simplicial refinement of F</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i69" onclick="swap_content('span:i69'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i69" onclick="swap_content('span:i69'); return false;"><b>wfan</b></a>
			(fan, homog)  

			
				→ WeightedComplex<div style="display: inline" id="span:i69" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function takes a fan::PolyhedralFan and converts it into a WeightedComplex object.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">fan::PolyhedralFan</a></td><td class="param">fan</td><td><div class="descr"><div>The fan to be converted</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">homog</td><td><div class="descr"><div>Whether the fan is to be interpreted in homogeneous coordinates</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A complex with no weights but with the same geometry as fan.</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i19" onclick="swap_content('span:i19'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i19" onclick="swap_content('span:i19'); return false;"><span class="category">Basic polyhedral properties</span></a></h3></div><div style="display: inline" id="span:i19" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i117" onclick="swap_content('span:i117'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i117" onclick="swap_content('span:i117'); return false;"><b>check_cycle_equality</b></a>
			(X, Y, check_weights)  

			
				→ Bool<div style="display: inline" id="span:i117" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes two pure-dimensional polyhedral complexes and checks if they are equal
 i.e. if they have the same lineality space, the same rays (modulo lineality space)
 and the same cones. Optionally, it can also check if the weights are equal</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A weighted complex</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">check_weights</td><td><div class="descr"><div>Whether the algorithm should check for equality of weights. 
 This parameter is optional and true by default</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td><div class="descr"><div>Whether the cycles are equal</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i115" onclick="swap_content('span:i115'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i115" onclick="swap_content('span:i115'); return false;"><b>contains_point</b></a>
			(complex, point)  

			
				→ bool<div style="display: inline" id="span:i115" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a weighted complex and a point and computed whether that point lies in 
 the complex</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>perl::Object</td><td class="param">complex</td><td><div class="descr"><div>A weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">point</td><td><div class="descr"><div>An arbitrary vector in the same ambient
 dimension as complex</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>bool</td><td><div class="descr"><div>Whether the point lies in the support of complex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i116" onclick="swap_content('span:i116'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i116" onclick="swap_content('span:i116'); return false;"><b>is_zerocycle</b></a>
			(c)  

			
				→ Bool<div style="display: inline" id="span:i116" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Checks whether a cycle is the 0 cycle</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">c</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td/></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i17" onclick="swap_content('span:i17'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i17" onclick="swap_content('span:i17'); return false;"><span class="category">Creation functions for specific varieties</span></a></h3></div><div style="display: inline" id="span:i17" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i104" onclick="swap_content('span:i104'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i104" onclick="swap_content('span:i104'); return false;"><b>halfspace_complex</b></a>
			(constant, equation)  

			
				→ WeightedComplex<div style="display: inline" id="span:i104" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the halfspace complex defined by an rational vector g and a rational b, i.e. the 
 complex consisting of the two maximal cones g &gt;= a and g &lt;= a</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">constant</td><td><div class="descr"><div>The constant translation a</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">equation</td><td><div class="descr"><div>The defining equation g</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The resultin halfspace complex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i106" onclick="swap_content('span:i106'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i106" onclick="swap_content('span:i106'); return false;"><b>linear_nspace</b></a>
			(n, weight)  

			
				→ WeightedComplex<div style="display: inline" id="span:i106" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the tropical variety R^n as a single lineality space with a given weight </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The dimension of the linear space</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i110" onclick="swap_content('span:i110'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i110" onclick="swap_content('span:i110'); return false;"><b>linear_space_by_basis</b></a>
			(rowspace, weight)  

			
				→ WeightedComplex<div style="display: inline" id="span:i110" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates a tropical variety consisting of a single linear space defined by a matrix and assigns it a given weight.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rowspace</td><td><div class="descr"><div>A matrix whose rows span the lineality space of the fan (hence the whole fan) and
 which are assumed to be linearly independent. If you are not sure of this, call linear_space_by_matrix instead.</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i107" onclick="swap_content('span:i107'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i107" onclick="swap_content('span:i107'); return false;"><b>linear_space_by_matrix</b></a>
			(rowspace, weight)  

			
				→ WeightedComplex<div style="display: inline" id="span:i107" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates a tropical variety consisting of a single linear space defined by a matrix and assigns it a given weight.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">rowspace</td><td><div class="descr"><div>A matrix whose rows span the lineality space of the fan (hence the whole fan).
 If you already know that your rows are a basis, call linear_space_by_basis instead.</div></div></td></tr><tr><td>int</td><td class="param">weight</td><td><div class="descr"><div>The weight of the linear space (1, if none is given)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The resulting fan object in non-homog. coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i109" onclick="swap_content('span:i109'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i109" onclick="swap_content('span:i109'); return false;"><b>point_variety</b></a>
			(point, weight)  

			
				→ WeightedComplex<div style="display: inline" id="span:i109" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes an n-dimensional vector and creates a subdivision of R^n into
 2^n orthants with apex at the given vertex. This can for example
 be used to subdivide a given variety such that a certain vertex is 
 contained in its polyhedral structure by intersecting with this variety</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">point</td><td><div class="descr"><div>The apex of the variety in homogeneous coordinates</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">weight</td><td><div class="descr"><div>Optional. All cones will have this weight. 1 by default</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The variety described above</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i103" onclick="swap_content('span:i103'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i103" onclick="swap_content('span:i103'); return false;"><b>tropical_cube</b></a>
			(n, k)  

			
				→ WeightedComplex<div style="display: inline" id="span:i103" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the tropical cube T^n_k, i.e. the tropical variety obtained by
 glueing together the 2^n L^n_k obtained by applying all possible sign 
 changes, such that the vertices form the k-skeleton of 
 the n-dimensional cube
 Alternatively, T^n_k can be seen as
 max(1,x_1,...,x_n,-x_1,...,-x_n)^(n-k) * R^n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The ambient dimension</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the cube</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>cube</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i108" onclick="swap_content('span:i108'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i108" onclick="swap_content('span:i108'); return false;"><b>tropical_lnk</b></a>
			(n, k)  

			
				→ WeightedComplex<div style="display: inline" id="span:i108" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the linear tropical space L^n_k. This tropical fan is defined in the following way: 
 As rays we take -e_i,i=1,...,n, where e_i is the i-th standard basis vector of R^n and 
 e_0 = e_1 + ... + e_n. As maximal cones we take the cones generated by rays {e_i, i in S}, where
 S runs over all k-element subsets of {0,..,n}.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The ambient dimension of the fan.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the fan (should be smaller equal n, otherwise an error is thrown).</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A PolyhedralFan object representing L^n_k</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i105" onclick="swap_content('span:i105'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i105" onclick="swap_content('span:i105'); return false;"><b>zero_cycle</b></a>
			()  

			<div style="display: inline" id="span:i105" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Returns the zero cycle, i.e. a weighted complex with no rays and cones</p></div></div><div class="level3"/></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i14" onclick="swap_content('span:i14'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i14" onclick="swap_content('span:i14'); return false;"><span class="category">Divisors</span></a></h3></div><div style="display: inline" id="span:i14" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i84" onclick="swap_content('span:i84'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i84" onclick="swap_content('span:i84'); return false;"><b>divisor</b></a>
			(X, f, k)  

			
				→ WeightedComplex<div style="display: inline" id="span:i84" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes as input a tropical variety X , a rational function f and optionally a positive integer k. It then computes
 f^k * X
Note that this still produces a meaningful result, if the WeightedComplex is not balanced: The "divisor" of a given function is computed by taking all codim-1-faces, at which the complex is balanced and computing weights there.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical variety</div></div></td></tr><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td class="param">f</td><td><div class="descr"><div>A rational function, defined on a DOMAIN that contains X. It can be in non-homog. coordinates, even if X is in homog. coordinates, but then its AMBIENT_DIM has to be equal to X-&gt;AMBIENT_DIM. Otherwise you will have to homogenize the function first.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>Optional. A positive integer. Note that if k is given, this will actually ignore f-&gt;POWER. Is 1 by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The divisor f^k * X. It is given in homogeneous coordinates, if and only if: Either X is already in homogeneous coordinates or f is defined on a homogeneous DOMAIN (which is implicitly the case for ALL MinMaxFunctions)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i83" onclick="swap_content('span:i83'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i83" onclick="swap_content('span:i83'); return false;"><b>divisorByPLF</b></a>
			(fan, A)  

			
				→ The<div style="display: inline" id="span:i83" class="foldit"><br/><div class="descr_func"><div class="descr"><p> NOTE: Deprecated. Use divisor(..) instead
 Computes the divisor of a MinMaxFunction on a given tropical variety. The result will be 
 in homogeneous coordinates, whether the tropical variety uses them or not. The function 
 should be given on the affine coordinates of the variety, NOT the homogeneous ones.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">fan</td><td><div class="descr"><div>A tropical variety, on which the divisor is computed</div></div></td></tr><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td class="param">A</td><td><div class="descr"><div>function whose DOMAIN should be equal to the affine coordinate 
 space of the variety, i.e. AMBIENT_DIM-1, if the variety uses homogeneous coordinates, 
 AMBIENT_DIM otherwise.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>corresponding divisor as a tropical variety in homogeneous coordinates.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i85" onclick="swap_content('span:i85'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i85" onclick="swap_content('span:i85'); return false;"><b>divisorByValue</b></a>
			(A, An)  

			
				→ The<div style="display: inline" id="span:i85" class="foldit"><br/><div class="descr_func"><div class="descr"><p> NOTE: Deprecated. Use divisor(...) or divisor_nr(...) in conjunction with function_value instead
 Takes as input a tropical fan / tropical variety and an array of rational values. The array length should coincide 
 with the number of <a href="#i192" onclick="unfold('span:i192');">CMPLX_RAYS</a> of the fan plus the dimension of the lineality space and will be interpreted as a rational
 function, where each value has been assigned to the rays given by $fan-&gt;CMPLX_RAYS and to the generators given by $fan-&gt;LINEALITY_SPACE  
 (in that order). Missing values will be filled up by 0's, superfluous ones will be ignored. 
 The function will then compute the corresponding Weil divisor and return it as a 
 tropical variety given as a fan. The fan uses homogeneous coordinates, if and only the input fan does.
 Note that this still produces a meaningful result, if the WeightedComplex is not balanced: The "divisor" of a given function is computed by taking all codim-1-faces, at which the complex is balanced and computing weights there.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">A</td><td><div class="descr"><div>tropical variety on which the divisor is computed. </div></div></td></tr><tr><td>values</td><td class="param">An</td><td><div class="descr"><div>array of rational values that define an integer affine map on the fan. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>divisor of the function defined by values on the given fan, as a tropical variety.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i82" onclick="swap_content('span:i82'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i82" onclick="swap_content('span:i82'); return false;"><b>divisor_list</b></a>
			(X, A)  

			
				→ WeightedComplex<div style="display: inline" id="span:i82" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes as input a tropical variety and an arbitrarily long list of RationalFunction objects f_1,...,f_k. It will 
# then compute f_k^a_k * ... * f_1^a_1 * X, where a_i = f_i-&gt;POWER</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td/></tr><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td class="param">A</td><td><div class="descr"><div>list f_1,...,f_k</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The divisor f_k * ... * f_1 * X. It is in homogeneous coordinates, if and only if: 
 Either X is in homog. coordinates or one of the functions is defined on a homogeneous domain.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i87" onclick="swap_content('span:i87'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i87" onclick="swap_content('span:i87'); return false;"><b>divisor_list_nr</b></a>
			()  

			<div style="display: inline" id="span:i87" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Works exactly as divisor_list(WeightedComplex;RationalFunction,...). Should be called
 ONLY,when every occuring function f is defined on a DOMAIN equal to X (in the sense that
 all properties like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not
 sufficient). In this case this function will in general be faster.</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i86" onclick="swap_content('span:i86'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i86" onclick="swap_content('span:i86'); return false;"><b>divisor_nr</b></a>
			()  

			<div style="display: inline" id="span:i86" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Works exactly as divisor(WeightedComplex, RationalFunction;Int). Should be called ONLY,
 when the function f is defined on a DOMAIN equal to X (in the sense that all properties
 like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not sufficient). In this
 case this function will in general be faster.</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i81" onclick="swap_content('span:i81'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i81" onclick="swap_content('span:i81'); return false;"><b>piecewise_divisor</b></a>
			(F, cones, coefficients)  

			
				→ WeightedComplex<div style="display: inline" id="span:i81" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a divisor of a linear sum of certain piecewise polynomials on a simplicial 
 fan. </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A simplicial fan without lineality space in non-homog.
 coordinates</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">IncidenceMatrix</a></td><td class="param">cones</td><td><div class="descr"><div>A list of cones of F (not maximal, but all of the same 
 dimension). Each cone t corresponds to a piecewise polynomial psi_t, defined by 
 subsequently applying the rational functions that are 1 one exactly one ray of t and 
 0 elsewhere. </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Integer&gt;</a></td><td class="param">coefficients</td><td><div class="descr"><div>A list of coefficients a_t corresponding to the 
 cones. </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The divisor sum_t a_t psi_t * F</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i80" onclick="swap_content('span:i80'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i80" onclick="swap_content('span:i80'); return false;"><b>random_minmaxdivisor</b></a>
			(ambient_dim, equations, upperCoeffBound, onlyFunction)  

			
				→ MinMaxFunction<div style="display: inline" id="span:i80" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates a divisor of a random MinMaxFunction (in n-space) with integer coefficients.
 By default this uses maximum for the function. It returns an array containing the function created and its
 divisor in n-space.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">ambient_dim</td><td><div class="descr"><div>The ambient dimension of the divisor</div></div></td></tr><tr><td>int</td><td class="param">equations</td><td><div class="descr"><div>The number of equations to use</div></div></td></tr><tr><td>int</td><td class="param">upperCoeffBound</td><td><div class="descr"><div>The upper bound  of the coefficients (All coefficients are positive)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">onlyFunction</td><td><div class="descr"><div>If true, only the random function is returned and the divisor is not computed. FALSE by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td><div class="descr"><div>if onlyFunction is TRUE, otherwise an array containing first the function, then the divisor</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i6" onclick="swap_content('span:i6'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i6" onclick="swap_content('span:i6'); return false;"><span class="category">External programs / Gfan</span></a></h3></div><div style="display: inline" id="span:i6" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i31" onclick="swap_content('span:i31'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i31" onclick="swap_content('span:i31'); return false;"><b>variety_from_gfan</b></a>
			(path)  

			
				→ WeightedComplex<div style="display: inline" id="span:i31" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes as input the name of a file containing a fan with or without multiplicities in standard gfan output 
 (not xml!) and parses it into a tropical variety</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">path</td><td><div class="descr"><div>The file name</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td/></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i15" onclick="swap_content('span:i15'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i15" onclick="swap_content('span:i15'); return false;"><span class="category">Hurwitz cycles</span></a></h3></div><div style="display: inline" id="span:i15" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i91" onclick="swap_content('span:i91'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i91" onclick="swap_content('span:i91'); return false;"><b>hurwitz_cycle</b></a>
			(k, degree, points)  

			
				→ WeightedComplex<div style="display: inline" id="span:i91" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function computes the Hurwitz cycle H_k(x), x = (x_1,...,x_n)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the Hurwitz cycle, i.e. the number of moving vertices</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Int&gt;</a></td><td class="param">degree</td><td><div class="descr"><div>The degree x. Should add up to 0</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">points</td><td><div class="descr"><div>Optional. Should have length n-3-k. Gives the images of 
 the fixed vertices (besides 0). If not given all fixed vertices are mapped to 0
 and the function computes the recession fan of H_k(x)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>H_k(x), in homogeneous coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i90" onclick="swap_content('span:i90'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i90" onclick="swap_content('span:i90'); return false;"><b>hurwitz_graph</b></a>
			(cycle)  

			
				→ graph::Graph<div style="display: inline" id="span:i90" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes as input a hurwitz curve and computes the corresponding graph object
 Directional rays are modeled as terminal vertices. Each vertex (including directional rays)
 is labeled with its combinatorial type</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">cycle</td><td><div class="descr"><div>A Hurwitz curve object (or any one-dimensional cycle in M_0,n)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">graph::Graph</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i88" onclick="swap_content('span:i88'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i88" onclick="swap_content('span:i88'); return false;"><b>hurwitz_pair</b></a>
			(k, degree, points, list)  

			
				→ An<div style="display: inline" id="span:i88" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function computes hurwitz_subdivision and hurwitz_cycle at the same time, 
 returning the result in an array</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the Hurwitz cycle, i.e. the number of moving vertices</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Int&gt;</a></td><td class="param">degree</td><td><div class="descr"><div>The degree x. Should add up to 0</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">points</td><td><div class="descr"><div>Optional. Should have length n-3-k. Gives the images of 
 the fixed vertices (besides 0). If not given all fixed vertices are mapped to 0
 and the function computes the subdivision of M_0,n containing the recession fan of H_k(x)</div></div></td></tr><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td class="param">list</td><td><div class="descr"><div>Optional. A list of rational curves. If given, the computation
 will be performed locally around the given combinatorial types.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>An</td><td><div class="descr"><div>array, containing first the subdivision of M_0,n, then the Hurwitz cycle</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i89" onclick="swap_content('span:i89'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i89" onclick="swap_content('span:i89'); return false;"><b>hurwitz_subdivision</b></a>
			(k, degree, points)  

			
				→ WeightedComplex<div style="display: inline" id="span:i89" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function computes a subdivision of M_0,n containing the Hurwitz cycle
 H_k(x), x = (x_1,...,x_n) as a subfan. If k = n-4, this subdivision is the unique
 coarsest subdivision fulfilling this property</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of the Hurwitz cycle, i.e. the number of moving vertices</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Int&gt;</a></td><td class="param">degree</td><td><div class="descr"><div>The degree x. Should add up to 0</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">points</td><td><div class="descr"><div>Optional. Should have length n-3-k. Gives the images of 
 the fixed vertices (besides 0). If not given all fixed vertices are mapped to 0
 and the function computes the subdivision of M_0,n containing the recession fan of H_k(x)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>A subdivision of M_0,n, in homogeneous coordinates</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i7" onclick="swap_content('span:i7'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i7" onclick="swap_content('span:i7'); return false;"><span class="category">Integer and lattice arithmetic</span></a></h3></div><div style="display: inline" id="span:i7" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i32" onclick="swap_content('span:i32'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i32" onclick="swap_content('span:i32'); return false;"><b>gcdext</b></a>
			(a, b, s, t)  

			
				→ Integer<div style="display: inline" id="span:i32" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the gcd of <em class="param">a</em> and <em class="param">b</em> and returns it. <em class="param">s</em> and <em class="param">t</em> are set
 such that gcd(<em class="param">a</em>,<em class="param">b</em>) = <em class="param">s</em> * <em class="param">a</em> + <em class="param">t</em> * <em class="param">b</em>
 This is just a wrapper for the GMP's gcdext method</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">a</td><td><div class="descr"><div>first argument of gcd(,)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">b</td><td><div class="descr"><div>second argument of gcd(,)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">s</td><td><div class="descr"><div>coefficient of a (will be set)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">t</td><td><div class="descr"><div>coefficient of b (will be set)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i36" onclick="swap_content('span:i36'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i36" onclick="swap_content('span:i36'); return false;"><b>latticeBasis</b></a>
			(cone)  

			
				→ Matrix&lt;Integer&gt;<div style="display: inline" id="span:i36" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a cone and computes a Z-basis of the vector space spanned by the cone,
 returned as row vectors of a matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">cone</td><td><div class="descr"><div>A cone for which a Z-basis is to be computed</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i34" onclick="swap_content('span:i34'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i34" onclick="swap_content('span:i34'); return false;"><b>latticeNormal</b></a>
			(taumatrix, sigmamatrix, additionalRay)  

			
				→ Vector<div style="display: inline" id="span:i34" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes two matrices whose rows define the dual of the linear span of cone tau and sigma.
 Assuming that tau is a codimension one face of sigma, computes a representative of
 the primitive lattice normal vector of sigma with respect to tau</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">taumatrix</td><td><div class="descr"><div>a codimension one face of sigma, given as a matrix defining its linear span</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">sigmamatrix</td><td><div class="descr"><div>an arbitrary cone, given as a matrix defining its linear span</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td class="param">additionalRay</td><td><div class="descr"><div>A ray that is contained in sigma, but not in tau. Used to calculate proper orientation of the normal vector.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i37" onclick="swap_content('span:i37'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i37" onclick="swap_content('span:i37'); return false;"><b>latticeNormalByCone</b></a>
			(tau, sigma)  

			
				→ Vector<div style="display: inline" id="span:i37" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Assuming that tau is a codimension one face of sigma, computes a representative of
 the primitive lattice normal vector of sigma with respect to tau</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">tau</td><td><div class="descr"><div>a codimension one face of tau, given as a cone</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">polytope::Cone</a></td><td class="param">sigma</td><td><div class="descr"><div>an arbitrary cone</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i33" onclick="swap_content('span:i33'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i33" onclick="swap_content('span:i33'); return false;"><b>lattice_index</b></a>
			(m)  

			
				→ Integer<div style="display: inline" id="span:i33" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This computes the index of a lattice of rank n in Z^n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">m</td><td><div class="descr"><div>A list of (row) generators of the lattice. The matrix must have
 full column rank, otherwise an error is thrown</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The index of the lattice</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i35" onclick="swap_content('span:i35'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i35" onclick="swap_content('span:i35'); return false;"><b>randomInteger</b></a>
			(max_arg, n)  

			
				→ Array&lt;Integer&gt;<div style="display: inline" id="span:i35" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Returns n random integers in the range 0.. (max_arg-1),inclusive
 Note that this algorithm is not optimal for real randomness:
 If you change the range parameter and then change it back, you will
 usually get the exact same sequence as the first time</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">max_arg</td><td><div class="descr"><div>The upper bound for the random integers</div></div></td></tr><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The number of integers to be created</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Array&lt;Integer&gt;</a></td><td/></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i11" onclick="swap_content('span:i11'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i11" onclick="swap_content('span:i11'); return false;"><span class="category">Intersection products</span></a></h3></div><div style="display: inline" id="span:i11" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i57" onclick="swap_content('span:i57'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i57" onclick="swap_content('span:i57'); return false;"><b>cycle_intersection</b></a>
			(X, Y)  

			
				→ WeightedComplex<div style="display: inline" id="span:i57" class="foldit"><br/><div class="descr_func"><div class="descr"><p> DEPRECATED: Use intersect(...)
 Computes the intersection product of two tropical cycles in a common ambient vector space V.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>The first tropical variety</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>The second tropical variety. Should have the same actual ambient 
 dimension (homogeneous coordinates don't count) as X.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>Z The intersection product X*Y in V = R^n, where n is the ambient 
 dimension of X and Y. The result has homogeneous coordinates in any case</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i54" onclick="swap_content('span:i54'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i54" onclick="swap_content('span:i54'); return false;"><b>degree</b></a>
			(complex)  

			
				→ Int<div style="display: inline" id="span:i54" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the property <a href="#i204" onclick="unfold('span:i204');">DEGREE</a> of a tropical variety as the degree of the 
 0-dimensional complex obtained when intersecting 
 the variety with an appropriate linear space L^n_k</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i58" onclick="swap_content('span:i58'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i58" onclick="swap_content('span:i58'); return false;"><b>diagonal_functions</b></a>
			(n)  

			
				→ MinMaxFunction<div style="display: inline" id="span:i58" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Produces the rational functions that cut out the diagonal of the vector space R^n.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>Them dimension of the vector space (and the number of functions)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td><div class="descr"><div>A list of the functions (max(x_n,y_n),...,max(x_1,y_1)), where x_i is the coordinate function in the first part of the cartesian product R^n x R^n and y_i is the coordinate function of the second half.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i60" onclick="swap_content('span:i60'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i60" onclick="swap_content('span:i60'); return false;"><b>diagonal_unk</b></a>
			(n, k)  

			
				→ RationalFunction<div style="display: inline" id="span:i60" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a list of functions that cut out the diagonal on B(U<sup>n</sup><sub>k</sub>) x B(U<sup>n</sup><sub>k</sub>) </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The ambient dimension of B(U<sup>n</sup><sub>k</sub>) </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension of B(U<sup>n</sup><sub>k</sub>) </div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td><div class="descr"><div>An array of functions, cutting out the diagonal</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i61" onclick="swap_content('span:i61'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i61" onclick="swap_content('span:i61'); return false;"><b>intersect</b></a>
			(X, Y)  

			
				→ WeightedComplex<div style="display: inline" id="span:i61" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the intersection product of two tropical cycles in R^n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical cycle</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A tropical cycle, living in the same space as X</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection product, always in homogeneous coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i56" onclick="swap_content('span:i56'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i56" onclick="swap_content('span:i56'); return false;"><b>intersect_by_diagonal</b></a>
			(X, Y, liftLineality, Verbose, A)  

			
				→ WeightedComplex<div style="display: inline" id="span:i56" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes an intersection product of two tropical cycles by using a prescribed list of rational functions
 cutting out the diagonal</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical cycle</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A tropical cycle living in the same ambient space</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">liftLineality</td><td><div class="descr"><div>Whether both cycles have to be lifted into an ambient space of one dimension higher. If true, this will be done by adding a zero at the end of RAYS and LINEALITY_SPACE and adding the LINEALITY generator (1,..,1). After computation, this additional coordinate will be stripped </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">Verbose</td><td><div class="descr"><div>Whether some progress information should be printed</div></div></td></tr><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td class="param">A</td><td><div class="descr"><div>list of rational functions cutting out the diagonal. These will be applied to the 
 cartesian product X x Y</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection product defined via the diagonal functions</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i63" onclick="swap_content('span:i63'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i63" onclick="swap_content('span:i63'); return false;"><b>intersect_in_lnk</b></a>
			(X, Y, k, verbose)  

			
				→ WeightedComplex<div style="display: inline" id="span:i63" class="foldit"><br/><div class="descr_func"><div class="descr"><p> # Computes the intersection product of two cycles living in L^n_k </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A cycle in L^n_k </div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A cycle in L^n_k </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension k of the ambient L^n_k </div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">verbose</td><td><div class="descr"><div>Optional. Whether some progress information should be printed</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection product of X and Y in B(U_n,k)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i55" onclick="swap_content('span:i55'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i55" onclick="swap_content('span:i55'); return false;"><b>intersect_in_unk</b></a>
			(X, Y, k, verbose)  

			
				→ WeightedComplex<div style="display: inline" id="span:i55" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the intersection product of two cycles living in the bergman fan B(U_n,k) (given in max-coordinates)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A cycle in B(U_n,k)</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A cycle in B(U_n,k)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The dimension k of the ambient U_n,k</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">verbose</td><td><div class="descr"><div>Optional. Whether some progress information should be printed</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection product of X and Y in B(U_n,k)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i59" onclick="swap_content('span:i59'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i59" onclick="swap_content('span:i59'); return false;"><b>matroid_intersection_by_flats</b></a>
			(X, Y, M)  

			
				→ WeightedComplex<div style="display: inline" id="span:i59" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes an intersection product of two cycles in a bergman fan</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical cycle living either in B(M) or B(M)/L, where
 L = &lt;(1,...,1)&gt;</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A tropical cycle, also living in B(M) or B(M)/L (both 
 cycles should of course live in the same fan)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">M</td><td><div class="descr"><div>A matroid. X and Y live in B(M) or B(M)/L. The function 
 will detect automatically which of these cases applies</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The intersection product of X and Y in B(M) (possibly mod L)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i62" onclick="swap_content('span:i62'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i62" onclick="swap_content('span:i62'); return false;"><b>minkowski_multiplicity</b></a>
			(X, Y)  

			
				→ Integer<div style="display: inline" id="span:i62" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the minkowski multiplicity of two fans, i.e. it computes the unique weight of 
 the Minkowski sum X + (-Y). More precisely: If the sum is not full-dimensional, it 
 returns 0. Otherwise it runs over all pairs of cones whose Minkowski sum contains
 a generic vector and adds the product of their weight times the lattice index of the sum
 of the lattices.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>A tropical fan</div></div></td></tr><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">Y</td><td><div class="descr"><div>A tropical fan, living in the same space as X</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The Minkowski multiplicity</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i10" onclick="swap_content('span:i10'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i10" onclick="swap_content('span:i10'); return false;"><span class="category">Inverse problems</span></a></h3></div><div style="display: inline" id="span:i10" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i52" onclick="swap_content('span:i52'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i52" onclick="swap_content('span:i52'); return false;"><b>cutting_functions</b></a>
			(F, weight_aim)  

			
				→ Matrix&lt;Rational&gt;<div style="display: inline" id="span:i52" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a weighted complex and a list of desired weights on its codimension one
 faces and computes all possible rational functions on (this subdivision of )
 the complex</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A tropical variety</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Integer&gt;</a></td><td class="param">weight_aim</td><td><div class="descr"><div>A list of weights, whose length should be equal to the number of CODIM_1_FACES. Gives the desired weight on each 
 codimension one face</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td><div class="descr"><div>The space of rational functions defined on this 
 particular subdivision. Each row is a generator. The columns correspond to 
 values on CMPLX_RAYS and LINEALITY_SPACE, except the last one, which is either 0 (then this 
 function cuts out zero and can be added to any solution) or non-zero (then 
 normalizing this entry to -1 gives a function cutting out the desired weights 
 on the codimension one skeleton</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i53" onclick="swap_content('span:i53'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i53" onclick="swap_content('span:i53'); return false;"><b>simplicial_diagonal_system</b></a>
			()  

			<div style="display: inline" id="span:i53" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function computes the inhomogeneous version of simplicial_piecewise_system
 in the sense that it computes the result of the above mentioned function (i.e. 
 which coefficients for the piecewise polynomials yield the zero divisor)
 and adds another column at the end where only the entries corresponding to the 
 diagonal cones are 1, the rest is zero. This can be seen as asking for a 
 solution to the system that cuts out the diagonal (all solutions whose last entry is 1)</p></div></div><div class="level3"/></div></div></div></li><li><div class="li"><div class="icon" id="icon:i50" onclick="swap_content('span:i50'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i50" onclick="swap_content('span:i50'); return false;"><b>simplicial_piecewise_system</b></a>
			(F)  

			
				→ Matrix&lt;Rational&gt;<div style="display: inline" id="span:i50" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function takes a d-dimensional simplicial fan F and computes the linear system 
 defined in the following way: For each d-dimensional cone t in the diagonal subdivision of FxF, let psi_t be the 
 piecewise polynomial defined by subsequently applying the rational functions that 
 are 1 one exactly one ray of t and 0 elsewhere. Now for which coefficients a_t
 is sum_t a_t psi_t * (FxF) = 0?</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A simplicial fan without lineality space in non-homog.
 coordinates</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td><div class="descr"><div>The above mentioned linear system. The rows 
 are equations, the columns correspond to d-dimensional cones of FxF in the order given 
 by skeleton_complex(simplicial_with_diagonal(F), d,TRUE)</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i51" onclick="swap_content('span:i51'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i51" onclick="swap_content('span:i51'); return false;"><b>simplicial_with_diagonal</b></a>
			(F)  

			
				→ WeightedComplex<div style="display: inline" id="span:i51" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function takes a simplicial fan F (in non-homog. coordinates and without 
 lineality space) and computes the coarsest subdivision of F x F containing all 
 diagonal rays (r,r)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">F</td><td><div class="descr"><div>A simplicial fan without lineality space in non-homog.
 coordinates</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The product complex FxF subdivided such that it contains 
 all diagonal rays</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i18" onclick="swap_content('span:i18'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i18" onclick="swap_content('span:i18'); return false;"><span class="category">Linear algebra</span></a></h3></div><div style="display: inline" id="span:i18" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i114" onclick="swap_content('span:i114'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i114" onclick="swap_content('span:i114'); return false;"><b>computeMatrixBases</b></a>
			(m)  

			
				→ IncidenceMatrix<div style="display: inline" id="span:i114" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a list of sets of column indices of a matrix such that
 the corresponding column sets form a basis of the column space</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">IncidenceMatrix</a></td><td><div class="descr"><div>A list of all bases</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i113" onclick="swap_content('span:i113'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i113" onclick="swap_content('span:i113'); return false;"><b>computeMatrixColoops</b></a>
			(m)  

			
				→ Set&lt;int&gt;<div style="display: inline" id="span:i113" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a list of column indices of a matrix such that 
 the corresponding columns are contained in every column basis of the
 column space. Is equal to the intersection of all sets returned by
 <a href="#i114" onclick="unfold('span:i114');">computeMatrixBases</a></p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Set&lt;int&gt;</a></td><td><div class="descr"><div>A list of all coloops </div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i112" onclick="swap_content('span:i112'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i112" onclick="swap_content('span:i112'); return false;"><b>linearRepresentation</b></a>
			(v, generators)  

			
				→ Vector<div style="display: inline" id="span:i112" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a vector v and a matrix with column dimension equal to the dimension of v. Assuming that 
 v is in the row span of the matrix, it computes one(!) possible representation of v in these 
 generators</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td class="param">v</td><td><div class="descr"><div>The vector supposed to be contained in the row span of the generators</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td class="param">generators</td><td><div class="descr"><div>A set of row vectors whose linear span should contain v</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector</a></td><td><div class="descr"><div>A vector (a1,..,an) such that v = (a1,...,an) * generators. It returns a vector of
 dimension 0, if v is not in the span of the generators. An error is thrown if the dimensions of 
 v and the generators mismatch</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i111" onclick="swap_content('span:i111'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i111" onclick="swap_content('span:i111'); return false;"><b>lllHNF</b></a>
			(matrix, tfmatrix, kdim)  

			
				→ Matrix<div style="display: inline" id="span:i111" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a row HNF of an integer matrix A
 It returns the normal form and stores the unimodular transformation matrix and the kernel dimension  of the transposed matrix in
 the last two parameters. The algorithm is the LLL-based HNF alg. by Havas, Majevski, Matthews</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">matrix</td><td><div class="descr"><div>the matrix for which the transformation is computed</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">tfmatrix</td><td><div class="descr"><div>The matrix that will contain the transformation matrix, i.e. tfmatrix * matrix = lllHNF(matrix)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td class="param">kdim</td><td><div class="descr"><div>This will be set to dim Ker(A^T) (This is needed for lattice normal computation)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i140" onclick="swap_content('span:i140'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i140" onclick="swap_content('span:i140'); return false;"><b>lllHNF</b></a>
			(matrix)  

			
				→ Matrix&lt;Integer&gt;<div style="display: inline" id="span:i140" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Does exactly the same as lllHNF(Matrix,Matrix,Integer), except 
 that it only returns the normal form and not the transformation matrix
 and the transposed kernel dimension</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td class="param">matrix</td><td><div class="descr"><div>the matrix for which the transformation is computed</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Integer&gt;</a></td><td><div class="descr"><div>The HNF of matrix</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i8" onclick="swap_content('span:i8'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i8" onclick="swap_content('span:i8'); return false;"><span class="category">Local computations</span></a></h3></div><div style="display: inline" id="span:i8" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i39" onclick="swap_content('span:i39'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i39" onclick="swap_content('span:i39'); return false;"><b>local_codim_1</b></a>
			(complex, face)  

			
				→ WeightedComplex<div style="display: inline" id="span:i39" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a weighted complex and an index of one of its codimension one faces
 (The index is in CODIM_1_FACES) and computes the complex locally restricted
 to that face</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>An arbitrary weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">face</td><td><div class="descr"><div>An index of a face in CODIM_1_FACES</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The complex locally restricted to the given face</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i38" onclick="swap_content('span:i38'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i38" onclick="swap_content('span:i38'); return false;"><b>local_point</b></a>
			(complex, v)  

			
				→ WeightedComplex<div style="display: inline" id="span:i38" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a weighted complex and an arbitrary vertex in homogeneous 
 coordinates that is supposed to lie in the support of the complex
 It then refines the complex such that the vertex is a cell in the polyhedral 
 structure and returns the complex localized at this vertex</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>An arbitrary weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vertex in homogeneous coordinates. It should lie
 in the support of the complex (otherwise an error is thrown)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The complex localized at the vertex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i41" onclick="swap_content('span:i41'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i41" onclick="swap_content('span:i41'); return false;"><b>local_ray</b></a>
			(complex, ray)  

			
				→ WeightedComplex<div style="display: inline" id="span:i41" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a weighted complex and an index of one of its vertices or rays 
 (the index is to be understood in RAYS)
 It then localizes the variety at this vertex / ray. The index should never
 correspond to a directional ray in a complex, since this would not be a cone</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>An arbitrary weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">ray</td><td><div class="descr"><div>The index of a ray/vertex in RAYS</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The complex locally restricted to the given ray/vertex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i40" onclick="swap_content('span:i40'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i40" onclick="swap_content('span:i40'); return false;"><b>local_restrict</b></a>
			(complex, cones)  

			
				→ WeightedComplex<div style="display: inline" id="span:i40" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a tropical variety and an IncidenceMatrix describing a set
 of cones (not necessarily maximal ones) of this variety. It will then
 create a variety that contains all compatible maximal cones and is
 locally restricted to the given cone set.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">complex</td><td><div class="descr"><div>An arbitrary weighted complex</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">IncidenceMatrix</a></td><td class="param">cones</td><td><div class="descr"><div>A set of cones, indices refer to RAYS</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The same complex, locally restricted to the given
 cones</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i13" onclick="swap_content('span:i13'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i13" onclick="swap_content('span:i13'); return false;"><span class="category">Matroid fans</span></a></h3></div><div style="display: inline" id="span:i13" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i77" onclick="swap_content('span:i77'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i77" onclick="swap_content('span:i77'); return false;"><b>bergman_fan</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<div style="display: inline" id="span:i77" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the bergman fan of a given matroid fan. This is much slower than using bergman_fan_matroid. The difference in output is that is always gives the coarsest polyhedral structure (but this can also be achieved by calling coarsen(...) after creation).</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">m</td><td><div class="descr"><div>A matroid</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i76" onclick="swap_content('span:i76'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i76" onclick="swap_content('span:i76'); return false;"><b>bergman_fan_flats</b></a>
			(m, mod_out_lineality)  

			
				→ WeightedComplex<div style="display: inline" id="span:i76" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function computes the bergman fan of a matroid in the subdivision associated
 to its lattice of flats. In particular it computes all flats and all maximal
 chains of flats. Hence this function is terribly inefficient and should only be used
 with comparatively small matroids</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">m</td><td><div class="descr"><div>Any matroid</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">mod_out_lineality</td><td><div class="descr"><div>Optional. Whether the lineality space should be 
 modded out by setting the last coordinate to be the sum of the remaining ones. True 
 by default.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The bergman fan B(m) in its flat subdivision</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i79" onclick="swap_content('span:i79'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i79" onclick="swap_content('span:i79'); return false;"><b>bergman_fan_linear</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<div style="display: inline" id="span:i79" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Uses an algorithm by Felipe Rincón to compute the bergman fan of a linear matroid</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">m</td><td><div class="descr"><div>Any rational matrix</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates. This algorithm uses min-convention, so the rays are given as (0,1) - vectors</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i78" onclick="swap_content('span:i78'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i78" onclick="swap_content('span:i78'); return false;"><b>bergman_fan_matroid</b></a>
			(m, modOutLineality, projectionCoordinate)  

			
				→ The<div style="display: inline" id="span:i78" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Uses an algorithm by Felipe Rincón to compute the bergman fan of a general matroid</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">matroid::Matroid</a></td><td class="param">m</td><td><div class="descr"><div>Any matroid</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">modOutLineality</td><td><div class="descr"><div>Optional argument. If set to TRUE, the lineality space is divided out before returning the 
 fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE</div></div></td></tr><tr><td>int</td><td class="param">projectionCoordinate</td><td><div class="descr"><div>Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is the last coordinate n-1.</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td>The</td><td><div class="descr"><div>bergman fan of the given matroid in non-homog. coordinates. This algorithm uses min-convention, so the rays are given as (0,1) - vectors</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i9" onclick="swap_content('span:i9'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i9" onclick="swap_content('span:i9'); return false;"><span class="category">Morphisms and functions</span></a></h3></div><div style="display: inline" id="span:i9" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i46" onclick="swap_content('span:i46'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i46" onclick="swap_content('span:i46'); return false;"><b>add_minmax_functions</b></a>
			(f, g, reduce)  

			
				→ MinMaxFunction<div style="display: inline" id="span:i46" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the sum of two MinMaxFunctions that are supposed to have the same USES_MIN property. You can also simply use '+' instead.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td class="param">f</td><td><div class="descr"><div>An arbitrary MinMaxFunction</div></div></td></tr><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td class="param">g</td><td><div class="descr"><div>A MinMaxFunction such that f-&gt;<a href="#i224" onclick="unfold('span:i224');">USES_MIN</a> == g-&gt;<a href="#i224" onclick="unfold('span:i224');">USES_MIN</a> and is
 defined on the same domain.</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Bool</a></td><td class="param">reduce</td><td><div class="descr"><div>Optional. False by default. If true, the function reduces the amount of terms by computing vertices of the newton polytope of the sum</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td><div class="descr"><div>The sum of both functions. </div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i47" onclick="swap_content('span:i47'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i47" onclick="swap_content('span:i47'); return false;"><b>function_mult</b></a>
			(f, r)  

			
				→ RationalFunction<div style="display: inline" id="span:i47" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This function takes a RationalFunction f and a Rational value r and returns the RationalFunction r*f. You can
 also use the operator '*' instead.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td class="param">f</td><td/></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">r</td><td/></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i44" onclick="swap_content('span:i44'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i44" onclick="swap_content('span:i44'); return false;"><b>function_point</b></a>
			(v)  

			
				→ RationalFunction<div style="display: inline" id="span:i44" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This takes a vector v = (v_1,...,v_n) in R^n and returns the list of rational functions
 max(v_i,x_i), i.e. a list of functions cutting out that point in R^n. </p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Rational&gt;</a></td><td class="param">v</td><td><div class="descr"><div>A vector in R^n (in non-homog. coordinates!)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td><div class="descr"><div>An array of RationalFunction objects, which cut out v in R^n</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i43" onclick="swap_content('span:i43'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i43" onclick="swap_content('span:i43'); return false;"><b>function_value</b></a>
			(X, An)  

			
				→ RationalFunction<div style="display: inline" id="span:i43" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This creates a RationalFunction by directly defining it on the rays of a given domain. It is simply shorter then
 using the constructor.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td class="param">X</td><td><div class="descr"><div>The DOMAIN of the function</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">An</td><td><div class="descr"><div>arbitrarily long list of function values. It should be as long as 
 X-&gt;CMPXL_RAYS-&gt;rows + X-&gt;LINEALITY_DIM. Function values will be assigned in this order 
 (first rays, then linspace). Superfluous values are ignored, missing ones replaced by zero</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td/></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i45" onclick="swap_content('span:i45'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i45" onclick="swap_content('span:i45'); return false;"><b>input_vars</b></a>
			(A)  

			
				→ Array&lt;String&gt;<div style="display: inline" id="span:i45" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates an array of strings from a single string by separating along commas. This can be used as a convenience 
 function to pass a list of variables to INPUT_VARS of MinMaxFunction, e.g. ...INPUT_VARS=&gt;input_vars("x,y,z")</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">A</td><td><div class="descr"><div>single string containing a comma-separated list of variable names</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Array&lt;String&gt;</a></td><td><div class="descr"><div>The list of the variable names</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i42" onclick="swap_content('span:i42'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i42" onclick="swap_content('span:i42'); return false;"><b>projection_map</b></a>
			(n, s)  

			
				→ Morphism<div style="display: inline" id="span:i42" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This creates a linear projection from R^n to a given set of coordinates</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The dimension of the domain</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Set&lt;Int&gt;</a></td><td class="param">s</td><td><div class="descr"><div>The set of coordinates to which this map should project (starting 
 the count at 0)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>The corresponding projection as a global linear map</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i138" onclick="swap_content('span:i138'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i138" onclick="swap_content('span:i138'); return false;"><b>projection_map</b></a>
			(n, m)  

			
				→ Morphism<div style="display: inline" id="span:i138" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This computes the projection from R^n to R^m (for m &lt; n) onto the first m coordinates</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>Dimension of domain</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">m</td><td><div class="descr"><div>Dimension of image</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>The corresponding projection</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i48" onclick="swap_content('span:i48'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i48" onclick="swap_content('span:i48'); return false;"><b>pull_back</b></a>
			(f, phi)  

			
				→ RationalFunction<div style="display: inline" id="span:i48" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the pull back of a rational function (or a list of rational functions) along a morphism</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td class="param">f</td><td><div class="descr"><div>A morphism</div></div></td></tr><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td class="param">phi</td><td><div class="descr"><div>A rational function whose <a href="#i176" onclick="unfold('span:i176');">DOMAIN</a> contains the image of f (or a list of RationalFunction objects)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i29" onclick="unfold('span:i29');">RationalFunction</a></td><td><div class="descr"><div>The pull-back of phi along f. If a list of functions was given, it returns an array with the pull-backs of these functions</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i49" onclick="swap_content('span:i49'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i49" onclick="swap_content('span:i49'); return false;"><b>scale_minmax_function</b></a>
			(f, a)  

			
				→ MinMaxFunction<div style="display: inline" id="span:i49" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Scales a MinMaxFunction by a given Rational a. You can also use the '*' operator.</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td class="param">f</td><td><div class="descr"><div>An arbitrary MinMaxFunction</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Rational</a></td><td class="param">a</td><td><div class="descr"><div>A scalar values</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i27" onclick="unfold('span:i27');">MinMaxFunction</a></td><td><div class="descr"><div>The scaled function</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li><li><div class="li"><div class="icon" id="icon:i20" onclick="swap_content('span:i20'); return false;">
		 	
	</div><div class="unfoldable"><div class="category"><h3><a href="#" class="javalink" id="i20" onclick="swap_content('span:i20'); return false;"><span class="category">The moduli space M_0,n</span></a></h3></div><div style="display: inline" id="span:i20" class="foldit"><ul><li><div class="li"><div class="icon" id="icon:i120" onclick="swap_content('span:i120'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i120" onclick="swap_content('span:i120'); return false;"><b>count_mn_cones</b></a>
			(n, k)  

			
				→ Integer<div style="display: inline" id="span:i120" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the number of k-dimensional cones of the tropical moduli space M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be &gt;= 3</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The number of bounded edges. This argument is optional and n-3 by default</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The number of k-dimensional cones of M_0,n</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i122" onclick="swap_content('span:i122'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i122" onclick="swap_content('span:i122'); return false;"><b>count_mn_rays</b></a>
			(n)  

			
				→ Integer<div style="display: inline" id="span:i122" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the number of rays of the tropical moduli space M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td>int</td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be &gt;= 3</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Integer</a></td><td><div class="descr"><div>The number of rays</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i126" onclick="swap_content('span:i126'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i126" onclick="swap_content('span:i126'); return false;"><b>evaluation_map</b></a>
			(n, r, Delta, i)  

			
				→ Morphism<div style="display: inline" id="span:i126" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This creates the i-th evaluation function on M_0,n^(lab)(R^r,Delta)
 (which is actually realized as M_0,(n+|Delta|) x R^r)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of marked points</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">r</td><td><div class="descr"><div>The dimension of the embedding space</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">Delta</td><td><div class="descr"><div>The directions of the unbounded edges (given as row vectors)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">i</td><td><div class="descr"><div>The index of the marked point that should be evaluated. i 
 should lie in between 1 and n
 Note that the i-th marked point is realized as the |Delta|+i-th leaf in M_0,(n+|Delta|)
 and that the R^r - coordinate is interpreted as the position of the n-th leaf. 
 In particular, ev_n is just the projection to the R^r-coordinates</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>ev_i. Its domain is the ambient space of the moduli space 
 in matroid coordinates cross R^r</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i142" onclick="swap_content('span:i142'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i142" onclick="swap_content('span:i142'); return false;"><b>evaluation_map</b></a>
			(n, r, d, i)  

			
				→ Morphism<div style="display: inline" id="span:i142" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This creates the i-th evaluation function on M_0,n^(lab)(R^r,d)
 (which is actually realized as M_0,(n+d(r+1)) x R^r)
 This is the same as calling the function
 evaluation_map(Int,Int,Matrix&lt;Rational&gt;,Int) with the standard d-fold
 degree as matrix</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of marked points</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">r</td><td><div class="descr"><div>The dimension of the embedding space</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">d</td><td><div class="descr"><div>The degree of the embedding. The direction matrix will be
 the standard d-fold directions (first (1,...1), then the inverted unit vectors)</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">i</td><td><div class="descr"><div>The index of the marked point that should be evaluated. i 
 should lie in between 1 and n</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>ev_i. Its domain is the ambient space of the moduli space 
 in matroid coordinates cross R^r</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i123" onclick="swap_content('span:i123'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i123" onclick="swap_content('span:i123'); return false;"><b>forgetful_map</b></a>
			(n, S)  

			
				→ Morphism<div style="display: inline" id="span:i123" class="foldit"><br/><div class="descr_func"><div class="descr"><p> This computes the forgetful map from the moduli space M_0,n to M_0,(n-|S|)</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves in the moduli space M_0,n</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Set&lt;Int&gt;</a></td><td class="param">S</td><td><div class="descr"><div>The set of leaves to be forgotten. Should be a subset of (1,..,n)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i25" onclick="unfold('span:i25');">Morphism</a></td><td><div class="descr"><div>The forgetful map. It will identify the remaining leaves 
 i_1,..,i_(n-|S|) with the leaves of M_0,(n-|S|) in order of their size</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i125" onclick="swap_content('span:i125'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i125" onclick="swap_content('span:i125'); return false;"><b>local_m0n</b></a>
			(A)  

			
				→ WeightedComplex<div style="display: inline" id="span:i125" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the moduli space M_0,n locally around a given list of combinatorial
 types. More precisely: It computes the weighted complex consisting of all
 maximal cones containing any of the given combinatorial types and localizes 
 at these types 
 This should only be used for curves of small codimension. What the function 
 actually does, is that it combinatorially computes the cartesian products 
 of M_0,v's, where v runs over the possible valences of vertices in the curves
 For max(v) &lt;= 8 this should terminate in a reasonable time (depending on the 
 number of curves)
 The coordinates are the same that would be produced by the function 
 tropical_m0n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#i28" onclick="unfold('span:i28');">RationalCurve</a></td><td class="param">A</td><td><div class="descr"><div>list of rational curves (preferrably in the same M_0,n)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The local complex</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i121" onclick="swap_content('span:i121'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i121" onclick="swap_content('span:i121'); return false;"><b>moduliDimensionFromLength</b></a>
			(k)  

			
				→ Int<div style="display: inline" id="span:i121" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Takes a positive length (of a vector) and assumes it is of the form (n over 2)
 It then computes n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">k</td><td><div class="descr"><div>The length = (n over 2)</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td><div class="descr"><div>n. If k is not of the form (n over 2), an error is thrown.</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i119" onclick="swap_content('span:i119'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i119" onclick="swap_content('span:i119'); return false;"><b>psi_class</b></a>
			(n, i)  

			
				→ WeightedComplex<div style="display: inline" id="span:i119" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes the i-th psi class in the moduli space of n-marked rational tropical curves
 M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves in M_0,n</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">i</td><td><div class="descr"><div>The leaf for which we want to compute the psi class ( in 1,..,n )</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The corresponding psi class</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i118" onclick="swap_content('span:i118'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i118" onclick="swap_content('span:i118'); return false;"><b>psi_product</b></a>
			(n, exponents)  

			
				→ WeightedComplex<div style="display: inline" id="span:i118" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Computes a product of psi classes psi_1^k_1 * ... * psi_n^k_n on the moduli space
 of rational n-marked tropical curves M_0,n</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves in M_0,n</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Vector&lt;Int&gt;</a></td><td class="param">exponents</td><td><div class="descr"><div>The exponents of the psi classes k_1,..,k_n. If the 
 vector does not have length n or if some entries are negative, an error is thrown</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The corresponding psi class divisor</div></div></td></tr></table></div></div></div></div></li><li><div class="li"><div class="icon" id="icon:i124" onclick="swap_content('span:i124'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i124" onclick="swap_content('span:i124'); return false;"><b>tropical_m0n</b></a>
			(n)  

			
				→ WeightedComplex<div style="display: inline" id="span:i124" class="foldit"><br/><div class="descr_func"><div class="descr"><p> Creates the moduli space of abstract rational n-marked curves. Its coordinates are
 given as the coordinates of the bergman fan of the matroid of the complete graph on 
 n-1 nodes (but not computed as such)
 The lineality space (1,..,1) is modded out by setting the last coordinate
 equal to minus the sum of the others
 The isomorphism to the space of curve metrics is obtained by choosing
 the last leaf as special leaf</p></div></div><div class="level3"><h5>Parameters</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Int</a></td><td class="param">n</td><td><div class="descr"><div>The number of leaves. Should be at least 4</div></div></td></tr></table><h5>Returns</h5><table class="args"><tr><td><a href="#i26" onclick="unfold('span:i26');">WeightedComplex</a></td><td><div class="descr"><div>The tropical moduli space M_0,n</div></div></td></tr></table></div></div></div></div></li></ul></div></div></div></li></ul></div><div class="level3"><h2>Common Option Lists</h2><ul class="unfoldable"><li><div class="li"><div class="icon" id="icon:i21" onclick="swap_content('span:i21'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i21" onclick="swap_content('span:i21'); return false;"><b>Visual::RationalCurve::decorations</b></a><div style="display: inline" id="span:i21" class="foldit"><div class="descr"><p> Visualization options for RationalCurve-&gt;VISUAL</p></div><b>imports from: </b>Visual::Graph::decorations, Visual::Wire::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">LengthLabels</td><td><div class="descr"><div>if set to "hidden", the labels indicating the lengths are hidden</div></div></td></tr></table></div></div></div></li><li><div class="li"><div class="icon" id="icon:i22" onclick="swap_content('span:i22'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i22" onclick="swap_content('span:i22'); return false;"><b>Visual::WeightedComplex::BoundingDecorations</b></a><div style="display: inline" id="span:i22" class="foldit"><div class="descr"><p> Some options only used for BB_VISUAL</p></div><b>imports from: </b>Visual::WeightedComplex::decorations, Visual::Polygons::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Flexible&lt;Rational&gt;</a></td><td class="param">BoundingDistance</td><td><div class="descr"><div>The distance of the border of the bounding box from the smallest box 
 containing the affine points of the complex. This is only relevant, if BoundingMode is "relative"</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Matrix&lt;Rational&gt;</a></td><td class="param">BoundingBox</td><td><div class="descr"><div>A fixed bounding box, determined by two row vectors that specify two of its vertices (on "on top" and one "at the bottom"). Is only relevant, if BoundingMode is "absolute"</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">BoundingMode</td><td><div class="descr"><div>If set to "relative", the function determines the smallest possible box containing all affine points of the complex and then enlarges the box by BoundingDistance to all sides. If set to "absolute", BoundingBox must be specified and the complex will be intersected with that box. By default this is set to "relative".</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">Array&lt;String&gt;</a></td><td class="param">ConeLabels</td><td><div class="descr"><div>A list of strings to be displayed as labels for the maximal cones. If this is empty, the weight labels (if present and not suppressed by WeightLabels=&gt;"hidden") are displayed</div></div></td></tr></table></div></div></div></li><li><div class="li"><div class="icon" id="icon:i23" onclick="swap_content('span:i23'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i23" onclick="swap_content('span:i23'); return false;"><b>Visual::WeightedComplex::decorations</b></a><div style="display: inline" id="span:i23" class="foldit"><div class="descr"><p> A set of basic visualization options for CMPLX_VISUAL and BB_VISUAL</p></div><b>imports from: </b>Visual::Polygons::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">Flexible&lt;Rational&gt;</a></td><td class="param">DirScale</td><td><div class="descr"><div>length of the directional rays</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">WeightLabels</td><td><div class="descr"><div>if set to "hidden", the labels indicating the weights are hidden</div></div></td></tr><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">CoordLabels</td><td><div class="descr"><div>If set to "show", the labels indicating the vertex coordinates are displayed, otherwise they are not. Note that this is expensive and significantly increases computation time.
 String If set to "pinned", the weight labels are computed in such a way that they stick to the cells when the complex
 is exploded. Note that this is expensive and significantly increases computation time.
 WeightLabelStyle =&gt; "",</div></div></td></tr></table></div></div></div></li><li><div class="li"><div class="icon" id="icon:i24" onclick="swap_content('span:i24'); return false;">
		 	
	</div><div class="unfoldable"><a href="#" class="javalink" id="i24" onclick="swap_content('span:i24'); return false;"><b>Visual::WeightedComplex::FunctionDecorations</b></a><div style="display: inline" id="span:i24" class="foldit"><div class="descr"><p>Visualization options for RationalFunction/Morphism-&gt;BB_VISUAL</p></div><b>imports from: </b>Visual::WeightedComplex::BoundingDecorations, Visual::WeightedComplex::decorations, Visual::Polygons::decorations, Visual::PointSet::decorations<br/><br/><h5>Options</h5><table class="args"><tr><td><a href="#" onclick="unfold('span:');">String</a></td><td class="param">FunctionLabels</td><td><div class="descr"><div>, if set to "show", the function labels indicatingt the affine linear representation of each function
 on each cone are computed</div></div></td></tr></table></div></div></div></li></ul></div></div>
