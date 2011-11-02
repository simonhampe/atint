/*
$Id: functions.js 9563 2010-03-14 15:04:28Z smoeser $

This file contains the JavaScript functions used in the documentation of 
polymake.

Author: Silke Möser, 2010
*/

// make content appear or disappear upon clicking
function swap_content( span_id ) {
	displayType = ( document.getElementById( span_id ).style.display == 'none' ) ? 'inline' : 'none';
	styleImage = ( document.getElementById( span_id ).style.display == 'none' ) ? 'url(images/minus.png)' : 'url(images/plus.png)';
	document.getElementById( span_id ).style.display = displayType;
	var id; id=span_id.split(':')[1];
	var icon; icon='icon:'+id;
	document.getElementById( icon ).style.backgroundImage = styleImage;
}

// unfold element (and also all its ancestors)
function unfold( id ) {
	var obj; obj=document.getElementById( id );
	while (obj) {
		if (obj.tagName=='SPAN') {
			obj.style.display = 'inline';
			id = obj.id;
			var icon; icon='icon:'+id.split(':')[1];
			document.getElementById( icon ).style.backgroundImage='url(images/minus.png)';
		}
		obj=obj.parentNode;
	}
}

// executed upon loading 
// fold everything (if JavaScript is activated)
function start() {
	var list; list=document.getElementsByTagName("span");
	for(var i=0; i< list.length; i++){
  		var temp;
  		temp=list[i];
  		if(temp.id){
			temp.style.display = 'none';
		}
	}

	// jumps to anchor if there is one
	if (location.href.split('#')[1]) {
		var anchor; anchor = location.href.split('#')[1];
		var span; span='span:'+anchor;
		unfold(span);
	}
}
